/**
 * @file albacore.c
 * @brief Main file for ALBACORE library.
 * @author - SCEC 
 * @version 1.0
 *
 * @section DESCRIPTION
 *
 * Delivers ALBACORE, Offshore Southern California Velocity Model
 *
 */

#include "albacore.h"

/**
 * Initializes the ALBACORE plugin model within the UCVM framework. In order to initialize
 * the model, we must provide the UCVM install path and optionally a place in memory
 * where the model already exists.
 *
 * @param dir The directory in which UCVM has been installed.
 * @param label A unique identifier for the velocity model.
 * @return Success or failure, if initialization was successful.
 */
int albacore_init(const char *dir, const char *label) {
	int tempVal = 0;
	char configbuf[512];

	// Initialize variables.
	albacore_configuration = calloc(1, sizeof(albacore_configuration_t));
	albacore_velocity_model = calloc(1, sizeof(albacore_model_t));

	// Configuration file location.
	sprintf(configbuf, "%s/model/%s/data/config", dir, label);

	// Read the configuration file.
	if (albacore_read_configuration(configbuf, albacore_configuration) != SUCCESS)
		tempVal = FAIL;

	// Set up the data directory.
	sprintf(albacore_data_directory, "%s/model/%s/data/%s/", dir, label, albacore_configuration->model_dir);

	// Can we allocate the model, or parts of it, to memory. If so, we do.
	tempVal = albacore_try_reading_model(albacore_velocity_model);

	if (tempVal == SUCCESS) {
		fprintf(stderr, "WARNING: Could not load model into memory. Reading the model from the\n");
		fprintf(stderr, "hard disk may result in slow performance.");
	} else if (tempVal == FAIL) {
		print_error("No model file was found to read from.");
		return FAIL;
	}

	// In order to simplify our calculations in the query, we want to rotate the box so that the bottom-left
	// corner is at (0m,0m). Our box's height is total_height_m and total_width_m. We then rotate the
	// point so that is is somewhere between (0,0) and (total_width_m, total_height_m). How far along
	// the X and Y axis determines which grid points we use for the interpolation routine.

	albacore_total_height_m = sqrt(pow(albacore_configuration->top_left_corner_n - albacore_configuration->bottom_left_corner_n, 2.0f) +
						  pow(albacore_configuration->top_left_corner_e - albacore_configuration->bottom_left_corner_e, 2.0f));
	albacore_total_width_m  = sqrt(pow(albacore_configuration->top_right_corner_n - albacore_configuration->top_left_corner_n, 2.0f) +
						  pow(albacore_configuration->top_right_corner_e - albacore_configuration->top_left_corner_e, 2.0f));

	// Let everyone know that we are initialized and ready for business.
	albacore_is_initialized = 1;

	return SUCCESS;
}

/**
 * Queries ALBACORE at the given points and returns the data that it finds.
 *
 * @param points The points at which the queries will be made.
 * @param data The data that will be returned (Vp, Vs, density, Qs, and/or Qp).
 * @param numpoints The total number of points to query.
 * @return SUCCESS or FAIL.
 */
int albacore_query(albacore_point_t *points, albacore_properties_t *data, int numpoints) {
	int i = 0;
	int load_x_coord = 0, load_y_coord = 0, load_z_coord = 0;
	double x_percent = 0, y_percent = 0, z_percent = 0;
	albacore_properties_t surrounding_points[8];
        double lon_e, lat_n;

	int zone = 11;
	int longlat2utm = 0;

        double delta_lon = (albacore_configuration->top_right_corner_e - albacore_configuration->bottom_left_corner_e)/(albacore_configuration->nx - 1);
        double delta_lat = (albacore_configuration->top_right_corner_n - albacore_configuration->bottom_left_corner_n)/(albacore_configuration->ny - 1);

	for (i = 0; i < numpoints; i++) {
		lon_e = points[i].longitude; 
		lat_n = points[i].latitude; 

		// Which point base point does that correspond to?
		load_y_coord = (int)(round((lat_n - albacore_configuration->bottom_left_corner_n) / delta_lat));
		load_x_coord = (int)(round((lon_e - albacore_configuration->bottom_left_corner_e) / delta_lon));
		load_z_coord = (int)((points[i].depth)/1000);
//fprintf(stderr, "XXX y %d x %d z %d XXX \n", load_y_coord, load_x_coord, load_z_coord);

		// Are we outside the model's X and Y boundaries?
		if (load_x_coord > albacore_configuration->nx - 1 || load_y_coord > albacore_configuration->ny - 1 || load_x_coord < 0 || load_y_coord < 0) {
			data[i].vp = -1;
			data[i].vs = -1;
			data[i].rho = -1;
			data[i].qp = -1;
			data[i].qs = -1;
			continue;
		}

		// Get the X, Y, and Z percentages for the bilinear or trilinear interpolation below.
		x_percent =fmod((lon_e - albacore_configuration->bottom_left_corner_e), delta_lon) /delta_lon;
		y_percent = fmod((lat_n - albacore_configuration->bottom_left_corner_n), delta_lat)/
delta_lat;
		z_percent = fmod(points[i].depth, albacore_configuration->depth_interval) / albacore_configuration->depth_interval;
//fprintf(stderr, "XXX y_p %lf x_p  %lf z_p %lf XXX \n", y_percent, x_percent, z_percent);

		if (load_z_coord == 0 && z_percent == 0) {
			// We're below the model boundaries. Bilinearly interpolate the bottom plane and use that value.
			load_z_coord = 0;
                   if(albacore_configuration->interpolation) {

			// Get the four properties.
			albacore_read_properties(load_x_coord,     load_y_coord,     load_z_coord,     &(surrounding_points[0]));	// Orgin.
			albacore_read_properties(load_x_coord + 1, load_y_coord,     load_z_coord,     &(surrounding_points[1]));	// Orgin + 1x
			albacore_read_properties(load_x_coord,     load_y_coord + 1, load_z_coord,     &(surrounding_points[2]));	// Orgin + 1y
			albacore_read_properties(load_x_coord + 1, load_y_coord + 1, load_z_coord,     &(surrounding_points[3]));	// Orgin + x + y, forms top plane.

			albacore_bilinear_interpolation(x_percent, y_percent, surrounding_points, &(data[i]));
                  } else {
			albacore_read_properties(load_x_coord,     load_y_coord,     load_z_coord,     &(data[i]));	// Orgin.
                  }

		} else {
		  if( albacore_configuration->interpolation) {
			// Read all the surrounding point properties.
			albacore_read_properties(load_x_coord,     load_y_coord,     load_z_coord,     &(surrounding_points[0]));	// Orgin.
			albacore_read_properties(load_x_coord + 1, load_y_coord,     load_z_coord,     &(surrounding_points[1]));	// Orgin + 1x
			albacore_read_properties(load_x_coord,     load_y_coord + 1, load_z_coord,     &(surrounding_points[2]));	// Orgin + 1y
			albacore_read_properties(load_x_coord + 1, load_y_coord + 1, load_z_coord,     &(surrounding_points[3]));	// Orgin + x + y, forms top plane.
			albacore_read_properties(load_x_coord,     load_y_coord,     load_z_coord - 1, &(surrounding_points[4]));	// Bottom plane origin
			albacore_read_properties(load_x_coord + 1, load_y_coord,     load_z_coord - 1, &(surrounding_points[5]));	// +1x
			albacore_read_properties(load_x_coord,     load_y_coord + 1, load_z_coord - 1, &(surrounding_points[6]));	// +1y
			albacore_read_properties(load_x_coord + 1, load_y_coord + 1, load_z_coord - 1, &(surrounding_points[7]));	// +x +y, forms bottom plane.

			albacore_trilinear_interpolation(x_percent, y_percent, z_percent, surrounding_points, &(data[i]));
                    } else {
                        // no interpolation, data as it is
			albacore_read_properties(load_x_coord,     load_y_coord,     load_z_coord,     &(data[i]));	// Orgin.
                    }
		}

		//??? Calculate density.
		//data[i].rho = albacore_calculate_density(data[i].vs);
		data[i].qp = -1;
		data[i].qs = -1;
	}

	return SUCCESS;
}

/**
 * Retrieves the material properties (whatever is available) for the given data point, expressed
 * in x, y, and z co-ordinates.
 *
 * @param x The x coordinate of the data point.
 * @param y The y coordinate of the data point.
 * @param z The z coordinate of the data point.
 * @param data The properties struct to which the material properties will be written.
 */
void albacore_read_properties(int x, int y, int z, albacore_properties_t *data) {

	// Set everything to -1 to indicate not found.
	data->vp = -1;
	data->vs = -1;
	data->rho = -1;
	data->qp = -1;
	data->qs = -1;

	float *ptr = NULL;
	FILE *fp = NULL;

	int location = z * (albacore_configuration->nx * albacore_configuration->ny) + (y * albacore_configuration->nx) + x;

	// Check our loaded components of the model.
	if (albacore_velocity_model->vs_status == 2) {
		// Read from memory.
		ptr = (float *)albacore_velocity_model->vs;
		data->vs = ptr[location];
	} else if (albacore_velocity_model->vs_status == 1) {
		// Read from file.
		fp = (FILE *)albacore_velocity_model->vs;
		fseek(fp, location * sizeof(float), SEEK_SET);
		fread(&(data->vs), sizeof(float), 1, fp);
	}

	// Check our loaded components of the model.
	if (albacore_velocity_model->vp_status == 2) {
		// Read from memory.
		ptr = (float *)albacore_velocity_model->vp;
		data->vp = ptr[location];
	} else if (albacore_velocity_model->vp_status == 1) {
		// Read from file.
		fseek(fp, location * sizeof(float), SEEK_SET);
		fread(&(data->vp), sizeof(float), 1, fp);
	}

	// Check our loaded components of the model.
	if (albacore_velocity_model->rho_status == 2) {
		// Read from memory.
		ptr = (float *)albacore_velocity_model->rho;
		data->rho = ptr[location];
	} else if (albacore_velocity_model->rho_status == 1) {
		// Read from file.
		fseek(fp, location * sizeof(float), SEEK_SET);
		fread(&(data->rho), sizeof(float), 1, fp);
	}
}

/**
 * Trilinearly interpolates given a x percentage, y percentage, z percentage and a cube of
 * data properties in top origin format (top plane first, bottom plane second).
 *
 * @param x_percent X percentage
 * @param y_percent Y percentage
 * @param z_percent Z percentage
 * @param eight_points Eight surrounding data properties
 * @param ret_properties Returned data properties
 */
void albacore_trilinear_interpolation(double x_percent, double y_percent, double z_percent,
							 albacore_properties_t *eight_points, albacore_properties_t *ret_properties) {
	albacore_properties_t *temp_array = calloc(2, sizeof(albacore_properties_t));
	albacore_properties_t *four_points = eight_points;

	albacore_bilinear_interpolation(x_percent, y_percent, four_points, &temp_array[0]);

	// Now advance the pointer four "albacore_properties_t" spaces.
	four_points += 4;

	// Another interpolation.
	albacore_bilinear_interpolation(x_percent, y_percent, four_points, &temp_array[1]);

	// Now linearly interpolate between the two.
	albacore_linear_interpolation(z_percent, &temp_array[0], &temp_array[1], ret_properties);

	free(temp_array);
}

/**
 * Bilinearly interpolates given a x percentage, y percentage, and a plane of data properties in
 * origin, bottom-right, top-left, top-right format.
 *
 * @param x_percent X percentage.
 * @param y_percent Y percentage.
 * @param four_points Data property plane.
 * @param ret_properties Returned data properties.
 */
void albacore_bilinear_interpolation(double x_percent, double y_percent, albacore_properties_t *four_points, albacore_properties_t *ret_properties) {

	albacore_properties_t *temp_array = calloc(2, sizeof(albacore_properties_t));

	albacore_linear_interpolation(x_percent, &four_points[0], &four_points[1], &temp_array[0]);
	albacore_linear_interpolation(x_percent, &four_points[2], &four_points[3], &temp_array[1]);
	albacore_linear_interpolation(y_percent, &temp_array[0], &temp_array[1], ret_properties);

	free(temp_array);
}

/**
 * Linearly interpolates given a percentage from x0 to x1, a data point at x0, and a data point at x1.
 *
 * @param percent Percent of the way from x0 to x1 (from 0 to 1 interval).
 * @param x0 Data point at x0.
 * @param x1 Data point at x1.
 * @param ret_properties Resulting data properties.
 */
void albacore_linear_interpolation(double percent, albacore_properties_t *x0, albacore_properties_t *x1, albacore_properties_t *ret_properties) {

	ret_properties->vp  = (1 - percent) * x0->vp  + percent * x1->vp;
	ret_properties->vs  = (1 - percent) * x0->vs  + percent * x1->vs;
	ret_properties->rho = (1 - percent) * x0->rho + percent * x1->rho;
}

/**
 * Called when the model is being discarded. Free all variables.
 *
 * @return SUCCESS
 */
int albacore_finalize() {
	if (albacore_velocity_model->vs) free(albacore_velocity_model->vs);
	if (albacore_velocity_model->vp) free(albacore_velocity_model->vp);
	if (albacore_velocity_model->rho) free(albacore_velocity_model->rho);

	free(albacore_configuration);

	return SUCCESS;
}

/**
 * Returns the version information.
 *
 * @param ver Version string to return.
 * @param len Maximum length of buffer.
 * @return Zero
 */
int albacore_version(char *ver, int len)
{
  int verlen;
  verlen = strlen(albacore_version_string);
  if (verlen > len - 1) {
    verlen = len - 1;
  }
  memset(ver, 0, len);
  strncpy(ver, albacore_version_string, verlen);
  return 0;
}

/**
 * Reads the configuration file describing the various properties of CVM-S5 and populates
 * the configuration struct. This assumes configuration has been "calloc'ed" and validates
 * that each value is not zero at the end.
 *
 * @param file The configuration file location on disk to read.
 * @param config The configuration struct to which the data should be written.
 * @return Success or failure, depending on if file was read successfully.
 */
int albacore_read_configuration(char *file, albacore_configuration_t *config) {
	FILE *fp = fopen(file, "r");
	char key[40];
	char value[80];
	char line_holder[128];

	// If our file pointer is null, an error has occurred. Return fail.
	if (fp == NULL) {
		print_error("Could not open the configuration file.");
		return FAIL;
	}

	// Read the lines in the configuration file.
	while (fgets(line_holder, sizeof(line_holder), fp) != NULL) {
		if (line_holder[0] != '#' && line_holder[0] != ' ' && line_holder[0] != '\n') {
			sscanf(line_holder, "%s = %s", key, value);

			// Which variable are we editing?
			if (strcmp(key, "utm_zone") == 0)
  				config->utm_zone = atoi(value);
			if (strcmp(key, "model_dir") == 0)
				sprintf(config->model_dir, "%s", value);
			if (strcmp(key, "nx") == 0)
  				config->nx = atoi(value);
			if (strcmp(key, "ny") == 0)
  			 	config->ny = atoi(value);
			if (strcmp(key, "nz") == 0)
  			 	config->nz = atoi(value);
			if (strcmp(key, "depth") == 0)
  			 	config->depth = atof(value);
			if (strcmp(key, "top_left_corner_e") == 0)
				config->top_left_corner_e = atof(value);
			if (strcmp(key, "top_left_corner_n") == 0)
		 		config->top_left_corner_n = atof(value);
			if (strcmp(key, "top_right_corner_e") == 0)
				config->top_right_corner_e = atof(value);
			if (strcmp(key, "top_right_corner_n") == 0)
				config->top_right_corner_n = atof(value);
			if (strcmp(key, "bottom_left_corner_e") == 0)
				config->bottom_left_corner_e = atof(value);
			if (strcmp(key, "bottom_left_corner_n") == 0)
				config->bottom_left_corner_n = atof(value);
			if (strcmp(key, "bottom_right_corner_e") == 0)
				config->bottom_right_corner_e = atof(value);
			if (strcmp(key, "bottom_right_corner_n") == 0)
				config->bottom_right_corner_n = atof(value);
			if (strcmp(key, "depth_interval") == 0)
				config->depth_interval = atof(value);
			if (strcmp(key, "p0") == 0)
				config->p0 = atof(value);
			if (strcmp(key, "p1") == 0)
				config->p1 = atof(value);
			if (strcmp(key, "p2") == 0)
				config->p2 = atof(value);
			if (strcmp(key, "p3") == 0)
				config->p3 = atof(value);
			if (strcmp(key, "p4") == 0)
				config->p4 = atof(value);
			if (strcmp(key, "p5") == 0)
				config->p5 = atof(value);
			if (strcmp(key, "interpolation") == 0) {
                                if (strcmp(value, "on") == 0) {
                                     config->interpolation = 1;
                                     } else {
                                          config->interpolation = 0;
                                }
                        };

		}
	}

	// Have we set up all configuration parameters?
	if (config->utm_zone == 0 || config->nx == 0 || config->ny == 0 || config->nz == 0 || config->model_dir[0] == '\0' ||
		config->top_left_corner_e == 0 || config->top_left_corner_n == 0 || config->top_right_corner_e == 0 ||
		config->top_right_corner_n == 0 || config->bottom_left_corner_e == 0 || config->bottom_left_corner_n == 0 ||
		config->bottom_right_corner_e == 0 || config->bottom_right_corner_n == 0 || config->depth == 0 ||
		config->depth_interval == 0 || config->p0 == 0 || config->p1 == 0 || config->p2 == 0 || config->p3 == 0 ||
		config->p4 == 0 || config->p5 == 0) {
		print_error("One configuration parameter not specified. Please check your configuration file.");
		return FAIL;
	}

	fclose(fp);

	return SUCCESS;
}

/**
 * Calculates the density based off of Vs. Based on Nafe-Drake scaling relationship.
 *
 * @param vs The Vs value off which to scale.
 * @return Density, in g/m^3.
 */
double albacore_calculate_density(double vs) {
	double retVal;
	vs = vs / 1000;
	retVal = albacore_configuration->p0 + albacore_configuration->p1 * vs + albacore_configuration->p2 * pow(vs, 2) +
			 albacore_configuration->p3 * pow(vs, 3) + albacore_configuration->p4 * pow(vs, 4) + albacore_configuration->p5 * pow(vs, 5);
	retVal = retVal * 1000;
	return retVal;
}

/**
 * Prints the error string provided.
 *
 * @param err The error string to print out to stderr.
 */
void print_error(char *err) {
	fprintf(stderr, "An error has occurred while executing CVM-S5. The error was:\n\n");
	fprintf(stderr, "%s", err);
	fprintf(stderr, "\n\nPlease contact software@scec.org and describe both the error and a bit\n");
	fprintf(stderr, "about the computer you are running CVM-S5 on (Linux, Mac, etc.).\n");
}

/**
 * Tries to read the model into memory.
 *
 * @param model The model parameter struct which will hold the pointers to the data either on disk or in memory.
 * @return 2 if all files are read to memory, SUCCESS if file is found but at least 1
 * is not in memory, FAIL if no file found.
 */
int albacore_try_reading_model(albacore_model_t *model) {
	double base_malloc = albacore_configuration->nx * albacore_configuration->ny * albacore_configuration->nz * sizeof(float);
	int file_count = 0;
	int all_read_to_memory = 1;
	char current_file[128];
	FILE *fp;

	// Let's see what data we actually have.
	sprintf(current_file, "%s/vp.dat", albacore_data_directory);
	if (access(current_file, R_OK) == 0) {
		model->vp = malloc(base_malloc);
		if (model->vp != NULL) {
			// Read the model in.
			fp = fopen(current_file, "rb");
			fread(model->vp, 1, base_malloc, fp);
			fclose(fp);
			model->vp_status = 2;
		} else {
			all_read_to_memory = 0;
			model->vp = fopen(current_file, "rb");
			model->vp_status = 1;
		}
		file_count++;
	}

	sprintf(current_file, "%s/vs.dat", albacore_data_directory);
	if (access(current_file, R_OK) == 0) {
		model->vs = malloc(base_malloc);
		if (model->vs != NULL) {
			// Read the model in.
			fp = fopen(current_file, "rb");
			fread(model->vs, 1, base_malloc, fp);
			fclose(fp);
			model->vs_status = 2;
		} else {
			all_read_to_memory = 0;
			model->vs = fopen(current_file, "rb");
			model->vs_status = 1;
		}
		file_count++;
	}

	sprintf(current_file, "%s/rho.dat", albacore_data_directory);
	if (access(current_file, R_OK) == 0) {
		model->rho = malloc(base_malloc);
		if (model->rho != NULL) {
			// Read the model in.
			fp = fopen(current_file, "rb");
			fread(model->rho, 1, base_malloc, fp);
			fclose(fp);
			model->rho_status = 2;
		} else {
			all_read_to_memory = 0;
			model->rho = fopen(current_file, "rb");
			model->rho_status = 1;
		}
		file_count++;
	}

	if (file_count == 0)
		return FAIL;
	else if (file_count > 0 && all_read_to_memory == 0)
		return SUCCESS;
	else
		return 2;
}

// The following functions are for dynamic library mode. If we are compiling
// a static library, these functions must be disabled to avoid conflicts.
#ifdef DYNAMIC_LIBRARY

/**
 * Init function loaded and called by the UCVM library. Calls albacore_init.
 *
 * @param dir The directory in which UCVM is installed.
 * @return Success or failure.
 */
int model_init(const char *dir, const char *label) {
	return albacore_init(dir, label);
}

/**
 * Query function loaded and called by the UCVM library. Calls albacore_query.
 *
 * @param points The basic_point_t array containing the points.
 * @param data The basic_properties_t array containing the material properties returned.
 * @param numpoints The number of points in the array.
 * @return Success or fail.
 */
int model_query(albacore_point_t *points, albacore_properties_t *data, int numpoints) {
	return albacore_query(points, data, numpoints);
}

/**
 * Finalize function loaded and called by the UCVM library. Calls albacore_finalize.
 *
 * @return Success
 */
int model_finalize() {
	return albacore_finalize();
}

/**
 * Version function loaded and called by the UCVM library. Calls albacore_version.
 *
 * @param ver Version string to return.
 * @param len Maximum length of buffer.
 * @return Zero
 */
int model_version(char *ver, int len) {
	return albacore_version(ver, len);
}

int (*get_model_init())(const char *, const char *) {
        return &albacore_init;
}
int (*get_model_query())(albacore_point_t *, albacore_properties_t *, int) {
         return &albacore_query;
}
int (*get_model_finalize())() {
         return &albacore_finalize;
}
int (*get_model_version())(char *, int) {
         return &albacore_version;
}



#endif

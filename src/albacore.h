/**
 * @file albacore.h
 * @brief Main header file for ALBACORE library.
 * @author - SCEC 
 * @version 1.0
 *
 * Delivers ALBACORE, Offshore Southern California Velocity Model
 *
 */

// Includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

// Constants
#ifndef M_PI
	/** Defines pi */
	#define M_PI 3.14159265358979323846
#endif

/** Defines a return value of success */
#define SUCCESS 0
/** Defines a return value of failure */
#define FAIL 1
/** Defines a return value of NA from model */
#define NA -1 

// Structures
/** Defines a point (latitude, longitude, and depth) in WGS84 format */
typedef struct albacore_point_t {
	/** Longitude member of the point */
	double longitude;
	/** Latitude member of the point */
	double latitude;
	/** Depth member of the point */
	double depth;
} albacore_point_t;

/** Defines the material properties this model will retrieve. */
typedef struct albacore_properties_t {
	/** P-wave velocity in meters per second */
	double vp;
	/** S-wave velocity in meters per second */
	double vs;
	/** Density in g/m^3 */
	double rho;
} albacore_properties_t;

/** The ALBACORE configuration structure. */
typedef struct albacore_configuration_t {
	/** The zone of UTM projection */
	int utm_zone;
	/** The model directory */
	char model_dir[128];
	/** Number of x points */
	int nx;
	/** Number of y points */
	int ny;
	/** Number of z points */
	int nz;
	/** Depth in meters */
	double depth;
	/** Top left corner easting */
	double top_left_corner_e;
	/** Top left corner northing */
	double top_left_corner_n;
	/** Top right corner easting */
	double top_right_corner_e;
	/** Top right corner northing */
	double top_right_corner_n;
	/** Bottom left corner easting */
	double bottom_left_corner_e;
	/** Bottom left corner northing */
	double bottom_left_corner_n;
	/** Bottom right corner easting */
	double bottom_right_corner_e;
	/** Bottom right corner northing */
	double bottom_right_corner_n;
	/** Z interval for the data */
	double depth_interval;
	/** Brocher 2005 scaling polynomial coefficient 10^0 */
	double p0;
	/** Brocher 2005 scaling polynomial coefficient 10^1 */
	double p1;
	/** Brocher 2005 scaling polynomial coefficient 10^2 */
	double p2;
	/** Brocher 2005 scaling polynomial coefficient 10^3 */
	double p3;
	/** Brocher 2005 scaling polynomial coefficient 10^4 */
	double p4;
	/** Brocher 2005 scaling polynomial coefficient 10^5 */
	double p5;
        /** Bilinear or Trilinear Interpolation on or off (1 or 0) */
        int interpolation;

} albacore_configuration_t;

/** The model structure which points to available portions of the model. */
typedef struct albacore_model_t {
	/** A pointer to the Vs data either in memory or disk. Null if does not exist. */
	void *vs;
	/** Vs status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
	int vs_status;
	/** A pointer to the Vp data either in memory or disk. Null if does not exist. */
	void *vp;
	/** Vp status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
	int vp_status;
	/** A pointer to the rho data either in memory or disk. Null if does not exist. */
	void *rho;
	/** Rho status: 0 = not found, 1 = found and not in memory, 2 = found and in memory */
	int rho_status;
	/** A pointer to the Qp data either in memory or disk. Null if does not exist. */
} albacore_model_t;

// Constants
/** The version of the model. */
const char *albacore_version_string = "ALBACORE";

// Variables
/** Set to 1 when the model is ready for query. */
int albacore_is_initialized = 0;

/** Location of the binary data files. */
char albacore_data_directory[128];

/** Configuration parameters. */
albacore_configuration_t *albacore_configuration;
/** Holds pointers to the velocity model data OR indicates it can be read from file. */
albacore_model_t *albacore_velocity_model;

/** The height of this model's region, in meters. */
double albacore_total_height_m = 0;
/** The width of this model's region, in meters. */
double albacore_total_width_m = 0;

// UCVM API Required Functions

#ifdef DYNAMIC_LIBRARY

/** Initializes the model */
int model_init(const char *dir, const char *label);
/** Cleans up the model (frees memory, etc.) */
int model_finalize();
/** Returns version information */
int model_version(char *ver, int len);
/** Queries the model */
int model_query(albacore_point_t *points, albacore_properties_t *data, int numpts);

#endif

// ALBACORE Related Functions

/** Initializes the model */
int albacore_init(const char *dir, const char *label);
/** Cleans up the model (frees memory, etc.) */
int albacore_finalize();
/** Returns version information */
int albacore_version(char *ver, int len);
/** Queries the model */
int albacore_query(albacore_point_t *points, albacore_properties_t *data, int numpts);

// Non-UCVM Helper Functions
/** Reads the configuration file. */
int albacore_read_configuration(char *file, albacore_configuration_t *config);
void print_error(char *err);
/** Retrieves the value at a specified grid point in the model. */
void albacore_read_properties(int x, int y, int z, albacore_properties_t *data);
/** Attempts to malloc the model size in memory and read it in. */
int albacore_try_reading_model(albacore_model_t *model);
/** Calculates density from Vs. */
double albacore_calculate_density(double vs);

// Interpolation Functions
/** Linearly interpolates two albacore_properties_t structures */
void albacore_linear_interpolation(double percent, albacore_properties_t *x0, albacore_properties_t *x1, albacore_properties_t *ret_properties);
/** Bilinearly interpolates the properties. */
void albacore_bilinear_interpolation(double x_percent, double y_percent, albacore_properties_t *four_points, albacore_properties_t *ret_properties);
/** Trilinearly interpolates the properties. */
void albacore_trilinear_interpolation(double x_percent, double y_percent, double z_percent, albacore_properties_t *eight_points,
							 albacore_properties_t *ret_properties);

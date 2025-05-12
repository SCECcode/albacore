/**
 * @file test.c
 * @brief Bootstraps the test framework for the ALBACORE library.
 * @author - SCEC
 * @version 1.0
 *
 * Tests the ALBACORE library by loading it and executing the code as
 * UCVM would do it.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "albacore.h"

/**
 * Initializes and runs the test program. Tests link against the
 * static version of the library to prevent any dynamic loading
 * issues.
 *
 * @param argc The number of arguments.
 * @param argv The argument strings.
 * @return A zero value indicating success.
 */
int main(int argc, const char* argv[]) {

	// Declare the structures.
	albacore_point_t pt;
	albacore_properties_t ret;

	// Initialize the model.
	char *envstr=getenv("UCVM_INSTALL_PATH");
        if(envstr != NULL) {
           assert(albacore_init(envstr, "albacore") == 0);
           } else {
             assert(albacore_init("..", "albacore") == 0);
        }

	printf("Loaded the model successfully.\n");

	// Query a point.
	pt.longitude = -118;
	pt.latitude = 34;
	pt.depth = 0;

	albacore_query(&pt, &ret, 1);

	assert(ret.vs > 0);
	assert(ret.vp > 0);
	assert(ret.rho > 0);

	printf("Query was successful.\n");

	// Close the model.
	assert(albacore_finalize() == 0);

	printf("Model closed successfully.\n");

	printf("\nALL ALBACORE TESTS PASSED");

	return 0;
}

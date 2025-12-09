#ifndef STRUCTS_NETCDF
#define STRUCTS_NETCDF

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <stdbool.h>
#include <pthread.h>

#include <netcdf.h>
#include <gsl/gsl_interp.h>

#define TYPE_VAR_NC_BYTE signed char
#define TYPE_VAR_NC_CHAR char
#define TYPE_VAR_NC_SHORT short
#define TYPE_VAR_NC_INT int
#define TYPE_VAR_NC_FLOAT float
#define TYPE_VAR_NC_DOUBLE double
#define TYPE_VAR_NC_UBYTE unsigned char
#define TYPE_VAR_NC_USHORT unsigned short
#define TYPE_VAR_NC_UINT unsigned int
#define TYPE_VAR_NC_INT64 long long
#define TYPE_VAR_NC_UINT64 unsigned long long
#define TYPE_VAR_NC_STRING char *

#define VALUE_ERR 9969209968386869046778552952102584320.0
#define VALUE_NAN NAN

#define GET_START gettimeofday(&begin, 0);

#define NODE_POOL_SIZE 1000

#define GET_END                                                                        \
    {                                                                                  \
        gettimeofday(&end, 0); /* Stop measuring time and calculate the elapsed time*/ \
        seconds = end.tv_sec - begin.tv_sec;                                           \
        microseconds = end.tv_usec - begin.tv_usec;                                    \
        elapsed = seconds + microseconds * 1e-6;                                       \
        /* printf("---> Time measured: %.3f seconds.\n", elapsed);*/                   \
        printf("%.3f,", elapsed);                                                      \
    }

/* Init - NetCDF data structure */
typedef struct
{
    char name[NC_MAX_NAME + 1];
    size_t len;
} Dimension;

typedef struct
{
    char name[NC_MAX_NAME + 1];
    nc_type type;
    int ndims;
    int id;
    int invalid_count;
    int num_win_valid_training;
    int num_win_valid_prediction;
    double invalid_percentage;
    double rmse;
    double bias;
    void *num_valid_window;
    void *data;
    void *created_data;
} Variable;

typedef struct
{
    int ncid_in;
    int ndims;
    int nvars;
    Dimension *dim;
    Variable *var;
} NetCDF;
/* End - NetCDF data structure */

/* Init - Data Segment structure */
typedef struct
{
    unsigned int num_thread;
    int start_prediction;
    int end_prediction;
    int start_training;
    int end_training;
    int k;
    int win_size;
    int win_size_interpolation;
    int num_Na;
    int win_Na;
    int argc;
    int indice_generic;
    int win_count;
    float current_best_distance;
    NetCDF *predicted_file;
    NetCDF *predictor_file;
} DataSegment;

typedef void (*process_func)(NetCDF *, DataSegment *);

/* End - Data Segment structure */

/* Init - Processing structure */
typedef struct
{
    unsigned int window_index;
    double distance;
} ClosestPoint;

/* End - Processing structure */

#endif
#ifndef PREPROCESS_NETCDF
#define PREPROCESS_NETCDF

#include "structs.h"

#define PRINT_LN_VALUES(TYPE, FMT)                           \
    case TYPE:                                               \
    {                                                        \
        printf(FMT "\t", ((TYPE_VAR_##TYPE *)var->data)[i]); \
    };                                                       \
    break;

#define VALIDATE_DATA(TYPE)                                                      \
    case TYPE:                                                                   \
    {                                                                            \
        for (int j = 0; j < file->dim->len; j++)                                 \
        {                                                                        \
            if (isnan((double)(((TYPE_VAR_##TYPE *)var->data)[j])))              \
            {                                                                    \
                invalid_count++;                                                 \
                continue;                                                        \
            }                                                                    \
            if (((TYPE_VAR_##TYPE *)var->data)[j] == (TYPE_VAR_##TYPE)VALUE_ERR) \
            {                                                                    \
                ((TYPE_VAR_##TYPE *)var->data)[j] = (TYPE_VAR_##TYPE)VALUE_NAN;  \
                invalid_count++;                                                 \
            }                                                                    \
        }                                                                        \
    };                                                                           \
    break;

#define IF_ISNAN(TYPE)                                                  \
    case TYPE:                                                          \
    {                                                                   \
        if (isnan((double)((TYPE_VAR_##TYPE *)var->data)[j]))           \
        {                                                               \
            if (!d_time[0] && j > 0)                                    \
            {                                                           \
                d_time[0] = ((TYPE_VAR_##TYPE *)var_zero->data)[j - 1]; \
                d_interp[0] = ((TYPE_VAR_##TYPE *)var->data)[j - 1];    \
                d_index[0] = j;                                         \
            }                                                           \
            count++;                                                    \
            continue;                                                   \
        }                                                               \
        if (d_time[0] && !d_time[1])                                    \
        {                                                               \
            if (count > ds->win_size_interpolation)                     \
            {                                                           \
                count = 0;                                              \
                d_time[0] = 0;                                          \
                d_interp[0] = 0;                                        \
                d_index[0] = 0;                                         \
                continue;                                               \
            };                                                          \
            d_time[1] = ((TYPE_VAR_##TYPE *)var_zero->data)[j];         \
            d_interp[1] = ((TYPE_VAR_##TYPE *)var->data)[j];            \
            d_index[1] = j;                                             \
        }                                                               \
    };                                                                  \
    break;

#define INTERP_DATA(TYPE)                                                   \
    case TYPE:                                                              \
    {                                                                       \
        gsl_interp_init(interp, d_time, d_interp, 2);                       \
        for (int l = d_index[0]; l <= d_index[1]; l++)                      \
        {                                                                   \
            if (isnan((double)((TYPE_VAR_##TYPE *)var->data)[l]))           \
            {                                                               \
                ((TYPE_VAR_##TYPE *)var->data)[l] =                         \
                    gsl_interp_eval(interp,                                 \
                                    d_time,                                 \
                                    d_interp,                               \
                                    ((TYPE_VAR_##TYPE *)var_zero->data)[l], \
                                    acc);                                   \
            }                                                               \
        }                                                                   \
        d_changes += count;                                                 \
        count = 0;                                                          \
        d_time[0] = 0;                                                      \
        d_time[1] = 0;                                                      \
        d_interp[0] = 0;                                                    \
        d_interp[1] = 0;                                                    \
        d_index[0] = 0;                                                     \
        d_index[1] = 0;                                                     \
    };                                                                      \
    break;

#define COUNT_VALID_WIN(TYPE)                                                        \
    case TYPE:                                                                       \
    {                                                                                \
        for (int j = start; j <= end; j++)                                           \
        {                                                                            \
            invalid_count = 0;                                                       \
            for (int i = 0; i < ds->win_size; i++)                                   \
                if (isnan((double)(((TYPE_VAR_##TYPE *)var->data)[(j - ds->k) + i]))) \
                {                                                                    \
                    invalid_count++;                                                 \
                    break;                                                           \
                }                                                                    \
            if (invalid_count == 0)                                                  \
                num_isvalid++;                                                       \
        }                                                                            \
    };                                                                               \
    break;

time_t convert_time(char *);
int binary_search(NetCDF *, int);
void analyze_data(NetCDF *, DataSegment *, process_func);
void print_data_values(NetCDF *, DataSegment *);
void count_invalid_values(NetCDF *, DataSegment *);
void count_valid_window(NetCDF *, DataSegment *);
void print_info_percentage(NetCDF *, DataSegment *);
void interpolation_values(NetCDF *, DataSegment *);

#endif
#ifndef RANDW_NETCDF
#define RANDW_NETCDF

#include "structs.h"

#define ALLOCATE_MEMORY(TYPE, LENGTH)                    \
    case TYPE:                                           \
        data = malloc(LENGTH * sizeof(TYPE_VAR_##TYPE)); \
        break;

void handle_error(int);
NetCDF *create_struct(DataSegment *, char *[]);
void *allocate_memory(char, size_t);
void deallocate_memory(NetCDF *, int);
void read_dimensions(NetCDF *);
void read_variables(NetCDF *);
void read_header_file(NetCDF *, DataSegment *);
void read_data_file(NetCDF *, DataSegment *);
void write_file(NetCDF *, char *);

#endif
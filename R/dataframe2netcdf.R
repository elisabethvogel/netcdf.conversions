#' Writing a data frame to a netcdf file.
#'
#' @param data_frame          Data frame that includes dimension names and at least one variable.
#' @param netcdf_file         File path for netcdf file to write to.
#' @param dim_names           Names of dimensions in data frame.
#' @param dim_units           Units for all dimensions.
#' @param var_names           Names of variables in data frame.
#' @param var_units           Units for all variables.
#' @param overwrite_existing  If file exists already, should it be overwritten?
#' @param fill_dimensions     If TRUE, dimensions will be extended to cover the whole world.
#' @param max_resolution_testing Maximum number of resolutions to be tested, for dimension filling.
#' @param verbose             Include print statements.
#' @return Returns file name of netcdf file, if saved successfully.
#' @export
dataframe2netcdf = function(data_frame,
                            netcdf_file,
                            dim_names = intersect(names(data_frame),
                                                  c("lon", "lat", "longitude", "latitude", "time")),
                            dim_units = NULL,
                            var_names = names(data_frame)[!names(data_frame) %in% dim_names],
                            var_units = rep("", length(var_names)),
                            overwrite_existing = FALSE,
                            fill_dimensions = FALSE,
                            max_resolution_testing = 5,
                            verbose = FALSE) {

  if (verbose) print("function: dataframe2netcdf")

  # Some tests
  stopifnot(is.data.frame(data_frame), nrow(data_frame) > 0)
  if (!is.null(dim_units)) stopifnot(length(dim_names) == length(dim_units))
  if (!is.null(var_units)) stopifnot(length(var_units) == length(var_names))
  stopifnot(is.logical(fill_dimensions))

  x = 1

  # Get units for each dimension
  if (is.null(dim_units)) {
    dim_units = sapply(dim_names, FUN = function(x) {
      if (grepl("lon", x)) return("longitude")
      if (grepl("lat", x)) return("latitude")
      if (x == "time") return(NA) # time units are automatically determined later in the code
    })
  }

  if (file.exists(netcdf_file) & overwrite_existing) {
    file.remove(netcdf_file)
  } else if (file.exists(netcdf_file) & !overwrite_existing) {
    stop("File exists already.")
  }

  # select only dimensions and variables that were passed to function
  data_frame = data_frame[, c(dim_names, var_names)]

  # for all dimensions, define a netcdf dimension
  dim_nc = list()
  dim_values = list()
  n_dims = length(dim_names)

  for (dim_idx in 1:n_dims) {

    # to do: create better implementation for time unit here

    dim_name = dim_names[dim_idx]

    if (dim_name == "time") {
      temp = time2nctime(data_frame$time)
      temp$time_vals = as.numeric(temp$time_vals)
      data_frame$time = temp$time_vals
      dim_unit = temp$time_unit
      dim_vals = sort(unique(data_frame[, dim_name, drop = TRUE]))

    } else {
      dim_unit = dim_units[dim_idx]
      dim_vals = sort(unique(data_frame[, dim_name, drop = TRUE]))

      if (length(dim_vals) > 2) {

        # check for evenly distributed values for lat/lon
        # (otherwise the netcdf file looks wrong)

        # round because of subtraction error
        dist = round(dim_vals[2:length(dim_vals)] - dim_vals[1:length(dim_vals) - 1], 10)

        res_temp = min(dist) # guess the resolution by taking the smallest difference

        # if there are unevenly distributed values, try to guess resolution and
        # fill up with NAs
        if (length(unique(dist)) != 1) { # unevenly distributed

          # test different resolutions
          for (k in 1:max_resolution_testing) {

            if (verbose) print(sprintf("Testing resolution: %.5f", res_temp/k))

            dim_vals_new = seq(min(dim_vals), max(dim_vals), by = res_temp/k)

            # test if all old dim_vals are completely included in dim_vals_new
            if (all(dim_vals %in% dim_vals_new)) {
              if (verbose) print("Ok, this worked...")
              dim_vals = dim_vals_new
              break # stop here
            }

            if (k == max_resolution_testing) {
              # if it didn't work until here
              stop(sprintf("Uneven values in %s dimension. Could not create netcdf file.", dim_name))
            }
          }
        }

        # fill dimensions (lon: 0 - 360, lat: -90 - 90)
        if (fill_dimensions && dim_name %in% c("lon", "longitude", "lat", "latitude")) {

          if (dim_name %in% c("lon", "longitude")) {
            if (any(dim_vals > 180)) {
              # longitudes from 0 to 360
              dim_vals_new = seq(0 + res_temp/2, 360 - res_temp/2, by = res_temp)
            } else {
              dim_vals_new = seq(-180 + res_temp/2, 180 - res_temp/2, by = res_temp)
            }

          } else if (dim_name %in% c("lat", "latitude")) {
            if (any(dim_vals > 90)) {
              dim_vals_new = seq(0 + res_temp/2, 180 - res_temp/2, by = res_temp)
            } else {
              dim_vals_new = seq(-90 + res_temp/2, 90 - res_temp/2, by = res_temp)
            }
          }

          # test if all old dim_vals are completely included in dim_vals_new
          if (all(dim_vals %in% dim_vals_new)) {
            dim_vals = dim_vals_new
          } else {
            stop(sprintf("Could not fill %s dimension. Stopped.", dim_name))
          }
        }
      }
    }

    dim_values[[dim_name]] = dim_vals
    dim_nc[[dim_name]] = ncdf4::ncdim_def(name = dim_name, units = dim_unit, vals = dim_vals)
  } #dims

  # create all possible combinations of dimensions
  dimensions_new = expand.grid(dim_values)
  stopifnot(nrow(dimensions_new) >= nrow(data_frame))

  if (nrow(dimensions_new) > nrow(data_frame)) {
    data_frame = dplyr::full_join(dimensions_new, data_frame)
  }

  # for all variables, define a netcdf variable
  var_nc = list()
  for (var_idx in 1:length(var_names)) {
    var_name = var_names[var_idx]
    if (!is.numeric(data_frame[, var_name, drop = TRUE])) next
    # to do: implement non-numeric variables

    var_unit = var_units[var_idx]
    var_nc[[var_name]] = ncdf4::ncvar_def(name = var_name,
                                          units = var_unit,
                                          dim = dim_nc,
                                          prec = "double")
  }

  # create netcdf or add variables to existing netcdf file
  if (!file.exists(netcdf_file)) {
    output_dir = dirname(netcdf_file)
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    netcdf = ncdf4::nc_create(filename = netcdf_file, vars = var_nc)
  } else {
    netcdf = ncdf4::nc_open(filename = netcdf_file, write = TRUE)
    for (var in var_nc) {
      netcdf = ncdf4::ncvar_add(netcdf, v = var)
    }
  }

  # put values into netcdf
  for (var_idx in 1:length(var_names)) {
    var_name = var_names[var_idx]

    if (!is.numeric(data_frame[, var_name, drop = TRUE])) next
    # to do: implement non-numeric variables

    var_vals = data_frame[, c(dim_names, var_name)]

    # create an n-dimensional array using acast
    formula = formula(paste(dim_names, collapse = " ~ "))
    var_vals = reshape2::acast(var_vals, formula = formula, value.var = var_name)
    ncdf4::ncvar_put(nc = netcdf, varid = var_name, vals = as.matrix(var_vals))
  }

  # close netcdf
  ncdf4::nc_close(netcdf)

  if (verbose) print(sprintf("Netcdf file created: %s", netcdf_file))
  return(netcdf_file)
}

#' Reading a netcdf file into a data frame.
#'
#' @param netcdf_file         Path to a netcdf file.
#' @param variables           Vector of variable names to read in.
#' @param remove_NA           Specify if empty observations should be removed.
#' @param time_format         Format vector. Specify how time variable should be returned.
#' @param start_idx           Specify the start index of the read in.
#' @param end_idx,count       Specify either the end index (end_idx) or the number of steps to read in (count).
#' @param grid_cells          Data frame with lat/lon columns that specify the coordinates to read in.
#' @param return_time_columns If TRUE, it returns columns containing time information (e.g. year, month, day, etc.)
#' @param years               Specify which years should be read in.
#' @param verbose             Include print statements.
#' @return A data frame with the dimensions and variables stored in the netcdf file, as specified.
#' @examples
#' # Download a file from the HadCRU dataset:
#' url = "https://crudata.uea.ac.uk/cru/data/temperature/HadSST.3.1.1.0.median.nc"
#' netcdf_file = "HadSST.3.1.1.0.median.nc"
#' download.file(url, netcdf_file)
#'
#' df = netcdf2dataframe(netcdf_file,
#'                       variables = "sst",
#'                       time_format = "%Y-%m-%d")
#' df = netcdf2dataframe(netcdf_file,
#'                       variables = "sst",
#'                       time_format = "%Y-%m-%d",
#'                       remove_NA = TRUE)
#' df = netcdf2dataframe(netcdf_file,
#'                       variables = "sst",
#'                       time_format = "%Y-%m-%d",
#'                       years = 2000)
#' df = netcdf2dataframe(netcdf_file,
#'                       variables = "sst",
#'                       time_format = "%Y-%m-%d",
#'                       years = 2000:2010,
#'                       grid_cells = data.frame(
#'                       lon = -17.5, lat = -72.5))
#' @export
netcdf2dataframe = function(netcdf_file, variables = "all", remove_NA = FALSE,
                            time_format = "%Y-%m-%d %H:%M:%S",
                            start_idx = NULL, end_idx = NULL,
                            count = NULL,
                            grid_cells = NULL,
                            years = NULL,
                            return_time_columns = FALSE,
                            verbose = FALSE) {

  if (verbose) print("function: netcdf2dataframe")

  # start timer
  if (verbose) start.time <- Sys.time()

  # some tests
  stopifnot(is.character(netcdf_file))
  stopifnot(is.character(variables))
  stopifnot(is.logical(remove_NA))
  stopifnot(is.null(time_format) || is.character(time_format))
  stopifnot(is.null(start_idx) || (is.numeric(start_idx) && all(start_idx >= 1)))
  stopifnot(is.null(end_idx) || (is.numeric(end_idx) && all(end_idx >= 1)))
  stopifnot(is.null(count) || is.numeric(count))
  # don't provide both count and end_idx (redundant)
  stopifnot(is.null(count) || is.null(end_idx))
  # check grid cells is a data frame with two columns: latitudes and longitudes
  if (!is.null(grid_cells)) {
    stopifnot(is.data.frame(grid_cells) &&
                length(intersect(names(grid_cells),
                                 c("lat", "lon", "latitude", "longitude")) == 2))
  }


  # preparations ========================================

  # option 1: both start and end idx given for all variables
  if (!is.null(start_idx) && !is.null(end_idx)) {
    stopifnot(length(start_idx) == length(end_idx))
    # create count vector from end_idx vector
    count = end_idx - start_idx + rep(1, length(end_idx))

    # option 2: only start vector, no end vector:
    # set count vector to -1 (read all values to the end)
  } else if (!is.null(start_idx) && is.null(end_idx)) {
    count = rep(-1, length(start_idx))

    # option 3: start and count vector given: use directly
  } else if (!is.null(start_idx) && !is.null(count)) {
    stopifnot(length(start_idx) == length(count))
  }

  # get dimension names and values
  dim_info = ncdf.tools::infoNcdfDims(netcdf_file)
  # get variable names and dimensions
  var_info = ncdf.tools::infoNcdfVars(netcdf_file)

  # if bnds is part of the dims or vars, remove it (from cdo averaging)
  dim_info = dim_info[dim_info$name != "bnds", ]
  var_info = var_info[var_info$name != "bnds", ]

  dim_names = dim_info$name
  n_dims = length(dim_names)
  dim_size = dim_info$length
  names(dim_size) = dim_names

  # coordinates = readNcdfCoordinates(netcdf_file) # this takes forever
  nc = ncdf4::nc_open(netcdf_file)
  coordinates = list()
  if (verbose) print("Reading in dimensions.")
  for (dim_name in dim_names) {
    if (verbose) print(sprintf("- %s", dim_name))

    # only read in dimensions that have numerical values (therefore try expression)
    try(
      expr = {coordinates[[dim_name]] = ncdf4::ncvar_get(nc, varid = dim_name)},
      silent = TRUE
    )

    if (is.null(coordinates[[dim_name]])) {
      coordinates[[dim_name]] = rep(NA, dim_size[dim_name])
    }
  }

  # check if latitude and longitude are in the netcdf file and which name they have
  lat_name = intersect(dim_names, c("lat", "latitude"))
  if (length(lat_name) == 0) stop("No latitude dimension in netcdf file.")
  lon_name = intersect(dim_names, c("lon", "longitude"))
  if (length(lon_name) == 0) stop("No longitude dimension in netcdf file.")
  stopifnot(length(lat_name) == 1, length(lon_name) == 1)

  # check if the netcdf file has a time dimension
  # to do: implement that time dimension can have any name, e.g. "year"
  if ("time" %in% dim_names) {
    time_df = nctime_get(netcdf_file)
    has_time = TRUE
  } else {
    has_time = FALSE
  }

  if (!has_time) time_format = NULL

  if (!has_time && (!is.null(years)))
    stop("This netcdf file does not have a time dimension.")

  if (!is.null(time_format) || !is.null(years)) {
    # coordinates$time = format(time_df$POSIXct, format = time_format)
    coordinates$time = format(time_df$POSIXct) # need the full time information here
    # the correct time format will be set later in the code (after reading in everything)
    stopifnot(length(coordinates$time) == length(unique(coordinates$time))) # no duplicates
  }

  # define start_idx and end_idx from lat/lon information in grid_cells if applicable
  x = 1
  if (!is.null(grid_cells)) {
    grid_cells$lat = grid_cells[, names(grid_cells) %in% c("lat", "latitude")]
    grid_cells$lon = grid_cells[, names(grid_cells) %in% c("lon", "longitude")]

    lat_idx_start = which(coordinates[[lat_name]] == min(grid_cells$lat))
    lat_idx_end = which(coordinates[[lat_name]] == max(grid_cells$lat))
    lon_idx_start = which(coordinates[[lon_name]] == min(grid_cells$lon))
    lon_idx_end = which(coordinates[[lon_name]] == max(grid_cells$lon))

    if (length(lat_idx_start) == 0 || length(lat_idx_end) == 0 ||
        length(lon_idx_start) == 0 || length(lon_idx_end) == 0) {
      stop("Provided grid cells do not have the same format as netcdf file and cannot be matched. Please check resolution, and the range of latitudes / longitudes in the netcdf file and provided grid cells.")
    }

    lat_count = abs(lat_idx_end - lat_idx_start) + 1
    lon_count = abs(lon_idx_end - lon_idx_start) + 1

    if (lat_idx_end < lat_idx_start) lat_idx_start = lat_idx_end
    if (lon_idx_end < lon_idx_end) lon_idx_start = lon_idx_end

    # create new start_idx and count vector
    start_idx = rep(1, n_dims) # 1: start at the beginning as default
    names(start_idx) = dim_names
    start_idx[lat_name] = lat_idx_start
    start_idx[lon_name] = lon_idx_start

    count = rep(-1, n_dims) # -1: all values are read in as default
    names(count) = dim_names
    count[lat_name] = lat_count
    count[lon_name] = lon_count
  }

  # adjust the time_start_idx and time_count to only read in the desired years
  if (!is.null(years)) {
    stopifnot("year" %in% names(time_df)) # check that there is indeed a year information in the data frame
    time_start_idx = min(which(time_df$year >= min(years)))
    time_end_idx = max(which(time_df$year <= max(years)))
    time_count = time_end_idx - time_start_idx + 1

    # set the start idx for time
    if (!is.null(start_idx)) {
      if (is.null(names(start_idx))) names(start_idx) = dim_names
      start_idx["time"] = time_start_idx
    } else {
      start_idx = rep(1, n_dims) # 1: start at the beginning as default
      names(start_idx) = dim_names
      start_idx[dim_names == "time"] = time_start_idx
    }

    # set the count idx for time
    if (!is.null(count)) {
      if (is.null(names(count))) names(count) = dim_names
      count[dim_names == "time"] = time_count
    } else {
      count = rep(-1, n_dims) # -1: all values are read in as default
      names(count) = dim_names
      count["time"] = time_count
    }
  }

  # test that start and count have the same length as number of dimensions
  stopifnot(is.null(start_idx) || length(start_idx) == n_dims)
  stopifnot(is.null(count) || length(count) == n_dims)

  # adjust coordinates if not whole netcdf is read in
  if (!is.null(start_idx)) {
    for (i in 1:n_dims) {
      if (count[i] == -1) { # if all values to the end are read in
        coordinates[[i]] = coordinates[[i]][ start_idx[i]:length(coordinates[[i]]) ]
      } else {
        coordinates[[i]] = coordinates[[i]][ start_idx[i]:(start_idx[i] + count[i] - 1) ]
      }
    }
    dim_size = sapply(coordinates, length)
  }

  # filter only variable names that are passed to the function
  if (!(variables[1] == "all")) {
    if (!(all(variables %in% var_info$name))) {
      print(sprintf("One of the variables not in the netcdf file: %s", paste(variables, collapse = ", ")))
      print(sprintf("netcdf file: %s", netcdf_file))
      stop()
    }
    var_info = var_info[var_info$name %in% variables, ]
  }
  var_names = var_info$name

  # check if all variables use the same dimensions (then merging not necessary, which can be slow)
  dim_columns = dplyr::select(var_info, dplyr::ends_with(".dim"))
  number_of_dims = apply(dim_columns, 2, function(x) length(unique(x)))

  if (all(number_of_dims == 1)) {
    consistent_dimensions = TRUE
  } else {
    consistent_dimensions = FALSE
  }

  # set start_idx, count to NA, if they are still NULL (that's what's needed in ncvar_get function)
  if (is.null(start_idx)) start_idx = NA
  if (is.null(count)) count = NA

  ################################################################################
  # read in data
  ################################################################################

  # loop through variables and create dataframe
  data_frame = NULL

  if (verbose) print("Reading in data.")
  for (var_number in 1:nrow(var_info)) {
    var_name = var_info[var_number, 2]
    if (verbose) print(sprintf("- %s", var_name))

    var_dimensions = var_info[var_number, ]
    var_dimensions = dplyr::select(var_dimensions, dplyr::ends_with(".dim"))
    var_dimensions = t(var_dimensions)
    var_dimensions = var_dimensions[!is.na(var_dimensions)]

    # test if all var_dimensions were read in into coordinates, if not skip variable
    if (!all(var_dimensions %in% names(coordinates))) {
      print(sprintf("Not all dimension data for variable %s could be read in - skipping this variable.", var_name))
      next
    }

    # extract only start_idx and count for this var
    # reverse as this is how it's done in ncdf4 package
    if (length(start_idx) == 1 && is.na(start_idx)) {
      start_idx_var = NA
    } else {
      start_idx_var = rev(start_idx[var_dimensions])
    }

    if (length(count) == 1 && is.na(count)) {
      count_var = NA
    } else {
      count_var = rev(count[var_dimensions])
    }

    var_vals = ncdf4::ncvar_get(nc = nc, varid = var_name, start = start_idx_var,
                                count = count_var, collapse_degen = FALSE)

    # create named array
    dim_size_var = rev(dim_size[var_dimensions])
    dimnames_var = rev(coordinates[var_dimensions])
    var_vals = array(data = var_vals, dim = dim_size_var, dimnames = dimnames_var)

    if (consistent_dimensions == FALSE) {

      # if dimensions are not consistent (i.e. some variables use other dimensions than other variables, then
      # expand variables --> make them all dependent on the same dimensions by replicating the variable for un-used
      # dimensions)

      # at the moment, however, it looks like as if R reads in the array with ALL dimensions, even if the specific
      # variable does not use them --> I'll check, when I actually have a netcdf with different dimensions
      # (at the moment not necessary to be implemented)
    }

    # melt array into dataframe
    var_data_frame = reshape2::melt(var_vals)
    names(var_data_frame)[names(var_data_frame) == "value"] = var_name

    # cbind variable with result data frame
    if (is.null(data_frame)) {
      data_frame = var_data_frame
    } else {
      data_frame = dplyr::full_join(data_frame, var_data_frame)
      names(data_frame)[ncol(data_frame)] = var_name
    }

  } # loop: var_number


  ################################################################################
  # post processing
  ################################################################################

  # remove grid cells that are not needed (so far, a rectangle has been cut out)
  if (!is.null(grid_cells)) {
    grid_cells[, lat_name] = grid_cells$lat
    grid_cells[, lon_name] = grid_cells$lon # change the names to the ones used in the netcdf file
    grid_cells = grid_cells[, c(lon_name, lat_name)]
    data_frame = dplyr::inner_join(data_frame, grid_cells)
  }

  # remove years that are not needed
  if (!is.null(years)) {
    data_frame$year = lubridate::year(data_frame$time)
    data_frame = data_frame[data_frame$year %in% years, ]
    data_frame$year = NULL
  }

  # if applicable, add time columns (year, month, day, ...) to data frame
  if (return_time_columns) {
    time_df$time = as.factor(time_df$POSIXct)
    data_frame = dplyr::left_join(data_frame, time_df)
  }

  # set the correct time format
  if (!is.null(time_format)) {
    temp_format = data.frame(time = unique(data_frame$time))
    temp_format$new_time = format(as.POSIXct(temp_format$time), format = time_format)
    data_frame$time = temp_format$new_time[match(data_frame$time, temp_format$time)]

  } else if (is.null(time_format) && !is.null(years)) {
    # if no format was given, set to original time values in the dataframe again
    warning("Specific years selected, but no time format for the time vector was provided. The time dimension values in the original netcdf file will be used. Please check if you want to provide a time format to return (e.g. %Y-%m-%d).")
    data_frame$time = time_df$vals[match(data_frame$time, as.character(time_df$POSIXct))]
  }

  # remove NAs (remove lat-lon-times, where no variable has any values)
  if (remove_NA == TRUE & length(var_names) == 1) {
    keep_lines = !is.na(data_frame[, var_names])
    data_frame = data_frame[keep_lines, ]
  } else if (remove_NA == TRUE & length(var_names > 1)) {
    keep_lines = apply(data_frame[, var_names], 1, function(x) any(!is.na(x)))
    data_frame = data_frame[keep_lines, ]
  }

  # order data
  if (has_time == TRUE && lat_name == "lat" && lon_name == "lon") {
    data_frame = dplyr::arrange_(data_frame, "lon", "lat", "time")
  } else if (has_time == TRUE && lat_name == "latitude" && lon_name == "longitude") {
    data_frame = dplyr::arrange_(data_frame, "longitude", "latitude", "time")
  } else if (has_time == FALSE && lat_name == "lat" && lon_name == "lon") {
    data_frame = dplyr::arrange_(data_frame, "lon", "lat", "time")
  } else if (has_time == FALSE && lat_name == "latitude" && lon_name == "longitude") {
    data_frame = dplyr::arrange_(data_frame, "longitude", "latitude")
  }

  if (verbose) print(sprintf("Data frame was created."))

  # close everything
  ncdf4::nc_close(nc) # netcdf4 close

  # end timer
  if (verbose) {
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(sprintf("Time elapsed: %.3f sec", time.taken))
  }

  # return data frame
  return(data_frame)
}

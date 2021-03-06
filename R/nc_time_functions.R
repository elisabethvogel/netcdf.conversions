nctime_get = function(netcdf_file,
                      return_variables = c("POSIXct", "year", "month", "day",
                                           "hour", "minute", "second", "origin",
                                           "unit", "vals"),
                      verbose = FALSE) {

    if (verbose)
        print("function: nctime_get")

    # some tests
    stopifnot(is.character(netcdf_file))
    stopifnot(is.character(return_variables))
    stopifnot(all(return_variables %in% c("POSIXct", "year", "month", "day", "hour", "minute", "second", "origin", "unit", "vals")))

    # test if netcdf has a time variable
    dims = ncdf.tools::infoNcdfDims(netcdf_file)
    if (!("time" %in% dims$name)) {
        warning("Netcdf file does not have a time dimension.")
        return(NULL)
    }

    # open netcdf file
    ncdf = ncdf4::nc_open(filename = netcdf_file)

    # read in time information
    time_vals = ncdf$dim$time$vals
    time_unit = ncdf$dim$time$units

    # with relative time axis (e.g. 'days since ...') :
    if (grepl("since", time_unit)) {

        time_origin = unlist(strsplit(x = time_unit, split = " since "))[2]
        time_unit = unlist(strsplit(x = time_unit, split = " since "))[1]

        # find the time origin with several attempts
        time_origin_date = lubridate::parse_date_time(x = time_origin, orders = c("ymd", "ymdHM", "ymdHMS", "y"))

        # if all values are integers, use lubridate::period instead of lubridate::duration, because the time steps will likely be
        # absolute (e.g. 1 year, 2 year etc. since time_origin)
        if (all(as.integer(time_vals) == time_vals)) {

          time_steps = lubridate::period(paste(time_vals, time_unit))
          time_df = data.frame(POSIXct = time_origin_date + time_steps)

            # if values are not integers, but fractional values, use lubridate::duration
        } else {
          time_steps = lubridate::duration(time_vals, time_unit)
          time_df = data.frame(POSIXct = time_origin_date + time_steps)

        }
        # for absolute time axis with years (other formats will be implemented when needed)
    } else {
        if (time_unit %in% c("year", "years")) {
            time_df = data.frame(POSIXct = lubridate::ymd(sprintf("%i-01-01", time_vals)))
        }
    }

    # add year, month, day, hour, minute, second, origin, unit and vals to time_df
    time_df$year = lubridate::year(time_df$POSIXct)
    time_df$month = lubridate::month(time_df$POSIXct)
    time_df$day = lubridate::day(time_df$POSIXct)
    time_df$hour = lubridate::hour(time_df$POSIXct)
    time_df$minute = lubridate::minute(time_df$POSIXct)
    time_df$second = lubridate::second(time_df$POSIXct)
    time_df$origin = time_origin
    time_df$unit = time_unit
    time_df$vals = time_vals

    # subset to only variables of interest
    time_df = time_df[, return_variables]

    # close netcdf file
    ncdf4::nc_close(ncdf)

    if (verbose)
        print("Time info returned.")
    return(time_df)
}




# This function returns a list with a time unit and time values for a relative
# time vector (days since ...) that can be used for storing in a netcdf file
# Reads in:
#   - time: vector of either characters or POSIXcts
# Returns:
#   - list with: time_unit, time_vals

time2nctime = function(time, orders = c("ymd", "ymdHM", "ymdHMS", "y", "ym")) {

  # some tests
  stopifnot(is.character(time) || lubridate::is.POSIXt(time)  ||
              is.factor(time) || is.numeric(time))

  if (is.factor(time)) {
    time = as.character(time)
    time = lubridate::parse_date_time(time, orders = orders)
  } else if (is.character(time)) {
    time = lubridate::parse_date_time(time, orders = c("ymd", "ymdHM", "ymdHMS", "y", "ym"))
  } else if (is.numeric(time) && all(time >= 1000)) { # arbitrary cut-off, but assumption is that values above 1000 are all years (as opposed to just time steps)
    time = lubridate::parse_date_time(time, orders = "y")
  } # otherwise, time stays a numeric time vector

  if (lubridate::is.POSIXct(time)) {
    origin = as.character(min(time))
    time_vals = time - min(time)
    unit = units(time_vals)

    if (unit == "secs") {
      time_vals = time_vals / (24*60*60)
      unit = "days"
    } else if (unit == "mins") {
      time_vals = time_vals / (24*60)
      unit = "days"
    }
    time_unit = sprintf("%s since %s", unit, origin)

  } else if (is.numeric(time)) {
    time_vals = time
    time_unit = "time steps"
  }

  return(list(time_vals = time_vals, time_unit = time_unit))
}

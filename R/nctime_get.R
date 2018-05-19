nctime_get = function(netcdf_file, return_variables = c("POSIXct", "year", "month", "day", "hour", "minute", "second", "origin", 
    "unit", "vals"), verbose = FALSE) {
    
    if (verbose) 
        print("function: nctime_get")
    
    # some tests
    stopifnot(is.character(netcdf_file), is.readable(netcdf_file))
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
            
            if (time_unit == "years") {
                time_df = data.frame(POSIXct = time_origin_date + lubridate::years(round(time_vals)))
            } else if (time_unit == "months") {
                time_df = data.frame(POSIXct = time_origin_date + months(time_vals))
            } else if (time_unit == "days") {
                time_df = data.frame(POSIXct = time_origin_date + days(time_vals))
            } else if (time_unit == "hours") {
                time_df = data.frame(POSIXct = time_origin_date + hours(time_vals))
            }
            # if values are not integers, but continuous values, use lubridate::duration
        } else {
            time_df = data.frame(POSIXct = time_origin_date + duration(time_vals, time_unit))
        }
        # for absolute time axis with years (other formats will be implemented when needed)
    } else {
        if (time_unit %in% c("year", "years")) {
            time_df = data.frame(POSIXct = ymd(sprintf("%i-01-01", time_vals)))
        }
    }
    
    # add year, month, day, hour, minute, second, origin, unit and vals to time_df
    time_df = mutate(time_df, year = year(POSIXct), month = month(POSIXct), day = day(POSIXct), hour = hour(POSIXct), minute = minute(POSIXct), 
        second = second(POSIXct), origin = time_origin, unit = time_unit, vals = time_vals)
    
    # subset to only variables of interest
    time_df = time_df[, return_variables]
    
    # close netcdf file
    nc_close(ncdf)
    
    if (verbose) 
        print("Time info returned.")
    return(time_df)
}

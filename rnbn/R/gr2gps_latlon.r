#' Covert grid reference to latitude and longitude
#' 
#' @export
#' @param gridref a vector of strings giving grid references to be changed 
#' @param precision gives the precision of the grid references. Must be of
#'        the same length as gridref. precision is given in meters as a numeric
#'        i.e. 10km square = 10000, 1km square = 1000. If NULL then the function
#'        tries to work out the precision itself.
#' @param centre If \code{TRUE} the coordinates for the centre of the cell are given
#'        else if \code{FALSE} the coordinates for the bottom left corner are given. 
#' @return A dataframe of results are returned
#' @keywords mapping latitude longitude grid-reference
#' @examples
#' 
#' gr2gps_latlon('SU616896')

gr2gps_latlon = function(gridref, precision = NULL, centre = TRUE){
    # NOTE FOR UTM30 gridrefs then no need to use helmert transformation as lat long are already in correct projection
    
    # Setup up variable to hold final output
    out_latlon = data.frame(LATITUDE = rep(NA, length(gridref)), LONGITUDE = rep(NA, length(gridref)))
    
    # Determine eastings & northings
    org_en = gr_let2num(gridref, centre = centre, gr_prec = precision, return_projection = TRUE)
    
    # Determine lat lon (original projection)
    out_latlon = OSGridstoLatLong(org_en$EASTING, org_en$NORTHING, org_en$PROJECTION)
    
    # Determine indices of gridrefs that are not UTM30 (which will need transformed)
    i_trans = which(!org_en$PROJECTION %in% c("UTM30","WGS84"))
    
	# Determine whether any supplied grid references need reprojection (i.e. are any grid references from UK, Ireland)
	if(length(i_trans) > 0){
		# Determine Cartesian (Original projection)
		org_cart = LatLong_Cartesian(out_latlon$LATITUDE[i_trans], out_latlon$LONGITUDE[i_trans], org_en$PROJECTION[i_trans])
		
		# Apply Helmert transformation to convert original projections to WGS84 (Cartesian in WGS84)
		helm_tran = helmert_trans(x =org_cart$x, y = org_cart$y, z = org_cart$z, trans = paste(org_en$PROJECTION[i_trans],"toWGS84", sep=""))
		
		# Convert Cartesian coordinates to Latitude/Longitude (Lat Lon in WGS84)
		out_latlon[i_trans,] = Cartesian_LatLong(helm_tran$x, helm_tran$y, helm_tran$z, "UTM30")
    }
    # Return output
    return(out_latlon)
}

#### ADDINTIONAL FUNCTIONS ####

gr_det_country <-
    function(gridref){
        # Create variable to store output
        cty_out = rep(NA, length(gridref))
                
        "(^[[:upper:]]{1,2}[[:digit:]]{2}([[:upper:]]?|[[:upper:]]{2})$)|(^[[:upper:]]{1,2}[[:digit:]]{2,}$)"
        
        # Find British Gridrefs
        cty_out[grepl('(^[[:upper:]]{2}[[:digit:]]{2}([[:upper:]]?|[[:upper:]]{2})$)|(^[[:upper:]]{2}[[:digit:]]{2,}$)',gridref) & !grepl('^(WA)|(WV)',gridref)] = "OSGB"
        
        # Find Irish Gridrefs
        cty_out[grepl('(^[[:upper:]]{1}[[:digit:]]{2}([[:upper:]]?|[[:upper:]]{2})$)|(^[[:upper:]]{1}[[:digit:]]{2,}$)',gridref)] = "OSNI"
        
        # Find Channel Islands Gridrefs
        cty_out[grepl('(^(WA)|(WV)[[:digit:]]{2}([[:upper:]]?|[[:upper:]]{2})$)|(^(WA)|(WV)[[:digit:]]{2,}$)',gridref)] = "UTM30"
        
        # Return output object
        return(cty_out)
    }


gr_let2num <-
    function(gridref, centre = FALSE, gr_prec = NULL, return_projection = FALSE){
        # Function required to calculate easting in Letter Grid
        spmod = function(x, mod){
            ret_obj = x %% mod
            ret_obj[ret_obj == 0] = mod
            return(ret_obj)
        }
        
        # Validate grid reference format (removing any spaces or hyphens)
        gr_comps = gr_components(gridref)
        gridref = gsub("[ -]","",toupper(gridref))
        
        # Setup variable to hold output
        len_grvec = nrow(gr_comps)
        if(return_projection){
            ret_obj = data.frame( EASTING = rep(NA,len_grvec), NORTHING = rep(NA, len_grvec), PROJECTION = rep(NA, len_grvec), row.names = NULL ) # row.names set to null to stop duplicate row names error
        } else {
            ret_obj = data.frame( EASTING = rep(NA,len_grvec), NORTHING = rep(NA, len_grvec), row.names = NULL ) # row.names set to null to stop duplicate row names error
        }
        
        # First find all British Gridrefs
        cty_inds = which(grepl('(^[[:upper:]]{2}[[:digit:]]{2}([[:upper:]]?|[[:upper:]]{2})$)|(^[[:upper:]]{2}[[:digit:]]{2,}$)',gr_comps$VALID_GR) & !grepl('^(WA)|(WV)',gr_comps$VALID_GR))
        
        # If British Gridrefs found then calc easting and northings
        if(length(cty_inds) > 0){
            # Get position of gridref letters in grid
            l1 = match(substr(gr_comps$CHARS[cty_inds],1,1), LETTERS[-9])
            l2 = match(substr(gr_comps$CHARS[cty_inds],2,2), LETTERS[-9])
            # Determine initial easting northing digits based on 500km square
            e = (spmod(l1,5) - 1)*5
            n = floor(abs(l1 - 25)/5)*5
            # Modify initial easting/northing digits based on 100km square
            e = e + (spmod(l2,5) - 1)
            n = n + floor(abs(l2 - 25)/5)
            # Recalulate for false origin (SV) of British Grid
            e = e - 10
            n = n - 5
            # Extend so easting and northing digits so that nchars of east_num/north_num = 5 right padded with zeros
            east_num = gsub(" ","0", format(gr_comps$DIGITS_EAST[cty_inds], width = 5))
            north_num = gsub(" ","0", format(gr_comps$DIGITS_NORTH[cty_inds], width = 5))
            # append numeric part of references to grid index
            e = paste(e,east_num, sep="")
            n = paste(n,north_num, sep="")
            # Overwrite placeholder values in ret_obj
            ret_obj[cty_inds,c("EASTING", "NORTHING")] = c(as.numeric(e), as.numeric(n))
            if(return_projection){
                ret_obj[cty_inds,"PROJECTION"] = "OSGB"
            }
        }
        
        # Find all Irish Gridrefs
        cty_inds = which(grepl('(^[[:upper:]]{1}[[:digit:]]{2}([[:upper:]]?|[[:upper:]]{2})$)|(^[[:upper:]]{1}[[:digit:]]{2,}$)',gr_comps$VALID_GR))
        # If Irish Gridrefs found the calc easting and northings
        if(length(cty_inds) > 0){
            # Get position of gridref letters in grid
            l2 = match(substr(gr_comps$CHARS[cty_inds],1,1), LETTERS[-9])
            # Determine initial easting/northing digits based on 100km square
            e = spmod(l2,5) - 1
            n = floor(abs(l2 - 25)/5)
            # Extend so easting and northing digits so that nchars of east_num/north_num = 5 right padded with zeros
            east_num = gsub(" ","0", format(gr_comps$DIGITS_EAST[cty_inds], width = 5))
            north_num = gsub(" ","0", format(gr_comps$DIGITS_NORTH[cty_inds], width = 5))
            # append numeric part of references to grid index
            e = paste(e,east_num, sep="")
            n = paste(n,north_num, sep="")
            # Overwrite placeholder values in ret_obj
            ret_obj[cty_inds,c("EASTING", "NORTHING")] = c(as.numeric(e), as.numeric(n))
            if(return_projection){
                ret_obj[cty_inds,"PROJECTION"] = "OSNI"
            }
        }
        
        # Find all Channel Islands Gridrefs
        cty_inds = which(grepl('(^(WA)|(WV)[[:digit:]]{2}([[:upper:]]?|[[:upper:]]{2})$)|(^(WA)|(WV)[[:digit:]]{2,}$)',gr_comps$VALID_GR))
        # If CI gridrefs found then calc easting and northings
        if(length(cty_inds) > 0){
            # Determine initial easting/northing based on letters
            e = rep(5, length(cty_inds))
            n = ifelse(grepl('^(WA)[[:digit:]]{2,}$',gr_comps$VALID_GR[cty_inds]),55,54)
            # Extend so easting and northing digits so that nchars of east_num/north_num = 5 right padded with zeros
            east_num = gsub(" ","0", format(gr_comps$DIGITS_EAST[cty_inds], width = 5))
            north_num = gsub(" ","0", format(gr_comps$DIGITS_NORTH[cty_inds], width = 5))
            # append numeric part of references to grid index
            e = paste(e,east_num, sep="")
            n = paste(n,north_num, sep="")
            # Overwrite placeholder values in ret_obj
            ret_obj[cty_inds,c("EASTING", "NORTHING")] = c(as.numeric(e), as.numeric(n))
            if(return_projection){
                ret_obj[cty_inds,"PROJECTION"] = "UTM30"
            }
        }
        
        # If any grid refs contained a tetrad or quadrant code then need to modify the easting and northing
        # Find all grid refs with a tetrad code
        cty_inds = which(!is.na(gr_comps$TETRAD))
        # For these grid refs modify the eastings and northings accordingly
        if(length(cty_inds) > 0){
            # Get list of tetrad codes (all letters but O)
            tet_codes = LETTERS[-15]
            # Determine the position in the tetrad codes vector for each code
            code_match = match(gr_comps$TETRAD[cty_inds], tet_codes)
            # Modify easting
            ret_obj[cty_inds, "EASTING"] = ret_obj[cty_inds,"EASTING"] + ((code_match-1) %/% 5)*2000
            # Modify northing
            ret_obj[cty_inds, "NORTHING"] = ret_obj[cty_inds,"NORTHING"] + ((code_match-1) %% 5)*2000
        }
        
        # Find all grid refs with a quadrant code
        cty_inds = which(!is.na(gr_comps$QUADRANT))
        # For these grid refs modify the eastings and northings accordingly
        if(length(cty_inds) > 0){
            # Get list of tetrad codes (all letters but O)
            quad_codes = c("SW","NW","SE","NE")
            # Determine the position in the tetrad codes vector for each code
            code_match = match(gr_comps$QUADRANT[cty_inds], quad_codes)
            # Modify easting
            ret_obj[cty_inds, "EASTING"] = ret_obj[cty_inds,"EASTING"] + ((code_match-1) %/% 2)*5000
            # Modify northing
            ret_obj[cty_inds, "NORTHING"] = ret_obj[cty_inds,"NORTHING"] + ((code_match-1) %% 2)*5000
        }
        
        # If centre is true the determine precision and give easting and northing for centre of gridref
        if(centre){
            # Determine precision (if gr_prec not supplied)
            if(is.null(gr_prec)){
                gr_prec = gr_comps$PRECISION
            }
            # Add half of precision to easting and northing
            ret_obj[,c("EASTING", "NORTHING")] = ret_obj[,c("EASTING", "NORTHING")] + (gr_prec/2)
        }
        # Return easting and Northings	
        return(ret_obj)
    }


OSGridstoLatLong <-
    function(Easting, Northing, Datum = "OSGB", datum_params = NULL, full_output = FALSE) {
        # If datum_params is null then defaults data.frame included with package will be used
        if(is.null(datum_params)){

            # Set datum_params to datum_vars
            datum_params = datum_vars
        }
        
        # Determine length of Easting & check same as Northing
        east_len = length(Easting)
        if(length(Northing) != east_len){
            stop("ERROR: 'Easting' & 'Northing' are of different lengths")
        }
        
        # Get list of Datums to be converted from Datum
        # If length of Datum input not 1 or same length as Easting/Northing then stop
        if(!length(Datum) %in% c(1,east_len)){
            stop("ERROR: Length of 'Datum' does not match length of Easting/Northing values")
        }
        datum_list = unique(Datum)
        # Check all Datum are in datum_params data frame
        miss_datum = which(!datum_list %in% datum_params$Datum)
        if(length(miss_datum) > 0){
            stop(paste("ERROR: Datum present in data for which parameters have not been given (",paste(sQuote(datum_list[miss_datum]), collapse=","), ")", sep=""))
        }
        
        # Setup object to hold output data
        if(full_output){
            ret_obj = data.frame(EASTING = Easting, NORTHING = Northing, DATUM = Datum, LATITUDE = rep(NA, east_len), LONGITUDE =rep(NA, east_len), stringsAsFactors = FALSE)
        } else {
            ret_obj = data.frame(LATITUDE = rep(NA, east_len), LONGITUDE =rep(NA, east_len), stringsAsFactors = FALSE)
        }
        
        # Loop through datum and extract params and then calculate lat long
        for(i_datum in 1:length(datum_list)){
            
            # Extract indices of easting/northing values corresponding to current datum
            if(length(Datum) == 1){
                dat_inds = 1:east_len
            } else {
                dat_inds = which(Datum == datum_list[i_datum])
            }
            
            # Extract datum params from datum_params (at same time covert lat0 and lon0 from degrees to radians
            par_ind = which(datum_params$Datum == datum_list[i_datum])
            # Check only one row of datum parameters data frame matches current datum
            if(length(par_ind) > 1){
                stop(paste("ERROR: More than one match for current Datum (",sQuote(datum_list[i_datum]),") found in datum_params data frame",sep=""))
            }
            a = datum_params$a[par_ind]
            b = datum_params$b[par_ind]
            F0 = datum_params$F0[par_ind]
            lat0 = datum_params$lat0[par_ind] * (pi/180) # Converting to radians
            lon0 = datum_params$lon0[par_ind] * (pi/180) # Converting to radians
            N0 = datum_params$N0[par_ind]
            E0 = datum_params$E0[par_ind]
            
            
            # Calculate other derived variables for projection/datum
            e2 = 1 - (b^2)/(a^2)		# eccentricity squared
            n = (a-b)/(a+b)
            n2 = n^2   					# Calculate n squared (used several times so store as variable)
            n3 = n^3					# Calculate n cubed (used several times so store as variable)
            
            # Iterate to estimate value of lat and M
            # Intial values
            # Lat
            lat = ( (Northing[dat_inds] - N0) / (a * F0) ) + lat0
            # Meridional Arc
            Ma = ( 1 + n + ((5/4)*n2) + ((5/4)*n3) ) * (lat-lat0)
            Mb = ( (3*n) + (3*n2) + ((21/8)*n3) ) * sin(lat-lat0) * cos(lat+lat0)
            Mc = ( ((15/8)*n2) + ((15/8)*n3) ) * sin(2*(lat-lat0)) * cos(2*(lat+lat0))
            Md = ( (35/24)*n3) * sin(3*(lat-lat0)) * cos(3*(lat+lat0))
            M = (b * F0) * (Ma - Mb + Mc - Md)			# Calculate Developed Meridional arc (M)
            # Setup index for lats that need further iteration
            iter_inds = which(abs( Northing[dat_inds] - N0 - M) >= 0.01)
            while(length(iter_inds) > 0){
                # Lat
                lat[iter_inds] = ( (Northing[dat_inds][iter_inds] - N0 - M[iter_inds]) / (a * F0) ) + lat[iter_inds]
                # Meridional Arc
                Ma[iter_inds] = ( 1 + n + ((5/4)*n2) + ((5/4)*n3) ) * (lat[iter_inds]-lat0)
                Mb[iter_inds] = ( (3*n) + (3*n2) + ((21/8)*n3) ) * sin(lat[iter_inds]-lat0) * cos(lat[iter_inds]+lat0)
                Mc[iter_inds] = ( ((15/8)*n2) + ((15/8)*n3) ) * sin(2*(lat[iter_inds]-lat0)) * cos(2*(lat[iter_inds]+lat0))
                Md[iter_inds] = ( (35/24)*n3) * sin(3*(lat[iter_inds]-lat0)) * cos(3*(lat[iter_inds]+lat0))
                M[iter_inds] = (b * F0) * (Ma[iter_inds] - Mb[iter_inds] + Mc[iter_inds] - Md[iter_inds])			# Calculate Developed Meridional arc (M)
                # Recalculate iteration index
                iter_inds = which(abs( Northing[dat_inds] - N0 - M) >= 0.01)
            }
            
            sinLat = sin(lat)
            tanLat = tan(lat)
            tan2Lat = tanLat^2
            tan4Lat = tanLat^4
            tan6Lat = tanLat^6
            secLat = 1 / cos(lat)
            
            elon = Easting[dat_inds] - E0
            
            nu = (a * F0)*(1-e2*sinLat^2)^-0.5              # transverse radius of curvature
            rho = (a * F0) * (1-e2)*(1-e2*sinLat^2)^-1.5   # meridional radius of curvature
            eta2 = (nu/rho)-1								# East-west component of the deviation of the vertical squared
            
            VII = tanLat / (2*rho*nu)
            VIII = ( tanLat / (24*rho*nu^3) ) * ( 5 + (3*tan2Lat) + eta2 - (9 * tan2Lat * eta2) )
            IX = (tanLat/(720*rho*nu^5)) * (61 + (90*tan2Lat) + (45*tan4Lat))
            X = secLat / nu
            XI = (secLat / (6 * nu^3)) * ( (nu/rho) + 2*tan2Lat)
            XII = (secLat / (120 * nu^5)) * (5 + (28*tan2Lat) + (24*tan4Lat))
            XIIA = (secLat / (5040 * nu^7)) * (61 + (662*tan2Lat) + (1320*tan4Lat) + (720*tan6Lat))
            
            Lat_rad = lat - (VII*elon^2) + (VIII*elon^4) - (IX*elon^6)
            Long_rad = lon0 + (X*elon) - (XI*elon^3) + (XII*elon^5) - (XIIA*elon^7)
            
            Latitude = Lat_rad * (180/pi)
            Longitude = Long_rad * (180/pi)
            
            ret_obj[dat_inds,c("LATITUDE","LONGITUDE")] = data.frame(Latitude, Longitude)
        }
        return (ret_obj)
    }


LatLong_Cartesian <-
    function(Latitude, Longitude, Datum = "OSGB", datum_params = NULL, H = NULL, full_output = FALSE){
        # If datum_params is null then defaults data.frame included with package will be used
        if(is.null(datum_params)){
            # Set datum_params to datum_vars
            datum_params = datum_vars
        }
        
        # Determine length of Latitude & check same as Longitude
        lat_len = length(Latitude)
        if(length(Longitude) != lat_len){
            stop("ERROR: 'Latitude' & 'Longitude' are of different lengths")
        }
        
        # Get list of Datums to be converted from Datum
        # If length of Datum input not 1 or same length as Latitude/Longitude then stop
        if(!length(Datum) %in% c(1,lat_len)){
            stop("ERROR: Length of 'Datum' does not match length of Latitude/Longitude values")
        }
        datum_list = unique(Datum)
        # Check all Datum are in datum_params data frame
        miss_datum = which(!datum_list %in% datum_params$Datum)
        if(length(miss_datum) > 0){
            stop(paste("ERROR: Datum present in data for which parameters have not been given (",paste(sQuote(datum_list[miss_datum]), collapse=","), ")", sep=""))
        }
        
        # Setup object to hold output data (Cartesian Coordinates)
        if(full_output){
            ret_obj = data.frame(LATITUDE = Latitude, LONGITUDE = Longitude, DATUM = Datum, x = rep(NA, lat_len), y = rep(NA, lat_len), z = rep(NA, lat_len),  stringsAsFactors = FALSE)
        } else {
            ret_obj = data.frame(x = rep(NA, lat_len), y = rep(NA, lat_len), z = rep(NA, lat_len), stringsAsFactors = FALSE)
        }
        
        # Loop through datum and extract params and then calculate Cartesian Coordinates
        for(i_datum in 1:length(datum_list)){
            
            # Extract indices of lat/long values corresponding to current datum
            if(length(Datum) == 1){
                dat_inds = 1:lat_len
            } else {
                dat_inds = which(Datum == datum_list[i_datum])
            }
            
            # Extract datum params from datum_params (at same time covert lat0 and lon0 from degrees to radians
            par_ind = which(datum_params$Datum == datum_list[i_datum])
            # Check only one row of datum parameters data frame matches current datum
            if(length(par_ind) > 1){
                stop(paste("ERROR: More than one match for current Datum (",sQuote(datum_list[i_datum]),") found in datum_params data frame",sep=""))
            }
            a = datum_params$a[par_ind]
            b = datum_params$b[par_ind]
            
            # Height (H) assummed to be 0 unless specified
            if(is.null(H)){
                H = 0
            } else {
                H = H
            }
            
            # Convert latitude/longitude values into radians
            lat = Latitude[dat_inds] * (pi/180)
            lon = Longitude[dat_inds] * (pi/180)
            
            # Calculate eccentricity squared (e2)
            e2 = (a^2 - b^2)/a^2
            
            # Calculate Cartezian coordinates 
            v = a / sqrt(1 - e2*sin(lat)^2)
            x = (v + H)*cos(lat)*cos(lon)
            y = (v + H)*cos(lat)*sin(lon)
            z = (((1-e2)*v) + H)*sin(lat)
            
            # Write values to ret_obj
            ret_obj[dat_inds,c("x","y","z")] = data.frame(x,y,z)
            
        }
        return(ret_obj)
    }


helmert_trans <-
    function(x, y, z, trans = "OSNItoOSGB", trans_params = NULL, full_output = FALSE){
        # If trans_params is null then check where helmert_trans_vars data.frame has been loaded from package
        if(is.null(trans_params)){

            # set trans_params = helmert_trans_vars
            trans_params = helmert_trans_vars
        }
        
        # Determine length of input vars x,y,z and check all 3 variables are same length
        len_x = length(x)
        if(length(y) != len_x | length(z) != len_x){
            stop("ERROR: input variables 'x', 'y' and 'z' are of different lengths")
        }
        
        # Get list of trans to be performed from trans
        # If length of trans input not 1 or same length as input vars then stop
        if(!length(trans) %in% c(1,len_x)){
            stop("ERROR: Length of 'trans' does not match length of input 'x','y' & 'z' values")
        }
        trans_list = unique(trans)
        # Check all Datum present in data are in datum_params data frame
        miss_trans = which(!trans_list %in% trans_params$TRANS)
        if(length(miss_trans) > 0){
            stop(paste("ERROR: Transformation requested for which parameters have not been given (",paste(sQuote(trans_list[miss_trans]), collapse=","), ")", sep=""))
        }
        
        # Setup object to hold output data (x,y,z)
        if(full_output){
            ret_obj = data.frame(org_x = x, org_y = y, org_z = z, TRANS = trans, x = rep(NA, len_x), y = rep(NA, len_x), z = rep(NA, len_x), stringsAsFactors = FALSE)
        } else {
            ret_obj = data.frame(x = rep(NA, len_x), y = rep(NA, len_x), z = rep(NA, len_x), stringsAsFactors = FALSE)
        }
        
        # Loop through transformations and extract params and then perform helmert transformation
        for(i_trans in 1:length(trans_list)){
            # Extract indices of x,y,z values corresponding to current transformation
            if(length(trans) == 1){
                dat_inds = 1:len_x
            } else {
                dat_inds = which(trans == trans_list[i_trans])
            }
            # Select relevant helmert transformation parameters
            tp = trans_params[trans_params$TRANS == trans_list[i_trans],]
            # Convert trans paras for rotation from degress to radians
            tp[,c("rx","ry","rz")] = (tp[,c("rx","ry","rz")]/3600) * (pi/180)
            # Normalise s to ppm
            tp[,"s"] = tp$s*1e-6
            # Apply transformation
            x_out = x[dat_inds] + (x[dat_inds]*tp$s) - (y[dat_inds]*tp$rz) + (z[dat_inds]*tp$ry) + tp$tx
            y_out = (x[dat_inds]*tp$rz) + y[dat_inds] + (y[dat_inds]*tp$s) - (z[dat_inds]*tp$rx) + tp$ty
            z_out = (-1*x[dat_inds]*tp$ry) + (y[dat_inds]*tp$rx) + z[dat_inds] + (z[dat_inds]*tp$s) + tp$tz
            
            # Insert transformed cartesian coordincates into output variable
            ret_obj[dat_inds,c("x","y","z")] = data.frame(x = x_out, y = y_out, z = z_out)
        }
        return(ret_obj)
    }


Cartesian_LatLong <-
    function(x,y,z, Datum = "OSGB", datum_params = NULL, full_output = FALSE ){
        # If datum_params is null then defaults data.frame included with package will be used
        if(is.null(datum_params)){
            # Determine if datum variables data.frame is already loaded from package

            # Set datum_params to datum_vars
            datum_params = datum_vars
        }
        
        # Determine length of input vars x,y,z and check all 3 variables are same length
        len_x = length(x)
        if(length(y) != len_x | length(z) != len_x){
            stop("ERROR: input variables 'x', 'y' and 'z' are of different lengths")
        }
        
        # Get list of Datums to be converted from Datum
        # If length of Datum input not 1 or same length as input vars then stop
        if(!length(Datum) %in% c(1,len_x)){
            stop("ERROR: Length of 'Datum' does not match length of input 'x','y' & 'z' values")
        }
        datum_list = unique(Datum)
        # Check all Datum present in data are in datum_params data frame
        miss_datum = which(!datum_list %in% datum_params$Datum)
        if(length(miss_datum) > 0){
            stop(paste("ERROR: Datum present in data for which parameters have not been given (",paste(sQuote(datum_list[miss_datum]), collapse=","), ")", sep=""))
        }
        
        # Setup object to hold output data (Latitude/Longitude)
        if(full_output){
            ret_obj = data.frame(x = x, y = y, z = z, DATUM = Datum, LATITUDE = rep(NA, len_x), LONGITUDE = rep(NA, len_x), HEIGHT = rep(NA, len_x), stringsAsFactors = FALSE)
        } else {
            ret_obj = data.frame(LATITUDE = rep(NA, len_x), LONGITUDE = rep(NA, len_x), stringsAsFactors = FALSE)
        }
        
        # Setup precision variable (used to determine when to stop iteration)
        prec = 1
        
        # Loop through datum and extract params and then calculate Latitude and Longitude
        for(i_datum in 1:length(datum_list)){
            
            # Extract indices of lat/long values corresponding to current datum
            if(length(Datum) == 1){
                dat_inds = 1:len_x
            } else {
                dat_inds = which(Datum == datum_list[i_datum])
            }
            
            # Extract datum params from datum_params (at same time covert lat0 and lon0 from degrees to radians
            par_ind = which(datum_params$Datum == datum_list[i_datum])
            # Check only one row of datum parameters data frame matches current datum
            if(length(par_ind) > 1){
                stop(paste("ERROR: More than one match for current Datum (",sQuote(datum_list[i_datum]),") found in datum_params data frame",sep=""))
            }
            a = datum_params$a[par_ind]
            b = datum_params$b[par_ind]
            
            
            # Calculate eccentricity squared (e2)
            e2 = (a^2 - b^2)/a^2
            
            # Determine longitude
            lon = atan(y[dat_inds]/x[dat_inds])
            # Iteratively determine latitude
            p = sqrt(x[dat_inds]^2 + y[dat_inds]^2)
            lat = atan( z[dat_inds] / (p * (1 - e2)) )
            v = a / sqrt(1 - e2*sin(lat)^2)
            lat2 = atan( (z[dat_inds] + e2*v*sin(lat))/p )
            d_lat = lat - lat2
            iter_inds = which(d_lat > prec)
            while(length(iter_inds) > 0){
                # Change lat to lat2 and recalculate lat2
                lat[iter_inds] = lat2[iter_inds]
                v[iter_inds] = a[dat_inds][iter_inds] / sqrt(1 - e2*sin(lat[iter_inds])^2)
                lat2[iter_inds] = atan( (z[dat_inds][iter_inds] + e2*v[iter_inds]*sin(lat[iter_inds]))/p[iter_inds] )
                d_lat[iter_inds] = lat[iter_inds] - lat2[iter_inds]
                # Recalculate which rows need further iteration
                iter_inds = which(d_lat > prec)
            }
            # Calculate Height (although in majority of cases not really needed/wanted)
            H = (p / cos(lat)) - v
            # Convert Lat / long back to degrees
            lat = lat * (180/pi)
            lon = lon * (180/pi)
            
            # Setup return object
            if(full_output){
                ret_obj[dat_inds,c("LATITUDE","LONGITUDE", "HEIGHT")] = data.frame(lat, lon, H)
            } else {
                ret_obj[dat_inds,c("LATITUDE","LONGITUDE")] = data.frame(lat, lon)
            }
            
        }
        # Return output object
        return(ret_obj)
    }

fmt_gridref = function(gridref, gr_fmt = NULL){
    # Setup object to hold output grid refs
    gr_out = rep(NA, length(gridref))
    # convert gridref string to upper case
    gridref = toupper(gridref)
    # Replace any spaces, punctuation or control characters
    gridref = gsub("[[:space:][:cntrl:][:punct:]]", "", gridref)
    # Check that gridref conforms to grid reference pattern after removals
    # Get indices of gridrefs which are in a valid format
    gr_inds = which(grepl("^[[:upper:]]{1,2}[[:digit:]]{2,}([[:upper:]]?|[[:upper:]]{2})$", gridref))
    # Copy valid gridrefs to output object
    gr_out[gr_inds] = gridref[gr_inds]
    # Extract components where gr_fmt is not NULL (1 = Whole gridref minus tet/quad codes, 2 = Inital letter(s), 3 = Digits only, 4 = Tetrad/Quad only, 5 = Tetrad only, 6 = Quadrant only)
    if(!is.null(gr_fmt)){
        gr_out = gsub("^(([[:upper:]]{1,2})([[:digit:]]{2,}))(([[:upper:]]?)|([[:upper:]]{2}))$", paste("\\",gr_fmt, sep=""), gr_out)
    }
    # Return formatted gridref
    return(gr_out)
}

det_gr_precision <- function(
    gridref
){
    # Convert letters to uppercase
    gridref = toupper(gridref)
    
    # Set up variable to store output
    prec_out = rep(NA,length(gridref))
    
    # Find valid gridrefs
    v_inds = which(grepl("(^[[:upper:]]{1,2}[[:digit:]]{2}([[:upper:]]?|[[:upper:]]{2})$)|(^[[:upper:]]{1,2}[[:digit:]]{2,}$)", gridref) & nchar(gsub("^(([[:upper:]]{1,2})([[:digit:]]{2,}))(([[:upper:]]?)|([[:upper:]]{2}))$", "\\3", gridref)) %% 2 == 0)
    
    # Split into components
    gr_char = gsub("^(([[:upper:]]{1,2})([[:digit:]]{2,}))(([[:upper:]]?)|([[:upper:]]{2}))$", "\\2", gridref[v_inds])
    gr_digits = gsub("^(([[:upper:]]{1,2})([[:digit:]]{2,}))(([[:upper:]]?)|([[:upper:]]{2}))$", "\\3", gridref[v_inds])
    gr_tet = gsub("^(([[:upper:]]{1,2})([[:digit:]]{2,}))(([[:upper:]]?)|([[:upper:]]{2}))$", "\\5", gridref[v_inds])
    gr_quad = gsub("^(([[:upper:]]{1,2})([[:digit:]]{2,}))(([[:upper:]]?)|([[:upper:]]{2}))$", "\\6", gridref[v_inds])
    
    # Determine number of digits pairs
    n_pairs = nchar(gr_digits)/2
    
    # Determine precison based on gr
    gr_prec = 10^5 / 10^n_pairs
    
    # If gr_tet contains valid letter then ignore precision based on gridref length and assign 2000
    gr_prec[gr_tet %in% LETTERS[-15]] = 2000
    # If not valid tetrad code in gr_tet then set to NA (i.e. O is not a valid tetrad code)
    gr_prec[!gr_tet %in% LETTERS[-15] & gr_tet != ""] = NA
    
    # If gr_quad contains valid letter then ignore precision based on gridref length and assign 5000
    gr_prec[gr_quad %in% c("NW","NE","SW","SE")] = 5000
    # If not valid quadrant code in gr_tet then set to NA
    gr_prec[!gr_quad %in% c("NW","NE","SW","SE") & gr_quad != ""] = NA
    
    # Write gr_prec values to output variable
    prec_out[v_inds] = gr_prec
    
    # Return output variable
    return(prec_out)
}
gr_components <- function(gridref, output_col = NULL){
    
    # Convert letters to uppercase
    gridref = toupper(gridref)
    
    # Set up variable to store output
    gr_comps = data.frame(GRIDREF = gridref, VALID_GR = NA, PRECISION = NA, CHARS = NA, DIGITS = NA, DIGITS_EAST = NA, DIGITS_NORTH = NA, TETRAD = NA, QUADRANT = NA)
    
    # Check values for output_col if supplied
    if(!is.null(output_col)){	
        if(!all(toupper(output_col) %in% names(gr_comps))){
            stop("Supplied output column name not recognised, valid values are:\n\t", paste(shQuote(names(gr_comps)), collapse = ", "))
        }
    }
    
    # Find valid gridrefs
    v_inds = which(grepl("(^[[:upper:]]{1,2}[[:digit:]]{2}([[:upper:]]?|[[:upper:]]{2})$)|(^[[:upper:]]{1,2}[[:digit:]]{2,}$)", gsub("[ -]","",gridref)) & nchar(gsub("^(([[:upper:]]{1,2})([[:digit:]]{2,}))(([[:upper:]]?)|([[:upper:]]{2}))$", "\\3", gsub("[ -]","",gridref))) %% 2 == 0)
    
    if(length(v_inds) > 0){
        
        # Update valid gridrefs to data.frame
        gr_comps[v_inds,"VALID_GR"] = gsub("[ -]","", gridref[v_inds])
        
        # Split into components
        # Characters
        gr_comps[v_inds,"CHARS"] = gsub("^(([[:upper:]]{1,2})([[:digit:]]{2,}))(([[:upper:]]?)|([[:upper:]]{2}))$", "\\2", gr_comps$VALID_GR[v_inds])
        # Digits
        gr_comps[v_inds,"DIGITS"] = gsub("^(([[:upper:]]{1,2})([[:digit:]]{2,}))(([[:upper:]]?)|([[:upper:]]{2}))$", "\\3", gr_comps$VALID_GR[v_inds])
        # split digits into east and north
        len_digit = nchar(gr_comps$DIGITS[v_inds])
        gr_comps[v_inds, "DIGITS_EAST"] = substr(gr_comps$DIGITS[v_inds], 1, len_digit/2)
        gr_comps[v_inds, "DIGITS_NORTH"] = substr(gr_comps$DIGITS[v_inds], (len_digit/2)+1, len_digit)
        # Determine precision based on digits
        gr_comps[v_inds,"PRECISION"] = 10^5 / 10^(len_digit[v_inds]/2)
        # Tetrad
        gr_comps[v_inds,"TETRAD"] = gsub("^(([[:upper:]]{1,2})([[:digit:]]{2,}))(([[:upper:]]?)|([[:upper:]]{2}))$", "\\5", gr_comps$VALID_GR[v_inds])
        # Find blank string values in gr_tet and replace with NA
        na_inds = which(!gr_comps$TETRAD %in% LETTERS[-15])
        if(length(na_inds) > 0){
            gr_comps[na_inds,"TETRAD"] = NA
        }
        # Quadrant
        gr_comps[v_inds,"QUADRANT"] = gsub("^(([[:upper:]]{1,2})([[:digit:]]{2,}))(([[:upper:]]?)|([[:upper:]]{2}))$", "\\6", gr_comps$VALID_GR[v_inds])
        # Find blank string values in gr_tet and replace with NA
        na_inds = which(!gr_comps$QUADRANT %in% c("NW","NE","SW","SE"))
        if(length(na_inds) > 0){
            gr_comps[na_inds,"QUADRANT"] = NA
        }
        
        # Modify precision if tetrad or quadrant is not null
        # Tetrad
        gr_comps[grepl("^[[:upper:]]{1,2}[[:digit:]]{2}[[:upper:]]{1,2}$", gr_comps$VALID_GR) & !is.na(gr_comps$TETRAD),"PRECISION"] = 2000
        # quadrant
        gr_comps[grepl("^[[:upper:]]{1,2}[[:digit:]]{2}[[:upper:]]{1,2}$", gr_comps$VALID_GR) & !is.na(gr_comps$QUADRANT),"PRECISION"] = 5000
        
        # Modify valid_gr to remove any gridrefs that had non-valid tetrad/quadrant codes
        na_inds = which(grepl("^[[:upper:]]{1,2}[[:digit:]]{2}[[:upper:]]{1,2}$", gr_comps$VALID_GR) & is.na(gr_comps$TETRAD) & is.na(gr_comps$QUADRANT))
        if(length(na_inds) > 0){
            gr_comps[na_inds,-1] = NA
        }
    }		
    
    # If output_col is null then output all columns otherwise only output requested columns
    if(is.null(output_col)){
        out_obj = gr_comps
    } else {
        out_obj = gr_comps[,toupper(output_col)]
    }
    
    # Return output object
    return(out_obj)
}

det_tet_quad <- function(
    gridref,
    precision = NULL,
    prec_out = NULL
){
    # Setup output object
    out_obj = data.frame(TETRAD_GR = rep(NA, length(gridref)), QUADRANT_GR = NA)
    
    # Split gridref into components (will also check gridref sytax is correct)
    gr_comps = gr_components(gridref)
    
    gr_inds = which(nchar(gr_comps$DIGITS) >= 2 & !is.na(gr_comps$VALID_GR))	
    sq10 = rep(NA, length(gridref))
    sq10[gr_inds] = paste(gr_comps$CHARS[gr_inds], substr(gr_comps$DIGITS_EAST[gr_inds],1,1), substr(gr_comps$DIGITS_NORTH[gr_inds],1,1), sep="")
    
    # All gridrefs that have valid tetrad/quadrant codes can be filled directly to out_obj
    # Tetrads
    gr_inds = which(!is.na(gr_comps$TETRAD))
    if(length(gr_inds) > 0){
        out_obj[gr_inds,"TETRAD_GR"] = paste(sq10[gr_inds], gr_comps$TETRAD[gr_inds], sep="")
    }
    #Quadrants
    gr_inds = which(!is.na(gr_comps$QUADRANT))
    if(length(gr_inds) > 0){
        out_obj[gr_inds,"QUADRANT_GR"] = paste(sq10[gr_inds], gr_comps$QUADRANT[gr_inds], sep="")
    }
    
    # If original gridref had tetrad code then where can be mapped directly to quadrant then do so
    # Find gridrefs with calculated precision of 2000 i.e. those that had a tetrad code in the grid ref
    gr_inds = which(gr_comps$PRECISION == 2000)
    if(length(gr_inds) > 0){
        # SW
        gr_inds = which(gr_comps$PRECISION == 2000 & gr_comps$TETRAD %in% c("A","B","F","G"))
        if(length(gr_inds) > 0){
            out_obj[gr_inds,"QUADRANT_GR"] = paste(sq10[gr_inds], "SW", sep="") 
        }
        # NW
        gr_inds = which(gr_comps$PRECISION == 2000 & gr_comps$TETRAD %in% c("D","E","J","I"))
        if(length(gr_inds) > 0){
            out_obj[gr_inds,"QUADRANT_GR"] = paste(sq10[gr_inds], "NW", sep="") 
        }
        # SE
        gr_inds = which(gr_comps$PRECISION == 2000 & gr_comps$TETRAD %in% c("Q","R","V","W"))
        if(length(gr_inds) > 0){
            out_obj[gr_inds,"QUADRANT_GR"] = paste(sq10[gr_inds], "SE", sep="") 
        }
        # NE
        gr_inds = which(gr_comps$PRECISION == 2000 & gr_comps$TETRAD %in% c("T","U","Y","Z"))
        if(length(gr_inds) > 0){
            out_obj[gr_inds,"QUADRANT_GR"] = paste(sq10[gr_inds], "NE", sep="") 
        }
    }		
    
    
    # Find all valid gridrefs where precision is high enough
    gr_inds = which(gr_comps$PRECISION <= 1000)
    
    if(length(gr_inds) > 0){
        
        # Extract relevant digits for tetrad/quadrant estimation
        code_e = as.numeric(substr(gr_comps$DIGITS_EAST[gr_inds],2,2))
        code_n = as.numeric(substr(gr_comps$DIGITS_NORTH[gr_inds],2,2))
        
        # Determine tetrad code from easting and northing digits (can perhaps
        # Get list of tetrad codes (no O)
        tet_codes = LETTERS[-15]
        # Get tetrad code for each gridref of suitable precision
        tets = tet_codes[(code_e%/%2)*5 + (code_n%/%2)+1]
        out_obj[gr_inds,"TETRAD_GR"] = paste(sq10[gr_inds],tets, sep="")
        
        # Determine quadrant codes from easting and northing digits
        # Build vector of quadrant codes (in correct order)
        quad_codes = c("SW","NW","SE","NE")
        # Get tetrad code for each gridref of suitable precision
        quads = quad_codes[(code_e%/%5)*2 + (code_n%/%5)+1]
        out_obj[gr_inds,"QUADRANT_GR"] = paste(sq10[gr_inds],quads, sep="")
        
    } 
    
    # If precision is supplied determine if it matches value estimated from gridref
    if(!is.null(precision)){
        # Check that precision is either a vector of length gridref or a single value
        if( !(length(precision) == 1 | length(precision) == length(gridref)) ){
            stop("precision does not match gridref length, supply either single value or vector of same length as gridref")
        }
        # Find non-matching precisions
        # Find non-matching non padded gridrefs (i.e. ones that shouldn't be wrong)
        gr_inds = which(gr_comps$PRECISION != precision & gr_comps$PRECISION != 1000 & (precision != 2000 | precision != 5000) )
        # Find all non-matched precision values
        if(length(gr_inds) > 0){
            stop("Precision supplied does not match determined precision for", length(gr_inds), " grid refs")
        }
        # If estimated precision = 1000 and supplied prec = 2000 or 5000 then assume tetrad/quadrant padded values (e.g. BRC type data)
        # Padded tetrads
        # Note if tetrad grid ref and tetrad code in (C,H,K,L,M,N,P,S,X) then quadrant cannot be uniquely determined
        # Find grid refs matching this criteria
        gr_inds = which(gr_comps$PRECISION != precision & gr_comps$PRECISION == 1000 & precision == 2000 & grepl("^[[:alpha:]]{1,2}[[:digit:]]{2}[CHKLMNPSX]$",out_obj$TETRAD_GR))
        # Set quadrant to NA
        if(length(gr_inds) > 0){
            out_obj[gr_inds,"QUADRANT_GR"] = NA
        }
        # Padded Quadrant
        # For all padded quadrant codes tetrad cannot be uniquely determined
        gr_inds = which(gr_comps$PRECISION != precision & gr_comps$PRECISION == 1000 & precision == 5000)
        # Set quadrant to NA
        if(length(gr_inds) > 0){
            out_obj[gr_inds,"TETRAD_GR"] = NA	
        }
        
    }		
    
    # Return output object based on prec_out value
    if(is.null(prec_out)){
        return(out_obj)
    } else if(prec_out == 2000){
        return(out_obj$TETRAD_GR)
    } else if(prec_out == 5000){
        return(out_obj$QUADRANT_GR)
    }
}
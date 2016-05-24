# Copyright (C) 2014 - 2015  James Balamuta, Stephane Guerrier, Roberto Molinari
#
# This file is part of GMWM R Methods Package
#
# The `gmwm` R package is free software: you can redistribute it and/or modify it
# under the terms of the Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
# included within the packages source as the LICENSE file.
#
# The `gmwm` R package is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the Attribution-NonCommercial-ShareAlike 4.0 International 
# (CC BY-NC-SA 4.0) along with `gmwm`.  If not, see <http://www.smac-group.com/licensing/>.


#' Latest Version of Package on CRAN
#' 
#' Determines the version number on cran via obtaining the packages page
#' 
#' @param pkg         Name of Package
#' @param cran_url    URL to CRAN Packages
#' @return A \code{vector} of \code{string}s that contain:
#' \itemize{
#' \item the verison of the package on cran
#' \item release date of the package on cran
#' }
#' @examples 
#' library(gmwm)
#' packageVersion("gmwm")
#' packageVersionCRAN("gmwm")
#' 
#' @author JJB
packageVersionCRAN = function(pkg, cran_url="http://cran.r-project.org/web/packages/")
{
  
  cran_pkg_loc = paste0(cran_url,pkg)
  
  suppressWarnings( conn <- try( url(cran_pkg_loc) , silent=TRUE ) )
  
  if ( all( class(conn) != "try-error") ) {
    suppressWarnings( cran_pkg_page <- try( readLines(conn) , silent=TRUE ) )
    close(conn)
  } else {
    return(NULL)
  }
  
  version_line = cran_pkg_page[grep("Version:",cran_pkg_page)+1]
  version_num_cran = gsub("<(td|\\/td)>","",version_line)
  
  publish_line = cran_pkg_page[grep("Published:",cran_pkg_page)+1]
  publish_date_cran = gsub("<(td|\\/td)>","",publish_line)
  
  c(version_num_cran,publish_date_cran)  
}


#' GM Conversion
#' 
#' Convert from AR1 to GM and vice-versa
#' @param theta        A \code{numeric vector} containing the theta values
#' @param process.desc A \code{character vector} containing the names of parameters.
#' @param freq         A \code{double} indicating the frequency of the data.
#' @keywords internal
#' @author JJB
#' @rdname gm_conv
conv.ar1.to.gm = function(theta, process.desc, freq){
  idx = process.desc %in% c("BETA","SIGMA2_GM")
  theta[idx] = ar1_to_gm(theta[idx],freq)
  
  theta
}

#' @rdname gm_conv
conv.gm.to.ar1 = function(theta, process.desc, freq){
  idx = process.desc %in% c("BETA","SIGMA2_GM")
  theta[idx] = gm_to_ar1(theta[idx],freq)
  
  theta
}


#' @title Print GMWM Data Object
#' @description 
#' Pretty formatting for \code{gts}, \code{imu}, and \code{lts} objects.
#' @param x         A \code{gts}, \code{imu}, \code{lts} object.
#' @param obs       A \code{integer} the specifies how many from the beginning and end of the data set to show.
#' @param row.names A \code{boolean} that indicates whether row names should be displayed or surpressed.
#' @param ...       Further arguments passed to or from other methods.
#' @return 
#' A \code{logical} value that indicates whether the object is of that class (TRUE) or not (FALSE).
#' @author JJB
#' @rdname print_data
#' @export
print.imu = function(x,
                     obs = 10L,
                     row.names = TRUE, ...)
{
  if(!is.null(attr(x,"name"))){
    cat("Data Name:",attr(x,"name"),"\n")
  }
  if(!is.null(attr(x,"stype"))){
    cat("Sensor:",attr(x,"stype"),"@",attr(x,"freq"),"Hz\n")
    cat("Obs:", nrow(x), " over ", round(nrow(x)/attr(x,"freq")/3600,2),"Hours \n")
  }else{
    cat("Freq:",attr(x,"freq"),"Hz\n")
  }
  outf(x, obs, row.names, ...)
}

#' @rdname print_data
#' @export
print.lts = function(x,
                     obs = 10L,
                     row.names = TRUE, ...)
{
  outf(x, obs, row.names, ...)
}

#' @rdname print_data
#' @export
print.gts = function(x,
                     obs = 10L,
                     row.names = TRUE, ...)
{
  outf(x, obs, row.names, ...)
}

#' @rdname print_data
outf = function(x, obs = 10L, row.names = TRUE){
  if(!is.numeric(obs)){ obs = 100L }
  if(!is.infinite(obs)){ obs = as.integer(obs) }
  
  if (obs*2 < nrow(x)) {
    print_lines = rbind(head(x,obs), tail(x, obs))
    rn = c(seq_len(obs), seq.int(to=nrow(x), length.out=obs))
    print_dashes = TRUE
  } else {
    print_lines = head(x,nrow(x))
    rn = seq_len(nrow(x))
    print_dashes = FALSE
  }
  
  if (isTRUE(row.names)){
    rownames(print_lines) = paste(format(rn,right=TRUE,scientific=FALSE),":",sep="")
  }else{
    rownames(print_lines) = rep("", nrow(print_lines))
  }
  if(is.null(colnames(x))){
    colnames(print_lines) = rep("NA", ncol(print_lines))
  }
  if(print_dashes) {
    print_lines = rbind(head(print_lines,obs),"---"="",tail(print_lines,obs))
    rownames(print_lines) = format(rownames(print_lines),justify="right")
  }
  
  print.default(print_lines,right=TRUE,quote=FALSE)
  return(invisible())
}

#' @title Is GMWM Object
#' @description 
#' Is the object a
#' \code{gts}, \code{imu}, \code{lts}, \code{wvar}, or \code{gmwm} object?
#' @param x  A \code{gts}, \code{imu}, \code{lts} object.
#' @return 
#' A \code{logical} value that indicates whether the object is of that class (TRUE) or not (FALSE).
#' @details
#'  Uses \code{\link[base]{inherits}} over \code{\link[methods]{is}} for speed. 
#' @author JJB
#' @rdname is_func
#' @export
is.gts = function(x){ inherits(x, "gts") }

#' @rdname is_func
#' @export
is.imu = function(x){ inherits(x, "imu") }

#' @rdname is_func
#' @export
is.lts = function(x){ inherits(x, "lts") }

#' @rdname is_func
#' @export
is.gmwm = function(x){ inherits(x, "gmwm") }

#' @rdname is_func
#' @export
is.wvar = function(x){ inherits(x, "wvar") }

#' @rdname is_func
#' @export
is.ts.model = function(x){ inherits(x, "ts.model") }


#' @title Obtain the value of an object's properties
#' @description 
#' Used to access different properties of the
#'  \code{gts}, \code{imu}, or \code{lts} object.
#' @param x      A \code{gts}, \code{imu}, or \code{lts} object.
#' @param type   A \code{string} indicating the field to be retrieved.
#' @return 
#' The method will return a single numeric or string result depending on the
#' slot being accessed.
#' @details 
#' To access information about \code{imu} properties use:
#' \describe{
#'  \item{\code{"accel"}}{Returns the number of accelerometers}
#'  \item{\code{"gyro"}}{Returns the number of gyroscopes}
#'  \item{\code{"sensors"}}{Returns total number of sensors}
#' }
#' @author JJB
#' @export
value = function(x, type){
  UseMethod("value")
}

#' @describeIn value Access \code{imu} object properties
#' @export
value.imu = function(x, type){
  switch(type,
         accel   = attr(x, 'num.sensor')[1],
         gyro    = attr(x, 'num.sensor')[2],
         sensors = sum(attr(x, 'num.sensor')),
         stop("The `type` specified is not an available slot")
  ) 
}

#' @title Obtain the value of an object's properties
#' @description 
#' Used to access different properties of the
#'  \code{gts}, \code{imu}, or \code{lts} object.
#' @param x      A \code{gts}, \code{imu}, or \code{lts} object.
#' @param type   A \code{string} indicating the field to be retrieved.
#' @return 
#' The method will return a single TRUE or FALSE response
#' @details 
#' To access information about \code{imu} properties use:
#' \describe{
#'  \item{\code{"accel"}}{Returns whether accelerometers have been specified}
#'  \item{\code{"gyro"}}{Returns whether accelerometers have been specified}
#'  \item{\code{"sensors"}}{Returns whether there exists both types of sensors}
#' }
#' @author JJB
#' @export
has = function(x, type){
  UseMethod("has")
}

#' @describeIn has Access \code{imu} object properties
#' @export
has.imu = function(x, type){
  switch(type,
         accel   = attr(x, 'num.sensor')[1] > 0,
         gyro    = attr(x, 'num.sensor')[2] > 0,
         sensors = attr(x, 'num.sensor')[1] > 0 & attr(x, 'num.sensor')[2] > 0,
         stop("The `type` specified is not an available slot")
  ) 
}
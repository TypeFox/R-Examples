# Purpose        : Extracts PROJ4 parameters and checks if they are compatible with ref_CRS
# Maintainer     : Pierre Roudier (pierre.roudier@landcare.nz);
# Contributions  : Tomislav Hengl (tom.hengl@wur.nl); Dylan Beaudette (debeaudette@ucdavis.edu);  
# Dev Status     : Alpha
# Note           : p4s_parameters list of proj4 parameter/value strings; Uses the string parsing functionality from the 'plyr' package;

extractProjValue <- function(p4s_parameters, param){
  
  # Locating the current PROJ4 parameter
  query <- ldply(p4s_parameters, str_locate, pattern = param)
  idx <- which(!is.na(query[, 1]) & !is.na(query[, 2]))
  
  # If the PROJ4 parameter is found we extract its value
  if (length(idx) > 0) {
    # Selecting the good param string
    param_value <- p4s_parameters[idx]
    # Extract the parameter value
    value <- strsplit(param_value, param)[[1]]
    value <- value[value != ""]
  }
  else { stop(paste("Proj4string does not contain", param, "parameter.\n Consider converting to the referent CRS", get("ref_CRS", envir = plotKML.opts),"manually."))
  }

  return(value)
}


## parse string:
parse_proj4 <- function(p4s, params){

  if(missing(params)) {
  ref_CRS = get("ref_CRS", envir = plotKML.opts)
  value <- strsplit(ref_CRS, "\\+")[[1]]
  value <- value[value != ""]
  param_names <- sapply(strsplit(value, "="), function(x){x[1]})
  params <- as.list(paste("\\+", sapply(strsplit(value, "="), function(x){x[1]}), "=", sep="")) 
  }

  # Splitting the whole PROJ4 string
  p4s_parameters <- str_split(p4s, " ")[[1]]
  # Extraction of the values of parameters specified above
  x <- laply(params, extractProjValue, p4s_parameters = p4s_parameters)
  # colnames for better looking result
  value <- sapply(sapply(params, strsplit, "\\+"), function(x){x[2]})
  param_names <- sapply(strsplit(value, "="), function(x){x[1]})
  names(x) <- param_names

  return(x)
}

## Get proj4string from an object
getCRS.Spatial <- function(obj) {
  CRSargs(CRS(proj4string(obj)))
}

setMethod("getCRS", "Spatial", getCRS.Spatial)

getCRS.Raster <- function(obj) {
  CRSargs(projection(obj, asText = FALSE))
}

setMethod("getCRS", "Raster", getCRS.Raster)


## check projection for Raster objects
setMethod("is.projected", signature(obj = "Raster"),
	function(obj) {
		p4str <- getCRS(obj)
		if (is.na(p4str) || nchar(p4str) == 0) 

			return(as.logical(NA))
		else {

			x <- grep("longlat", p4str, fixed=TRUE)

			if (length(x) == 0)
				return(TRUE)
			else
				return(FALSE)
		}
	}
)

## check proj4string
check_projection <- function(obj, control = TRUE, ref_CRS = get("ref_CRS", envir = plotKML.opts)){
  
  if(is.na(proj4string(obj))){
    stop("Proj4 string missing")
  } 

  #  First, check if it is in the metric system or unprojected:
  if(ref_CRS=="+proj=longlat +datum=WGS84"&is.projected(obj)){
    ret = FALSE
  }
  
  else {

  # Using PROJ.4 to get the PROJ4 string
  p4s <- getCRS(obj)

  # Parsing the PROJ4 string for proj and datum values
  params <- parse_proj4(p4s)

  # the default target proj4 string:
  value <- strsplit(ref_CRS, "\\+")[[1]]
  value <- value[value != ""]
  target_params <- stringr::str_trim(sapply(strsplit(value, "="), function(x){x[2]}))
  names(target_params) <- sapply(strsplit(value, "="), function(x){x[1]})

  # if already projection type is missing the string is invalid
  if(params["proj"] != ""){  

  # If test fails
  if (sum(is.na(match(params, target_params)))>0) {
    if (control==FALSE)
      stop(paste("'", ref_CRS, "' coordinate system required"))
    else
      ret <- FALSE
  }
  # If test succeed
  else {
      ret <- TRUE
    }
  }
  
  else {
    stop("A valid proj4string required. See 'CRS-methods' for more details.")
    }
  }

  return(ret)
}

# end of script;

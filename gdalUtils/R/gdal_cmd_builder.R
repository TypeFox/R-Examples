#' gdal_cmd_builder
#' 
#' Helper function for building GDAL commands.
#' 
#' @param executable Character. The GDAL command to use (e.g. "gdal_translate")
#' @param parameter_variables List. A list of parameter names, organized by type.
#' @param parameter_values List. A list of the parameters names/values.
#' @param parameter_order Character. The order of the parameters for the GDAL command.
#' @param parameter_noflags Character. Parameters which do not have a flag.
#' @param parameter_doubledash Character. Parameters which should have a double dash "--".
#' @param parameter_noquotes Character. Parameters which should not be wrapped in quotes (vector parameters only, at present).
#' @param gdal_installation_id Numeric. The ID of the GDAL installation to use.  Defaults to 1.
#' 
#' @return Formatted GDAL command for use with system() calls. 
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net})
#' 
#' @details This function takes the executable name (e.g. "gdal_translate"),
#' a list of parameter names organized by logical, vector,
#' scalar, character, repeatable, a list of values of these parameters, 
#' the order they should be used in the GDAL command, and a list of
#' parameters that should not have a flag, and returns a properly
#' formatted GDAL command (with the full path-to-executable) that
#' should work with a system() call.
#' 
#' Sometimes, a user may not want to use the most recent GDAL install
#' (gdal_installation_id=1), so the gdal_installation_id can be used
#' to set a different install.  This is often used with gdal_chooseInstallation
#' if, for instance, the particular GDAL installation required needs
#' a specific driver that may not be available in all installations.
#' 
#' In general, an end user shouldn't need to use this function -- it
#' is used by many of the GDAL wrappers within gdalUtils.
#'
#' @references \url{http://www.gdal.org/gdal_translate.html}
#' @examples \dontrun{ 
#' # This builds a gdal_translate command.
#' executable <- "gdal_translate"
#' 
#' parameter_variables <- list(
#' 			logical = list(
#' 					varnames <- c("strict","unscale","epo",
#' 					"eco","q","sds","stats")),
#' 			vector = list(
#' 					varnames <- c("outsize","scale","srcwin",
#' 					"projwin","a_ullr","gcp")),
#' 			scalar = list(
#' 					varnames <- c("a_nodata")),
#' 			character = list(
#' 					varnames <- c("ot","of","mask","expand","a_srs",
#' 					"src_dataset","dst_dataset")),
#' 			repeatable = list(
#' 					varnames <- c("b","mo","co")))
#' 
#' parameter_order <- c(
#' 			"strict","unscale","epo","eco","q","sds","stats",
#' 			"outsize","scale","srcwin","projwin","a_ullr","gcp",
#' 			"a_nodata",
#' 			"ot","of","mask","expand","a_srs",
#' 			"b","mo","co",
#' 			"src_dataset","dst_dataset")
#' 
#' parameter_noflags <- c("src_dataset","dst_dataset")
#' 
#' # Now assign some parameters:
#' parameter_values = list(
#' 	src_dataset = "input.tif",
#' 	dst_dataset = "output.envi",
#' 	of = "ENVI",
#' 	strict = TRUE
#' )
#' 
#' cmd <- gdal_cmd_builder(
#' 			executable=executable,
#' 			parameter_variables=parameter_variables,
#' 			parameter_values=parameter_values,
#' 			parameter_order=parameter_order,
#' 			parameter_noflags=parameter_noflags)
#' 
#' cmd
#' system(cmd,intern=TRUE) 
#' }
#' @export

#TODO: additional commands
#TODO: work without parameters (executable only)
#TODO: command aliases (e.g. for commands with a hyphen, 
#	since R doesn't allow that in variable naming).

gdal_cmd_builder <- function(executable,parameter_variables=c(),
		parameter_values=c(),parameter_order=c(),parameter_noflags=c(),
		parameter_doubledash=c(),
		parameter_noquotes=c(),
		gdal_installation_id=1)
{
	# path to executable check in here?
	
	gdal_setInstallation()
	if(is.null(getOption("gdalUtils_gdalPath"))) return()
	
	executable <- normalizePath(list.files(
					getOption("gdalUtils_gdalPath")[[gdal_installation_id]]$path,
					executable,full.names=TRUE))
	
	parameter_variables_types <- names(parameter_variables)
	defined_variables <- names(parameter_values)[sapply(parameter_values,function(X) class(X) != "name")]
	
	if(any("logical" %in% parameter_variables_types))
	{
		parameter_variables_logical <- parameter_variables$logical[[1]]
		parameter_variables_logical_defined <- defined_variables[defined_variables %in% parameter_variables_logical]
		# Only set the flag if TRUE
		if(length(parameter_variables_logical_defined)>0)
		{
			parameter_variables_logical_defined_true <- sapply(parameter_variables_logical_defined,
					function(X,parameter_values)
					{
						return(parameter_values[[which(names(parameter_values)==X)]])
					},parameter_values=parameter_values)
			
			parameter_variables_logical_strings <- sapply(parameter_variables_logical_defined,
					function(X,parameter_doubledash)
					{
						if(X %in% parameter_noflags)
						{
							flag=NULL
						} else
						{
							if(X %in% parameter_doubledash)
							{
								flag=paste("--",X," ",sep="")	
							} else
							{
								flag=paste("-",X," ",sep="")
							}
						}
						return(flag)
					},parameter_doubledash=parameter_doubledash)	
			names(parameter_variables_logical_strings) <- names(parameter_variables_logical_defined_true)
		} else
		{
			parameter_variables_logical_strings <- NULL
		}
		
	}
	
	if(any("vector" %in% parameter_variables_types))
	{
		parameter_variables_vector <- parameter_variables$vector[[1]]
		parameter_variables_vector_defined <- defined_variables[defined_variables %in% parameter_variables_vector]
		if(length(parameter_variables_vector_defined)>0)
		{
			parameter_variables_vector_strings <- sapply(parameter_variables_vector_defined,
					function(X,parameter_values,parameter_doubledash)
					{
						if(X %in% parameter_noflags)
						{
							flag=NULL
						} else
						{
							if(X %in% parameter_doubledash)
							{
								flag=paste("--",X," ",sep="")	
							} else
							{
								flag=paste("-",X," ",sep="")
							}
						}
						
						if(X %in% parameter_noquotes)
						{
							parameter_variables_vector_string <- paste(flag,
									paste(parameter_values[[which(names(parameter_values)==X)]],collapse=" "),
									sep="")
						} else
						{						
							parameter_variables_vector_string <- paste(flag,
									qm(paste(parameter_values[[which(names(parameter_values)==X)]],collapse=" ")),
									sep="")
						}
						return(parameter_variables_vector_string)
					},parameter_values=parameter_values,parameter_doubledash=parameter_doubledash)			
		} else
		{
			parameter_variables_vector_strings <- NULL
		}
	}
	
	if(any("scalar" %in% parameter_variables_types))
	{
		parameter_variables_scalar <- parameter_variables$scalar[[1]]
		parameter_variables_scalar_defined <- defined_variables[defined_variables %in% parameter_variables_scalar]
		if(length(parameter_variables_scalar_defined)>0)
		{
			parameter_variables_scalar_strings <- sapply(parameter_variables_scalar_defined,
					function(X,parameter_values,parameter_doubledash)
					{
						if(X %in% parameter_noflags)
						{
							flag=NULL
						} else
						{
							if(X %in% parameter_doubledash)
							{
								flag=paste("--",X," ",sep="")	
							} else
							{
								flag=paste("-",X," ",sep="")
							}
						}
						parameter_variables_scalar_string <- paste(flag,
								qm(parameter_values[[which(names(parameter_values)==X)]]),
								sep="")
						return(parameter_variables_scalar_string)
					},parameter_values=parameter_values,parameter_doubledash=parameter_doubledash)			
		} else
		{
			parameter_variables_scalar_strings <- NULL
		}
	}
	
	if(any("character" %in% parameter_variables_types))
	{
		# Do we need to embed quotes in the command?
		parameter_variables_character <- parameter_variables$character[[1]]
		parameter_variables_character_defined <- defined_variables[defined_variables %in% parameter_variables_character]
		if(length(parameter_variables_character_defined)>0)
		{
			parameter_variables_character_strings <- sapply(parameter_variables_character_defined,
					function(X,parameter_values,parameter_noflags,parameter_doubledash)
					{
						if(X %in% parameter_noflags)
						{
							flag=NULL
						} else
						{
							if(X %in% parameter_doubledash)
							{
								flag=paste("--",X," ",sep="")	
							} else
							{
								flag=paste("-",X," ",sep="")
							}
						}
						parameter_variables_character_string <- paste(flag,
								qm(parameter_values[[which(names(parameter_values)==X)]]),
								sep="")
						return(parameter_variables_character_string)
					},parameter_values=parameter_values,parameter_noflags=parameter_noflags,parameter_doubledash=parameter_doubledash)			
		} else
		{
			parameter_variables_character_strings <- NULL
		}
	}
	
	if(any("repeatable" %in% parameter_variables_types))
	{
		parameter_variables_repeatable <- parameter_variables$repeatable[[1]]
		parameter_variables_repeatable_defined <- defined_variables[defined_variables %in% parameter_variables_repeatable]
		if(length(parameter_variables_repeatable_defined)>0)
		{
			parameter_variables_repeatable_strings <- sapply(parameter_variables_repeatable_defined,
					function(X,parameter_values,parameter_doubledash)
					{
						if(X %in% parameter_noflags)
						{
							flag=NULL
						} else
						{
							if(X %in% parameter_doubledash)
							{
								flag=paste("--",X," ",sep="")	
							} else
							{
								flag=paste("-",X," ",sep="")
							}
						}
						parameter_variables_repeatable_string <- paste(
								paste(flag,
										qm(parameter_values[[which(names(parameter_values)==X)]]),
										sep=""),
								collapse=" ")
						return(parameter_variables_repeatable_string)
					},parameter_values=parameter_values,parameter_doubledash=parameter_doubledash)			
		} else
		{
			parameter_variables_repeatable_strings <- NULL
		}
	}
	
	if(!is.null(parameter_noflags))
	{
#		parameter_variables_noflag <- parameter_variables$noflag[[1]]
#		parameter_variables_noflag_defined <- defined_variables[defined_variables %in% parameter_variables_noflag]
#		if(length(parameter_variables_noflag_defined)>0)
#		{
		parameter_variables_noflag_strings <- sapply(parameter_noflags,
				function(X,parameter_values)
				{
					parameter_variables_noflag_string <- paste(
							parameter_values[[which(names(parameter_values)==X)]],
							sep="")
					return(parameter_variables_noflag_string)
				},parameter_values=parameter_values)			
#		} else
#		{
#			parameter_variables_noflag_strings <- NULL
#		}
	}
	
	parameter_vector <- c(
			parameter_variables_logical_strings,
			parameter_variables_vector_strings,
			parameter_variables_scalar_strings,
			parameter_variables_character_strings,
			parameter_variables_repeatable_strings,
			parameter_variables_noflag_strings
	)
	
	# Reorder the parameters if neccessary
	if(!missing(parameter_order))
	{
		parameter_order_defined <- parameter_order[which(parameter_order %in% names(parameter_vector))]
		parameter_vector <- parameter_vector[parameter_order_defined]
	}
	
	# Collapse multiple parameter entries:
	parameter_vector <- sapply(parameter_vector,function(x) paste(x,collapse=" "))
	
	cmd <- paste(c(qm(executable),parameter_vector),collapse=" ")
	
	return(cmd)
	
}
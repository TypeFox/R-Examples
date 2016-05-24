#' gdal_setInstallation
#' 
#' Sets local GDAL installation options
#' 
#' @param search_path Character. Force a search in a specified directory.  This directory should contain the gdalinfo(.exe) executable.  If a valid GDAL install is found in this path, this will force gdalUtils to use this installation.  Remember to set rescan=TRUE if you have already set an install.
#' @param rescan Logical. Force a rescan if neccessary (e.g. if you updated your GDAL install).
#' @param ignore.full_scan Logical. If FALSE, perform a brute-force scan if other installs are not found.  Default is TRUE.
#' @param verbose Logical. Enable verbose execution? Default is FALSE.  

#' @return Sets an option "gdalUtils_gdalPath" with GDAL installation information.
#' @author Jonathan A. Greenberg (\email{gdalUtils@@estarcion.net}) and Matteo Mattiuzzi
#' 
#' @details This function searches the local system for valid installations of
#' GDAL, and returns a list, one item per valid GDAL install, containing 
#' the path to the installation, the version, the release date, available drivers,
#' and available python utilities.  The list will be sorted by release date, so
#' in general the first entry is the one that is used by the various GDAL utilities.
#' Note that this will automatically run every time a GDAL wrapper function is called,
#' so the user does not have to explicitly run it.
#' 
#' gdal_setInstallation is designed to invoke consecutively more 
#' rigorous searches in able to find a valid GDAL install.  Understanding
#' the search routine may help debug problems on your system.  The order
#' of the searches is as follows, noting that as soon as a valid install
#' is found (determined by running gdalinfo --version and getting the
#' correct output), gdal_setInstallation stops further searching:
#' \enumerate{
#' \item Checks a pre-determined location given by the search_path parameter.
#' \item Checks using Sys.which().  This is typically defined
#' 		in the system's PATH, so will override any other install.
#' \item Checks in common install locations (OS specific).  
#' \item (optional, if ignore.full_scan=FALSE) Finally, if it can't find a valid GDAL install anywhere else,
#' 		it will brute-force search the entire local system (which may
#' 		take a long time).  
#' }
#' @references \url{http://www.gdal.org/gdal_translate.html}
#' @examples \dontrun{ 
#' # Assumes you have GDAL installed on your local machine.
#' getOption("gdalUtils_gdalPath")
#' gdal_setInstallation()
#' getOption("gdalUtils_gdalPath")
#' # If there is more than one installation of GDAL, this is the 
#' # most recent installation:
#' getOption("gdalUtils_gdalPath")[[1]]
#' # The version number:
#' getOption("gdalUtils_gdalPath")[[1]]$version
#' }
#' @importFrom R.utils listDirectory
## @import utils
## @importFrom utils shortPathName
#' @export

# TODO: interface with gdal_chooseInstallation to remove some installs.
# TODO: set permanently in e.g. .Rprofile (Matteo)
# TODO: force a re-scan in gdal_setInstallation
# TODO: allow the user to specify a single installation path in gdal_setInstallation
# TODO: if nothing is found, give suggestions on where to download GDAL
# TODO: check if the user has permission to execute the commands

gdal_setInstallation <- function(search_path=NULL,rescan=FALSE,
		ignore.full_scan=TRUE,
		verbose=FALSE)
{
	
# Returns the available GDAL python utilities
	gdal_python_utilities <- function(path)
	{
		if(missing(path)) { path <- gdal_path() }
		sapply(path,
				listDirectory,
#				list.files,
				pattern="\\.py")
	}
	
# Returns the available GDAL drivers
	gdal_drivers <- function(path, verbose=FALSE)   
	{
		if(missing(path)) path <- gdal_path(checkValidity=TRUE)
		
		cmd <- file.path(path, "gdalinfo")
		cmd <- paste0('"',cmd,'"'," --formats")
		
		drivers_raw <- lapply(cmd,system,intern=TRUE)
		
		result <- vector(mode='list',length(path))
		names(result) <- path
		for(i in seq_along(drivers_raw))
		{
			drivers_raw[[i]] <- drivers_raw[[i]][-1]
			drivers=strsplit(drivers_raw[[i]],":")
			driver_names=gsub("^ ","",sapply(drivers,function(x) { x[2] })) # Need to remove spaces
			driver_codes_perm=strsplit(sapply(drivers,function(x) { x[1] }),"\\(")
			driver_codes=gsub(" ","",sapply(driver_codes_perm,function(x) { x[1] }),fixed=TRUE)
			driver_perm=gsub("\\)","",sapply(driver_codes_perm,function(x) { x[2] }))
			
			r <- w <- u <- v <- s <- rep(FALSE,length(driver_perm))
			r[grep(driver_perm, pattern="r")]   <- TRUE
			w[grep(driver_perm, pattern="w")]   <- TRUE
			u[grep(driver_perm, pattern="\\+")] <- TRUE
			v[grep(driver_perm, pattern="v")]   <- TRUE
			s[grep(driver_perm, pattern="s")]   <- TRUE	
			
			result[[i]] <- data.frame(format_code=driver_codes,read=r,write=w,update=u,virtualIO=v,subdatasets=s,format_name=driver_names)
		}
		return(result[[1]])
	}
	
# Returns the GDAL version
	gdal_version <- function(path, newerThan=NULL, verbose=FALSE)
	{
		if(missing(path)) { path <- gdal_path() }
		
		cmd <- normalizePath(
				listDirectory(
#				list.files(
						path, "^gdalinfo$|^gdalinfo\\.exe$",
						fullNames=TRUE)
#						full.names=TRUE)
		)
		cmd <- paste0('"',cmd,'"'," --version")
		
		result <- lapply(cmd,system,intern=TRUE)
		
		# It seams that system works also in Windows (MAC?), do you confirm?
		# Does shell work here?
		# if (.Platform$OS=="unix") 
		#{
		# gdal_version <- system(cmd,intern=TRUE) 
		# else 
		#{
		#gdal_version <- shell(cmd,intern=TRUE)
		#}
		
		res <- sapply(result,length)
		
		if(sum(res)!=length(result))
		{
			message("Probably broken install of gdal at '",paste0(path[which(res!=1)],collapse= "' and '"),"'")
		}
		result <- result[res==1]
		
		date <- version <- vector(mode = "list", length = length(result))
		
		for(i in seq_along(result))
		{
			ingd         <- strsplit(result[[i]],",")[[1]]
			version[[i]] <- gsub(ingd[1],pattern="GDAL ",replacement="")
			ind          <- grep(ingd,pattern="releas") # should this be: glob2rx("GDAL*")?
			date[[i]]    <- as.character(as.Date(gsub(ingd[ind],pattern=" released ",replacement=""),format="%Y/%m/%d"))
		}
		
		if(!is.null(newerThan))
		{
			test <- try(as.Date(newerThan),silent=TRUE)
			if(!inherits(test,"try-error"))
			{
				datein <- lapply(date,as.Date)
				res    <- sapply(datein,">=",as.Date(newerThan)) 
			} else
			{
				version   <- gsub(tolower(version),pattern="[a-z]",replacement="")
				res       <- sapply(version,strsplit,"\\.") 
				newerThan <- strsplit(newerThan,"\\.")[[1]]
				
				for(i in seq_along(res))
				{
					difs <- as.numeric(res[[i]]) - as.numeric(newerThan)
					difs <- sign(difs)
					
					if(sum(difs==-1)==0)
					{
						res[[i]] <- TRUE
					} else
					{
						if(difs[1]<0)
						{
							res[[i]] <- FALSE
						} else if(difs[1]>0)
						{
							res[[i]] <- TRUE
						} else if(difs[1]==0)
						{
							if(difs[2]<0)
							{
								res[[i]] <- FALSE
							} else if(difs[2]>0)
							{
								res[[i]] <- FALSE
							} else
							{  
								if(difs[3]>=0)
								{
									res[[i]] <- TRUE                  
								} else if (difs[3]<0)
								{
									res[[i]] <- FALSE
								}
							}
						}
					}
				}
			}
			names(res) <- path
			return(res)
		}
		result <- as.data.frame(cbind(path=path[res==1],version=version,date=date), stringsAsFactors=FALSE)
		return(result)
	}
	
	correctPath <- function(x)
	{
		if(!is.null(x))
		{
			if (.Platform$OS.type=="windows")
			{
				x <- utils::shortPathName(x)
			} else
			{
				x <- path.expand(x)
			}
			x      <- gsub(x,pattern="\\\\",replacement="/") # switch "\\" to "/
			ind    <- substr(x,nchar(x),nchar(x))!="/"       # some x without "/" at the end?
			x[ind] <- paste0(x[ind],"/")                     # add "/" at the end
		}
		return(x)
	}
	
# Checks if GDAL is functional
	gdal_check_validity <- function(path,utility="^gdalinfo$|^gdalinfo\\.exe$")
	{
#		full_utility_list <- c(
#				"^gdalinfo$|^gdalinfo\\.exe$",
#				"^gdal_rasterize$|^gdal_rasterize\\.exe$",
#				"^gdal_translate$|^gdal_translate\\.exe$",
#				"^gdaladdo$|^gdaladdo\\.exe$",
#				"^gdalbuildvrt$|^gdalbuildvrt\\.exe$",
#				"^gdaldem$|^gdaldem\\.exe$",
#				"^gdalsrsinfo$|^gdaldem\\.exe$",
#				
#				)
	 	
		
		checkValidity <- sapply(path,
				function(x)
				{
					cmd <- normalizePath(
							listDirectory(
#							list.files(
#									path=x,pattern="^gdalinfo$|^gdalinfo\\.exe$",
									path=x,pattern=utility,
									fullNames=TRUE)
#									full.names=TRUE)
					)
					
					if(length(cmd)==0)
					{
						return(FALSE)
					} else
					{
						cmd <- paste0('"',cmd,'"'," --version")
						validity = length(try(gdal <- system(cmd,intern=TRUE),silent=TRUE))
						
						return(as.logical(validity))
					}
				}
		)
	}
	
# Determines the path to GDAL installations
	gdal_path <- function(
			search_path=NULL,
			ignore.options=FALSE,
			ignore.which=FALSE,
			ignore.common=FALSE,
			ignore.full_scan=FALSE,
#			force_full_scan = FALSE, 
			checkValidity, 
			search_path_recursive=FALSE,
			verbose = FALSE)
	{
#		browser()
		owarn <- getOption("warn")
		options(warn=-2)
		on.exit(options(warn=owarn))
		
		if(missing(checkValidity))
		{
			if(is.null(getOption("gdalUtils_gdalPath"))) checkValidity=TRUE else checkValidity=FALSE
		}
		
		path <- NULL
		# Rescan will override everything.
		if(ignore.full_scan)
		{
#			if(verbose) message("No forced full scan...")
			# Check options first.
			
			if(!ignore.options)
			{
				if(verbose) message("Checking the gdalUtils_gdalPath option...")
				option_paths <- unlist(
						sapply(getOption("gdalUtils_gdalPath"),function(x) return(x$path)))
				if(!is.null(option_paths) && checkValidity)
				{
					option_paths_check <- gdal_check_validity(option_paths)
					option_paths <- option_paths[option_paths_check]
				}
				path <- c(path,option_paths)
			}
			
			# Next, try scanning the search path
#			if(!missing(search_path) && length(path)==0)
			if(!is.null(search_path) && length(path)==0)
			
			{
				if(verbose) message("Checking the search path...")
				if (.Platform$OS=="unix")
				{
					search_paths <- 
							listDirectory(
									path=search_path,pattern="^gdalinfo$|^gdalinfo\\.exe$",
									recursive=search_path_recursive,
									fullNames=TRUE)
				} else
				{
					search_paths <- 
							list.files(
									path=search_path,pattern="^gdalinfo$|^gdalinfo\\.exe$",
									recursive=search_path_recursive,
									full.names=TRUE)
				}
				if(length(search_paths)==0) search_paths <- NULL else 
					search_paths <- normalizePath(dirname(search_paths))
				if(!is.null(search_paths) && checkValidity)
				{
					search_paths_check <- gdal_check_validity(search_paths)
					search_paths <- search_paths[search_paths_check]
				}
				path <- c(path,search_paths)
			}
			
			# Next try Sys.which unless ignored:
			if(!ignore.which && length(path)==0)
			{
				if(verbose) message("Checking Sys.which...")
				Sys.which_path <- dirname(Sys.which("gdalinfo"))
				if(Sys.which_path=="") Sys.which_path <- NULL
				if(!is.null(Sys.which_path) && checkValidity)
				{
					Sys.which_path_check <- gdal_check_validity(Sys.which_path)
					Sys.which_path <- Sys.which_path[Sys.which_path_check]
				}
				path <- c(path,Sys.which_path)
			}
			
			# If nothing is still found, look in common locations
			if(!ignore.common && length(path)==0)
			{
				if(verbose) message("Checking common locations...")
				if (.Platform$OS=="unix")
				{
					common_locations <- c(
							### UNIX systems
							"/usr/bin",
							"/usr/local/bin",
							### Mac
							# Kyngchaos frameworks:
							"/Library/Frameworks/GDAL.framework",
							# MacPorts:
							"/opt/local/bin",
							# Homebrew:
							"/usr/local/Cellar/gdal/"
					)
				}
				
				if (.Platform$OS=="windows")
				{
					common_locations <- c(
							"C:\\Program Files",
							"C:\\Program Files (x86)",
							"C:\\OSGeo4W",
							"C:\\OSGeo4W64"
					)
				}
				
				
				if(length(common_locations != 0))
				{
					common_paths <- unlist(sapply(common_locations,
									function(x)
									{
										if (.Platform$OS=="unix")
										{
											search_common_paths <- 
													#normalizePath(dirname(
													listDirectory(
#														list.files(
															path=x,pattern="^gdalinfo$|^gdalinfo\\.exe$",recursive=TRUE,
#																full.names=TRUE)
															fullNames=TRUE)
											#	))
										} else
										{
											search_common_paths <- 
													#normalizePath(dirname(
													list.files(
#														list.files(
															path=x,pattern="^gdalinfo$|^gdalinfo\\.exe$",recursive=TRUE,
#																full.names=TRUE)
															full.names=TRUE)
										}
										if(length(search_common_paths)==0)
											return(search_common_paths)
										else
											return(normalizePath(dirname(search_common_paths)))
									}))
					if(length(common_paths)==0) common_paths <- NULL
					if(!is.null(common_paths) && checkValidity)
					{
						common_paths_check <- gdal_check_validity(common_paths)
						common_paths <- common_paths[common_paths_check]
					}
					path <- c(path,common_paths)
				}
			}
			if(ignore.full_scan && length(path)==0)
			{
				ignore.full_scan=FALSE
			}
		}
		
		if(!ignore.full_scan)
		{
			if(verbose) message("Scanning your root-dir for available GDAL installations,... This could take some time...")
			if (.Platform$OS=="unix")
			{
				root_dir <- "/"	
			}
			
			if (.Platform$OS=="windows")
			{
				root_dir <- "C:\\"
			}
			if (.Platform$OS=="unix")
			{
				search_full_path <- 
						# normalizePath(dirname(
						listDirectory(
#							list.files(
								path=root_dir,pattern="^gdalinfo$|^gdalinfo\\.exe$",
								recursive=TRUE,
								fullNames=TRUE)
#									full.names=TRUE)
				# ))
			} else
			{
				search_full_path <- 
						# normalizePath(dirname(
						list.files(
#							list.files(
								path=root_dir,pattern="^gdalinfo$|^gdalinfo\\.exe$",
								recursive=TRUE,
								full.names=TRUE)
			}
			if(length(search_full_path)==0) search_full_path <- NULL
			else search_full_path <- normalizePath(dirname(search_full_path))
			if(!is.null(search_full_path) && checkValidity)
			{
				search_full_path_check <- gdal_check_validity(search_full_path)
				search_full_path <- search_full_path[search_full_path_check]
			}
			path <- c(path,search_full_path)
		}
		
		if(length(path)==0)
		{
			#add QGIS?
			return(NULL)
		} else
		{	
			return(correctPath(unique(path)))
		}
	}
	
# Returns the full GDAL installation status
	gdal_installation <- function(
			return_versions=TRUE,
			return_drivers=TRUE,
			return_python_utilities=TRUE,
			sort_most_current=TRUE,
			rescan=FALSE,
			search_path=NULL,
			ignore.full_scan=FALSE,
			verbose=FALSE
	)
	{
		if(verbose) message("Scanning for GDAL installations...")
		path <- gdal_path(ignore.options=rescan,search_path=search_path,ignore.full_scan=ignore.full_scan,
				verbose=verbose)
		#	browser()
		if(is.null(path)) return(NULL)
		
		gdal_installation_results <- lapply(path,
				function(x,return_drivers,return_python_utilities,return_versions)
				{
					#		browser()
					result <- list(path=x)
					
					if(return_versions)
					{
						version <- gdal_version(x)
						result$version <- unlist(version$version)
						result$date <- unlist(version$date)
					}
					
					if(return_drivers)
					{
						result$drivers <- gdal_drivers(x)    
					}
					
					if(return_python_utilities)
					{
						result$python_utilities <- gdal_python_utilities(x)    
					}
					return(result)
				},return_drivers=return_drivers,
				return_python_utilities=return_python_utilities,return_versions=return_versions)
		if(sort_most_current)
		{
			versions <- unlist(sapply(gdal_installation_results,function(x) return(x$date)))
			gdal_installation_results <- gdal_installation_results[
					order(as.Date(unlist(versions)),decreasing=TRUE)]
		}
		return(gdal_installation_results)    
	}
	
# Sets the installation for this session.
	if(is.null(getOption("gdalUtils_gdalPath")))
	{
		rescan=TRUE	
	}
	gdal_installation_out <- gdal_installation(search_path=search_path,rescan=rescan,ignore.full_scan=ignore.full_scan,
			verbose=verbose)
	options(gdalUtils_gdalPath=gdal_installation_out)
	if(is.null(getOption("gdalUtils_gdalPath")))
	{
		warning("No GDAL installation found. Please install 'gdal' before continuing:\n\t- www.gdal.org (no HDF4 support!)\n\t- www.trac.osgeo.org/osgeo4w/ (with HDF4 support RECOMMENDED)\n\t- www.fwtools.maptools.org (with HDF4 support)\n") # why not stop?
		if(ignore.full_scan) warning("If you think GDAL is installed, please run:\ngdal_setInstallation(ignore.full_scan=FALSE)")
	}
	else
	{
		if(verbose) message("GDAL version ",unlist(getOption("gdalUtils_gdalPath")[[1]]$version))
	}
}
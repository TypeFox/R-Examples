# Utils for R package projects
# 
# Author: Renaud Gaujoux
# Created: May 1, 2013
###############################################################################

find_devpackage <- function (x) 
{
	
	is_package_path <- function(x, check=FALSE) {
		if (is.null(x)) return(FALSE)
		x <- normalizePath(x, mustWork = FALSE)
		x <- gsub("\\\\$", "", x)
		desc_path <- file.path(x, "DESCRIPTION")
		if( !check ){
			file.exists(x) && file.exists(desc_path)
		}else{
			if (!file.exists(x)) {
				stop("Can't find directory ", x, call. = FALSE)
			}
			if (!file.info(x)$isdir) {
				stop(x, " is not a directory", call. = FALSE)
			}
			desc_path <- file.path(x, "DESCRIPTION")
			if (!file.exists(desc_path)) {
				stop("No DESCRIPTION file found in ", x, call. = FALSE)
			}
			TRUE
		}
	}
	
	if (is_package_path(x)) {
		return(x)
	}
	
	config_path <- "~/.Rpackages"
	if (!file.exists(config_path)) {
		return(NULL)
	}
	config_path <- path.expand(config_path)
	lookup <- source(config_path)$value
	if (is_package_path(lookup[[x]])) {
		return(lookup[[x]])
	}
	if (!is.null(lookup$default)) {
		default_loc <- lookup$default(x)
		if (is_package_path(default_loc, check=TRUE)) {
			return(default_loc)
		}
	}
	NULL
}


is_Mac <- function(check.gui=FALSE){
	is.mac <- (length(grep("darwin", R.version$platform)) > 0)
	# return TRUE is running on Mac (adn optionally through GUI)
	is.mac && (!check.gui || .Platform$GUI == 'AQUA')
}

R_OS <- function(){
	if( is_Mac() ) 'MacOS'
	else .Platform$OS.type
}

packageMakefile <- function(package=NULL, template=NULL, temp = FALSE, print = TRUE){
	
	capture.output(suppressMessages({
		library(pkgmaker)
#		library(methods)
		library(devtools)					
	}))
#	defMakeVar <- pkgmaker::defMakeVar
#	subMakeVar <- pkgmaker::subMakeVar
	
	project_path <- getwd()
	project_name <- basename(project_path)
	subproject_path_part <- ''
	if( is.null(package) || isString(package) ){
		if( isString(package) && !nzchar(package) ) package <- NULL
		lookup_dir <- c('pkg', '.')
		if( !is.null(package) ){
			lookup_dir <- c(lookup_dir, file.path('pkg', package))
			subproject_path_part <- file.path(package, '')
		}
		pdir <- file.path(lookup_dir, 'DESCRIPTION')
		if( !length(sd <- which(is.file(pdir))) ){
			stop("Could not detect package source directory")
		}
		package <- pdir[sd[1L]]
	}
	package <- normalizePath(package)
	p <- pkg <- as.package(dirname(package));
	pdir <- package_dir <- p[['path']];
	
	## create makefile from template
	# load template makefile
	if( is.null(template) ){
		template <- packagePath('package.mk', package='pkgmaker')
	}
	l <- paste(readLines(template), collapse="\n")
	
	# user
	cuser <- Sys.info()["user"]
	l <- defMakeVar('AUTHOR_USER', cuser, l)
	l <- defMakeVar('R_PACKAGE', pkg$package, l)
	# R_PACKAGE_PATH
	l <- defMakeVar('R_PACKAGE_PATH', package_dir, l)
	# R_PACKAGE_PROJECT
	l <- defMakeVar('R_PACKAGE_PROJECT', project_name, l)
	# R_PACKAGE_PROJECT_PATH
	l <- defMakeVar('R_PACKAGE_PROJECT_PATH', project_path, l)
	l <- defMakeVar('R_PACKAGE_SUBPROJECT_PATH_PART', subproject_path_part, l)
	# R_BIN
	l <- subMakeVar('R_BIN', R.home('bin'), l)
	# R_PACKAGE_TAR_GZ
	pkg_targz <- paste0(p[['package']], '_', p[['version']], '.tar.gz')
	l <- defMakeVar('R_PACKAGE_TAR_GZ', pkg_targz, l)
	# R_PACKAGE_TYPE	
	l <- defMakeVar('R_PACKAGE_OS', R_OS(), l)
	#

    # auto-conf variables
    init_var <- list(version = pkg$version)
    if( is.dir(file.path(package_dir, 'vignettes')) ) 
        init_var <- c(init_var, has_vignettes=TRUE)
    # dump variables
    if( length(init_var) ){
        init_var <- setNames(init_var, paste0('R_PACKAGE_', toupper(names(init_var))))
        init_var_str <- str_out(init_var, Inf, use.names = TRUE, sep = "\n")
        l <- subMakeVar('INIT_CHECKS', init_var_str, l)
    }
	
	# R_CMD_CHECK
	rlibs <- ''
    if( is.dir(devlib <- file.path(dirname(pdir), 'lib')) ){
		rlibs <- paste0("R_LIBS=", devlib, ' ')
	}
    l <- subMakeVar('R_LIBS', rlibs, l)
	#

	# create makefile
	mk <- if( temp ) tempfile('package_', tmpdir='.', fileext='.mk') else 'package.mk'
	cat(l, file=mk)
	if ( print ){
		cat(mk)
	}
	invisible(l)
}
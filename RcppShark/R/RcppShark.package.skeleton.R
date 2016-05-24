## RcppShark -- An R interface to the Shark machine learning library
##
## Copyright (C)  2014  Qiang Kou
## This file was part of RcppMLPACK.
##
## Copyright (C) 2015  Aydin Demircioglu
##
## This file is part of the RcppShark library for GNU R.
## It is made available under the terms of the GNU General Public
## License, version 3, or at your option, any later version,
## incorporated herein by reference.
##
## This program is distributed in the hope that it will be
## useful, but WITHOUT ANY WARRANTY; without even the implied
## warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
## PURPOSE.  See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public
## License along with this program; if not, write to the Free
## Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
## MA 02111-1307, USA



#' create a skeleton for a package using RcppShark
#'
#' @param name    character string: the package name and directory name for your package.
#' @param list	character vector naming the R objects to put in the package. 
#'				Usually, at most one of "list", "environment", or "code_files" will be supplied.  See "Details" of package.skeleton.
#' @param environment 	an environment where objects are looked for.  See "Details" of package.skeleton.
#' @param path		 path to put the package directory in.
#' @param force	 If "FALSE" will not overwrite an existing directory.
#' @param code_files	 a character vector with the paths to R code files to build the package around.  See "Details" of package.skeleton.
#' @param example_code	add example code to package?
#'
#' @export
RcppShark.package.skeleton <- function(name="anRpackage", list=character(),
										environment=.GlobalEnv,
										path=".", force=FALSE, 
										code_files=character(), 
										example_code=TRUE) 
{

	env <- parent.frame(1)

	if (! length(list)) {
		fake <- TRUE
		assign( "Rcpp.fake.fun", function(){}, envir = env )
	} else {
		fake <- FALSE
	}

	## first let the traditional version do its business
	call <- match.call()
	call[[1]] <- as.name("package.skeleton")
	if ("example_code" %in% names(call)){
		## remove the example_code argument
		call[["example_code"]] <- NULL
	}
	if (fake) {
		call[["list"]] <- "Rcpp.fake.fun"
	}

	tryCatch(eval(call, envir=env),
			error = function(e) {
				stop("error while calling `package.skeleton`")
			})

	message("\nAdding RcppShark settings")

	## now pick things up 
	root <- file.path(path, name)

	## Add Rcpp to the DESCRIPTION
	DESCRIPTION <- file.path(root, "DESCRIPTION")
	if (file.exists(DESCRIPTION)) {
		x <- cbind(read.dcf( DESCRIPTION), 
				"Imports" = sprintf("Rcpp (>= %s)", 
									utils::packageDescription("Rcpp")[["Version"]]), 
				"LinkingTo" = "Rcpp, RcppShark, BH")
		write.dcf(x, file=DESCRIPTION)
		message(" >> added Imports: Rcpp")
		message(" >> added LinkingTo: Rcpp, RcppShark, BH")
	}

	## add a useDynLib to NAMESPACE, 
	NAMESPACE <- file.path( root, "NAMESPACE")
	lines <- readLines( NAMESPACE )
	if (! grepl("useDynLib", lines)) {
		lines <- c(sprintf("useDynLib(%s)", name),
				"importFrom(Rcpp, evalCpp)",        ## ensures Rcpp instantiation
				lines)
		writeLines(lines, con = NAMESPACE)
		message( " >> added useDynLib and importFrom directives to NAMESPACE")
	}

	## lay things out in the src directory
	src <- file.path(root, "src")
	if (!file.exists(src)) {
		dir.create(src)
	}
	skeleton <- system.file("skeleton", package="RcppShark")
	Makevars <- file.path(src, "Makevars")
	if (!file.exists(Makevars)) {
		file.copy(file.path(skeleton, "Makevars"), Makevars)
		message(" >> added Makevars from ", file.path(skeleton, "Makevars"), " to ", Makevars)
	}

	Makevars.win <- file.path(src, "Makevars.win")
	if (! file.exists( Makevars.win)) {
		file.copy(file.path(skeleton, "Makevars.win"), Makevars.win)
		message(" >> added Makevars.win file with RcppShark settings")
	}

	if (example_code) {
		file.copy(file.path(skeleton, "KMeansTutorial.cpp"), src)
		file.copy(file.path(skeleton, "utils.h"), src)
		file.copy(file.path(skeleton, "utils.cpp"), src)
		message(" >> added example src file using RcppShark classes")
		Rcpp::compileAttributes(root)
		message(" >> invoked Rcpp::compileAttributes to create wrappers")
	}

	if (fake) {
		rm("Rcpp.fake.fun", envir=env)
		unlink(file.path(root, "R"  , "Rcpp.fake.fun.R"))
		unlink(file.path(root, "man", "Rcpp.fake.fun.Rd"))
	}

	invisible(NULL)
}

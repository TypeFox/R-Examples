#' Create Documentation file from r-source
#' 
#' Create Documentation file (especially section arguments) from r-source
#' 
#' @details This assumes the following structure of source code:\cr 
#'    \code{MyFun <- function(}\cr 
#'    \code{arg1, # Explanation of this item}\cr 
#'    \code{arg2=TRUE, # Ditto, with default}\cr 
#'    \code{arg3)}\cr 
#'    \code{'{'}\cr
#'    \code{computations'}'}\cr
#'    The opening bracket line may ONLY contain the curly brace
#' 
#' @note This might be deprecated in favor of Roxygen
#' @return None. Cats documentation for \code{fun} in \code{path}/man. Only usage, arguments and author section are filled, the rest is empty (but the frame is there).
#' @section Warning: This is highly specific to my way of working, don't rely blindly on it.\cr 
#'     If a file already exists, it is not overwritten, instead a new file \code{fun_2.Rd} or \code{fun_3.Rd} (up to 99), is created.\cr 
#'     Empty (or space-only) lines are silently ignored.\cr A line with two arguments will throw a warning, 
#'     as they can't be listed in the argument section. They should be written normally into the usage section.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, June + Dec 2014
#' @seealso \code{\link{package.skeleton}}, \code{\link{prompt}}, \code{\link{scan}}, 
#'          \code{\link{cat}}, Roxygen and Roxygen2 \url{http://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html}
#' @keywords documentation
#' @export
#' @examples
#' 
#' #createDoc("textField")
#' 
#' @param fun Character string or unquoted name. Function (== filename) with structure described in 'Details' in source code.
#' @param path Path to package in development containing folders 'R' and 'man'. DEFAULT: "S:/Dropbox/Public/berryFunctions"
#' 
createDoc <- function(
  fun,
  path="S:/Dropbox/Public/berryFunctions"
  )
{
warning("createDoc will be deprecated in late 2016. Use createFun() instead.")
fun <- deparse(substitute(fun))
fun <- gsub("\"", "", fun, fixed=TRUE)
if(length(fun) >1) stop("'fun' must be a single function name.")
if(length(path)>1) stop("'path' must be a single character string.")
# work PC path change:
if(!file.exists(path)) substr(path, 1,1) <- "D"
# laptop linux path change:
if(!file.exists(path)) path <- gsub("D:", "~", path)
# new work PC path change:
if(!file.exists(path)) path <- gsub("~", "C:/Users/boessenkool", path)
# path control
if(!file.exists(path)) stop("path does not exist. ", path)
owd <- setwd(path)
# set wd back to old working directory:
on.exit(setwd(owd))
#
                            rfilename <- paste0("R/",fun,".r")
if(!file.exists(rfilename)) rfilename <- paste0("R/",fun,".R") # capital R
if(!file.exists(rfilename)) stop("File ", path, "/", rfilename, " does not exist")
# case control Windows:
if(!any(dir("R")==paste0(fun,".r")|dir("R")==paste0(fun,".R"))) stop("'", fun,
   "' does not match capitalization of files in ", path, "/R")
# read file
rfile <- readLines(rfilename)
rfile <- rfile[removeSpace(rfile)!=""]
anf <- grep("<- function", rfile)[1]      # begin line of argument section
end <- which(removeSpace(rfile)=="{")[1]  # end line
if (end < anf) stop("Argument section (begin=",anf,", end=",end,") was not correctly identified!")
#
rdfile <- paste0("man/",fun,".Rd")
# File name check:
Newfilecreated <- FALSE   ; file_nr <- 1
while(file.exists(rdfile))
    {
    rdfile <- paste0("man/",fun,"_", file_nr,".Rd")
    file_nr <- file_nr + 1
    Newfilecreated <- TRUE
    }
if(Newfilecreated) warning("File already existed. Created the file ", path, "/", rdfile)
# HEADER
cat(paste0("\\name{", fun, "}\n\\alias{", fun, "}\n"), file=rdfile)
cat(paste0("\\title{}\n\\description{}\n\\usage{", fun , "("), file=rdfile, append=TRUE)
# USAGE
# ignore comments:
usage <- sapply(strsplit(rfile[(anf+1):(end-1) ], "#"), "[", 1)
# Remove leading and trailing white spaces:
usage <- removeSpace(usage)
# ignore empty lines:
usage <- usage[usage!=""]
# double the backslashes:
usage <- gsub("\\n", "\\\\n", usage, fixed=TRUE)
usage <- gsub("\\t", "\\\\t", usage, fixed=TRUE)
# write section:
cat(paste(usage, collapse=" "), file=rdfile, append=TRUE)
# ARGUMENTS
cat(paste0("}\n\\arguments{\n"), file=rdfile, append=TRUE)
for(i in (anf+1):(end-1) )
  {
  # Split argument and explanation:
  arg_expl <- strsplit(rfile[i], "#")[[1]]
  # Remove leading and trailing white spaces:
  arg_expl <- removeSpace(arg_expl)
  # double the backslashes:
  arg_expl <- gsub("\\n", "\\\\n", arg_expl, fixed=TRUE)
  arg_expl <- gsub("\\t", "\\\\t", arg_expl, fixed=TRUE)
  #arg_expl <- gsub("\\", "\\\\", arg_expl, fixed=TRUE)
  # remove trailing comma in Argument:
  if(grepl("[,]$", arg_expl[1]))
     arg_expl[1] <- substring(arg_expl[1], 1, nchar(arg_expl[1])-1)
  # split arg name and default value:
  if(grepl("=", arg_expl[1]))
     arg_expl[c(1,3)] <- strsplit(
       sub("=", "firstEqualSplit", arg_expl[1]), "firstEqualSplit")[[1]]
  # Dots:
  if(grepl("...", arg_expl[1], fixed=TRUE)) arg_expl[1] <- "\\dots"
  # Write to Rd-File:
  if( arg_expl[1] != ")" )
    {
    # write arg name and explanation
    cat(paste0("  \\item{",arg_expl[1],"}{",arg_expl[2]), file=rdfile, append=TRUE)
    # write default value:
    if(length(arg_expl)==3 & !grepl("default:", arg_expl[2], ignore.case=TRUE) )
      cat(paste0(if(!grepl("[.?!]$", arg_expl[2]))".", " DEFAULT: ", arg_expl[3]),
          file=rdfile, append=TRUE)
    #
    cat("}\n", file=rdfile, append=TRUE)
    }
  } # End Loop
# FOOTER
cat(paste0("}
\\details{}
\\value{}
\\section{Warning}{}
\\note{}
\\author{Berry Boessenkool, \\email{berry-b@gmx.de}, ",format(Sys.Date(), "%b %Y"),"}
\\references{}
\\seealso{\\code{\\link{help}} }
\\examples{

}
\\keyword{}
\\keyword{}
"), file=rdfile, append=TRUE)
# warning if there are less argument lines than arguments themselves:
if(removeSpace(rfile[end-1]) ==")") end <- end-1
source(rfilename, local=TRUE)
n_missing <- length(formals(fun))  -  (end-anf-1)
if(n_missing != 0) warning(n_missing, " items are missing in the arguments section.
There possibly are several arguments on one line in ", path, "/", rfilename,
"\nlikely, the 'DEFAULT: ' at the line end is corrupted as well.
beginline=",anf,", endline=",end,", nlines=",end-anf-1,", length formals=",length(formals(fun)))
# file location:
if(!Newfilecreated) message("Created the file ", path, "/", rdfile)
} # End of Function

#' @title Find matlab path
#'
#' @description This tries to find matlab's path using a system which
#' command, and then, if not found, looks at \code{getOption("matlab.path")}.  If not path is found, it fails.
#' @export
#' @return Character of command for matlab
get_matlab = function(){
  # find.matlab <- system("which matlab", ignore.stdout=TRUE)
  mat = paste0("matlab", 
               ifelse(.Platform$OS.type %in% "windows", ".exe", "")
  )
  find.matlab = as.numeric(Sys.which(mat) == "")
  matcmd <- paste0(mat, ' -nodesktop -nosplash -nodisplay -r ')
  if (find.matlab != 0){
    mpath = getOption("matlab.path")
    stopifnot(!is.null(mpath))
    stopifnot(file.exists(mpath))
    mpath = shQuote(mpath)
    matcmd <- file.path(mpath, matcmd)
  }
  return(matcmd)
}

#' @title Logical check if MALTAB is accessible
#' @description Uses \code{get_matlab} to check if MATLAB is accessible 
#' or the option
#' \code{matlab.path} is set and returns logical
#' @return Logical TRUE is MALTAB is accessible, FALSE if not
#' @export
#' @examples
#' have_matlab()
have_matlab = function(){
  x = suppressWarnings(try(get_matlab(), silent=TRUE))
  return(!inherits(x, "try-error"))
}


#' @title Run matlab script
#'
#' @description This function runs a matlab script, and 
#' returns exit statuses
#' @param fname Filename of matlab script (.m file)
#' @param ... Options passed to \code{\link{system}}
#' @export
#' @return Exit status of matlab code
run_matlab_script = function(fname, ...){
  stopifnot(file.exists(fname))
  matcmd = get_matlab()
  cmd = paste0(' "', "try, run('", fname, "'); ",
               "catch err, disp(err.message); ", 
               "exit(1); end; exit(0);", '"')  
  cmd = paste0(matcmd, cmd)
  x <- system(cmd, ...)
  return(x)
}

#' @title Runs matlab code
#'
#' @description This function takes in matlab code, where
#' the last line must end with a ;, and returns the exit
#' status
#' @param code Character vector of code. 
#' @param endlines Logical of whether the semicolon (;) should be
#' pasted to each element of the vector.
#' @param verbose Print out filename to run
#' @param add_clear_all Add \code{clear all;} to the beginning of code
#' @param ... Options passed to \code{\link{run_matlab_script}}
#' @export
#' @return Exit status of matlab code 
#' @examples 
#' if (have_matlab()){
#'    run_matlab_code("disp(version)")
#'    run_matlab_code(c("disp('The version of the matlab is:')", "disp(version)"))
#'    run_matlab_code(c("x = 5", "disp(['The value of x is ', num2str(x)])"))
#' }
run_matlab_code = function(code, endlines = TRUE, verbose = TRUE,
                           add_clear_all = FALSE,
                           ...){
  matcmd = get_matlab()
  code = c(ifelse(add_clear_all, "clear all;", ""), 
           paste0("cd('", getwd(), "');"), code)
  sep = ifelse(endlines, ";", " ")
  code = paste0(code, sep = sep, collapse= "\n")
  code = gsub(";;", ";", code)
#   cmd <- paste(' "try \n')
#   cmd <- paste(cmd, code)
#   cmd <- paste(cmd, "\n catch err \n disp(err.message); \n exit(1); \n")
#   cmd <- paste0(cmd, 'end; \n exit(0);"')
#   cmd = gsub("\n", ";", cmd)
#   cmd = paste0(matcmd, cmd)
  cmd = code
  fname = tempfile(fileext = ".m")
  cat(cmd, file = fname)
  if (verbose){
    cat(paste0(fname, "\n"))
  }
  x = run_matlab_script(fname, ...)
  return(x)
}



#' @title Convert R vector to matlab cell mat
#'
#' @description This function takes in an R vector then turns it into 
#' a cell list
#' @param x Character vector of values
#' @param matname Object in matlab to be assigned
#' @export
#' @return Character scalar of matlab code
rvec_to_matlabclist = function(x, matname = NULL){
  x = paste0("{'", x, "'};")
  x = paste(x, collapse= " ")
  x = paste0('[', x, '];')
  if (!is.null(matname)) x = paste0(matname, " = ", x)
  x
}



#' @title Convert R vector to matlab cell mat
#'
#' @description This function takes in an R numeric and returns a
#' status
#' @param x Numeric vector of values
#' @param row Create row vector instead of column vector
#' @param matname Object in matlab to be assigned
#' @export
#' @return Character scalar of matlab code
#' @import stringr
rvec_to_matlab = function(x, row = FALSE,
                          matname = NULL){
  x = paste0(x, ifelse(row, ",", ";"))
  x = paste(x, collapse= " ")
  x = str_trim(x)
  x = gsub(paste0(ifelse(row, ",", ";"), "$"), "", x)
  x = paste0("[", x, "];")
  if (!is.null(matname)) x = paste0(matname, " = ", x)
  x
}

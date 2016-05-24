utils::globalVariables(c(".home"))
utils::suppressForeignCheck(".home")

#' Enables detection of the folder a script resides in with certain accuracy
#'
#' R does not have a default way to return and use, the location of a specific script.
#'
#' @export
#'
#' @importFrom utils tail
#'
#' @source https://github.com/molgenis/molgenis-pipelines/wiki/How-to-source-another_file.R-from-within-your-R-script
get_script_path <- function(){
  this.file = NULL
  # This file may be 'sourced'
  for (i in -(1:sys.nframe())) {
    if (identical(sys.function(i), base::source)) this.file = (normalizePath(sys.frame(i)$ofile))
  }

  if (!is.null(this.file)) return(dirname(this.file))

  # But it may also be called from the command line
  cmd.args = commandArgs(trailingOnly = FALSE)
  cmd.args.trailing = commandArgs(trailingOnly = TRUE)
  cmd.args = cmd.args[seq.int(from=1, length.out=length(cmd.args) - length(cmd.args.trailing))]
  res = gsub("^(?:--file=(.*)|.*)$", "\\1", cmd.args)

  # If multiple --file arguments are given, R uses the last one
  res = tail(res[res != ""], 1)


  if (length(res) > 0){
    res = tools::file_path_as_absolute(dirname(res))
    #assign(.home, res, envir = .GlobalEnv)
    #.home <<- res

    return(res)
  }


  # Both are not the case. Maybe we are in an R GUI?
  return(NULL)
}

#' @rdname  get_script_path
#' @export
detect_home=get_script_path

# # this dynamically detects path to the script
# ultraseq_home = detect_home()
# message("ultraseq root dir: ", ultraseq_home)

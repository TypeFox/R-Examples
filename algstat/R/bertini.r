#' Evaluate Bertini Code
#'
#' Write a Bertini file, evaluate it through a back-end connection to Bertini, and bring the output back into R.
#' 
#' @param code Bertini code as either a character string or function; see examples
#' @param dir directory to place the files in, without an ending /
#' @param opts options for bertini
#' @param quiet show bertini output
#' @return an object of class bertini
#' @export bertini
#' @examples
#' \dontrun{
#' 
#' # where does the circle intersect the line y = x?
#' code <- "
#' INPUT
#' 
#' variable_group x, y;
#' function f, g;
#' 
#' f = x^2 + y^2 - 1;
#' g = y - x;
#' 
#' END;
#' "
#' bertini(code)
#' 
#' c(sqrt(2)/2, sqrt(2)/2)
#'
#' 
#' 
#' 
#' # where do the surfaces
#' #  x^2 - y^2 - z^2 - 1/2
#' #  x^2 + y^2 + z^2 - 9
#' #  x^2/4 + y^2/4 - z^2
#' # intersect?
#' # 
#' code <- "
#' INPUT
#' 
#' variable_group x, y, z;
#' function f, g, h;
#' 
#' f = x^2 - y^2 - z^2 - 1/2;
#' g = x^2 + y^2 + z^2 - 9;
#' h = x^2/4 + y^2/4 - z^2;
#' 
#' END;
#' "
#' bertini(code)
#'
#' # algebraic solution :
#' c(sqrt(19)/2, 7/(2*sqrt(5)), 3/sqrt(5)) # +/- each ordinate
#' 
#' 
#'
#'
#' # example from bertini manual
#' code <- "
#' INPUT
#' 
#' variable_group x, y;
#' function f, g;
#' 
#' f = x^2 - 1;
#' g = x + y - 1;
#' 
#' END;
#' "
#' out <- bertini(code)
#' summary(out)
#' 
#' 
#' 
#' 
#' 
#' # non zero-dimensional example 
#' code <- "
#' CONFIG
#'   TRACKTYPE: 1;
#' END;
#' 
#' INPUT
#'   variable_group x, y, z;
#'   function f1, f2;
#'   f1 = x^2-y;
#'   f2 = x^3-z;
#' END;
#' "
#' out <- bertini(code)
#' # bertini(code, quiet = FALSE) # print broken here
#' 
#'
#'
#'
#' 
#'
#' 
#'
#' 
#'
#'
#' 
#'
#' 
#'
#' 
#' }
#' 
bertini <- function(code, dir = tempdir(), opts = "", quiet = TRUE){
	
  ## make dir to put bertini files in (within the tempdir) timestamped
  timeStamp <- as.character(Sys.time())
  timeStamp <- str_replace_all(timeStamp, "-", "_")
  timeStamp <- str_replace_all(timeStamp, " ", "_")  
  timeStamp <- str_replace_all(timeStamp, ":", "_")  
  dir2 <- paste(dir, timeStamp, sep = "/")
  suppressWarnings(dir.create(dir2))
	
  
  ## define a function to write the code to a file
  write_bert <- function(code, codeFile = "bertiniCode"){

    # pull code from function body
    if(is.function(code)) code <- as.character(body(code))[-1]
    if(is.character(code)){
      if(str_sub(code, 1, 1) != "\n") code <- paste("\n", code, sep = "")
      if(str_sub(code, -1, -1) != "\n") code <- paste(code, "\n", sep = "")      
      code <- strsplit(code, "\\n")[[1]][-1]
    }
  
    # write code file
    writeLines(code, con = paste(dir2, codeFile, sep = "/"))
    invisible(code)
  }	
  
  
  ## make bertini file
  write_bert(code)


  ## switch to temporary directory, run bertini
  oldWd <- getwd()
  setwd(dir2)
  
  
  ## run bertini  
  system2(
    paste(getOption("bertiniPath"), "bertini", sep = "/"),
    paste(opts, paste(dir2, "bertiniCode", sep = "/")),
    stdout = "bertiniOut"
  )

  
  ## print bertini output, if requested
  if(!quiet) cat(readLines("bertiniOut"), sep = "\n")
  
  
  ## figure out what files to keep them, and make bertini object
  files <- list.files()
  rawOutput <- as.list(vector(length = length(files)))
  names(rawOutput) <- files
  for(k in seq_along(files)) rawOutput[[k]] <- readLines(files[k])  



  ## convert the raw output into interesting output
  out <- rawOutput
  if("finite_solutions" %in% files) out$finite_solutions <- parse_bertini_finite_solutions(out)
  if("nonsingular_solutions" %in% files) out$nonsingular_solutions <- parse_bertini_nonsingular_solutions(out)  
  if("singular_solutions" %in% files) out$singular_solutions <- parse_bertini_singular_solutions(out)  
  if("real_finite_solutions" %in% files) out$real_finite_solutions <- parse_bertini_real_finite_solutions(out)
  if("raw_solutions" %in% files) out$raw_solutions <- parse_bertini_raw_solutions(out)      
  if("midpath_data" %in% files) out$midpath_data <- parse_bertini_midpath_data(out)      
  if("start" %in% files) out$start <- parse_bertini_start(out) 
  if("failed_paths" %in% files) out$failed_paths <- parse_bertini_failed_paths(out)            

  
  # add code and directory
  out$raw_output <- rawOutput
  out$bertiniCode <- code     
  out$dir <- dir2  
  
  ## migrate back to original working directory
  setwd(oldWd)
  

  
  # class and out
  class(out) <- "bertini"  
  out
}

























































parse_bertini_finite_solutions <- function(rawOutput){

  ## check for no finite solutions
  if(
    length(rawOutput$finite_solutions) == 1 && 
    rawOutput$finite_solutions == ""
  ) return(FALSE)  
  
  ## get variables
  vars <- str_replace(rawOutput$main_data[2], "Variables:  ", "")
  vars <- str_split(vars, " ")[[1]]
  p <- length(vars)  

  ## grab output, format and return
  fSolns <- rawOutput$finite_solutions
  fSolns <- fSolns[-c(1,2)]
  fSolns <- fSolns[-length(fSolns)]
  
  fSolns <- strsplit(fSolns, " ")[nchar(fSolns) > 0]
  fSolns <- sapply(fSolns, function(x){
    x <- zapsmall(as.numeric(x))
    complex(1, x[1], x[2])
  })
  fSolns <- matrix(fSolns, ncol = p, byrow = TRUE)
  colnames(fSolns) <- vars

  fSolns
}







parse_bertini_nonsingular_solutions <- function(rawOutput){

  ## check for no finite solutions
  if(str_sub(rawOutput$nonsingular_solutions[1], 1, 1) == "0") return(FALSE)  
  
  ## get variables
  vars <- str_replace(rawOutput$main_data[2], "Variables:  ", "")
  vars <- str_split(vars, " ")[[1]]
  p <- length(vars)  

  ## grab output, format and return
  nsSolns <- rawOutput$nonsingular_solutions
  nsSolns <- nsSolns[-c(1,2)]
  nsSolns <- nsSolns[-length(nsSolns)]
  
  nsSolns <- strsplit(nsSolns, " ")[nchar(nsSolns) > 0]
  nsSolns <- sapply(nsSolns, function(x){
    x <- zapsmall(as.numeric(x))
    complex(1, x[1], x[2])
  })
  nsSolns <- matrix(nsSolns, ncol = p, byrow = TRUE)
  colnames(nsSolns) <- vars

  nsSolns
}







parse_bertini_singular_solutions <- function(rawOutput){

  ## check for no finite solutions
  if(str_sub(rawOutput$singular_solutions[1], 1, 1) == "0") return(FALSE)  
  
  ## get variables
  vars <- str_replace(rawOutput$main_data[2], "Variables:  ", "")
  vars <- str_split(vars, " ")[[1]]
  p <- length(vars)  

  ## grab output, format and return
  sSolns <- rawOutput$singular_solutions
  sSolns <- sSolns[-c(1,2)]
  sSolns <- sSolns[-length(sSolns)]
  
  sSolns <- strsplit(sSolns, " ")[nchar(sSolns) > 0]
  sSolns <- sapply(sSolns, function(x){
    x <- zapsmall(as.numeric(x))
    complex(1, x[1], x[2])
  })
  sSolns <- matrix(sSolns, ncol = p, byrow = TRUE)
  colnames(sSolns) <- vars

  sSolns
}









parse_bertini_real_finite_solutions <- function(rawOutput){

  ## check for no finite solutions
  if(str_sub(rawOutput$real_finite_solutions[1], 1, 1) == "0") return(FALSE)  
  
  ## get variables
  vars <- str_replace(rawOutput$main_data[2], "Variables:  ", "")
  vars <- str_split(vars, " ")[[1]]
  p <- length(vars)  

  ## grab output, format and return
  rfSolns <- rawOutput$real_finite_solutions
  rfSolns <- rfSolns[-c(1,2)]
  rfSolns <- rfSolns[-length(rfSolns)]
  
  rfSolns <- strsplit(rfSolns, " ")[nchar(rfSolns) > 0]
  rfSolns <- sapply(rfSolns, function(x) zapsmall(as.numeric(x[1])))
  rfSolns <- matrix(rfSolns, ncol = p, byrow = TRUE)
  colnames(rfSolns) <- vars

  rfSolns
}










parse_bertini_raw_solutions <- function(rawOutput){

  ## check for no finite solutions
  if(str_sub(rawOutput$raw_solutions[1], 1, 1) == "0") return(FALSE)  
  
  ## get variables
  vars <- str_replace(rawOutput$main_data[2], "Variables:  ", "")
  vars <- str_split(vars, " ")[[1]]
  p <- length(vars)  

  ## grab output, format and return
  rawSolns <- rawOutput$raw_solutions
  rawSolns <- rawSolns[-c(1,2)]
  rawSolns <- rawSolns[-length(rawSolns)]
  rawSolns <- rawSolns[str_detect(rawSolns, " ")]
  
  rawSolns <- strsplit(rawSolns, " ")
  rawSolns <- sapply(rawSolns, function(x){
    x <- zapsmall(as.numeric(x))
    complex(1, x[1], x[2])
  })
  rawSolns <- matrix(rawSolns, ncol = p, byrow = TRUE)
  colnames(rawSolns) <- vars

  rawSolns
}










parse_bertini_midpath_data <- function(rawOutput){

  ## check for no finite solutions
  if(
    length(rawOutput$midpath_data) == 1 && 
    rawOutput$midpath_data == ""
  ) return(FALSE)  
  
  ## get variables
  vars <- str_replace(rawOutput$main_data[2], "Variables:  ", "")
  vars <- paste(vars, "homog")
  vars <- str_split(vars, " ")[[1]]
  p <- length(vars)  

  ## grab output, format and return
  mdpthPts <- rawOutput$midpath_data
  mdpthPts <- mdpthPts[-length(mdpthPts)]
  mdpthPts <- mdpthPts[nchar(mdpthPts) > 0]
  
  pthMarkerNdx <- which(unname(sapply(mdpthPts, nchar)) == 1)
  mdpthPts <- mdpthPts[-pthMarkerNdx]
  mdpthPts <- strsplit(mdpthPts, " ")
  mdpthPts <- sapply(mdpthPts, function(x){
    x <- zapsmall(as.numeric(x))
    complex(1, x[1], x[2])
  })  
  
  mdpthPts <- matrix(mdpthPts, ncol = p, byrow = TRUE)
  colnames(mdpthPts) <- vars

  mdpthPts
}











parse_bertini_start <- function(rawOutput){

  ## check for no finite solutions
  if(
    length(rawOutput$start) == 1 && 
    rawOutput$start == ""
  ) return(FALSE)  
  
  ## get variables
  vars <- str_replace(rawOutput$main_data[2], "Variables:  ", "")
  vars <- paste(vars, "homog")
  vars <- str_split(vars, " ")[[1]]
  p <- length(vars)  

  ## grab output, format and return
  startPts <- rawOutput$start
  startPts <- startPts[-length(startPts)]
  startPts <- startPts[str_detect(startPts, " ")]
  startPts <- str_replace_all(startPts, ";", "")
  
  startPts <- strsplit(startPts, " ")
  startPts <- sapply(startPts, function(x){
    x <- zapsmall(as.numeric(x))
    complex(1, x[1], x[2])
  })  
  
  startPts <- matrix(startPts, ncol = p, byrow = TRUE)
  colnames(startPts) <- vars

  startPts
}






parse_bertini_failed_paths <- function(rawOutput){

  ## check for no finite solutions
  if(
    length(rawOutput$failed_paths) == 1 && 
    rawOutput$failed_paths == ""
  ) return(FALSE)  
  
  rawOutput$failed_paths
}

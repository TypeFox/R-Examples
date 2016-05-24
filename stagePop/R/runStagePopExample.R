#' runStagePopExample
#'
#' This function is similar to the demo() function but requires less interaction
#' It is used to run the canned examples from the stagePop package.
#'
#' @param name Name of the example to run. If Name is NULL the list of examples will be printed. 
#'
#' @export
runStagePopExample <- function(name=NULL) {

  if (is.null(name)) {
    all_examples <- grep("functions",
                         sapply(
                                Sys.glob(paste(system.file("DemoFiles",package="stagePop"),"/*.R",sep="")),
                                function(filename) {
                                  return(sub("^([^.]*).*", "\\1", basename(filename)))
                                }),
                         invert=TRUE,
                         value=TRUE)

    cat("List of stagePop examples:\n\n")
    for (example in all_examples) {
      cat(paste(example,"\n"))
    }
    
  } else {
    example_file=system.file("DemoFiles",paste(name,".R",sep=""),package="stagePop")

    if (file.exists(example_file)) {
      print(paste("*** RUNNING EXAMPLE FILE",example_file))
      source(example_file)
    } else {
      stop(paste("No example named '",name,"' exists"))
    }   
  }
}

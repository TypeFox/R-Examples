#' Start MiTCR directly from the package.
#' 
#' @description
#' Start the MiTCR tools directly from the package with given settings.
#' 
#' @param .input,.output Input and output files.
#' @param .file.path Path prepending to \code{.input} and \code{.output}.
#'        If input and output is empty, but .file.path is specified, 
#'        than process all files from the folder \code{.file.path}
#' @param ... Specify input and output files and arguments 
#'            of the MITCR without first '-' to run it.
#' @param .mitcr.path Path to MiTCR .jar file.
#' @param .mem Volume of memory available to MiTCR.
#' 
#' @details
#' Don't use spaces in paths!
#' You should have insalled JDK 1.7 to make it works.
#' 
#' @examples
#' \dontrun{
#' # Equal to
#' # java -Xmx8g -jar ~/programs/mitcr.jar -pset flex 
#' #      -level 2 ~/data/raw/TwA1_B.fastq.gz ~/data/mitcr/TwA1_B.txt
#' startmitcr('raw/TwA1_B.fastq.gz', 'mitcr/TwA1_B.txt', .file.path = '~/data/',
#'   pset = 'flex', level = 1, 'debug', .mitcr.path = '~/programs/', .mem = '8g')
#' }
startmitcr <- function (.input = '', .output = '', ..., .file.path = '', .mitcr.path = '~/programs/', .mem = '4g') {  
  if (.input == '') {
    startmitcr(list.files(.file.path), paste0('/mitcr/', list.files(.file.path), '.txt'), ..., .file.path = .file.path, .mitcr.path = .mitcr.path, .mem = .mem)
  }
  else if (length(.input) > 1) {
    for (i in 1:length(.input)) {
      startmitcr(.input[i], .output[i], ..., .file.path = .file.path, .mitcr.path = .mitcr.path, .mem = .mem)
    }
  } else {
    if (.file.path[length(.file.path)] != '/' && .file.path != '') {
      .file.path <- paste(.file.path, '/', sep = '')
    }
    .input <- paste(.file.path, .input, sep = '')
    .output <- paste(.file.path, .output, sep = '')
    
    args <- list(...)
    args.names <- names(args)
    if (is.null(names(args))) {
      args.names <- rep('', length(args))
    }
    args.str <- ''
    for (i in 1:length(args)) {
      if (args.names[i] == '') {
        args.str <- paste(args.str, paste('-', args[[i]], sep = ''))
      } else {
        args.str <- paste(args.str, paste('-', args.names[i], sep = ''), args[[i]])
      }
    }
    args.str <- paste(args.str, .input, .output)
    
    system(paste('java',
                 paste('-Xmx', .mem, sep = ''),
                 '-jar',
                 paste(.mitcr.path, 'mitcr.jar', sep = ''),
                 args.str))
  }
}


# Make EPS more optimal, and useable by PostScriptTrace
# .ps2ps <- function (file.in, file.out) {
#   gsexe <- Sys.getenv("R_GSCMD")
#   if (.Platform$OS.type == "windows") {
#     if (!nzchar(gsexe)) 
#       gsexe <- Sys.getenv("GSC")
#     if (is.null(gsexe) || !nzchar(gsexe)) {
#       poss <- Sys.which(c("gswin64c.exe", "gswin32c.exe"))
#       poss <- poss[nzchar(poss)]
#       gsexe <- if (length(poss)) 
#         poss
#       else "gswin32c.exe"
#     }
#     else if (grepl(" ", gsexe, fixed = TRUE)) 
#       gsexe <- shortPathName(gsexe)
#     
#   }
#   else {
#     if (is.null(gsexe) || !nzchar(gsexe)) {
#       gsexe <- "gs"
#       rc <- system(paste(shQuote(gsexe), "-help > /dev/null"))
#       if (rc != 0) 
#         stop("sorry, 'gs' cannot be found")
#     }
#   }
#   
#   file.out = suppressWarnings( normalizePath(file.out) )
#   cmd <- paste(gsexe, " -q -sDEVICE=ps2write -sstdout=%stderr -sOutputFile=", 
#                file.out, " -dNOPAUSE -dBATCH -P- -dSAFER ", normalizePath(file.in), 
#                sep = "")
#   ret <- switch(.Platform$OS.type, unix = system(cmd), windows = system(cmd, invisible = TRUE))
#   if (ret != 0) stop(gettextf("status %d in running command '%s'", ret, cmd), domain = NA)
#   
#   return(cmd)
# }

# plotlogo <- function(file){
#   if(!file.exists(file)){
#     writeLines(sprintf('File "%s" does not exist!', file))
#     return(1)
#   }
#   data = readJPEG(file)
#   plot.new()
#   lim = par()$usr
#   rasterImage(data, lim[1], lim[3], lim[2], lim[4])
# }

# Display a sequence logo as an R plot
#
# This function will display a sequence logo of an already generated eps sequence logo. 
#
# param file path to the eps formatted file of the sequence logo
# 
# export
#  
# examples
# # Get path to an example eps sequence logo
# fpath = system.file("extdata", "example_logo.eps", package="RWebLogo")
# # Plot it!
# plotlogo(fpath)
# plotlogo <- function(file){
#   if(!file.exists(file)){
#     writeLines(sprintf('File "%s" does not exist!', file))
#     return(1)
#   }
#   if(!grepl('.eps$', basename(file))){
#     writeLines(sprintf('The file "%s" must be eps format to plot. To save your logo as in eps format, set the "format" option to "eps" in the weblogo call', file))
#     return(1)
#   }
#   
#   filenew = tempfile(pattern = 'filenew')
#   newps = .ps2ps(file, filenew)
#   
#   outfile = tempfile(pattern = 'outfile')
#   PostScriptTrace(file = filenew, outfilename = outfile)
#   pic = readPicture(outfile)
#   grid.newpage()
#   grid.picture(pic)
# }

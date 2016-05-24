check.utility <- function(x = c("muscle", "dssp", "stride", "mustang", "makeup"), 
   quiet = TRUE) {
  
  utilities <- match.arg(x, several.ok = TRUE)

  ##- Check on missing utility programs
  missing.util <- nchar(Sys.which(utilities)) == 0
  if( any(missing.util) ) {
    if(!quiet) {
       warning(paste0("  Checking for external utility programs failed\n",
         "    Please make sure '", paste(names(missing.util[missing.util]), collapse="', '"),
         "' is in your search path, see:\n",
         "    http://thegrantlab.org/bio3d/tutorials/installing-bio3d#utilities"))
    }
    pass = FALSE
  } else {
    if(!quiet) cat("External utility programs found\n")
    pass = TRUE
  }
  invisible(pass)
}  

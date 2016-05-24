service <- Sys.getenv("_R_CHECK_HAVE_PERLCSVXS_")

Sys.info()

if(identical(as.logical(service), TRUE)) {
   require("TSmisc") 
   if( ! any(c("XLS","XLSX") %in% gdata::xlsFormats())) stop("Perl libraries are not available." )
 }else {
   cat("Perl libraries not available. Skipping tests.\n")
   cat("_R_CHECK_HAVE_PERLCSVXS_ setting ", service, "\n")
   }

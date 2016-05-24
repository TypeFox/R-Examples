"boa.importBUGS" <-
function(prefix, path = NULL)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   link <- NULL
   filenames <- paste(path, paste(prefix, c(".ind", ".out"), sep = ""), sep = "/")
   import <- try({
      if(all(file.exists(filenames))) {
         ind <- scan(filenames[1], list(pnames = "", first = 0, last = 0))
         out <- scan(filenames[2], list(iter = 0, parms = 0))
         iter.first <- max(out$iter[ind$first])
         iter.last <- min(out$iter[ind$last])
         if(iter.first <= iter.last) {
            idx <- match(c(iter.first, iter.last), out$iter)
            iter <- out$iter[idx[1]:idx[2]]
            link <- matrix(out$parms[(iter.first <= out$iter) &
                                    (out$iter <= iter.last)],
                           nrow = length(iter), ncol = length(ind$pnames),
                           dimnames = list(iter, ind$pnames))
         } else {
            cat("Warning: import failed\n",
                "No common iterations to import.\n")
         }
      } else {
         cat("Warning: import failed\n",
             " Could not find '", filenames[1], "' or '", filenames[2], "'.\n",
             sep = "")
      }
   }, TRUE)
   
   if(inherits(import, "try-error"))
      cat("Warning: import failed\n",
          " File formats do not conform to the CODA standard.  Confirm that\n",
          " '", filenames[1], "' and '", filenames[2], "'\n",
          " are saved as ASCII text with the correct extensions.\n", sep = "")

   return(link)
}

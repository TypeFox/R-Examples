
showMaxControl <- function(object) {
   cat("A 'MaxControl' object with slots:\n")
   for(s in slotNames(object)) {
      if(s == "sann_cand") {
         ## This is a function or NULL, handle with care:
         if(is.null(slot(object, s))) {
            cat("sann_cand = <default Gaussian Markov kernel>\n")
         }
         else {
            cat("sann_cand =\n")
            print(str(slot(object, s)))
         }            
      }
      else {
         ## Just print
         cat(s, "=", slot(object, s), "\n")
      }
   }
}

setMethod("show", "MaxControl", showMaxControl)

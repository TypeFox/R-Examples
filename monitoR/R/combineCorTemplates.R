# For combining correlation template lists together
# Modified: 2013 JUNE 13

combineCorTemplates <-
function(...) {

   temps <- list(...)
   x <- NULL

   for(i in 1:length(temps)) {
      x <- c(x, temps[[i]]@templates)
   }

   templates <- new('corTemplateList', templates=x)

   while(any(dups <- duplicated(names(templates@templates)))) {
     warning("Found identical template names: ", paste(names(templates@templates)[dups], collapse=", "), ". Adding suffix.")   
     names(templates@templates)[dups] <- paste0(names(templates@templates)[dups], ".2")
   }

   return(templates)

}


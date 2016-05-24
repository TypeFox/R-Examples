acshgral <- function(dotdata, species = as.integer(species)) {
   if (is.null(class(dotdata)) | class(dotdata) != "dotdata") {
       cat("Argument is not of class 'dotdata' \n")
       return(invisible())
    }
   stopifnot(length(species) > 0)
   species <- unique(species)
   eluq <- setdiff(species, 1:length(dotdata$Label))
   stopifnot(length(eluq)==0) 
   pts <- w <- c() 
   for(i in species) {
    pts <- c(pts, dotdata$occupancy[[i]]) 
    w <- c(w, dotdata$MSTsp[[i]]$wght) 
   }
   aux <- sum(w*apply(dotdata$PtSpDist[pts, species, drop = FALSE], 1, max))
   return(1/length(species)*aux)
}


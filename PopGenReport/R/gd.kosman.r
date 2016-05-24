gd.kosman <- function(population){
  # check to see if the function was passed a genind object
  if (class(population) != "genind") {
    message("You did not provide a valid genind object! Script stopped!")
    return
  }
  
  # getting the number of loci
  n <- length(locNames(population)) 
  
  # getting the ploidy of an animal....
  uniqueploidy<-unique(population@ploidy)
  if(length(uniqueploidy)==1){
    ploidy<-uniqueploidy
  } else if(uniqueploidy>1){
    message("Your data set has multiple ploidies, please separate loci by ploidy")
    message("Script stopped!")
    return
  } else if(uniqueploidy<=0){
    message("Your dataset has an invalid ploidy (ploidy<=0). Script stopped!")
    return
  }
  
  # this calculates the manhattan distance between each individual and adjusts for ploidy..
  matrices <- lapply(seploc(population), function(l) as.matrix(dist(l@tab, "manhattan")/(2*ploidy)))
  
  # if the value is missing, mark it with a 1, if it is real, mark it 0
  missing <- lapply(matrices, function(m) ifelse(is.na(m), 1, 0))
  
  # This is going to replace missing numbers with 0. 
  replaced <- lapply(matrices, function(m) ifelse(is.na(m),0,m))
  
  
  loci.used<-(n-Reduce("+", missing))
  colnames(loci.used) <- indNames(population)
  rownames(loci.used) <- indNames(population)
  
  # This sums the values across lists and then divides by the number of loci compared less loci with missing numbers  
  d.fast<-(Reduce("+", replaced)/loci.used)
  colnames(d.fast) <- indNames(population)
  rownames(d.fast) <- indNames(population)
  
  # clean up matrices for export
  d.fast[upper.tri(d.fast, diag = FALSE)] <- NA
  loci.used[upper.tri(loci.used, diag = FALSE)] <- NA
  kosman.out <- list(geneticdist = d.fast, loci_used = loci.used)
  return(kosman.out)
}

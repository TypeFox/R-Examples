siarelicit <-
function(siardata) {

if(siardata$SHOULDRUN==FALSE || siardata$GRAPHSONLY ==TRUE) {
    cat("You must load in some data first (via option 1) in order to use this feature of the program. \n")
    cat("Press <Enter> to continue")
    readline()
    invisible()
    return(siardata)
}

cat("This function allows you to elicit parameters for the Dirichlet \n")
cat("distribution based on estimates for the mean proportions and the \n")
cat("standard deviations from one of them. \n\n")

sourcenames <-  as.character(siardata$sources[,1])
cat(paste("There are",length(sourcenames),"current sources. They are: \n"))
cat(sourcenames,"\n \n")

cat("===========================================================================\n")
cat("Mean proportion estimates. \n")
cat("===========================================================================\n\n")
cat("Please first enter your mean estimate of the proportions of each of these \n")
cat(" sources separated by a space. (Note: these should add up to 1) \n")
propsum <- 0
while(propsum!=1) {
  meanprops <- as.numeric(scan(what="list",nlines=1,quiet=TRUE,sep=" "))
  propsum <- sum(meanprops)
  if(length(meanprops)!=length(sourcenames)) {
    cat(paste("Only",length(sourcenames),"allowed. \n"))
    propsum <- 0
  }
  if(propsum!=1) cat("Does not sum to 1, try again. \n")
}

cat("Thank you. \n \n")

cat("===========================================================================\n")
cat("Standard deviation estimate. \n")
cat("===========================================================================\n\n")

qu <- "Which source would you like to enter the standard deviation for:"
sdans <- menu(sourcenames,title = qu)

cat(paste("Please now enter the standard deviation for source:",sourcenames[sdans],"\n"))
badsd <- TRUE
while(badsd==TRUE) {
  sourcesd <- as.numeric(scan(what="",nlines=1,quiet=TRUE))
  if(sourcesd < 0) {
    cat("Bad standard deviation entered. Try again. \n")
  } else if((meanprops[sdans]*(1-meanprops[sdans])/sourcesd^2)<1) {
    cat("Standard deviation does not satisfy Dirichlet constraints. Try again. \n")
  } else {
    badsd <- FALSE
  }
}

cat("Thank you. \n \n")

Q <- (meanprops[sdans]*(1-meanprops[sdans])/sourcesd^2)-1
alphapars <- meanprops*Q

cat("Dirichlet prior parameters are: \n")
cat(alphapars)
cat("\n")
cat("Please enter these as arguments in the siarmcmcdirichletv4 function. \n")

}


dmm.array <-
function(mdf,fixform,components,cohortform,posdef,gls,glsopt,dmeopt,ncomp.pcr,relmat,dmekeep,dmekeepfit)
# dmm.array()   -   do the traits pairwise and return a dmmarray object
{
  if(is.null(mdf$rel)) {
    df <- mdf
  } else {
    df <- mdf$df
  }
  if(!exists("Ymat",df)) {
    stop("dmm: dataframe must contain 'Ymat' for traitspairwise option\n")
  }
  traits <- colnames(df$Ymat)
  l <- length(traits)
  if(l <= 0) {
    stop("dmm: there must be at least one trait in 'Ymat'\n")
  }
  pdum <- matrix(diag(l),l,l,dimnames=list(traits,traits))
  fit <- array(make.dmmobj(pdum,,pdum),c(l,l))
  dimnames(fit) <- list(traits,traits)
  for(i in traits) {
    for(j in traits) {
      ymat <- cbind(df[,i],df[,j])
      dimnames(ymat) <- list(NULL,c(i,j))
      if(is.null(mdf$rel)) {
         mdf$Ymat <- ymat # put local ymat into mdf so passed to dmm()
       }
       else {
         mdf$df$Ymat <- ymat # put local ymat into mdf$df so passed to dmm()
       }
      fit[[i,j]] <- dmesolve(mdf,fixform,components,cohortform,posdef,gls,glsopt,dmeopt,ncomp.pcr,relmat,dmekeep,dmekeepfit)
      class(fit[[i,j]]) <- "dmm"
    }
  }
  class(fit) <- "dmmarray"
  return(fit)
}

dmm.blockarray <-
function(mdf,fixform,components,cohortform,posdef,gls,glsopt,dmeopt,ncomp.pcr,relmat,dmekeep,dmekeepfit,...)
#  dmm.blockarray()   - do traits in pairs of blocks of traits
#                    - return dmmblockarray object
{
  if(is.null(mdf$rel)) {
  df <- mdf
  }
  else {
  df <- mdf$df
  }
  if(!exists("Ymat",df)) {
    stop("dmm: dataframe must contain 'Ymat' for traitsblockwise option\n")
  }
  all.blocks <- list(...)
  blocknames <- names(all.blocks)
  b <- length(all.blocks)
  if(b <= 0) {
    stop("dmm: no blocks specified\n")
  }
  pdum <- matrix(diag(b),b,b,dimnames=list(blocknames,blocknames))
  fit <- array(make.dmmobj(pdum,,pdum),c(b,b))
  dimnames(fit) <- list(blocknames,blocknames)
  ntraits <- rep(0,b)
  for ( i in blocknames) {
    ntraits[i] <- length(all.blocks[[i]])
    for(j in blocknames) {
       subi <- all.blocks[[i]]
       subj <- all.blocks[[j]]
       ymat <- as.matrix(cbind(df[,subi],df[,subj]))
       dimnames(ymat) <- list(NULL,c(subi,subj))
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
  retobj <- list(array=fit, blocks=all.blocks)
  class(retobj) <- "dmmblockarray"
  return(retobj)
}

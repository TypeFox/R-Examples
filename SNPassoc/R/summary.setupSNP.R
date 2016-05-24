`summary.setupSNP` <-
function(object, ...)
 {
  if (!inherits(object, "setupSNP")) 
        stop("object must be an object of class 'setupSNP'")
  colSNPs<-attr(object,"colSNPs")
  if(length(colSNPs)>0){
  temp <- mclapply(object[,colSNPs, drop=FALSE], expandsetupSNP)

  ans<-temp[[1]] 
  i<-2
  while (i <= length(temp))
   { 
    ans<-rbind(ans,temp[[i]])
    i<-i+1 
   } 
  dimnames(ans)[[1]] <- attr(object,"label.SNPs")
  out<-as.matrix(ans)
  dimnames(out)[[2]][4] <- "missing (%)"
  print(out, quote=FALSE, na.print="-")
  } else {
    class(object)<-"data.frame"
    ans<-summary(object)
    print(ans)
  }
  invisible(ans)
 }


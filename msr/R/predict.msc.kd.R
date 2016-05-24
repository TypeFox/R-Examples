predict.msc.kd <- function (object, newdata, addExtrema=TRUE, ...) 
{
   
    xe <- model.matrix(object, newdata)
    x <- object$x
    nc <- ncol(x)
    nrxe <- nrow(xe)
   
    msLevel <- object$level[[object$predictLevel]]

    #create kernel density estimators
    d <- matrix(ncol=length(msLevel$mins), nrow=nrxe)
    xet <- t(xe);
    
    for(i in 1:length(msLevel$mins)){
      index = msc.level.ind(msLevel, i, addExtrema);
      xtmp <- t(x[index, ])
      nct <- ncol(xtmp)
      p <- .C("gkde", as.integer(nc), as.integer(nct), as.double(xtmp), 
              as.double(object$bw), as.integer(nrxe), as.double(xet), 
              p = double(nrxe))$p
      d[,i] <- p;
    }
    d <- d/rowSums(d)
    d
}


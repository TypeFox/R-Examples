model.matrix.msr <- function (object, newdata, ...) 
{
    ms <- object$ms;
    if( ! is.null(ms$nLevels) ){
      ms <- ms$mscl[[object$predictLevel]]
    }
    nc <- ncol(ms$x)
    nr <- nrow(newdata)    

    x <- model.matrix(ms, newdata)

    #compute weights for each crystals    
    d <- predict(ms, x)
    xnames <-colnames(ms$x)
    cxnames <- c()
    x <- as.matrix(x)

    #fit model to each crystals 
    mins <- ms$mins
    mAll = matrix(nrow=nr, ncol=(1+nc)*length(mins))
    j <- 1
    for(i in 1:length(mins)){
      mAll[,j] <- d[, i]
      mAll[,(j+1):(j+nc)] <- x * d[,i]
      j <- j+nc+1
      s = paste("c", as.character(i), sep="")
      cxnames <- append(cxnames, c(paste(s,"intercerpt", sep="."), paste(s,xnames, sep=".")))
    }
    colnames(mAll) <- cxnames
    mAll
}

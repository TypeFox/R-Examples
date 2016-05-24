"faces2" <-
function(mat, which=1:ncol(mat), labels=rownames(mat),
                  nrows=ceiling(nrow(mat)/ncols),
                  ncols=ceiling(sqrt(nrow(mat))),
                  byrow=TRUE,
                  scale=c("columns","all","center","none"),
                  fill=c(.5,.5,1,.5,.5,.3,.5,.5,.5,.5,.5,.5,.5,.5,
                    .5,.5,1,.5),
                  ...) {

  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))

  if(byrow){
    par(mfrow=c(nrows,ncols))
  } else {
    par(mfcol=c(nrows,ncols))
  }

  par(mar=rep(0,4))
            
  mat <- as.matrix(mat)
  scale <- match.arg(scale)

  if(scale=="columns") {
    mat <- sweep(mat,2, apply(mat,2,min,na.rm=TRUE), '-')
    mat <- sweep(mat,2, apply(mat,2,max,na.rm=TRUE), '/')
  } else if(scale=="all") {
    mat <- mat - min(mat,na.rm=TRUE)
    mat <- mat / max(mat,na.rm=TRUE)
  } else if(scale=="center"){
    mat <- sweep(mat, 2, apply(mat,2,mean,na.rm=TRUE), '-')
    mat <- sweep(mat, 2, apply(abs(mat),2,max,na.rm=TRUE), '/')
    mat <- (mat+1)/2
  }

  if(ncol(mat) > 18){
    warning("using only first 18 columns of input")
    mat <- mat[,1:18]
  }

  mat2 <- matrix(fill, ncol=18, nrow=nrow(mat), byrow=TRUE)

  mat2[,which] <- mat

  lo <- c(rep(0.2, 5), 0.1, 0.2, 0, 0.2, 0.1, 0.1, 0.3, 0.1, 0.3, rep(0.1, 4))
  hi <- c(0.8, 0.8, 1, 0.8, 0.8, 0.4, 0.8, 1, 0.8, 0.7, 0.9, 0.7, rep(0.9, 4),
          1, 0.9)
  df <- hi-lo



  mat2 <- sweep(mat2, 2, df, '*')
  mat2 <- sweep(mat2, 2, lo, '+')

  ## special handeling for column 8
  mat2[,8] <- (2*mat2[,8]-1)*mat2[,9]

  if(length(labels != nrow(mat2))){
    labels=rep(labels,nrow(mat2))[1:nrow(mat2)]
  }
  
     
  for (i in 1:nrow(mat2)){
    face2.plot(mat2[i,])
    text(0,-500,labels[i],...)
  }
  invisible()
}


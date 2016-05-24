tmean<-function(x, alpha=0.2, plotting=FALSE, new=TRUE, cols=c(1,4,8),...)
{
n <- nrow(x);
I <- MBD(x,plotting=FALSE)$ordering;
if (is.null(alpha)){stop("alpha must be a number greater than or equal to 0 and smaller than 1...")}

if (plotting){
  if (length(alpha)==1){
    N <- n-floor(alpha*n)
    m1 <- x[I[1:N],] 
    tm <- apply(m1,2,mean) 
    Gene.Expression<-t(x)
    matplot(Gene.Expression, type="l",col=cols[3],...)
    matlines(t(m1),col=cols[2],...)
    lines(tm,col=cols[1],...)
    legend("top",col=cols, lty=1, legend=c(paste(alpha,"-trimmed mean",sep=""), paste((1-alpha)*100,"% remaining samples",sep=""),paste(alpha*100,"% trimmed out samples", sep="")))
    return(list(tm=tm, tm.x=m1))
    }  ## end plotting one trimmed-mean
  else {
    alpha <- unique(sort(alpha))
    j<-0; tm<-list(); m1<-list()
    Gene.Expression <- c()
    for (a in alpha) {
      if ((a<0) | (a>=1)) {stop("alpha must contain numbers greater than or equal to 0 and smaller than 1...")}
      j <- j+1
      N <- n-floor(a*n)
      m1[[j]] <- x[I[1:N],]
      tm[[j]] <- apply(m1[[j]],2,mean) 
      Gene.Expression <- cbind(Gene.Expression, tm[[j]])
      }  ## end FOR a
    if (new) {matplot(Gene.Expression,type="l",col=colorRampPalette(c(cols[1],'gray'))(length(alpha)),...)} 
    else {matlines(Gene.Expression, col=colorRampPalette(c(cols[1],'gray'))(length(alpha)),...)}   
    return(list(tm=tm, tm.x=m1))
    }  ## end several trimmed means
  }  ## end IF plotting
else { ## not plotting
  if (length(alpha)>1) {
    alpha <- unique(sort(alpha))
    j<-0; tm<-list(); m1<-list()
    for (a in alpha) {
      if ((a<0) | (a>1)) {stop("alpha must contain numbers greater than or equal to 0 and smaller than 1...")}
      j <- j+1
      N <- n-floor(a*n)
      m1[[j]] <- x[I[1:N],]
      tm[[j]] <- apply(m1[[j]],2,mean) 
      }  ## end FOR a
    }  ##  end several tmeans
  else {
    if ((alpha<0) | (alpha>1)) {stop("alpha must be a number greater than or equal to 0 and smaller than 1...")}
    N <- n-floor(alpha*n)
    m1 <- x[I[1:N],] 
    tm <- apply(m1,2,mean) 
    return(list(tm=tm, tm.x=m1))
    }  ## end ELSE
  }  ##  end ELSE
}


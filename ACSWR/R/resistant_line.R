resistant_line <-
function(x,y,iterations)  {
  three_medians <- function(x,y) {
    n <- length(x)
    k <- n %% 3
    dix <- sort(x,index.return=TRUE)$ix
    x <- x[dix]; y <- y[dix]
    if(k==0)  {
      t <- n/3
      xleft <- x[1:t]; xmid <- x[(t+1):(2*t)]; xright <- x[(2*t+1):n]
      yleft <- y[1:t]; ymid <- y[(t+1):(2*t)]; yright <- y[(2*t+1):n]
    }
    if(k==1)  {
      t <- (n-1)/3
      xleft <- x[1:t]; xmid <- x[(t+1):(2*t+1)]; xright <- x[(2*t+2):n]
      yleft <- y[1:t]; ymid <- y[(t+1):(2*t+1)]; yright <- y[(2*t+2):n]
    }
    if(k==2)  {
      t <- (n-2)/3
      xleft <- x[1:t+1]; xmid <- x[(t+2):(2*t+1)]; xright <- x[(2*t+2):n]
      yleft <- y[1:t+1]; ymid <- y[(t+2):(2*t+1)]; yright <- y[(2*t+2):n]
    }
    xlm <- median(xleft); xmm <- median(xmid); xrm <- median(xright)
    ylm <- median(yleft); ymm <- median(ymid); yrm <- median(yright)
    xmed = c(xlm,xmm,xrm); ymed = c(ylm,ymm,yrm)
    b1 <- (yrm-ylm)/(xrm-xlm)
    # b0 <- ((ylm+ymm+yrm)-b1*(xlm+xmm+xrm))/3
    b0 <- mean(ymed-b1*(xmed-xmed[2]))
    bl <- (ymm-ylm)/(xmm-xlm); br <- (yrm-ymm)/(xrm-xmm)
    the_medians = data.frame(xmed=xmed,ymed=ymed)
    return(list(the_medians=the_medians,coeffs = c(b0,b1),bl=bl,br=br,xCenter=xmm))
  }  
  tms <- three_medians(x,y)$the_medians
  b0=b1=0
  bl=br=0
  resid <- y
  for(i in 1:iterations)  {
    currmodel <- three_medians(x,resid)
    b0 <- b0+currmodel$coeffs[1]
    b1 <- b1+currmodel$coeffs[2]
    if(i==1) {
      bl <- currmodel$bl
      br <- currmodel$br
    }
    resid <- y - (b0+b1*(x-currmodel$xCenter))  
  }
  return(list(coeffs=c(b0,b1),tms=tms,xCenter = currmodel$xCenter, hsr = br/bl))
}

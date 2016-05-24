MBD<- function(x, xRef=NULL, plotting=TRUE, grayscale=FALSE, band=FALSE, band.limits=NULL, 
      lty=1, lwd=2, col=NULL, cold=NULL, colRef=NULL, ylim=NULL, cex=1,...)
{
n <- nrow(x); d <- ncol(x) # n: number of observations (samples);  d: dimension of the data
x <- as.matrix(x)

if (length(xRef)==0) {  ## MBD with respect to the same sample

  ## depth computation
  if (ncol(x) == 1) {x <- t(x)}
  depth <- matrix(0,1,n)
  ordered.matrix <- x
  if (n>1) {
     for (columns in 1:d) {
        ordered.matrix[,columns] <- sort(x[,columns])
        for (element in 1:n) {
             index.1 <- length(which(ordered.matrix[,columns] < x[element,columns]))
             index.2 <- length(which(ordered.matrix[,columns] <= x[element,columns]))
             multiplicity <- index.2 - index.1
             depth[element] <- depth[element] + index.1 * (n - (index.2)) + multiplicity * (n - index.2 + index.1) + choose(multiplicity,2)
             }   ### end FOR element
        }  ### end FOR columns
     depth <- depth / (d * choose(n,2) )
     } ## end IF
  if (n==1) {deepest <- x; depth <- 0}
  ordering<-order(depth,decreasing=TRUE)

  ## plot

  if (plotting) {

     par(mar=c(4,5,3,3),xpd=FALSE)  	
     Gene.Expression<-t(x[ordering[n:2],])
     lwdd<-lwd[1]                     ## 'lwd' for the deepest point 
     if (is.null(cold)) {cold<-2}     ## 'col' for the deepest point
     if (is.null(ylim)) {ylim<-range(x)}

     if (band){
        lty<-lty[1]                   ## 'lty' for the deepest point
        if (is.null(band.limits)) {band.limits <- c(0.5,1)}
        else {band.limits<-unique(sort(band.limits))}
        if (floor(n*(band.limits[1])) < 2) {stop("Check the limits. The band must contain at least 2 curves...")}
        no.poly<-length(band.limits)
        if (is.null(col)) {
            if (grayscale) {color<-rev(gray((1:no.poly)/(no.poly+1)))}
            else {color<-3:(2+no.poly)}
            }  ## end IF 
        else {  #  overwrites grayscale
            if (length(col)<no.poly) {color<-rep(col,length.out=no.poly)[1:no.poly]} 
            else {color<-col[1:no.poly]}
            } ## end ELSE
        par(mar=c(4,5,5,3),xpd=FALSE)  	
        matplot(Gene.Expression, ylim=ylim, type="l", lty=0,...)

        # find the polygon limits
        for (poly in no.poly:1){
          limit<-band.limits[poly]
          no.points <- floor(limit*n)
          Gene.Expression<-t(x[ordering[no.points:1],])
          upper<-apply(Gene.Expression,1,max)
          lower<-apply(Gene.Expression,1,min)
          polygon(c(1:d,d:1),c(upper,rev(lower)),col=color[poly])
          }  # end FOR poly
        }  # end IF band

      else  {  ##  not band plots
        if (is.null(col)) {
            if (grayscale) {color<-rev(gray( (1:n)/(n+1) ))}
            else {color<-rep(8,n)}
            }  ## end IF 
        else {  #  overwrites grayscale
            if (length(col)<n) {color<-rep(col,length.out=n)[n:1]} 
            else {color<-col[n:1]}
            } ## end ELSE
        matplot(Gene.Expression, type="l", ylim=ylim, lty=lty, lwd=lwd, col=color[n:2],...)
        }   # end ELSE (no bands)

        lines(x[ordering[1],], lty=lty, lwd=lwdd, col=cold)
        par(xpd=TRUE)
        legend("top", legend="deepest sample",col=cold, lty=lty, lwd=lwdd,cex=cex)
        if (band) {legend("top", inset=-0.1*cex, horiz=TRUE, legend=band.limits, col=color, pch=15, title="Proportion of central samples",cex=cex,bty="n")}
      }  ## end IF plotting
  } ## end IF no reference sample

else {
  xRef <- as.matrix(xRef)
  if (ncol(xRef)!=d) {stop("Dimensions of x and xRef do not match")}
  n0 <- nrow(xRef)

  ## depth computations
  if (ncol(x) == 1) {x <- t(x)}
  depth <- matrix(0,1,n)
  ordered.matrix <- xRef
  if (n0>1) {
     for (columns in 1:d) {
        ordered.matrix[,columns] <- sort(xRef[,columns])
        for (element in 1:n) {
             index.1 <- length(which(ordered.matrix[,columns] < x[element,columns]))
             index.2 <- length(which(ordered.matrix[,columns] <= x[element,columns]))
             multiplicity <- index.2 - index.1
             depth[element] <- depth[element] + (index.1 + multiplicity ) * (n0 - index.1 - multiplicity) + multiplicity * ( index.1 + (multiplicity-1)/2)
             }   ### end FOR element
        }   ### end FOR columns
     depth <- depth / (d * choose(n0,2) )
     } ## end IF
  if (n==1) {deepest <- x; depth <- 0}
  ordering<-order(depth,decreasing=TRUE)

  ## plot
  if (plotting)  {
    par(mar=c(4,5,3,3),xpd=FALSE)
    Gene.Expression <- t(xRef)
    if (is.null(colRef)) {colRef<-4}
    else {colRef<-colRef[1]}  ###only one color for the reference data set
    if (is.null(cold)) {cold<-2}    ## 'col' for the deepest point
    if (is.null(ylim)) {ylim<-range(x,xRef)}
    lwdd<-lwd[1]                    ## 'lwd' for the deepest point

    if (band){
        lty<-lty[1]                  ## 'lty' for the deepest point
        if (is.null(band.limits)) {band.limits <- c(0.5,1)}
        band.limits<-unique(sort(band.limits))
        if (floor(n*(band.limits[1])) < 2) {stop("Check the limits. The band must contain at least 2 curves...")}
        no.poly<-length(band.limits)
        par(mar=c(4,5,5,3),xpd=FALSE)  	
        matplot(Gene.Expression, type="l", ylim=ylim, lty=lty, lwd=lwd/2, col=colRef,...)
        if (is.null(col)) {
            if (grayscale) {color<-rev(gray((1:no.poly)/(no.poly+1)))}
            else {color<-5:(no.poly+4)}
            }  ## end IF 
        else {  #  overrides grayscale
            if (length(col)<no.poly) {color<-rep(col,length.out=no.poly)[1:no.poly]} 
            else {color<-col[1:no.poly]}
            } ## end ELSE
        # find the polygon limits
        for (poly in no.poly:1){
          limit<-band.limits[poly]
          no.points <- floor(limit*n)
          Gene.Expression<-t(x[ordering[no.points:1],])
          upper<-apply(Gene.Expression,1,max)
          lower<-apply(Gene.Expression,1,min)
          polygon(c(1:d,d:1),c(upper,rev(lower)),col=color[poly])
          }  # end FOR poly
        lines(x[ordering[1],], lty=lty, lwd=lwdd, col=cold)
        }  # end IF band

      else  {  ##  no band plots
        matplot(Gene.Expression, type="l", ylim=ylim, lty=lty, lwd=lwd/2, col=colRef,...)
        if (is.null(col)) {
            if (grayscale) {color<-rev(gray(1/(n+1):(n/(n+1))))}
            else {color<-rep(8,n)}
            }  ## end IF 
        else {  #  overwrites grayscale
            if (length(col)<n) {color<-rep(col,length.out=n)[n:1]} 
            else {color<-col[n:1]}
            } ## end ELSE
        matlines(t(x[ordering[n:2],]), lwd=lwd, col=color[n:2], lty=lty)
        lines(x[ordering==1,],lwd=lwdd, lty=lty[1],col=cold)
        }   # end ELSE (no bands)
      par(xpd=TRUE)
      legend("top",legend=c("deepest sample","reference set"),col=c(cold,colRef), 
            lty=lty,lwd=c(lwdd,lwd/2),cex=cex)
      if (band) {legend("top", inset=-0.1*cex, horiz=TRUE, legend=band.limits, col=color, pch=15, title="Proportion of central samples",cex=cex,bty="n")}
      }  ## end IF plotting
  }  ## end ELSE
return(list(ordering=ordering,MBD=depth))
}

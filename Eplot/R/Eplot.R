#library(devtools)
#library(roxygen2)


# Nchars1 = unlist(regexec("x",getwd())) + 1
# file0 = paste(substr(getwd(),start = 1, stop = Nchars1),"myprofile.r", sep = "" )
# FilesToLoad = c(file0)
# sapply(FilesToLoad,source,.GlobalEnv)

#------------ package main -----------------

#setwd('C:/Users/MP-User/Dropbox/E2/Webfolder/posts')


# coll <- c(wes.palette(4, "Moonrise2")[4],wes.palette(4, "Moonrise2")[1:3],
#           wes.palette(4, "GrandBudapest")[1:4]) 


# @keywords K
# @seealso 
# @aliases 

coll <- c("#29211F","#798E87", "#C27D38", "#CCC591", "#F1BB7B","#FD6467",
          "#5B1A18","#D67236") # fancy colors from wesanderson package
lwd1 <- fontt <- 2
cexx=1.3
# parr <- par(mfrow = c(2,2), bg=wes.palette(4, "Chevalier")[3], bty ="n", fg = 'black' ,
#             font.lab = fontt, 
#             font.axis = fontt, font.main = fontt, col.axis = 'black',col.lab = 'black',
#             cex.axis = 1.4, pch = 21,  tck = -0.02, cex.lab = 1.6 ,cex.main = 2,las=1,
#             mar = c(5, 4, 4, 2) + 0.4) 


#' plott
#' 
#' Plotting longitudinal series.
#' 
#' Sets default parameters to get a nicer figure. If y is given then a scatter
#' plot is created. y can also be of class Date. If add=TRUE, add the series to
#' an existing device. Make a sensible choice as to which series should be
#' plotted first so that the ylim and xlim cover are sufficient.
#' 
#' @param x series to be plotted.
#' @param y possible second series, if provided a scatter plot is created.
#' @param add if \code{add = TRUE} the series is added to existing active device. The active device graphical
#' parameters must match, meaning it must be created using \code{return.to.default=F}.
#' @param pch same as in \code{\link{par}}. 
#' @param xlab a label for the x axis, same as in \code{\link{plot.default}}.
#' @param main main title of the plot, same as in \code{\link{plot.default}}.
#' @param col Color of the series, same as in \code{\link{par}}.
#' @param ty character indicating the type of plotting, any of the types as in \code{\link{plot.default}}.
#' @param return.to.default for reverting back to previous par settings.
#' Default is \code{return.to.default=TRUE}.
#' @param ... more graphical parameters can be given as arguments.
#' @return Called for its side effect.
#' @seealso \code{\link{par}}, \code{\link{plot.default}}
#' @export plott
plott = function(x,y=c(1:length(x)),add=FALSE,pch=19,xlab="",col=1,
           main=NULL,ty="b",return.to.default=TRUE,...){
 lwd1 <- fontt <- 2
if(add){
  lines(x~y,xlab="",ylab = "",ty=ty,lwd=lwd1,col=col,pch=pch)
}
if(!add){
  pardefault <- par(no.readonly = TRUE)
  parr = par(bg="white", bty ="n", fg = 'black' , font.lab = fontt, 
             font.axis = fontt, font.main = fontt, col.axis = 'black',col.lab = 'black',
             cex.axis = 1.4, pch = 21,  tck = -0.02, cex.lab = 1.6 ,cex.main = 2,las=1,
             mar = c(5, 4, 4, 2) + 0.4) 
    plot(x~y,xlab=xlab, ylab = "",ty=ty,lwd=2,main=main,pch=pch,col=col,...)
    grid(col = "grey")
  if(return.to.default) { par(pardefault)}
}  }

#' tsideplot
#' 
#' Create a plot of two series vertical axes on both left and right side.
#' 
#' Create a plot of two series with y-axes on both left and right side.
#' \code{Set return.to.default=TRUE} to keep the new settings, otherwise
#' default to revert to previous par values. xaxis parameter is the optional
#' xaxis, if not provided then \code{if(is.null(xaxis)) {xaxis=
#' c(1:length(series1))}} is used.
#' 
#' @param series1,series2 First and second series to be plotted.
#' @param main main title of the plot, same as in \code{\link{plot.default}}.
#' @param return.to.default for reverting back to previous par settings.
#' @param xaxis Optional, the xaxis to be used, see details.
#' @param col Color of the second series, same as in \code{\link{par}}.
#' @param ... more graphical parameters can be given as arguments.
#' @return Called for its side effect.
#' @seealso \code{\link{par}}, \code{\link{plot.default}}.
#' @export tsideplot
tsideplot <- function(series1,series2,main="", return.to.default=T,
                      xaxis=NULL, col='red',...){
  if (NROW(series1) != NROW(series2) ) { stop("Problem is that NROW(series1) != NROW(series2)") }
  lwd1 <- fontt <- 2
    pardefault <- par(no.readonly = TRUE)
  par(mfrow = c(1,1),bg="white", bty ="n", fg = 'black' , font.lab = fontt, 
      font.axis = fontt, font.main = fontt, col.axis = 'black',col.lab = 'black',
      cex.axis = 1.4, pch = 21,  tck = -0.02, cex.lab = 1.6 ,cex.main = 2,las=1,
      mar = c(5, 4, 4, 4) + 0.1) 
  if(is.null(xaxis)) {xaxis= c(1:length(series1))}
  plot(series1~xaxis,ty="b",pch=19,xlab="",ylab="",col = 1,lwd=lwd1,...)
  par(new=TRUE) 
  plot(series2~xaxis,main=main,bty="o",ty="b",pch=19,xlab="",ylab="",
       axes=F,col = col,lwd=lwd1,...)
  axis(4,at=pretty(series2),col=col,col.ticks=col,col.axis=col)
  grid(col = "grey")
  if(return.to.default){par(pardefault)}
}

#' mplott
#' 
#' Multivariate plot.
#' 
#' Multivariate plot. Limited to 5 series. Legend is added automatically using
#' \code{colnames(x)}.
#' 
#' @param x a matrix to be plotted
#' @param wherelegend where to place the legend
#' @param textlegend what should the legend read (see details).
#' @param main main title of the plot, same as in \code{\link{plot.default}}.
#' @param return.to.default for reverting back to previous par settings.
#' Default is \code{return.to.default=TRUE}
#' @param ... more graphical parameters can be given as arguments.
#' @return called for its side effect.
#' @seealso \code{\link{par}}, \code{\link{plot.default}}.
#' @export mplott
mplott <- function(x,wherelegend='bottomleft',textlegend=colnames(x),main="",
                   return.to.default=T,...){
  lwd1 <- fontt <- 2
  pardefault <- par(no.readonly = TRUE)
  l = NCOL(x)
  pchh <- c(19,1:4)
  if (l>5) { 
    l=5
    textlegend=colnames(x)[1:5]
    print("Only columns 1:5 are plotted") 
  }
  par(bg="white", bty ="n", fg = 'black' , font.lab = fontt, 
      font.axis = fontt, font.main = fontt, col.axis = 'black',col.lab = 'black',
      cex.axis = 1.4, pch = 21,  tck = -0.02, cex.lab = 1.6 ,cex.main = 2,las=1,
      mar = c(5, 4, 4, 4) + 0.1) 
  x <- as.matrix(x)
    plot(x[,1],ty = "b",ylim = c(min(na.omit(x)),max(na.omit(x))),xlab = "",
       ylab = "",lwd=2,pch=pchh[1],main=main,...)
  for (i in 2:l){
    if(l==1) {next}
    lines(x[,i],col=i,ty = "b",lty=i,lwd=2,pch=pchh[i],...)
  }
  if(is.null(textlegend)) { textlegend = paste("Series", 1:l)  }
  legend(wherelegend,textlegend,col=1:l,lty=1:l,cex=1.3,text.font=2,
         ncol=2,bty="n",text.col=1:l,lwd = rep(2,l),pch=pchh[1:l])
grid(col='grey')
if(return.to.default){par(pardefault)}
}
#mplott(x,wherelegend='topleft')
#x <- cbind(rnorm(100),rnorm(100),rnorm(100))

# mplott(cbind(x,y,z))

#' lagmat
#' 
#' Creates a lagged matrix with the desired number of lags.
#' 
#' @param x the series to be lagged
#' @param lags number of lags desired
#' @return matrix with dimension \code{[NROW(x),length(lags)]}
#' @examples
#' 
#' x = rnorm(100)
#' lx <- lagmat(x,2)
#' tail(lx)
#' tail(x)
#' 
#' @export lagmat
lagmat <- function(x,lags){
  TT <- NROW(x) ; colnam <- NULL
  lagmatt <- matrix(nrow = TT, ncol = lags) 
  for (i in 1:lags){
    lagmatt[(i+1):(TT),i]  <- x[1:(TT-i)]
    colnam[i] <- paste('Lag',i,sep="")
  }
  lagmatt[1:lags,]  <- NA
  colnames(lagmatt) <- colnam
  lagmatt
}

#' linpred
#' 
#' Provides linear regression based predictions from a \code{y~x} type model
#' using recursive or rolling regression.
#' 
#' The training is done using the direct method: \eqn{y_{1 : (t+h-1)} = \beta
#' x_{1:(t-1)} + \varepsilon_{1:(t+h-1)} } and the forecast is made at time
#' (t+h) as \eqn{\widehat{y}_{t+h} = \widehat{\beta} x_t}.
#' 
#' @param y a series to be predicted
#' @param x a numeric or matrix of explanatory variables
#' @param h The horizon for which you would like to have the prediction for
#' (see details)
#' @param wind the size of the rolling window or the initial training period if
#' recursive is used
#' @param rr recursive or rolling window? Possible values are
#' \code{c("Rec","Rol")}
#' @return vector of prediction values with the same dimension as the original
#' series. The first \code{wind} values are NA's
#' @examples
#' 
#' x = rnorm(100)
#' lx <- lagmat(x,2)
#' tail(lx)
#' tail(x)
#' out <- linpred(x,lx)
#' plott(x, return.to.default=FALSE)
#' plott(out,add=TRUE,col=2)
#' 
#' @export linpred
linpred <- function(y,x,h=1,wind=NULL,rr=c("Rec")){
  if (rr != "Rol" &  rr !="Rec") { stop("rr is not properly defined") }
  if (NROW(y) != NROW(x) ) { stop("Problem is that NROW(y) != NROW(x)") }
  if(is.null(wind)){wind=0.25*length(y)}
  predd <- NULL
  TT <- length(y)
  x=as.matrix(x) # if x is a numeric make it a matrix
  if(rr=='Rec'){
    for (j in (1+wind):TT){
      lm0 <- lm(y[1:(j-h)]~x[1:(j-h),]) # train the coefficietns, stop before i+wind-h.
      predd[j] <- lm0$coef%*%c(1,x[(j-h),])
    }}
  if(rr=='Rol'){
    for (j in 1:(TT-wind)){
      lm0 <- lm(y[j:(wind+j-h)]~x[j:(wind+j-h),]) # train the coefficietns, stop before i+wind-h.
      predd[(j+wind)] <- lm0$coef%*%c(1,x[(j+wind-h),])
    } }
  predd
}

#' FCIplot
#' 
#' Estimate and plot prediction standard deviation. Given the series, the
#' function estimate point prediction based on AR(1) model and, using the
#' resdiuals from this simple model, estimate an ARCH model to estimate the
#' prediction standard deviation. If \code{plott=TRUE}, a plot of the most
#' recent \code{k} values is created.
#' 
#' Estimate and plot prediction confidence intervals based on AR-ARCH model.
#' 
#' @param series series to be plotted.
#' @param plott should a plot be created? default is \code{plott=TRUE}.
#' @param wind1 window size for the AR component (see details).
#' @param wind2 window size for the ARCH component (see details).
#' @param k if \code{plott=TRUE}, \code{tail(series,k)} will be plotted.
#' @param rrr1 will the AR model be estimated using Recursive ("Rec") or
#' Rolling ("Rol") window?
#' @param rrr2 will the ARCH model be estimated using Recursive ("Rec") or
#' Rolling ("Rol") window?
#' @param main main title of the plot, same as in \code{\link{plot.default}}.
#' @return vector of prediction's standard deviation.
#' @examples
#' 
#' par(mfrow = c(2,1))
#' out <- FCIplot(rnorm(100),plott=TRUE,k=30)
#' plott(out,main="The out-of-sample standard deviation")
#' 
#' @export FCIplot
FCIplot <- function(series,plott=TRUE,wind1=24,wind2=60,k=60,rrr1="Rec",rrr2="Rec",main="series"){
  if (wind2< wind1) { stop("wind2 has to be > than wind1") }
  TT <- length(series)
  temp <- linpred(series[2:TT],series[1:(TT-1)],h=1, wind=wind1,rr=rrr1)
  ertemp <- (series[2:TT]-temp)^2 
  TT2 <- length(ertemp)
  sdestimate <- linpred(ertemp[2:TT2], ertemp[1:(TT2-1)], 
                        wind = wind2, rr = rrr2)
  # We lost 1 observation to construct the AR forecast and 1 to construct the ARCH
  up1 = series[(wind2+2):(TT-1)] + 1*sqrt(abs(na.omit(sdestimate)))
  down1 = series[(wind2+2):(TT-1)] - 1*sqrt(abs(na.omit(sdestimate)))
  up2 = series[(wind2+2):(TT-1)] + 2*sqrt(abs(na.omit(sdestimate)))
  down2 = series[(wind2+2):(TT-1)] - 2*sqrt(abs(na.omit(sdestimate)))
  if(k>length(down2)){k=length(down2) ; warning( "k was set to a lower value" )}
  if(plott){
plot(tail(series,k),ylim=c(min(tail(down2,k)),max(tail(up2,k))),ty="b",pch=19,xlab="",
     main=main,ylab="")
grid()
    polygon(x = c(rev(1:k),(1:k)), y = c(rev(tail(down2,k)),tail(up2,k)),
            col = rgb(0,0,0,0.2),border ='red')
    polygon(x = c(rev(1:k),(1:k)), y = c(rev(tail(down1,k)),tail(up1,k)),
            col = rgb(0,0,0,0.2),border ='blue')
  }
  return(sdestimate^.5)
}



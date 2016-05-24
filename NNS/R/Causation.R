#' Uni-Directional Causation
#'
#' Returns the uni-directional causality from observational data between two variables.  Causation nets out the univariate effect.
#'
#' @param x Variable
#' @param y Variable
#' @param tau Number of lagged observations to consider
#' @keywords causation
#' @author Fred Viole, OVVO Financial Systems
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{Uni.caus(x,y,3)}
#' @export


Uni.caus <- function(x,y,tau){


  min.length = min(length(x),length(y))

  x.vectors = list()
  y.vectors = list()

  ## Normalize x to x.tau
  for (i in 0:tau){

      x.vectors[[paste('x.tau.',i,sep="")]] <- numeric(0L)
  }


  for (i in 0:tau){
      start = tau-i+1
      end = min.length-i
      x.vectors[[i+1]] = x[start:end]
    }

  x.vectors.tau = do.call(cbind.data.frame, x.vectors)

  x.norm.tau <- VN.norm.caus(x.vectors.tau)[,1]


  ## Normalize y to y.tau
  for (i in 0:tau){

    y.vectors[[paste('x.tau.',i,sep="")]] <- numeric(0L)

  }
  for (i in 0:tau){
    start = tau-i+1
    end = min.length-i
    y.vectors[[i+1]] = y[start:end]
  }

  y.vectors.tau = do.call(cbind.data.frame, y.vectors)

  y.norm.tau <- VN.norm.caus(y.vectors.tau)[,1]


  ## Normalize x.norm.tau to y.norm.tau
  x.tau.y.tau = cbind(x.norm.tau,y.norm.tau)
  x.norm.to.y = VN.norm.caus(x.tau.y.tau)[,1]
  y.norm.to.x = VN.norm.caus(x.tau.y.tau)[,2]


  ## Conditional Probability from Normalized Variables P(x.norm.to.y |y.norm.to.x)
  P.x.given.y = UPM(0,min(x.norm.to.y),y.norm.to.x) - UPM(0,max(x.norm.to.y),y.norm.to.x)


  ## Correlation of Normalized Variables
  rho.x.y = VN.cor(x.norm.to.y,y.norm.to.x,
                   degree= ifelse(length(x)<100,0,1))

  #rho.x.y = cor(x.norm.to.y,y.norm.to.x)
  #print(rho.x.y)
  Causation.x.given.y= P.x.given.y*rho.x.y



  par(mfrow=c(3,1))

  ## Raw Variable Plot
  ymin = min(c(min(x),min(y)))
  ymax = max(c(max(x),max(y)))
  par(mar=c(2, 4, 0, 1))
  plot(y,type = 'l', ylim=c(ymin, ymax),ylab='RAW',col='red',lwd = 3)
  lines(x, col = 'steelblue',lwd = 3)
  legend('top',c("X","Y"), lty = 1,lwd=c(3,3),
         col=c('steelblue','red'),ncol=2,bty ="n")



  ## Time Normalized Variables Plot
  ymin = min(c(min(x.norm.tau),min(y.norm.tau)))
  ymax = max(c(max(x.norm.tau),max(y.norm.tau)))
  par(mar=c(2, 4, 0, 1))
  plot(y.norm.tau,type = 'l', ylim=c(ymin, ymax),ylab='TIME NORMALIZED',col='red',lwd = 3)
  lines(x.norm.tau, col='steelblue',lwd=3)
  legend('top',c("X","Y"), lty = 1,lwd=c(3,3),
         col=c('steelblue','red'),ncol=2,bty ="n")

  ## Time Normalized Variables Normalized to each other Plot
  ymin = min(c(min(x.norm.to.y),min(y.norm.to.x)))
  ymax = max(c(max(x.norm.to.y),max(y.norm.to.x)))
  par(mar=c(2, 4, 0, 1))
  plot(y.norm.to.x,type = 'l', ylim=c(ymin, ymax),ylab='X & Y NORMALIZED',col='red',lwd = 3)
  lines(x.norm.to.y,col='steelblue',lwd=3)
  legend('top',c("X","Y"), lty = 1,lwd=c(3,3),
         col=c('steelblue','red'),ncol=2,bty ="n")

  return(Causation.x.given.y)

}

#' VN Causation
#'
#' Returns the causality from observational data between two variables
#'
#' @param x Variable
#' @param y Variable
#' @param tau Number of lagged observations to consider
#' @keywords causation
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{VN.caus(x,y,3)}
#' @export

VN.caus <- function(x,y,tau){

  if(abs(Uni.caus(x,y,tau))<abs(Uni.caus(y,x,tau))){
    return(c(Causation.x.given.y = Uni.caus(x,y,tau),
         Causation.y.given.x = Uni.caus(y,x,tau),
 "C(y--->x)" = Uni.caus(y,x,tau)-Uni.caus(x,y,tau)))}
  else{
    return(c(Causation.x.given.y = Uni.caus(x,y,tau),
             Causation.y.given.x = Uni.caus(y,x,tau),
             "C(x--->y)" = -Uni.caus(y,x,tau)+Uni.caus(x,y,tau)))
  }


}

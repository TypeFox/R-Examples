#' Henon map
#' @description
#' Generates a 2-dimensional time series using the Henon map.
#' @details
#' The Henon map is defined as follows:
#' \deqn{ x_n = 1 - a \cdot x_{n - 1}^2 + y_{n - 1}}{x[n] = 1 - a*x[n - 1]^2 + y[n - 1]}
#' \deqn{ y_n = b \cdot x_{n - 1}}{y[n] = b*x[n - 1].}
#' The default selection for both \emph{a} and \emph{b} parameters (\emph{a}=1.4 and \emph{b}=0.3) is known to
#' produce a deterministic chaotic time series. 
#' @param start A 2-dimensional vector indicating the starting values for the x and y Henon coordinates. 
#' If the starting point is not specified, it is generated randomly.
#' @param a The \emph{a} parameter. Default: 1.4.
#' @param b The \emph{b} parameter. Default: 0.3.
#' @param n.sample Length of the generated time series. Default: 5000 samples.
#' @param n.transient Number of transient samples that will be discarded. Default: 500 samples.
#' @param do.plot Logical value. If TRUE (default value), a plot of the generated Henon system is shown.
#' @return A list with two vectors named \emph{x} and \emph{y} containing the x-components and the 
#' y-components of the Henon map, respectively.
#' @note Some initial values may lead to an unstable system that will tend to infinity.
#' @references Strogatz, S.: Nonlinear dynamics and chaos: with applications to physics, biology, chemistry and engineering (Studies in Nonlinearity)
#' @author Constantino A. Garcia
#' @seealso \code{\link{logisticMap}, \link{lorenz}, \link{rossler}, \link{ikedaMap}, \link{cliffordMap}, \link{sinaiMap}, \link{gaussMap}}
#' @examples
#' \dontrun{
#' henon.map=henon(n.sample = 1000, n.transient=10,do.plot=TRUE,
#'                start=c(-0.006423277,-0.473545134))
#' # accessing the x coordinate and plotting it
#' plot(ts(henon.map$x))
#' }
#' @export henon
henon=function(start=runif(min=-0.5,max=0.5,n=2), a=1.4, b=0.3, n.sample = 5000, n.transient=500, do.plot=TRUE){
  
  ##########################################################################
  n.sample = n.sample + n.transient
  y = x = vector(mode = "numeric", length = n.sample)
  x[[1]] = start[[1]]
  y[[1]] = start[[2]]
  sampleVector = seq(2, n.sample)
  for (i in sampleVector) {
    x[[i]] = y[[i-1]] + 1 -a*x[[i-1]]^2
    y[[i]] = b*x[[i - 1]]
  }
  transientVector = seq(n.transient)
  x = x[-transientVector]
  y = y[-transientVector]
  
  # plotting
  if (do.plot){
    title=paste("Henon map\n","a = ",a," b = ",b)
    plot(x,y,xlab="x[n]",ylab="y[n]",main=title,type="p")
  }
  
  # return values
  return(list(x = x, y = y))
}

################################# logistic map ################################
#' Logistic map
#' @description
#' Generates a time series using the logistic map.
#' @details
#' The logistic map is defined as follows:
#' \deqn{x_n = r  \cdot  x_{n-1}   \cdot  (1 - x_{n-1})}{x[n] = r * x[n-1]  * (1 - x[n-1])}
#' The default selection for the \emph{r} parameter is known to
#' produce a deterministic chaotic time series.
#' @param start A numeric value indicating the starting value for the time series.
#' If the starting point is not specified, it is generated randomly.
#' @param r The \emph{r} parameter. Default: 4
#' @param n.sample Length of the generated time series. Default: 5000 samples.
#' @param n.transient Number of transient samples that will be discarded. Default: 500 samples.
#' @param do.plot Logical value. If TRUE (default value), a plot of the generated logistic system is shown.
#' @return A vector containing the values of the time series that has been generated.
#' @note Some initial values may lead to an unstable system that will tend to infinity.
#' @references Strogatz, S.: Nonlinear dynamics and chaos: with applications to physics, biology, chemistry and engineering (Studies in Nonlinearity)
#' @author Constantino A. Garcia
#' @seealso \code{\link{henon}, \link{lorenz}, \link{rossler}, \link{ikedaMap}, \link{cliffordMap}, \link{sinaiMap}, \link{gaussMap}}
#' @examples
#' \dontrun{
#' log.map=logisticMap(n.sample = 1000, n.transient=10,do.plot=TRUE)
#' }
#' @export logisticMap
logisticMap = function( r=4, start=runif(n=1,min=0,max=1), n.sample=5000,n.transient=500,do.plot=TRUE){
  ############################################################
  n.sample = n.sample + n.transient
  x = vector(mode = "numeric", length = n.sample)
  x[[1]] = start
  sampleVector = seq(2, n.sample)
  for(i in sampleVector){
    x[[i]] = r * x[[i-1]]  * (1 - x[[i-1]])
  }
  # plotting
  if (do.plot){
    title=paste("logistic map\n","r = ",r)
    plot(c(1,sampleVector),x,xlab="n",ylab="x[n]",main=title,type="l")
  }
  # eliminate transients
  transientVector = seq(n.transient)
  x = x[-transientVector]
  # return values
  return(x)
}

################################################################################

bifurcationDiagram=function(mappingFunctionName, parameterVector, n.transient = 500){
  FUN <- match.fun(mappingFunctionName) 
    
  x=FUN(parameterVector[[1]])
  total.length = length(x)
  effective.length = total.length - n.transient + 1
  # use plot the first time
  plot(rep(parameterVector[[1]],effective.length ),
       x[n.transient:total.length],'p',xlim=range(parameterVector),ylim=c(0,1),cex=0.1)
  parameterVector=parameterVector[-1]
  # use lines 
  for (r in parameterVector){
    x=FUN(r)
    lines(rep(r,effective.length ),x[n.transient:total.length],'p',cex=0.1)
  }
  
}

################################################################################
#' Lorenz system
#' @description
#' Generates a 3-dimensional time series using the Lorenz equations.
#' @details
#' The Lorenz system is a system of ordinary differential equations defined as:
#' \deqn{\dot{x} = \sigma(y-x)}{dx/dt = sigma*( y - x )}
#' \deqn{\dot{y} = \rho x-y-xz}{dy/dt = rho*x - y - xz}
#' \deqn{\dot{z} = -\beta z + xy}{dz/dt = -beta*z + xy}
#' The default selection for the system parameters (\eqn{\sigma=10, \rho=28, \beta=8/3}{sigma=10, rho=28, beta=8/3}) is known to
#' produce a deterministic chaotic time series.
#' @param start A 3-dimensional numeric vector indicating the starting point for the time series.
#' Default: c(-13, -14, 47).
#' @param sigma The \eqn{\sigma}{sigma} parameter. Default: 10.
#' @param rho The \eqn{\rho}{rho} parameter. Default: 28.
#' @param beta The \eqn{\beta}{beta} parameter. Default: 8/3.
#' @param time The temporal interval at which the system will be generated. Default: time=seq(0,50,by = 0.01).
#' @param do.plot Logical value. If TRUE (default value), a plot of the generated Lorenz system is shown.
#' @return A list with four vectors named \emph{time}, \emph{x}, \emph{y} and \emph{z} containing the time, the x-components, the 
#' y-components and the z-components of the Lorenz system, respectively.
#' @note Some initial values may lead to an unstable system that will tend to infinity.
#' @references Strogatz, S.: Nonlinear dynamics and chaos: with applications to physics, biology, chemistry and engineering (Studies in Nonlinearity)
#' @author Constantino A. Garcia
#' @seealso \code{\link{henon}, \link{logisticMap}, \link{rossler}, \link{ikedaMap}, \link{cliffordMap}, \link{sinaiMap}, \link{gaussMap}}
#' @examples 
#' \dontrun{
#' lor=lorenz(time=seq(0,30,by = 0.01))
#' # plotting the x-component 
#' plot(lor$time,lor$x,type="l")
#' }
#' @export lorenz
#' @import rgl
lorenz=function(sigma = 10, beta = 8/3, rho = 28, start=c(-13, -14, 47),
                time=seq(0,50,by = 0.01),do.plot=TRUE){
  
  params=c(sigma,beta,rho)
  lorenzEquations=function(coord,t,params){
    x=coord[[1]];y=coord[[2]];z=coord[[3]];sigma=params[[1]];beta=params[[2]];rho=params[[3]]
    return (c(sigma*(y-x), rho*x-y-x*z, x*y - beta*z))
  }
  l=rungeKutta(lorenzEquations,start,time,params)
  
  if (do.plot){
    plot3d(l[,1],l[,2],l[,3],pch=1,cex=1,xlab="x(t)",ylab="y(t)",zlab="z(t)")
  }
  
  
  return(list(time=time,x=l[,1],y=l[,2],z=l[,3]))
}
################################################################################
################################################################################
#' Rossler system
#' @description
#' Generates a 3-dimensional time series using the Rossler equations.
#' @details
#' The Rossler system is a system of ordinary differential equations defined as:
#' \deqn{\dot{x} = -(y + z)}{dx/dt = -(y + z)}
#' \deqn{\dot{y} = x+a \cdot y}{dy/dt = x + a*y}
#' \deqn{\dot{z} = b + z*(x-w)}{dz/dt = b + z*(x-w)}
#' The default selection for the system parameters (\emph{a} = 0.2, \emph{b} = 0.2, \emph{w} = 5.7) is known to
#' produce a deterministic chaotic time series.
#' @param start A 3-dimensional numeric vector indicating the starting point for the time series.
#' Default: c(-2, -10, 0.2).
#' @param a The \emph{a} parameter. Default:0.2.
#' @param b The \emph{b} parameter. Default: 0.2.
#' @param w The \emph{w} parameter. Default: 5.7.
#' @param time The temporal interval at which the system will be generated. Default: time=seq(0,50,by = 0.01).
#' @param do.plot Logical value. If TRUE (default value), a plot of the generated Lorenz system is shown.
#' @return A list with four vectors named \emph{time}, \emph{x}, \emph{y} and \emph{z} containing the time, the x-components, the 
#' y-components and the z-components of the Rossler system, respectively.
#' @note Some initial values may lead to an unstable system that will tend to infinity.
#' @references Strogatz, S.: Nonlinear dynamics and chaos: with applications to physics, biology, chemistry and engineering (Studies in Nonlinearity)
#' @author Constantino A. Garcia
#' @seealso \code{\link{henon}, \link{logisticMap}, \link{rossler}, \link{ikedaMap}, \link{cliffordMap}, \link{sinaiMap}, \link{gaussMap}}
#' @examples 
#' \dontrun{
#' r.ts = rossler(time=seq(0,30,by = 0.01))
#' }
#' @export rossler
#' @import rgl
rossler=function(a = 0.2, b = 0.2, w = 5.7, start=c(-2, -10, 0.2)
                ,time=seq(0,50,by = 0.01), do.plot=TRUE){
  
  params=c(a,b,w)
  rosslerEquations=function(coord,t,params){
    x=coord[[1]];y=coord[[2]];z=coord[[3]];a=params[[1]];b=params[[2]];w=params[[3]]
    return (c(-y-z, x+a*y, b + z*(x-w)))
  }
  r=rungeKutta(rosslerEquations,start,time,params)
  
  if (do.plot){
    plot3d(r[,1],r[,2],r[,3],col='black',pch=1,cex=1,xlab="x(t)",ylab="y(t)",zlab="z(t)")
  }
  
  
  return(list(time=time,x=r[,1],y=r[,2],z=r[,3]))
}
################################################################################
################################# ikeda map ################################
#' Ikeda map
#' @description
#' Generates a time series using the Ikeda map
#' @details
#' The Ikeda map is defined as follows:
#' \deqn{z_{n+1} = a + b \cdot z_n \cdot exp( ik-\frac{ic}{( 1+ |z_{n-1}|^2  )} )}{z[n+1] = a + b*z[n]*exp(ik-ic/( 1+ |z_[n-1]|^2  ) )}
#' The default selection for the \emph{a}, \emph{b}, \emph{c} and \emph{k} parameters is known to
#' produce a deterministic chaotic time series.
#' @param start a 2-dimensional numeric vector indicating the starting value for the time series.
#' If the starting point is not specified, it is generated randomly.
#' @param a The \emph{a} parameter. Default: 0.85.
#' @param b The \emph{b} parameter. Default: 0.9.
#' @param cc The \emph{c} parameter. Default: 7.7.
#' @param k The \emph{k} parameter. Default: 0.4.
#' @param n.sample Length of the generated time series. Default: 5000 samples.
#' @param n.transient Number of transient samples that will be discarded. Default: 500 samples.
#' @param do.plot Logical value. If TRUE (default value), a plot of the generated ikeda system is shown.
#' @return a list with 2 vectors named \emph{x} and \emph{y} the x-components and the 
#' y-components of the Ikeda map, respectively.
#' @note Some initial values may lead to an unstable system that will tend to infinity.
#' @references Strogatz, S.: Nonlinear dynamics and chaos: with applications to physics, biology, chemistry and engineering (Studies in Nonlinearity)
#' @author Constantino A. Garcia
#' @seealso \code{\link{henon}, \link{logisticMap}, \link{lorenz}, \link{rossler}, \link{cliffordMap}, \link{sinaiMap}, \link{gaussMap}}
#' @examples
#' \dontrun{
#' ikeda.map=ikedaMap(n.sample = 1000, n.transient=10, do.plot=TRUE)
#' }
#' @export ikedaMap

# \code{\link{http://www.pessoal.utfpr.edu.br/msergio/cap3/6Ikeda/index.html}}
ikedaMap=function(a = 0.85, b = 0.9, cc= 7.7, k = 0.4, start=runif(2), n.sample = 5000, n.transient=500, do.plot=TRUE){
  n.sample = n.sample + n.transient
  z = vector(mode = "complex", length = n.sample)
  z[[1]] =  complex(real=start[[1]],imaginary=start[[2]]) 
  sampleVector = seq(2, n.sample)
  a=as.complex(a)
  b=as.complex(b)
  cc=as.complex(cc)
  k = as.complex(k)
  for (n in sampleVector) {
    z[[n]] = a + b*z[[n-1]]*exp(1i*( k-cc/( abs(z[[n-1]])^2 +1 ) ))
  }
  
  z = z[-(1:n.transient)]
  x = Re(z); y = Im(z)
  # plotting
  if (do.plot){
    title=paste("Ikeda map\n","a = ",Re(a)," b = ",Re(b), " = ",Re(cc), " k = ",Re(k),"\n")
    plot(x,y,xlab="Re(z[n])",ylab="Im(z[n])",cex=0.3,main=title,type="p")
  }
  
  # return values
  return(list(x=x,y=y))
}
################################################################################
#' Clifford map
#' @description
#' Generates a 2-dimensional time series using the Clifford map.
#' @details
#' The Clifford map is defined as follows:
#' \deqn{x_{n+1} = sin(a \cdot y_n) + c \cdot cos(a \cdot x_n)}{x[n+1] = sin(a*y[n]) + c*cos(a*x[n])}
#' \deqn{y_{n+1} = sin(b \cdot x_n) + d \cdot cos(b \cdot y_n)}{y[n+1] = sin(b*x[n] + d*cos(b*y[n])}
#' The default selection for the \emph{a} \emph{b} \emph{c} and \emph{d} parameters is known to
#' produce a deterministic chaotic time series.
#' @param start a 2-dimensional vector indicating the starting values for the x and y Clifford coordinates. 
#' If the starting point is not specified, it is generated randomly.
#' @param a The \emph{a} parameter. Default: -1.4
#' @param b The \emph{b} parameter. Default: 1.6
#' @param cc The \emph{c} parameter. Default: 1.0
#' @param d The \emph{d} parameter. Default: 0.7
#' @param n.sample Length of the generated time series. Default: 5000 samples.
#' @param n.transient Number of transient samples that will be discarded. Default: 500 samples.
#' @param do.plot Logical value. If TRUE (default value), a plot of the generated Clifford system is shown.
#' @return A list with two vectors named x and y containing the x-components and the 
#' y-components of the Clifford map, respectively.
#' @note Some initial values may lead to an unstable system that will tend to infinity.
#' @author Constantino A. Garcia
#'  @seealso \code{\link{henon}, \link{logisticMap}, \link{lorenz}, \link{rossler}, \link{ikedaMap}, \link{sinaiMap}, \link{gaussMap}}
#' @examples
#' \dontrun{
#' clifford.map=cliffordMap(n.sample = 1000, n.transient=10,do.plot=TRUE)
#' # accessing the x coordinate and plotting it
#' plot(ts(clifford.map$x))}
#' @export cliffordMap 
cliffordMap=function(a = -1.4, b = 1.6, cc = 1.0, d = 0.7,start=runif(2),
                      n.sample = 5000, n.transient=500, do.plot=TRUE){
  n.sample = n.sample + n.transient
  x = y = vector(mode = "numeric", length = n.sample)
  x[[1]] = start[[1]]
  y[[1]] = start[[2]]
  sampleVector = seq(2, n.sample)
  for (n in sampleVector) {
    x[[n]] = sin(a * y[[n-1]]) + cc * cos(a * x[[n-1]])
    y[[n]] = sin(b * x[[n-1]]) + d * cos(b * y[[n-1]])
  }
  transient=1:n.transient
  x=x[-(transient)]
  y=y[-(transient)]
  # plotting
  if (do.plot){
    title=paste("clifford map\n","a = ",a," b = ",b, " c = ", cc, " d = ",d,"\n")
    plot(x,y,xlab="x[n]",ylab="y[n]",cex=0.3,main=title,type="p")
  }
  
  # return values
  return(list(x=x,y=y))
}

#######################   sinai map ########################################
#' Sinai map
#' @description
#' Generates a 2-dimensional time series using the Sinai map.
#' @details
#' The Sinai map is defined as follows:
#' \deqn{x_{n+1} = (x_{n} + y_{n} + a \cdot cos(2 \cdot pi \cdot y_{n}) )mod 1}{x[n+1] = (x[n] + y[n] + a*cos(2*pi*y[n]) )mod 1}
#' \deqn{y_{n+1} = (x_{n} + 2 \cdot y_n)mod 1}{y[n+1] = (x[n] + 2*y[n])mod 1}
#' The default selection for the \emph{a} parameter is known to produce a deterministic chaotic time series.
#' @param start A 2-dimensional vector indicating the starting values for the x and y Sinai coordinates. 
#' If the starting point is not specified, it is generated randomly.
#' @param a The \emph{a} parameter. Default: 0.1
#' @param n.sample Length of the generated time series. Default: 5000 samples.
#' @param n.transient Number of transient samples that will be discarded. Default: 500 samples.
#' @param do.plot Logical value. If TRUE (default value), a plot of the generated Sinai system is shown.
#' @return A list with two vectors named x and y containing the x-components and the 
#' y-components of the Sinai map, respectively.
#' @note Some initial values may lead to an unstable system that will tend to infinity.
#' @author Constantino A. Garcia
#' @seealso \code{\link{henon}, \link{logisticMap}, \link{lorenz}, \link{rossler}, \link{ikedaMap}, \link{cliffordMap}, \link{gaussMap}}
#' @references
#' Mcsharry, P. E. and P. R. Ruffino (2003). Asymptotic angular stability 
#' in nonlinear systems: rotation numbers and winding numbers. Dynamical Systems 18(3), 191-200.
#' @examples
#' \dontrun{
#' sinai.map = sinaiMap(n.sample = 1000, n.transient=10,do.plot=TRUE)
#' # accessing the x coordinate and plotting it
#' plot(ts(sinai.map$x))
#' }
#' @export sinaiMap
sinaiMap=function(a = 0.1,start=runif(2),
                   n.sample = 5000, n.transient=500, do.plot=TRUE){
  n.sample = n.sample + n.transient
  x = y = vector(mode = "numeric", length = n.sample)
  x[[1]] = start[[1]]
  y[[1]] = start[[2]]
  sampleVector = seq(2, n.sample)
  for (n in sampleVector) {
    x[[n]] = (x[[n-1]] + y[[n-1]] + a*cos(2*pi*y[[n-1]]) )%%1 #x mod 1 = fractional part of x
    y[[n]] = (x[[n-1]] + 2*y[[n-1]])%%1
  }
  transient=1:n.transient
  x=x[-(transient)]
  y=y[-(transient)]
  # plotting
  if (do.plot){
    title=paste("Sinai map\n","a = ",a,"\n")
    plot(x,y,xlab="x[n]",ylab="y[n]",cex=0.3,main=title,type="p")
  }
  
  # return values
  return(list(x=x,y=y))
}

#######################   Gauss map ########################################
#' Gauss map
#' @description
#' Generates a 1-dimensional time series using the Gauss map
#' @details
#' The Gauss map is defined as follows:
#' \deqn{x_{n+1}= exp(-a \cdot (x_n)^2) + b}{x[n+1] = exp(-a*(x[n])^2) + b}
#' The default selection for both \emph{a} and \emph{b} parameters is known to produce a deterministic chaotic time series.
#' @param start A numeric value indicating the starting value for the time series.
#' If the starting point is not specified, it is generated randomly.
#' @param a The \emph{a} parameter. Default: 4.9
#' @param b The \emph{b} parameter. Default: -0.58
#' @param n.sample Length of the generated time series. Default: 5000 samples.
#' @param n.transient Number of transient samples that will be discarded. Default: 500 samples.
#' @param do.plot Logical value. If TRUE (default value), a plot of the generated Gauss system is shown.
#' @return A vector containing the values of the time series that has been generated.
#' @note Some initial values may lead to an unstable system that will tend to infinity.
#' @author Constantino A. Garcia
#'  @seealso \code{\link{henon}, \link{logisticMap}, \link{lorenz}, \link{rossler}, \link{ikedaMap}, \link{cliffordMap}, \link{sinaiMap}}
#' @references
#'  Chaos and nonlinear dynamics: an introduction for scientists and engineers, by Robert C. Hilborn, 2nd Ed., Oxford, Univ. Press, New York, 2004.
#' @export gaussMap
#http://en.wikipedia.org/wiki/Gauss_iterated_map
gaussMap=function(a = 4.9, b = -0.58, start=runif(1,min=-0.5,max=0.5),
                        n.sample = 5000, n.transient=500, do.plot=TRUE){
  n.sample = n.sample + n.transient
  x = vector(mode = "numeric", length = n.sample)
  x[[1]] = start[[1]]
  sampleVector = seq(2, n.sample)
  for (n in sampleVector) {
    x[[n]] = exp(-a*(x[[n-1]])^2) + b
  }
  transient=1:n.transient
  x=x[-(transient)]
  # plotting
  if (do.plot){
    title=paste("Gauss map\n","a = ",a," b = ",b,"\n")
    length.series = length(x)
    plot(1:length.series,x,xlab="n",ylab="x[n]",cex=0.3,main=title,type="p")
  }
  # return values
  return(x)
}

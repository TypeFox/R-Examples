# MS or MM linear scaling, return transformation
simple.scale = function(sp, np, mm=FALSE) {
  A = if (mm) mean(np[,1])/mean(sp[,1]) else sd(sp[,2])/sd(np[,2])
  B = mean(sp[,2]) - A*mean(np[,2])
  return(list(A=A,B=B))
}

# function optimised in Lord-Stocking scaling
sl = function (x, sp, np, qp, qw) {
  A = x[1]
  B = x[2]
  np[,1] = np[,1]/A
  np[,2] = np[,2]*A + B
  dif = trf(ip=sp,x=qp)$f - trf(ip=np,x=qp)$f
  return(sum(dif*dif*qw))
}

# function optimised in Haebara scaling
hb = function (x, sp, np, qp, qw) {
  A = x[1]
  B = x[2]
  np[,1] = np[,1]/A
  np[,2] = np[,2]*A + B
  dif = irf(ip=sp,x=qp)$f - irf(ip=np,x=qp)$f
  return(sum(dif*dif*qw))
}

# versions with correction for back-equating
sl2 <- function (x, sp, np, qp, qw){
    A21 = x[1]
    K21 = x[2]
    A12 = 1/A21
    K12 = -K21/A21
    s = length(qp)/2
    np21 = np
    sp12 = sp
    np21[, 1] = np21[, 1]/A21
    np21[, 2] = np21[, 2]*A21 + K21
    Q1 <- trf(ip = sp, x = qp[1:s])$f -
          trf(ip = np21, x = qp[1:s])$f
    sp12[, 1] = sp12[, 1]/A12
    sp12[, 2] = sp12[, 2]*A12 + K12
    Q2 <- trf(ip = sp12, x = qp[(s+1):(2*s)])$f -
          trf(ip = np, x = qp[(s+1):(2*s)])$f
    dif <- c(Q1, Q2)
    return(sum(dif * dif * qw))
}

hb2 <- function (x, sp, np, qp, qw){
    A21 = x[1]
    K21 = x[2]
    A12 = 1/A21
    K12 = -K21/A21
    s = length(qp)/2
    np21 = np
    sp12 = sp
    np21[, 1] = np21[, 1]/A21
    np21[, 2] = np21[, 2]*A21 + K21
    Q1 <- irf(ip = sp, x = qp[1:s])$f -
          irf(ip = np21, x = qp[1:s])$f
    sp12[, 1] = sp12[, 1]/A12
    sp12[, 2] = sp12[, 2]*A12 + K12
    Q2 <- irf(ip = sp12, x = qp[(s+1):(2*s)])$f -
          irf(ip = np, x = qp[(s+1):(2*s)])$f
    dif <- rbind(Q1, Q2)
    return(sum(dif * dif * qw))
}

# do Lord-Stocking or Haebara scaling, return transformation
adv.scale = function(sp,np,sq=NULL,nq=NULL,haeb=FALSE,bec=FALSE) {
  if (is.null(sq)) stop("no quadrature for characteristic curve method")
  if (is.null(nq) && haeb) stop("Haebara method needs both old and new quadrature")
  qp = if (haeb) c(sq$quad.points,  nq$quad.points)  else sq$quad.points
  qw = if (haeb) c(sq$quad.weights, nq$quad.weights) else sq$quad.weights
  if (bec) {
    r  = if (haeb) optim(c(1,0),hb2,method="BFGS",sp=sp,np=np,qp=qp,qw=qw) else
      optim(c(1,0),sl2,method="BFGS",sp=sp,np=np,qp=qp,qw=qw)    
  } else {
    r  = if (haeb) optim(c(1,0),hb,method="BFGS",sp=sp,np=np,qp=qp,qw=qw) else
      optim(c(1,0),sl,method="BFGS",sp=sp,np=np,qp=qp,qw=qw)    
  }	  
  return(list(A=r$par[1],B=r$par[2])) 
}


#' Linear transformation of the IRT scale
#' 
#' Linearly transform a set of IRT parameters to bring them to the scale of
#' another set of parameters. Four methods are implemented: Mean/Mean,
#' Mean/Sigma, Lord-Stocking, and Haebara.
#' 
#' 
#' @param old.ip A set of parameters that are already on the desired scale
#' @param new.ip A set of parameters that must be placed on the same scale as
#' \code{old.ip}
#' @param old.items A vector of indexes pointing to those items in
#' \code{old.ip} that are common to both sets of parameters
#' @param new.items The indexes of the same items in \code{new.ip}
#' @param old.qu A quadrature object for \code{old.ip}, typically produced by
#' the same program that estimated \code{old.ip}. Only needed if
#' \code{method="LS"} or \code{method="HB"}
#' @param new.qu A quadrature object for \code{new.ip}, typically produced by
#' the same program that estimated \code{new.ip}. Only needed if
#' \code{method="HB"}
#' @param method One of "MM" (Mean/Mean), "MS" (Mean/Sigma), "SL"
#' (Stocking-Lord), or "HB" (Haebara). Default is "MS"
#' @param bec Use back-equating correction? When TRUE, the Stocking-Lord or
#' Hebaera procedures will be adjusted for back-equating (see Hebaera, 1980).
#' Ignored when method is MM or MS. Default is FALSE.
#' @return A list of: \item{slope}{The slope of the linear transformation}
#' \item{intercept}{The intercept of the linear transformation}
#' \item{scaled.ip}{The parameters in \code{new.ip} tranformed to a scale that
#' is compatible with \code{old.ip}}
#' @author Ivailo Partchev and Tamaki Hattori
#' @references Kolen, M.J. & R.L. Brennan (1995) Test Equating: Methods and
#' Practices. Springer.
#'
#' Haebara, T. (1980) Equating logistic ability scales by
#' a weighted lest squares method. Japanese Psychological Research, 22,
#' p.144--149
#' @keywords models
#' @export
#' @examples
#' 
#' \dontrun{
#' # a small simulation to demonstrate transformation to a common scale
#' # fake 50 2PL items
#' pa <- cbind(runif(50,.8,2), runif(50,-2.4,2.4), rep(0,50))
#' # simulate responses with  two samples of different ability levels
#' r.1 <- sim(ip=pa[1:30,],  x=rnorm(1000,-.5))
#' r.2 <- sim(ip=pa[21:50,], x=rnorm(1000,.5))
#' # estimate item parameters
#' p.1 <- est(r.1, engine="ltm")
#' p.2 <- est(r.2, engine="ltm")
#' # plot difficulties to show difference in scale
#' plot(c(-3,3), c(-3,3), ty="n", xlab="True",ylab="Estimated",
#'    main="Achieving common scale")
#' text(pa[1:30,2],  p.1$est[,2], 1:30)
#' text(pa[21:50,2], p.2$est[,2], 21:50, co=2)
#' # scale with the default Mean/Sigma method
#' pa.sc = sca(old.ip=p.1$est, new.ip=p.2$est, old.items=21:30, new.items=1:10)
#' # add difficulties of scaled items to plot
#' text(pa[21:50,2], pa.sc$scaled.ip[,2], 21:50, co=3)
#' }
#' 
sca = function(old.ip, new.ip, old.items, new.items,
  old.qu=NULL, new.qu=NULL, method="MS", bec=FALSE) {
  if (length(old.items)  != length(new.items)) stop("no of common items does not match")
  if (!all(old.items %in% 1:nrow(old.ip))) stop("bad index for some scaled item")
  if (!all(new.items %in% 1:nrow(new.ip))) stop("bad index for some new item")
  sp = old.ip[old.items, ]
  np = new.ip[new.items, ]
  r = switch(method,
    "MS"= simple.scale(sp,np,mm=FALSE),
    "MM"= simple.scale(sp,np,mm=TRUE),
    "HB"= adv.scale(sp,np,old.qu,new.qu,haeb=TRUE,bec=bec),
    "SL"= adv.scale(sp,np,old.qu,haeb=FALSE,bec=bec),
    stop(paste("unknown scaling method",method))
  )
  new.ip[,1] = new.ip[,1] / r$A
  new.ip[,2] = new.ip[,2] * r$A + r$B   
  return(list(slope=r$A, intercept=r$B, scaled.ip=new.ip))
}

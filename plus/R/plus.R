
########################################################
#####   the R code for plus  v. 1.4
#####   June 2009
########################################################

########################################################
# plus(), the main plus function; documented 
###
plus <- function(x,y, method = c("lasso", "mc+", "scad", "general"), m=2, gamma,v,t, 
  monitor=FALSE, normalize = TRUE, intercept = TRUE, 
  Gram, use.Gram = FALSE, eps=1e-15, max.steps=500, lam)
{ 
 ## fill missing values
 if (missing(gamma)) gamma <- 0
 if (missing(v)) v <- 0
 if (missing(t)) t <- 0
 if (missing(Gram)) Gram <- 0
 if (missing(lam)) lam <- -1
 lam.min <- min(lam)
 ## check data
 data.dim.error <- FALSE
 if (! is.matrix(x)) data.dim.error <- TRUE
 else {
   n <- dim(x)[1]
   if ( (n<=1) | (n!=length(y)) ) data.dim.error <- TRUE
   p <- dim(x)[2]
   if (p<= 1) data.dim.error <- TRUE
 }
 if (data.dim.error == TRUE){ 
   cat("Data dimention error.", "\n")
   return(list(data.dim.error = TRUE))
 }
 ## standardize x 
 orig.x <- x
 if (intercept)  x <- t(t(x)-apply(x,2,mean)) # center to mean zero
 if (normalize) {
   normalize.factor <- sqrt(apply(x^2,2,mean))
   x <- t(t(x)/normalize.factor)
 }
 ## compute Gram matrix if necessary
 ## if ((!use.Gram) & (p<=200)) use.Gram <- TRUE
 if ((use.Gram) & (p > 1000)) use.Gram <- FALSE
 if (use.Gram)
   if (! is.matrix(Gram)) Gram <- t(x)%*%x
   else if ((dim(Gram)[1]!=p) || (dim(Gram)[2]!=p)) Gram <- t(x)%*%x
 ## determine method and compute penalty        
 if (length(method)==1)
   if (method == "lasso") m <- 1
   else if (method == "mc+") m <- 2
   else if (method == "scad") m <- 3
 if (m==1) method <- "LASSO"
 else if (m==2) method <- "MC+"
 else if (m==3) method <- "SCAD"
 else method <- "PLUS" 
 if (monitor == TRUE) {
   cat("\n")
   cat("Sequence of", method, "moves:", "\n","\n")
 }
 if (m<4) { # provide default penalty for lasso, MCP and SCAD 
   if ((m>1) & (use.Gram) & (gamma == 0)) {
     max.corr <- max(abs(Gram[ row(Gram) != col(Gram) ]))/n
     if ((1 - max.corr) >= (1e-5) )   gamma <- 2/(1- max.corr)
   }
   if ((m==3) & (gamma <=1)) gamma <- 0  
   tmp <- plus.penalty(m,a=gamma)
   m <- tmp$m
   gamma <- tmp$gamma
   v <- tmp$v
   t <- tmp$t
 }
 ## declare variables 
 names <- dimnames(x)
 b <- matrix(0, max.steps+1,p)
 etatil <- matrix(0,2,p+1)
 tau <- rep(0, max.steps+1)
 my.env <-environment()
 my.env$x <-x
 my.env$y <- y
 my.env$tx <- t(x)/n
 my.env$eta <- b
 my.env$g.prime <- rep(0,p)
 my.env$grad <- rep(0,p)
 my.env$etatil2 <- etatil
 my.env$etatil2[1,] <- 1:(p+1)
 if (use.Gram) my.env$Sigma <- Gram/n
 else { 
   my.env$Sigma <- matrix(0,p,100)  # e$a.set.var will match a.set
   my.env$var.list <- 0
 }
 my.env$use.Gram <- use.Gram
 my.env$ties <- rep(FALSE,p)
 my.env$t.fun <- plus.knots(t,m)
 z <- my.env$tx %*%y 
 my.env$tau1 <- 1/max(abs(z))
 ## initialization: for the initial step k=0, tau^{(0)} is needed, b[1,] <- 0 done already
 tau[1] <- my.env$tau1   
 etatil[1,] <- 1:(p+1)   # the first etatil
 etatil[2,1:p] <- sign(z)
 etatil <- etatil[ , c(abs(z) >= max(abs(z)) - eps,TRUE) ]
 k <-1 
 exit.while <- FALSE
 last.valid.s <- rep(0,p)
 ## the main loop
 while (exit.while==FALSE) { # segment k: model eta[k+1,] begins with b[k,] and ends with b[k+1,]
   tmp <- plus.single.step(k, etatil, z, b[k,], tau[k], m, v, t, eps=eps,e=my.env)
   ## tmp holds return(new.eta, etatil, s, tau, b, singular.Q, forced.stop, full.path) 
   etatil <- tmp$etatil     # etatil = 0 if forced.stop, length(etatil)=2 if Ise.termination
   exit.while <- tmp$forced.stop # cannot find a valid new.eta
   full.path <- tmp$full.path 
   if (exit.while == FALSE) { # save data from tmp if not forced to stop 
     eta[k+1,] <- tmp$new.eta
     if (monitor == TRUE) { # print info to monitor the iterations
       var.change <- sign(abs(sign(eta[k+1,]))-abs(sign(eta[k,])))
       if (sum(abs(var.change))==1) {
         if(sum(var.change)==1)
           action <-"added "
         else
           action <-"removed " 
           var_name <- (order(-abs(var.change))[1])
         if (!is.null(names))
      var_name <- names[[2]][var_name]      
         cat("Step",k,": ",action,var_name,"\n")
       }  
     }
     last.valid.s <- tmp$s  # save the s for the last valid eta 
     tau[k+1]<- tmp$tau
     b[k+1,] <- tmp$b  # b[k+1] is at the end of segment eta[k+1,] 
     b[k+1,eta[k+1,]==0] <- 0 # more accurate b=0, b[k,eta[k+1,]==0] <- 0 is also ok
     k <- k+1
     if ((k > max.steps)|(full.path)|(tau[k]*lam.min > 1)) { # stop with a valid eta
       exit.while <- TRUE
       if ((full.path) & (m>1)) k <- k-1
     } 
   } # if (exit.while == FALSE) 
 } # end the main loop
 ## tmp holds return(new.eta, etatil, s, tau, b, singular.Q, forced.stop, full.path) 
 singular.Q <- tmp$singular.Q
 forced.stop <- tmp$forced.stop
 lam.path <- 1/tau[1:k]
 beta.path <- lam.path * b[1:k,]
 if (normalize) {
    beta.path <- t(t(beta.path)/normalize.factor)
 }
 lam.path[k] <- max(abs(my.env$tx %*%(y-x%*%beta.path[k,])))
 if (full.path) lam.path[k] <- 0
 eta <- eta[1:k,]
 total.steps <- k-1
 ## calculate the estimator beta at the given lam or the set of values of lam.path
 x <- orig.x
 if (sum(lam) == -1) lam <- sort(lam.path, decreasing=TRUE)
 tmp <- plus.hit.points(lam,lam.path) 
 total.hits <- tmp$total.hits
 if (total.hits > 0) {
   beta <- tmp$w1*beta.path[tmp$k1,] + (1-tmp$w1)*beta.path[tmp$k2,]
   lam <- lam[1:total.hits]
   dim <- apply(abs(beta)>0, 1, sum)
   if (total.hits == 1)
     r.square <- 1- sum((y - x%*%beta)^2)/sum(y^2)
   else
     r.square <- 1- apply((x%*%t(beta) - y)^2, 2, sum)/sum(y^2)   
     obj<-list(x=x,y=y,eta=eta, beta.path=beta.path, lam.path=lam.path, 
    beta=beta, lam=lam, dim=dim, r.square = r.square, total.hits=total.hits, 
  method = method, gamma = gamma, total.steps=total.steps, max.steps=max.steps, 
  full.path=full.path, forced.stop=forced.stop, singular.Q=singular.Q)
 }
 else 
   obj<-list(x=x,y=y,eta=eta, beta.path=beta.path, lam.path=lam.path, total.hits=total.hits, 
  method = method, gamma = gamma, total.steps=total.steps, max.steps=max.steps, 
  full.path=full.path, forced.stop=forced.stop, singular.Q=singular.Q)  
 class(obj) <- "plus"
 return(obj)  
}

########################################################
# print.plus(), to print plus() moves; documented 
###
print.plus <- function(x, print.moves = 20, ...) { 
  cat("\n")
  cat("Sequence of", x$method, "moves:", "\n","\n")
  print.moves <- round(print.moves)
  if (print.moves < 1) print.moves <- 20
  names <- dimnames(x$x)
  m <- 1
  k <- 1
  exit.while <- FALSE
  while (exit.while == FALSE) {
    var.change <- sign(abs(sign(x$eta[k+1,]))-abs(sign(x$eta[k,])))
    if (sum(abs(var.change))==1) {
      if(sum(var.change)==1)
        action <-"added "
      else
        action <-"removed "   
      var_name <- (order(-abs(var.change))[1])
      if (!is.null(names))
        var_name <- names[[2]][var_name]
      cat("Step",k,": ",action,var_name,"\n")
      m <- m+1
      if (m>print.moves) exit.while <- TRUE
    }
    k <- k+1
    if (k > x$total.steps) exit.while <- TRUE
  } # end while
} 

########################################################
# predict.plus(), to extract coefficients and preditions; documented 
###
predict.plus <- function(object,lam,newx, ...) { 
  if (missing(newx)) {
    flag <- FALSE 
    cat("Warning message: ")
    cat("no newx argument; newy not produced", "\n") 
  } 
  else flag <- TRUE
  if (missing(lam)) {
    lam <- sort(object$lam.path, decreasing=TRUE)
    cat("Warning message: ")
    cat("no lam argument; object$lam.path used", "\n") 
  }
  else if (max(lam) < min(object$lam.path)) {
    lam <- sort(object$lam.path, decreasing=TRUE)
    cat("Warning message: ")
    cat("lam not reached by the plus path; object$lam.path used", "\n") 
  }
  cat("\n")
  tmp <- plus.hit.points(lam,object$lam.path)
  total.hits <- tmp$total.hits
  beta <- tmp$w1* object$beta.path[tmp$k1,] + (1-tmp$w1)* object$beta.path[tmp$k2,]
  lam <- lam[1:total.hits]
  dim <- apply(abs(beta)>0, 1, sum)
  if (total.hits == 1) {
    r.square <- 1- sum((object$y - object$x%*%beta)^2)/sum(object$y^2)
    if (flag) {
      if (is.matrix(newx)) newy <- newx%*% beta
      else newy <- sum(beta * newx)
    }
  }
  else {
    r.square <- 1- apply((object$x%*%t(beta) - object$y)^2, 2, sum)/sum(object$y^2)
    if (flag) {
      if (is.matrix(newx)) newy <- beta %*% t(newx)
      else newy <- beta %*% newx
    }
  }
  if (flag) 
    return(list(lambda=lam, coefficients=beta, dimension=dim, r.square = r.square, step = tmp$k2-1,  
  method = object$method, newy = newy))
  else
    return(list(lambda=lam, coefficients=beta, dimension=dim, r.square = r.square, step = tmp$k2-1,  
  method = object$method))
}

########################################################
# coef.plus(), to extract coefficients; documented 
###
coef.plus <- function(object,lam, ...){ 
  if (missing(lam)) {
    lam <- sort(object$lam.path, decreasing=TRUE)
    cat("Warning message: ")
    cat("no lam argument; object$lam.path used", "\n") 
  }
  else if (max(lam) < min(object$lam.path)) {
    lam <- sort(object$lam.path, decreasing=TRUE)
    cat("Warning message: ")
    cat("lam not reached by the plus path; x$lam.path used", "\n") 
  }
  cat("\n")
  tmp <- plus.hit.points(lam,object$lam.path)
  beta <- tmp$w1* object$beta.path[tmp$k1,] + (1-tmp$w1)* object$beta.path[tmp$k2,]
  return(beta)
}

########################################################
# plot.plus(), to plot coefficients, preditions, penalty level, dimension or R-square 
# againt plus step or penalty level; documented 
###
plot.plus <- function(x, xvar=c("lam","step"), yvar=c("coef","newy","lam","dim","R-sq"), 
        newx, step.interval, lam.interval, predictors, ...) {
  # determin xvar in the plot
  if (length(unique(c(xvar,"lam")))==length(unique(c(xvar)))) xvar <- "lambda"
  else if (length(unique(c(xvar,"step")))==length(unique(c(xvar)))) xvar <- "step"
  else xvar <- "lambda"
  # determin yvar in the plot
  if (length(unique(c(yvar,"coef")))==length(unique(c(yvar)))) yvar <- "coefficients"
  else if (length(unique(c(yvar,"newy")))==length(unique(c(yvar)))) yvar <- "newy"
  else if (length(unique(c(yvar,"lam")))==length(unique(c(yvar)))) yvar <- "lambda"
  else if (length(unique(c(yvar,"dim")))==length(unique(c(yvar)))) yvar <- "dimension"
  else if (length(unique(c(yvar,"R-sq")))==length(unique(c(yvar)))) yvar <- "R-square"
  else yvar <- "coefficients"
  if (yvar != "newy") newx <- matrix(0,2,dim(x$x)[2])
  if ((yvar == "newy") & missing(newx)) {
    cat("Warning message: ")
    cat("no newx argument", "\n") 
    return()
  }
  if (missing(predictors)) predictors <- 1:(dim(x$x)[2])
  ## plot 
  if (xvar == "step") {
    xtmp <-  0:(x$total.steps)
    plot.set <- rep(TRUE, length(xtmp))
    if (! missing(step.interval)) {
      plot.set[xtmp > max(step.interval)] <- FALSE
      plot.set[xtmp < min(step.interval)-1] <- FALSE
      xtmp <- xtmp[plot.set]
    }
    if (yvar == "coefficients") ytmp <- x$beta.path[plot.set, ]
    if ((yvar == "newy") & is.matrix(newx)) ytmp <- x$beta.path[plot.set, ]%*%t(newx)
    if ((yvar == "newy") & (! is.matrix(newx))) ytmp <- x$beta.path[plot.set, ]%*%newx
    if (yvar == "lambda") ytmp <- x$lam.path[plot.set]
    if (yvar == "dimension") ytmp <- apply(abs(sign(x$eta[plot.set,])),1, sum)
    if (yvar == "R-square") 
      ytmp <- 1- apply((x$x%*%t(x$beta.path[plot.set, ]) - x$y)^2, 2, sum)/sum(x$y^2)
  }    
  if (xvar == "lambda") {
    xtmp <- sort(x$lam.path,decreasing=TRUE)
    tmp <- predict(x,xtmp,newx)
    plot.set <- rep(TRUE, length(xtmp))
    if (! missing(lam.interval)) {
      plot.set[xtmp > max(lam.interval)] <- FALSE
      plot.set[xtmp < min(lam.interval)] <- FALSE
      xtmp <- - xtmp[plot.set]
    }
    else xtmp <- -xtmp
    if (yvar == "coefficients") ytmp <- tmp$coefficients[plot.set, ]      
    if ((yvar == "newy") & is.matrix(newx)) ytmp <- tmp$newy[plot.set,]
    if ((yvar == "newy") & (! is.matrix(newx))) ytmp <- tmp$newy[plot.set]
    if (yvar == "lambda") ytmp <- tmp$lambda[plot.set]    
    if (yvar == "dimension") ytmp <- tmp$dimension[plot.set]
    if (yvar == "R-square") ytmp <- tmp$r.square[plot.set]
  }
  xmax <- max(xtmp)
  xmin <- min(xtmp)
  ymax <- max(ytmp)
  ymin <- min(ytmp)
  if (xvar == "lambda") {
    plot(c(xmin,xmax),c(ymin,ymax), xaxt = "n", xlab = xvar, ylab=yvar, type="l",lty=0, main=x$method)
    axis(1, at = axTicks(1), labels = abs(axTicks(1)) )
  } 
  else plot(c(xmin,xmax),c(ymin,ymax), xlab = xvar, ylab=yvar, type="l",lty=0, main=x$method)
  if ( ((yvar=="coefficients")|(yvar=="newy")) & is.matrix(ytmp)) {
     index <- 1: (dim(ytmp)[2])
     if (yvar=="coefficients") {
       f.set <- rep(FALSE, max(length(index),max(predictors)))
       f.set[predictors] <- TRUE
       f.set <- f.set[1:length(index)]
     }
     else f.set <- rep(TRUE, length(index))
     axis(4, at = ytmp[dim(ytmp)[1],f.set], labels = index[f.set])
     for (j in 1: (dim(ytmp)[2])) 
       if (f.set[j]) lines(xtmp, ytmp[,j]) 
  }
  else lines(xtmp, ytmp)
}

########################################################
# internal functions used by plus()  
###

########################################################
# one single step in the iterative algorithm 
###
plus.single.step <- function(k, etatil,  z, b, tau, m, v, t, eps=1e-15, e) { 
## input 
# k:  number of parallelepiped indicators already found, including 0, k-1 steps completed 
# etatil, 2 X (|C|+1): label and new indicator for max crossing with dummy variable p+1
# NOTE: |C| > 0 required here
# z: t(x)y/n; b, tau: old turning point b at old tau; m, v, t: penalty; eps: numerical zero, 
## out put:   return(new.eta, etatil, s, tau,b, singular.Q, forced.stop, full.path) 
#   including info at the end of segment new.eta, e.g. etatil for the eta beyond new.eta 
 n <- dim(e$x)[1] # no dummy
 p <- dim(e$x)[2]
 # this is in step k, eta^{(0)},...,eta^{(k-1)} have already been found 
 sides <- length(etatil)/2-1 # number of boundaries for possible crossing
 singular.Q <- FALSE # this and next values are needed in case plus.check.eta not called 
 forced.stop <- TRUE 
 if (sides == 1) { # one-at-a-time scenario 
   new.eta <- plus.new.eta(e$eta[k,],etatil,rep(TRUE,2),p) # T for crossing the boundary and dummy 
   ## always a new eta in theory, only check when k is a multiplier of 50 to prevent numerical loop    
   if (k != round(k/50)*50) check.new.eta <- TRUE
   else check.new.eta <- plus.not.in.loop(k, new.eta, e=e) 
   if (check.new.eta==TRUE) { # do not run plus.check.eta if eta is not new, try xi = 1 first 
     tmp <- plus.check.eta(e$eta[k,], new.eta, etatil, 1, z, v, eps=eps,e=e) 
     ## tmp holds return(xi, singular.Q, s, invalid.cross)      
     singular.Q <- tmp$singular.Q  # save if invalid cross is due to singular Q, an almost nonevent
     if ( (m>1) & ((tmp$invalid.cross) == TRUE) ) { # try xi = -1
       tmp <- plus.check.eta(e$eta[k,], new.eta, etatil, -1, z, v, eps=eps,e=e)
       if (singular.Q == FALSE) singular.Q <- tmp$singular.Q
     }
     forced.stop <- tmp$invalid.cross # due to singular Q or not 
   } # if (check.new.eta==T)
 }  # end of the one-at-a-time case 
 else { # sides > 1, first check ties
   etatil <- etatil[, c((! e$ties[ etatil[1,1:sides] ]),TRUE) ] # remove known ties 
   sides <- length(etatil)/2-1
   if (sides > 1) { # remove new ties 
     if (e$use.Gram) Sgm.B <- e$Sigma[etatil[1,1:sides], etatil[1,1:sides]]
     else Sgm.B <- e$tx[etatil[1,1:sides],] %*% e$x[,etatil[1,1:sides]]        
     for (i in  1:(sides-1)) for (j in (i+1):sides) 
       if (Sgm.B[i,j] > 1- eps) {
         e$ties[j] <- TRUE
        }
     etatil <- etatil[, c(! e$ties[ etatil[1,1:sides] ],TRUE) ] # remove known ties 
     sides <- length(etatil)/2-1
   }
   if (sides == 1) {# go back to the one-at-time case
     new.eta <- plus.new.eta(e$eta[k,],etatil,rep(TRUE,2),p) # T for crossing the boundary and dummy 
     ## always a new eta in theory, only check when k is a multiplier of 50 to prevent numerical loop    
     if (k != round(k/50)*50) check.new.eta <- TRUE
     else check.new.eta <- plus.not.in.loop(k, new.eta, e=e) 
     if (check.new.eta==TRUE) { # do not run plus.check.eta if eta is not new, try xi = 1 first 
       tmp <- plus.check.eta(e$eta[k,], new.eta, etatil, 1, z, v, eps=eps,e=e) 
       ## tmp holds return(xi, singular.Q, s, invalid.cross)      
       singular.Q <- tmp$singular.Q  # save if invalid cross is due to singular Q, an almost nonevent
       if ( (m>1) & ((tmp$invalid.cross) == TRUE) ) { # try xi = -1
         tmp <- plus.check.eta(e$eta[k,], new.eta, etatil, -1, z, v, eps=eps,e=e)
         if (singular.Q == FALSE) singular.Q <- tmp$singular.Q
       }
       forced.stop <- tmp$invalid.cross # due to singular Q or not 
     } # if (check.new.eta==T)
   }  # end of the one-at-a-time case 
   else if (sides > 1) { # general case, (more than) two-at-a-time 
     M <- matrix(FALSE, 2^sides, sides + 1) # generate all possible candidates for new eta
     M[, sides+1] <- TRUE # last T for dummy
     M[ (2^(sides-1)+1):(2^sides), 1 ] <- TRUE # T for cross the indicated side
     for (i in (2:sides) ) M[,i] <- M[(2^(sides-i)+1):(3*2^(sides-i)),i-1]
     M <- M[order(apply(M,1,sum)),] # closer eta (smaller number of crossings) has higher priority  
     # initial while loop look for new eta 
     continue.while <- TRUE 
     j <- 1
     while (continue.while == TRUE) {
       new.eta <- plus.new.eta(e$eta[k,],etatil,M[j+1,],p) # skip the first M[1,]=F...FT
       new.is.new <- plus.not.in.loop(k, new.eta, e=e) 
       if (new.is.new == TRUE) { # try xi =1
         tmp <- plus.check.eta(e$eta[k,], new.eta, etatil, 1,  z, v, eps=eps, e=e)
         ## tmp holds return(xi, singular.Q, s, invalid.cross)      
         ## invalid.cross == T iff  either (singular.Q == T) or (xi and s fail to cross to new.eta)
         if (singular.Q == FALSE) singular.Q <- tmp$singular.Q
         if ( (m>1) & ((tmp$invalid.cross) == TRUE)) { # try xi = -1
           tmp <- plus.check.eta(e$eta[k,], new.eta, etatil, -1, z, v, eps=eps, e=e)
           if (singular.Q == FALSE) singular.Q <- tmp$singular.Q  # singular.Q should be T if happens once 
         }
         forced.stop <- tmp$invalid.cross # forced.stop is true if all "not new" or invalid
         continue.while <- tmp$invalid.cross # exit while loop for valid cross to new eta 
       } # end if (new.is.new == T)
       j <- j+1 # continue.while does not change when new.is.new == F
       if (j==2^sides) continue.while <- FALSE
     } # end while, if forced.stop == F, then tmp holds info for a valid new eta 
   } # end of (more than) two-at-a-time
   else if (sides < 1) new.eta <- rep(m+2,p) 
 }
 full.path <- FALSE
 if (forced.stop == TRUE) { # terminate the program, need value for the 9 out of 10 var
   ## dummy output for return(new.eta, etatil, s, tau,b, singular.Q, forced.stop, full.path)
   ##  values have already been assigned to new.eta, singular.Q, forced.stop, full.path
   etatil <- c(p+1,m+2)  # etatil <- as.matrix(c(p+1,m+2)) 
   s <- rep(0,p)
   tau <- 0
   b <- rep(0,p)
 }
if (forced.stop==FALSE) { # valid new eta has been found as new.eta
  ## tmp holds return(xi, singular.Q, s, invalid.cross), singular.Q = invalid.cross = F      
  xi <- tmp$xi
  s <- tmp$s
  ## need values for return(etatil, tau,b, full.path); new.eta, singular.Q, forced.stop done 
  a.set <- (new.eta!=0)
  if ( k - 50*round(k/50) == 1) # directly calculate the gradient to prevent accumulation of error 
  # e$grad <- tau*z - as.matrix(e$Sigma[,e$a.set.var])%*%b[a.set]    
  if (sum(a.set)<2) 
     e$grad <- tau*z - (e$Sigma[,e$a.set.var])*b[a.set]          # gradient
  else
     e$grad <- tau*z - (e$Sigma[,e$a.set.var])%*%b[a.set]     
  ## compute DeltaJ and new etatil   
  DeltaJ <- rep(-1,p) # -1 represents infinity     
  new.etatil <- rep(m+2,p) # m+2 is the dummy value for etatil calculation 
  ind <- (xi*s>=eps) & (new.eta != 0) & (new.eta < m)
  if (sum(ind)>0) { # move towards higher eta value
    DeltaJ[ind] <- xi*(e$t.fun[new.eta[ind]+1+m]-b[ind])*plus.recip(s[ind],eps=eps)
    new.etatil[ind] <- new.eta[ind] + 1
  }
  ind <- (xi*s<= -eps) & (new.eta != 0) & (new.eta > -m)
  if (sum(ind)>0) { # move towards lower eta value
    DeltaJ[ind] <- xi*(e$t.fun[new.eta[ind]+m]-b[ind])*plus.recip(s[ind],eps=eps)
    new.etatil[ind] <- new.eta[ind] - 1
  }
  ind <- (xi*e$g.prime>= eps) & (new.eta == 0)
  if (sum(ind)>0) { # new variable to with positive b
    DeltaJ[ind] <- xi*(1-e$grad[ind])*plus.recip(e$g.prime[ind],eps=eps)
    new.etatil[ind] <- 1
  }   
  ind <- (xi*e$g.prime<= -eps) & (new.eta == 0)
  if (sum(ind)>0) { # new variable to with negative b 
    DeltaJ[ind] <- xi*(-1-e$grad[ind])*plus.recip(e$g.prime[ind],eps=eps)
    new.etatil[ind] <- -1
  }
  # compute Delta, new.etatil, full.path
  old.etatil <- rep(m+2,p+1) # again m+2 is the dummy value
  old.etatil[etatil[1,]] <- etatil[2,]
  old.etatil <- old.etatil[1:p]  # old etatil in long format 
  if ( all(DeltaJ == (-1)) ) { # all infinity
    Delta <- -1                              # lse has attained for m > 1 and will be the next for lasso  
    full.path <- TRUE
    index.new.etatil <- rep(FALSE,p) # no new sides, all done 
  }
  else { # some stopping time is finite
    Delta <- min(DeltaJ[DeltaJ != (-1)])
    # index.new.etatil <- (DeltaJ < (Delta + eps)) & (DeltaJ !=(-1)) # does not work well 
    index.new.etatil <- DeltaJ==Delta # sides hit, numerically delicate 
    ind <- (abs(s) < eps) & (new.eta != 0) & (old.etatil !=m+2)  
    if (sum(ind)>0) { # sides not moving 
      new.etatil[ind] <- old.etatil[ind]
      index.new.etatil[ind] <- TRUE
    }
    ind <- (abs(e$g.prime) < eps) & (new.eta == 0) & (old.etatil !=m+2) 
    if (sum(ind)>0) { # sides not moving  
      new.etatil[ind] <- old.etatil[ind]
      index.new.etatil[ind] <- TRUE
    }
  } # end of else 
  # calculate the new etatil
  index.new.etatil <- c(index.new.etatil,TRUE) # add dummy
  e$etatil2[2,1:p] <- new.etatil  # other entries of e$etatil2 never change value 
  etatil <- e$etatil2[,index.new.etatil] 
  # etatil is m+2 at dummy column p+1 when full.path == T 
  # still need to compute return(tau,b) 
  if (Delta != -1) {
    tau <- tau + xi*Delta
    b <- b + xi*Delta*s  # update b, b[k+1,] <- b[k,] + ...
    e$grad <- e$grad + xi*Delta*e$g.prime 
    if (tau * eps > 1/5) {  # end at perfect fit with zero gradient numerically
      full.path <- TRUE
      index.new.etatil <- rep(FALSE,p) # 
    } 
    if (tau < e$tau1 - eps) { # shall not return to the beginning
      forced.stop <- TRUE
    }
  }
  else {
    tau <- xi
    b <- s
  }
} # end if forced.stop == F 
return(list(new.eta=new.eta, etatil=etatil, s=s, tau=tau, b=b, singular.Q=singular.Q, 
                   forced.stop = forced.stop, full.path = full.path)) 
}

########################################################
# check the validity of a specific parallelepiped 
###
plus.check.eta <- function(old.eta, new.eta, etatil, xi, z, v, eps=1e-15, e) { 
### test a possible new eta and xi combination
## input
# old.eta, 1 X p: old parallelepiped indicator; new.eta: the new one to be tried 
# etatil, 2 X (|C|+1): label and new indicator for max crossing, with a dummy variable p+1
# xi: the sign of new.tau - old.tau
# z: t(x)y/n; v: 2nd derivative of penalty; eps: numerical zero; e: environment including x, y
## out put:  return(xi, s, singular.Q, invalid.cross)
# singular.Q: error for solve(Q, z[a.set])
# invalid.cross: either singular.Q or fail to cross to new.eta
 p <- dim(e$x)[2]
 a.set <- rep(TRUE,p)
 a.set[new.eta==0] <- FALSE
 ##  update e$Sigma, e$var.list and e$a.set.var
 ##  e$Sigma[,e$a.set.var] will equal to t(x)%*%x[,a.set]/n when done 
 if (e$use.Gram) e$a.set.var <- a.set # e$Sigma is complete
 else if ( sum(e$var.list) == 0 ) { # the initial case of e$Sigma = 0 
   e$var.list <- order(!a.set)[1:sum(a.set)]
   e$a.set.var <- 1:sum(a.set)
   e$Sigma[, e$a.set.var] <- e$tx %*% e$x[,a.set]
 } 
 else { # incrementally add new columns to e$Sigma if necessary 
   new.var <- a.set
   new.var[e$var.list] <- FALSE # new.var are those in a.set but nor in e$var.list 
   if (sum(new.var)>0) { # need to add 
     e$var.list <- c(e$var.list, order(!new.var)[1:sum(new.var)])
     if (length(e$var.list) > (dim(e$Sigma)[2])) { # need to make e$Sigma larger 
       tmp.Sigma <- e$Sigma
       e$Sigma <- matrix(0,p,dim(tmp.Sigma)[2]+100)
       e$Sigma[,1:(dim(tmp.Sigma)[2])] <- tmp.Sigma
     }
     e$Sigma[,(length(e$var.list) - sum(new.var)+1):length(e$var.list)] <- e$tx %*% e$x[,new.var]
   }
   e$a.set.var <- (1:length(e$var.list))[a.set[e$var.list]] # the equivalent labels of a.set
   e$a.set.var <- e$a.set.var[order(e$var.list[ e$a.set.var])]  
 }
 ## calculate for return(xi, singular.Q, s, invalid.cross), xi is given 
 if (sum(a.set)==1) Q <- e$Sigma[a.set,e$a.set.var] - v[ abs(new.eta[a.set]) ]
 else Q <- e$Sigma[a.set,e$a.set.var] - diag( v[ abs(new.eta[a.set]) ] )
 singular.Q <- FALSE
 # LINPACK=T generated an incorrect calculation, due to the use of single precision
 # s1 <- try(solve(Q, z[a.set], LINPACK = (length(Q)>100)), silent=T)
 # s1 <- try(solve(Q, z[a.set], symmetry=T), silent=T)  # option has no effect 
 s1 <- try(solve(Q, z[a.set]), silent=TRUE)
 if (inherits(s1,"try-error")==TRUE) {
   s1 <- rep(0, sum(a.set))
   singular.Q <- TRUE
 }
 s <- rep(0,p)
 if (singular.Q==FALSE) {
   s[ a.set ] <- s1
   ## check the new s
   etatil2 <- c(old.eta,0) # etatil in full length
   etatil2[ etatil[1,] ] <- etatil[2,]
   etatil2 <- etatil2[1:p] # remove dummy
   #  e$g.prime <- z - as.matrix((e$Sigma[,e$a.set.var]))%*%s1  # Eq. (2.13)
   if (length(s1)==1)  e$g.prime <- z - (e$Sigma[,e$a.set.var])*s1  # Eq. (2.13)
   else  e$g.prime <- z - (e$Sigma[,e$a.set.var])%*%s1 
   flag <- rep(TRUE,p)  # T means violation 
   flag[old.eta==etatil2] <- FALSE # noncritical j are fine
   # check according to (2.15) with eps = numerical zero
   flag[ (old.eta != new.eta) & (new.eta != 0) & (xi*(new.eta-old.eta)*s > -eps) ] <- FALSE
   flag[ (old.eta == new.eta) & (new.eta != 0) & (etatil2 != old.eta) & (xi*(etatil2-old.eta)*s < eps) ] <- FALSE 
   flag[ (old.eta != 0) & (new.eta == 0) & (xi*old.eta*e$g.prime < eps)] <- FALSE
   flag[ (old.eta == 0) & (new.eta == 0) & (etatil2 != 0) & (xi*etatil2*e$g.prime < eps)] <- FALSE
   invalid.cross <- sum(flag) > 0 # T here iff the new eta and xi pair is invalid and singular.Q==FALSE
 }
 else invalid.cross <- TRUE # invalid.cross is T if singular.Q == T
 return(list(xi=xi, singular.Q=singular.Q, s=s, invalid.cross=invalid.cross))
} 

########################################################
# find the segment and weight in the path for the hitting point
###
plus.hit.points <- function (lam, lam.path) {
	if (is.null(lam) || !is.numeric(lam) || is.nan(sum(lam)))
	{        
		lam <- sort(lam.path, decreasing = TRUE)
		cat("Warning message: ")
		cat("invalid lam; lam.path used\n")
	}
	if (max(lam)< min(lam.path))
	{        
		lam <- sort(lam.path, decreasing = TRUE)
		cat("Warning message: ")
		cat("lam not reached by the plus path; lam.path used\n")
	}
	lam <- mapply(function(x) min(x,max(lam.path)),lam)
	if (min(lam) < min(lam.path))
	{
		lam <- lam[lam>=min(lam.path)]
		cat("Warning message: ")
		cat("some lam not reached by the plus path and dropped\n")
	}
	lam.lesser.path <- mapply(function(x) x<=lam.path,lam)
	lam.greater.path <- mapply(function(x) x>=lam.path,lam)
	hit <- (lam.lesser.path[1:(length(lam.path)-1),] & 
			lam.greater.path[2:length(lam.path),])
	k1 <- apply(as.matrix(hit),2,function(bool.vec) which(bool.vec)[1])
	k2 <- k1+1
	w1 <- (lam - lam.path[k2])/(lam.path[k1]-lam.path[k2])
	return(list("k1"=k1,"k2"=k2,"w1"=w1,"total.hits"=length(lam)))
}


########################################################
# create a new eta indicating a new parallelepiped 
###
plus.new.eta <- function(eta, etatil, cross,p) {
  eta <- c(eta,0) # put the dummy in the zero interval
  eta[ etatil[1,cross] ] <- etatil[2, cross]
  eta <- eta[1:p] # remove dummy
  return(eta)
}

########################################################
# check if a newly created eta is indeed new 
###
plus.not.in.loop <- function( k, new.eta,e) {
## input: eta, (k+1)Xp; new.eta, 1Xp
## output is T iff new.eta is not any of eta[i,], i\le k
  new.is.new <- TRUE
  for (i in 1:k) 
    if( all(e$eta[i,]==new.eta)) return(FALSE)
  return(new.is.new)
}

########################################################
# compute the knots and second derivative of a penalty  
###
plus.penalty <- function(m,a=0) {
### lasso for m=1, MCP for m=2, SCAD for m=3; default gamma = a = 3.7
### vectors v and t with dummy are used for lasso to maintain vector attribute 
m <- round(m)
if (m < 1) m <- 2
if (m > 3) m <- 2
if (m == 2) {
  if (a==0) a <- 3.7
  v <- c(1/a,0)
  t <- c(0,a)
}
if (m == 1) {
  a <- 0
  v <- c(0,0) 
  t <- c(0,0)
}
if (m==3) {
  if (a==0) a <- 3.7
  v <- c(0,1/(a-1),0) 
  t <- c(0,1,a)
}
return(list(m=m, gamma=a, v=v,t=t))
}

########################################################
# stablized reciprocals 
###
plus.recip <- function(Q, eps=1e-15){
### compute reciprocal of reals    
Q[ (Q>=0) & (Q <eps) ] <- eps
Q[ (Q<0) & (Q> - eps) ] <- - eps
return(1/Q)
}

########################################################
# generate knots for both sides of zero
###
plus.knots <- function(t=0,m=0) {
### compute the function t() as in (2.6) without infinity and - infinity
### note that t_1=t(1)=t(0) in (2.6)
if (m==0) { # default MC penalty
  m <- 2 
  t <- c(0,3)
}  
t2 <- rep(0,2*m)
t2[1:m] <- - t[m:1]
t2[(m+1):(2*m)] <- t[1:m]
return(c(t2,0))
} 

profile.logistf <-
function(fitted,  which, variable, steps=100, pitch = 0.05, limits,
                    alpha = 0.05,  firth = TRUE,
                    legends = TRUE,  control, plcontrol, plot=FALSE, ...){

# by MP, 06.02.01
# adapted and renamed by GH, 10.03.11 (incredible! 10 years!)
# again adapted and renamed by GH, 13.05.13
# which  ... righthand formula des zu plottenden Term (z.B. ~B oder ~A:D)
# pitch  ... distances between points in std's
# limits ... vector of MIN & MAX in std's, default=extremes of both CI's
#            +- 0.5 std. of beta
#

formula<-fitted$formula
data<-fitted$data
if(is.null(data)) stop("Call logistf with dataout=TRUE.\n")

# Next line added by Harry Southworth, 22/10/02.
 if (missing(which) & missing(variable)) stop("You must specify a variable: either by which (a one-sided formula) or by variable.")
 if (missing(control)) control<-logistf.control()
 if (missing(plcontrol)) plcontrol<-logistpl.control()
 
   call <- match.call()

    
    mf<-model.frame(fitted$formula, data=fitted$data)
##    mf <- match.call(expand.dots =FALSE)
 ##   m <- match(c("formula", "data","weights", "na.action", 
##        "offset"), names(mf), 0L)
## #   mf<-model.frame(formula, data=data, weights=weights)
##    mf <- mf[c(1, m)]
##    mf$drop.unused.levels <- TRUE
 ##   mf[[1L]] <- as.name("model.frame")
 ##   mf <- eval(mf, parent.frame())
#    y <- model.response(mf)
    y <- fitted$y
    n <- length(y)
    x <- model.matrix(formula, data = data) ## Model-Matrix 
    cov.name <- labels(x)[[2]]
   # weight <- as.vector(model.weights(mf)  )
    offset <- model.offset(mf)   
    weight<-model.weights(mf)
    if (is.null(offset)) offset<-rep(0,n)
    else offset<-as.vector(offset)
    if (is.null(weight)) weight<-rep(1,n)

    
    cov.name <- labels(x)[[2]]
    k <- ncol(x)
    if (dimnames(x)[[2]][1] == "(Intercept)")  {
        int <- 1
        coltotest <- 2:k
    }

    else {
        int <- 0
        coltotest <-1:k
    }
  if(!missing(which)) cov.name2 <- labels(model.matrix(which, data = data))[[2]] ## Label des Test-Fakt.
  else cov.name2 <- variable
  pos <- match(cov.name2, cov.name) ## Position des Testfakors
  fit<-logistf.fit(x, y, weight=weight, offset=offset, firth=firth, control=control) 
  std.pos <- diag(fit$var)[pos]^0.5
 
 coefs <- fit$beta ## "normale" Koeffizienten
 covs <- fit$var ## Varianzen
# n <- nrow(data)
 n <- nrow(x)
 cov.name <- labels(x)[[2]]
 if(missing(limits)) {
  lim.pl<-numeric(0)
  LL.0 <- fit$loglik - qchisq(1 - alpha, 1)/2
  lower.fit<-logistpl(x, y, init=fit$beta, weight=weight, offset=offset, firth=firth, LL.0=LL.0, which=-1, i=pos, plcontrol=plcontrol)
  lim.pl[1]<-lower.fit$beta
  upper.fit<-logistpl(x, y, init=fit$beta, weight=weight, offset=offset, firth=firth, LL.0=LL.0, which=1, i=pos, plcontrol=plcontrol)
  lim.pl[2]<-upper.fit$beta
  lim.pl <- (lim.pl - coefs[pos])/std.pos
  limits <- c(min(qnorm(alpha/2), lim.pl[1]) - 0.5, max(qnorm(1 - alpha/2), lim.pl[2]) + 0.5)
 }

 limits <- c(floor(limits[1]/pitch) * pitch, ceiling(limits[2]/pitch) * pitch)

 knots <- seq(limits[1], limits[2], diff(limits)/steps)
 nn <- length(knots)
 res <- matrix(knots, nn, 3) #initialisiere Werte
 dimnames(res) <- list(1:nn, c("std", cov.name2, "log-likelihood"))
 for(i in 1:nn) {
  res[i, 2] <- coefs[pos] + covs[pos, pos]^0.5 * knots[i]
  if(i == 1){
     init<-lower.fit$betahist[nrow(lower.fit$betahist),]
     init[pos]<-res[i,2]
     xx <- logistf.fit(x, y, weight=weight, offset=offset, firth=firth, col.fit<-(1:k)[-pos], init=init,
                   control=control) 
  }     
  else {
     init<-xx$beta
     init[pos]<-res[i,2]
     xx <- logistf.fit(x, y, weight=weight, offset=offset, firth=firth, col.fit<-(1:k)[-pos], init=init,
                   control=control) # use solution from last step
  }
  res[i, 3] <- xx$loglik
 }

 #### Graphischer Output:

 if(plot==TRUE){
   my.par <- act.par <- par()
   my.par$mai[3] <- 1.65 * act.par$mai[3]
## if(legends) my.par$mai[1] <- 2 * act.par$mai[1]
   par(mai = my.par$mai)
   ind <- (1:nn)[round(4 * res[, 1]) == round(4 * res[, 1], 10)]
   if(length(ind) == 0) ind <- 1:nn
   pp <- max(res[, 3]) - 0.5 * res[, 1]^2

   plot(res[, -1], type = "l", xlab=expression(beta)) ##Profile likelihood

 #lines(res[,2], pp, lty=4)  #<<<Wald approximative profile lik. >>>

   points(res[res[, 1] == 0, 2], max(res[, 3])) ##Maximum of likelihood

   segments(min(res[, 2]), max(res[, 3]) - 0.5 * qchisq(1 - alpha, 1),
                max(res[, 2]), max(res[, 3]) - 0.5 * qchisq(1 - alpha, 1), lty = 3) ##refer.line

   yy <- par("usr")[4] - (par("usr")[4] - par("usr")[3]) * c(0.9, 0.95)

   segments(fit$beta[pos] - qnorm(alpha/2) * std.pos, yy[1], fit$beta[pos] - qnorm(1 - alpha/2) *
            std.pos, yy[1], lty = 6) ##Wald-CI
   segments(lower.fit$beta, yy[2], upper.fit$beta, yy[2], lty = 8) ##prof.pen.lik.-CI

   axis(side = 3, at = res[ind, 2], labels = res[ind, 1])

   mtext(expression(paste("distance from ", hat(beta)," in multiples of ", hat(sigma))), side = 3, line = 3)
  ## mtext(expression(paste(beta, " of ", cov.name2)), side=1, line = 3)
   par(mai = act.par$mai)
 
   if (legends)
    {
     legend(x=fit$beta[pos],
            y=min((min(res[,3])+max(res[,3]))/2,(max(res[, 3]) - 0.5 * qchisq(1 - alpha, 1))),
        legend=c("Profile penalized likelihood",
                 paste(100 * (1 - alpha),"%-reference line"),
                 "Wald confidence interval",
                 "Profile likelihood confidence interval"),
        lty=c(1,3,6,8), 
        text.col=c("black","black","black","black"), ncol=1, bty="n", xjust=0.5)
    }


    title(paste("Profile of penalized likelihood for Variable",cov.name2))
 }
 signed.root<-sqrt(2*(-res[,3]+max(res[,3])))*sign(res[,1])
 cdf<-pnorm(signed.root)
 
 results<-list(beta=res[,2], stdbeta=res[,1], profile=2*(res[,3]-max(res[,3])), loglike=res[,3], signed.root=signed.root, cdf=cdf)
 attr(results,"class")<-"logistf.profile"
 results
}


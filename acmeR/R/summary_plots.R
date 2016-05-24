#'@import graphics grDevices
#'
#    Reports and Summaries:
#       scav.plot()   Kaplan-Meyer Survival Plot w/overlay
#       srch.plot()   Empirical and Diminishing Proficiency Plots
#       scav.cont()   Contour plot of nllh.scav for Weibull removal
#       lscav.cont()  Contour plot of nllh.scav, on log scale
#       srch.cont()   Contour plot of nllh.srch, with bleed-though fixed
##################################################################
#
scav.plot <- function(rd,  spec="", verb=TRUE, add=FALSE, name="",
             ps="", pub=FALSE) {
  # Kaplan-Meyer Survival Plot, with optional overlay
  # of best-fit from Weibull and Exponential models
  
  spec <- names(sort(table(rd$scav$Species),decreasing=TRUE));
  rd <- subspec(rd, spec);
  spec <- names(sort(table(rd$scav$Species),decreasing=TRUE));
  if(missing(name)) {
    if(length(spec)>4) {
      name <- paste(c(spec[1:4],"..."),sep=",",collapse="+");
    } else {
      name <- spec2name(spec);
    }
  }
  scav <- rd$scav;
  spec <- class2spec(spec);
  ok   <- scav$Species %in% spec;
  scav <- scav[ok,];
  lo   <- as.numeric(difftime(scav$Lo, scav$Placed, units="days"));
  hi   <- as.numeric(difftime(scav$Hi, scav$Placed, units="days"));
  hi[hi>1e4] <- Inf;  # More than 27 years?  Carcass NEVER found.
  n    <- sum(ok)
  top  <- max(lo, hi[is.finite(hi)]);
  nfin <- sum(is.finite(hi))
  t.lo <- c(0,sort(lo));
  t.hi <- c(0,sort(hi));
  y    <- (n:0)/n;
  rv   <- list()
##################################################################
  mle.w <- mle.wei(rd, spec=spec, v=FALSE);
# Choose Sur.w evenly, not t, for smoother plots:
  Sur.w <- seq(1,exp(-(mle.w$rho * top)^mle.w$alp),,100);
  t     <- (- log(Sur.w))^(1/mle.w$alp)/mle.w$rho;
  msg.w <- bquote(paste("Weibull ",
              P(T>t), " = exp", group("{",-( rho * t )^ alpha, "}"), ": ",
              hat(alpha) == .(round(mle.w$alp,4)), ", ",
              hat(rho)   == .(round(mle.w$rho,4)), ", E[T] = ",
              .(signif(mle.w$tij,4)), "d"));
  rv$wei <- mle.w;
  mle.e  <- mle.exp(rd, spec, v=FALSE);
  Sur.e  <- exp(- (mle.e$rho * t));
#  msg.e  <- bquote(paste("Exponential: E[T] = ",
#                        .(signif(mle.e$tij,4)), "d"));
  msg.e <- bquote(paste("Exponential ",
                  P(T>t), " = exp", group("{",- rho * t, "}"), ": ",
                  hat(rho)  == .(round(mle.e$rho,4)), ", E[T] = ",
                        .(signif(mle.e$tij,4)), "d"));
  rv$exp <- mle.e;
##################################################################
  if(!add) {
    if(nchar(ps)>0) {
      if(pub) {  # Prepare for publication?
          postscript(file=ps, paper="special", height=5,
                     width=12, horizontal=FALSE);
      } else {
          postscript(file=ps, horizontal=TRUE);
      }
    }
    opar <- par(no.readonly=TRUE);
    par(mar=opar$mar+c(1,1,0,0));
  }
  
  plot(t.lo, 100*y, type="S", ylim=c(0,100),
       cex.axis=1.1, cex.lab=1.25,
       xlab="Carcass Age t (days)",
       ylab="Carcasses Remaining (%)");
  title("Carcass Removal",cex.main=1.75)
  lines(c(t.hi[1:nfin],top),
        100*c(y[1:nfin],1-nfin/n), type="s");
  if(TRUE) {
    text(top/2, 95, msg.w, cex=1.25, col="blue");
    text(top/2, 83, msg.e, cex=1.25, col="red");
    text(top/2, 71, paste("Carcass type:", name), cex=1.25);
  }
  lines(t, 100*Sur.w, col="blue");         # Best fit with Weibull model
  lines(t, 100*Sur.e, col="red", lty=2);   # Best fit with Exponential
##################################################################
  if(!add) {
    par(opar);
    if(nchar(ps)>0) {
      dev.off();
      print(paste("Plot stored as postscript file",ps));
    }
  }
  invisible(rv);
}  
##################################################################
#
srch.plot <- function(rd, spec="", width=5, dx=1.5, ab=numeric(0), bt1=FALSE,
              ttop=Inf, ytop=0, ps="", add=FALSE, verb=TRUE, pub=FALSE) {
  # Plots of Empirical & Expo Dim search success prob vs. time
  
  rd      <- subspec(rd, spec);
  if(!length(ab)) {
    if(bt1) {
      srch <- mle.bt1(rd, v=FALSE);
      ab  <- c(a=srch$a, b=srch$b);
    } else {
      srch <- mle.srch(rd, v=FALSE);
      ab  <- c(a=srch$a.hat, b=srch$b.hat, bto=srch$bt.hat/(1-srch$bt.hat));
    }
  }
  scav   <- rd$scav;  n.carc <- dim(scav)[1];
  srch   <- rd$srch;  n.srch <- dim(srch)[1];
# Get search ages
  ages   <- numeric(n.srch);            # Ages of carcasses at search
  for(i in 1:n.carc) {
    ok       <- srch$Id==scav$Id[i];    # Entries for i'th placed carcass
    ages[ok] <- as.numeric(difftime(srch$Date[ok],
                                    scav$Placed[i], units="days"));
  }
  ttop  <- min(ttop, max(ages));
  t    <- seq(0,ttop,,101);
  phat <- exp(-cbind(1,t) %*% ab[1:2]);
  if(length(ab)>2) {
    bt <- ab[3]/(1+ab[3]);
    phat <- phat * bt^floor(t/rd$Ik["mu"]);
  }
  if(!add) {
    if(nchar(ps)>0) {
      if(pub) {  # Prepare for publication?
        postscript(file=ps, paper="special",
             height=5, width=12, horizontal=FALSE);
      } else {
        postscript(file=ps, horizontal=TRUE);
      }
    }
    opar <- par(no.readonly=TRUE);
    par(mar=opar$mar+c(1,1,0,0));
  }
  yl <- c(0, max(ceiling(1.1*exp(-ab[1])*10)/10, min(1,ytop)));

  plot(t, phat, ylim=yl, xlab="Carcass Age t (days)",
       ylab=expression(paste("Estimated Search Proficiency ", hat(S)~(t) )),
       type="n", cex.axis=1.1, cex.lab=1.25);
  lines(t, phat, lwd=2, col="blue");                  #  Model fit
  succ <- rd$srch$Found;
# Empirical:
  points(ages+(runif(n.srch)-0.5)*dx, succ*yl[2], pch="|"); # Add uniform jitter
  M <- exp(-(1/width)*outer(t, ages, function(x,y) {abs(x-y)}));
  y <- (M %*% succ) / (M %*% rep(1,length(succ)));
  lines(t, y, col="black", lty=2, lwd=2);             # Empirical MA fit
  lines(range(t), rep(mean(succ[ages<=ttop]),2),       # Constant proficiency
                  col="red", lty=3, lwd=2);
  if(verb) {
    lab.s <- paste("Succ (", sum( succ[ages<=ttop]), "):",sep="");
    lab.f <- paste("Fail (", sum(!succ[ages<=ttop]), "):",sep="");
    summ  <- rbind(summary(ages[ages<=ttop & succ]),
                   summary(ages[ages<=ttop & !succ]));
    dimnames(summ)[[1]] <- c(lab.s,lab.f);
    print(summ);
    if(length(ab)==2) {
      print(paste("nllh.bt1:",  round(nllh.bt1(ab, rd),4)));
    } else {
      print(paste("nllh.srch:", round(nllh.srch(ab, rd),4)));
    }
  }
  if(bt1) {
    model <- expression(paste("Model: ",
        S(t) == e^{- a - b * t}));
    mod2 <- bquote(paste(
          hat(a) == .(round(ab[1],4)), ", ",
          hat(b) == .(round(ab[2],4)), "."));
  } else {
    model <- expression(paste("Model: ", S(t) ==
         B^k * e^{- a - b * t}));                      # B was: theta
    mod2 <- bquote(paste(
          hat(a) == .(round(ab[1],4)), ", ",
          hat(b) == .(round(ab[2],4)), ", ",
          hat(B) == .(round(ab[3]/(1+ab[3]),4)), "."));  # Was: theta
  }
  rv <- legend(x=quantile(t,0.3), y=yl[2], legend=c(
        model, "Constant", "Empirical"), bty="n",
        lty=c(1, 3, 2),  col=c("blue", "red", "black"),
        text.col=c("blue", "red", "black"),
        lwd=c(2,2,2), cex=1.25);
  text(x=rv$rect$left+1.1*rv$rect$w, pos=4,
       y=rv$rect$top-0.30*rv$rect$h, 
        mod2,cex=1.2,col="blue");
  spec <- names(sort(table(rd$scav$Species),decreasing=TRUE));
  if(length(spec)>4) {
    spec <- paste(c(spec[1:4],"..."),sep=",",collapse="+");
  } else {
    spec <- spec2name(spec);
  }
#  spec <- spec2name(spec);
  title(paste("Search Proficiency for: ", paste(spec, collapse=", ")),
        cex.main=1.75);
  if(!add) {
    par(opar);
    if(nchar(ps)>0) {
      dev.off();
      print(paste("Plot stored as postscript file",ps));
    }
  }
  invisible(list(mu=mean(succ), ab=ab));
}

##################################################################
#
scav.cont <- function(rd, al, rl, sd=1.5, nxy=41) {
  # Simple contour plot of Weibull nllh in (alpha,rho)
  if(missing(rd)) rd <- read.data(spec="");
  scav <- rd$scav;
  mle  <- mle.wei(rd,v=FALSE);
  if(missing(al))  al <- mle$alp * exp(c(-sd,sd)*mle$alp.se/mle$alp);
  if(missing(rl))  rl <- mle$rho * exp(c(-sd,sd)*mle$rho.se/mle$rho);
  alp <- seq(al[1],al[2],,nxy);
  rho <- seq(rl[1],rl[2],,nxy);
  z   <- matrix(NA,nrow=nxy,ncol=nxy);
  for(i in 1:nxy)
    for(j in 1:nxy)
      z[i,j] <- nllh.scav(c(alp[i],rho[j]),scav);
  contour(alp, rho, z,
          xlab=expression(alpha), ylab=expression(rho));
  text(mle$alp, mle$rho, "*", col="red", cex=2);
}  
##################################################################
#
lscav.cont <- function(rd, al=c(0.5,0.7), rl=c(0.04,0.10), sd=3, nxy=41) {
  # Contour plot (as above), on log scale
  if(missing(rd)) rd <- read.data(spec="");
  scav <- rd$scav;
  mle  <- mle.wei(rd,v=FALSE);
  if(missing(al))  al <- mle$alp * exp(c(-sd,sd)*mle$alp.se/mle$alp);
  if(missing(rl))  rl <- mle$rho * exp(c(-sd,sd)*mle$rho.se/mle$rho);
  alp <- seq(log(al[1]),log(al[2]),,nxy);
  rho <- seq(log(rl[1]),log(rl[2]),,nxy);
  z   <- matrix(NA,nrow=nxy,ncol=nxy);
  for(i in 1:nxy)
    for(j in 1:nxy)
      z[i,j] <- nllh.scav(c(exp(alp[i]),exp(rho[j])),scav);
  contour(alp, rho, z,
          xlab=expression(log (alpha)), ylab=expression(log (rho)));
  text(log(mle$alp), log(mle$rho), "*", col="red", cex=2);
}  
##################################################################
#
srch.cont <- function(rd, spec="", bt=1, al, bl, sd=1.0, nxy=41) {
  # Simple contour plot of Acme nllh in (a,b), fixed theta
  if(!missing(spec) || !nchar(spec)) rd <- subspec(rd, spec);
  mle  <- mle.srch(rd,v=FALSE);
  if(missing(al))  al <- mle$a.hat *
       exp(c(-sd,sd)*sqrt(mle$Sig[1,1])/mle$a.hat);
  if(missing(bl))  bl <- mle$b.hat *
       exp(c(-sd,sd)*sqrt(mle$Sig[2,2])/mle$b.hat);
  if(missing(bt))  bt <- mle$bt.hat;
  a <- seq(al[1],al[2],,nxy);
  b <- seq(bl[1],bl[2],,nxy);
  abt <- cbind(rep(a,nxy),rep(b,rep(nxy,nxy)),bt/(1-bt));
  z <- matrix(nllh.srch(abt, rd),nrow=nxy);
  contour(a,b,z, xlab="a", ylab="b", cex.lab=1.25);
  text(mle$a, mle$b, "*", col="red", cex=2);
}  


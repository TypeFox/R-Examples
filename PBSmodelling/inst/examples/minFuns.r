#------------------------
# Minimization Functions
#------------------------

calcMin <- function(pvec,func,tol=0) { # minimizes func given pvec
   Sfun <- function(S,pvec,func) { # function of surrogate parameters
      P <- restorePar(S,pvec);
      Uval <- func(P);
      #print(paste("S =",paste(round(S,3),collapse=","),"  ","P =",paste(round(P,3),collapse=",") )); #<-- for debugging
      return(Uval); };
   Sval <- scalePar(pvec); nS <- length(Sval);
   Fout <- nlm(f=Sfun,p=Sval,typsize=rep(0.005,nS), pvec=pvec, func=func)
   Pest <- Fout$estimate
   Pfin <- restorePar(Pest,pvec);
   Pmat <- cbind(Pfin,pvec[,2:4]);
   Pout <- list(start=pvec[,1], end=Pfin, surrogates=scalePar(Pmat), check=Pest,
           grad=Fout$grad, code=Fout$code, iters=Fout$iterations, fmin=Fout$minimum );
   cat("\n\n")
   return(Pout); };

scalePar <- function(pvec) { # Convert true parameters to surrogates
   Pval <- pvec[,1]; Pmin <- pvec[,2]; Pmax <- pvec[,3]; idx <- pvec[,4];
   Sval <- (Pval[idx]-Pmin[idx]) / (Pmax[idx]-Pmin[idx]);
   Sval <- pmax(Sval,0); Sval <- pmin(Sval,1);  # enforces the range
   S    <- (2/pi) * asin(sqrt(Sval));
   names(S) <- dimnames(pvec)[[1]][idx];
   return(S); }

restorePar <- function(S,pvec) { # Convert surrogates to true parameters
   Pval <- pvec[,1]; Pmin <- pvec[,2]; Pmax <- pvec[,3]; idx <- pvec[,4];
   if (sum(idx) != length(S)) stop("Warning: S & P not consistent/n");
   Pcon <- Pmin[idx] + (Pmax[idx]-Pmin[idx])*sin(pi*S/2)^2;
   P <- Pval; P[idx] <- Pcon; names(P) <- dimnames(pvec)[[1]];
   return(P); };

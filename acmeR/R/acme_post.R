#' Posterior Calculation of Mortality
#' 
#' Calculates and plots the posterior distribution of mortality count.
#' 
#' Assuming a Gamma(xi, lam) on the average daily mortality rate m, this model
#' treats the mortality M for the current period as Poisson-distributed with
#' mean m*I. The carcass count C will include "new" carcasses with a Bi(M,T)
#' distribution as well as "old" carcasses (if bt > 0). For derivation of 
#' resulting conditional pdf see Wolpert (2015).
#' 
#'
#'@param C Observed mortality count. Non-negative integer.
#'@param Rstar ACME inverse-inflation factor R*, reported by acme.summary() 
#'  as "Rstar."
#'@param T The first term in recursive calculation of Rstar, from acme.summary.
#'@param gam Values for highest posterior density credible interval.
#'@param I Interval length, days.
#'@param xlim 2-element vector of plotting ranges. Default first element of 0, 
#'second element of 2 greater than maximum calculated for larger hpd.
#'@param Mmax Maximimum value for which posterior probability is calculated.
#'@param xi First parameter of gamma prior. Default is 1/2 for Objective prior.
#'@param lam Second parameter of gamma prior. Default is 0 for Objective prior.
#'@param ps Postscript message. Default empty string suppresses output.
#'@param plotit Boolean to determine if plot should be created. Default is TRUE.
#'  
#'
#'@return The function invisibly returns a vector with input C, ACME estimate, posterior
#' mean, and credible interval ranges. If plotit = TRUE, it also plots the posterior 
#' probabilities for values in the range of xlim, and prints a short summery including
#' the true coverage probabilities.
#'
#'The parameter \code{plotit} should almost never be set to FALSE - if the user
#'desires the vector that is inivisibly returned, it is suggested to use the 
#'wrapper function \code{acme.table}.

#'@export
#'@import graphics stats grDevices
#'
#'@examples
#'acme.post(C=5, Rstar = .25, T = .2, gam = c(.9,.95), I = 5, xi = .5,lam = 0) 


###################################################################
acme.post <- function(C=0, Rstar=0.2496, T=0.1740, gam=c(0.5, 0.9), I=7,
                      xlim, Mmax = 200, xi=1/2, lam=0, ps="", plotit=TRUE) {
  if( T<0 | Rstar<T | 1<Rstar ) stop("ObjPost: Must have 0 < T <= Rstar < 1");
  if( xi+C<=0 | lam<0 | I<=0 ) stop("ObjPost: Check xi, C, lam, I");
  M <- seq(0, Mmax);
  
  #Table: for each C value - Posterior Mean, M estimate, credible ints from gam
  C_table <- matrix(numeric((3 + 2*length(gam))),ncol = (3 + 2*length(gam)));
  C_tab_names <- character(length(C_table))
  C_tab_names[1:3] <- c("C","M_hat","Post_Mean");
  for(i in 1:length(gam)){
    C_tab_names[2*(i+1)] <- paste(gam[i],"low",sep="_")
    C_tab_names[2*(i+1) + 1] <- paste(gam[i],"hi",sep="_")
  }
  colnames(C_table) <- C_tab_names;
  
  C_table[1] <- round(C);
  
  
  #Calculate posterior pdf
  if(Rstar>T) {                        # Positive bleed-through B>0
    lcon <- lgamma(xi+C+M) - lgamma(xi+C) - lgamma(M+1) +
      (xi+C)*log(Rstar+lam/I) - (xi+C+M)*log1p(Rstar-T+lam/I) +
      C*log1p(-T/Rstar) + M*log1p(-T);
    zed <- T * (1+Rstar-T+lam/I) /((1-T)*(Rstar-T));
    pmf <- exp(lcon) * my2F1(-C, -M, 1-xi-C-M, -zed);
  } else {                         # No bleed-through, B==0
    pmf <- dnbinom(M-C, xi+C,  (lam+Rstar*I)/(lam+I));
  }
  if(sum(pmf)+1e-5 < 1)
    warning(paste("Try increasing Mmax (now ",Mmax,"); sum(pmf)=",
                  sum(pmf),".",sep=""));
  mu <- sum(M*pmf);
  
  C_table[2] <- round(C/Rstar,2);
  C_table[3] <- round(mu,2);
  
  
  #Calculate hpd
  if((ng <-length(gam)) > 0) {
    lohi <- cbind(gam=gam, cvg=NA, lo=NA, hi=NA);
    o   <- order(-pmf);
    mycol = c("red", "blue", "black"); mycex=c(3, 2, 1); mypch=c(15,19,1);
    for(i in 1:ng) {
      hpd  <- o[1:min(which(cumsum(pmf[o]) >= gam[i]))];
      lohi[i,"cvg"] <- sum(pmf[hpd]);
      lohi[i,"lo"]  <- min(hpd);        # INDEX not value
      lohi[i,"hi"]  <- max(hpd);
      
    }
    Mtop <- max(lohi[,"hi"])+2;
  } else {
    Mtop <- min(which(cumsum(pmf) >= 0.90 ));
  }
  if(missing(xlim)) xlim <- c(0, M[Mtop]) else Mtop <- xlim[2];
  #
  if(nchar(ps)>0)
  { postscript(ps, paper="special", height=8, width=8, horizontal=F); }
  
  #get hpd intervals
  o    <- order(pmf,decreasing = TRUE);
  lg <- length(gam);
  gam.order <- sort(gam,decreasing=FALSE);
  for(i in 1:lg){
    max <- min(which((cumsum(pmf[o])>gam.order[i])))
    C_table[2*(i+1)] <- round(min(M[o[1:max]]))
    C_table[(2*(i+1) + 1)] <- round(max(M[o[1:max]]))
  }
  
  ###
  # Plotting
  if(plotit){
    opar <- par(no.readonly=TRUE);
    par(mar=(opar$mar+c(0,1,0,0)));  # Add some space on left for big ylab
    plot(M, pmf, xlab="Mortality m", cex.lab=1.5, cex.axis=1.5,type="n",
         ylim=range(c(0,pmf)), xlim=xlim,
         ylab=bquote(paste("P[ ", M[i] == m, " | ",
                           C[i] == .(C), " ]"))); 
    arrows(C/Rstar, max(pmf)/3,,0, col="black", lwd=2);
    text  (C/Rstar, max(pmf)/2.5, expression(hat(M[i])),cex=2); 
    arrows(mu,  max(pmf)/4,,0, col="black", lwd=2);
    text  (mu,  max(pmf)/3.5, expression(bar(M[i])),cex=2);
    if(ng > 0) {
      mycol <- c("red", "blue", "black");
      mycex <- c( 3,  2, 1);
      mypch <- c(15, 19, 1);
      ok <- rep(TRUE,Mtop);
      for(i in 1:ng) {
        hpd  <- seq(lohi[i,"lo"],lohi[i,"hi"]);
        yes  <- ok & (1:Mtop) %in% hpd;
        points(M[yes], pmf[yes], col=mycol[i], pch=mypch[i], cex=mycex[i]);
        ok[hpd] <- F;
      }
      points(M[ok],pmf[ok]); # The rest of them
    }
    legend(1.01*xlim[2], 1.10*max(pmf), xjust=1, yjust=1,
           c(expression(""),
             bquote(C[i] == .(C)),
             bquote(hat(M[i]) == .(round(C/Rstar,2))),
             bquote(bar(M[i]) == .(round(mu, 2))),
             bquote(I[.(100*gam[1])] == group("[",
                                              list(.(M[lohi[1,"lo"]]),.(M[lohi[1,"hi"]])),"]")),
             bquote(I[.(100*gam[2])] == group("[",
                                              list(.(M[lohi[2,"lo"]]),.(M[lohi[2,"hi"]])),"]"))),
           pch=c(NA,NA,NA,NA, mypch[1], mypch[2]),
           col=c(1,1,1,1,mycol[1],mycol[2]),
           text.col=c(1,1,1,1,mycol[1],mycol[2]),
           cex=1.75, bty="n");
    if(nchar(ps)>0) {
      dev.off();
      print(paste("Figure saved in file: `", ps, "'.",sep=""));
    }
    par(opar);
    lohi[,3:4] <- M[lohi[,3:4]];  # Value not index
    print(lohi)
  }
  invisible(C_table);
}
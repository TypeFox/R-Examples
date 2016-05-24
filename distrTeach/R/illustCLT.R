illustrateCLT <- function(Distr, len, sleep = 0){
     Distrname <- deparse(substitute(Distr))
     if(is.na(E(Distr)) || is.na(var(Distr)))
        stop(gettextf("Distribution %s does not have a variance/expectation.", 
           Distrname))
       graphics.off()
       Sn <- 0
       for(k in 1:len){
            o.warn <- getOption("warn"); options(warn = -1)
            on.exit(options(warn=o.warn))
            Sn <- Sn + Distr
            options(warn = o.warn)
            Tn <- make01(Sn)
            plotCLT(Tn,k, summands = Distrname)
            Sys.sleep(sleep)
           } 
}

illustrateCLT.tcl <- function(Distr, k, Distrname){
  if(is.na(E(Distr)) || is.na(var(Distr)))
          stop(gettextf("Distribution %s does not have a variance/expectation.", 
               Distrname ))
  if(is(Distr,"LatticeDistribution")||is(Distr,"AbscontDistribution"))
     Sn <- convpow(Distr,k)
  else {
     Sn <- 0
     o.warn <- getOption("warn"); options(warn = -1)
     on.exit(options(warn=o.warn))
     for(j in 1:k)
         Sn <- Sn + Distr
     options(warn=o.warn)
  }
  Tn <- make01(Sn)
  plotCLT(Tn,k, summands = Distrname)
  }


              ## graphical output

setMethod("plotCLT","DiscreteDistribution", function(Tn, k, summands = "") {
                N <- Norm()
                supp <- support(Tn)
                supp <- supp[supp >= -5]
                supp <- supp[supp <= 5]
                x <- seq(-5,5,0.01)
                dTn <- d(Tn)(supp)
                ymax <- max(1/sqrt(2*pi), dTn)
                opar <- par(no.readonly = TRUE)
    #            opar$cin <- opar$cra <- opar$csi <- opar$cxy <-  opar$din <- NULL
                on.exit(par(opar))
                dw <- min(diff(supp)) 
                facD <- min(dw*2,1)
                thin <- FALSE
                supp2 <- supp
                dTn2 <- dTn
                if (length(supp)>50){ 
                      nn <- seq(1,length(supp), by = length(supp) %/% 50)
                      thin <- TRUE
                      supp2 <- supp[nn]
                      dTn2 <- dTn[nn]
                }else nn <- seq(length(supp))
                par(mfrow = c(1,2), mar = c(5.1,4.1,4,2.1))
                plot(supp2, dTn2 / max(dTn2) * ymax, ylim = c(0, ymax), 
                     type = "h", ylab = gettext("densities/probabilities"), 
                     main = "", lwd = 4 * facD, xlab = "x")
                points(supp2, dTn2 / max(dTn2) * ymax, pch =16, 
                       cex = 2 * facD)      
                title(if(summands=="") gettextf("Sample size %i", k) else 
                      gettextf("Summands: %s", summands), line = 0.6)
                if (thin) text(supp[1],0.3, adj=0, "[grid thinned\n out]", 
                               cex = 0.7) 
                lines(x, d(N)(x), col = "orange", lwd = 2)
                mtext(gettext("Distribution of Standardized Centered Sums"),
                      cex=2,adj=.5,outer=TRUE,padj=1.4)
                pn <- p(Tn)(supp)
                plot(stepfun(supp, c(0,pn)), ylim = c(0, 1), 
                     ylab = gettext("cdfs"),
                     main = "", lwd = 4, cex.points = 2 * facD, pch = 16)
                points(supp, c(0,pn[-length(supp)]), pch = 21, cex = 2 * facD)      
                title(paste("Sample size", k),line=0.6)     
                lines(x, p(N)(x), col = "orange", lwd = 2)
                kd <- round(max(abs(pn-p(N)(supp))),4)
                text(1, 0.2, gettextf("Kolmogoroff-\nDistance:\n%1.4f",kd), 
                     adj = 0, cex = 0.8)
                legend("topleft", # x = supp[1], y = 1., 
                      legend = c(expression(italic(L)(frac(S[n]-E(S[n]), 
                                 sd(S[n])))), expression(italic(N)(0,1))), 
                       cex = .8, bty = "n", col = c("black", "orange"), 
                       lwd = c(4,2))
       })

setMethod("plotCLT","AbscontDistribution", function(Tn,k, summands = "") {
                N <- Norm()
                x <- seq(-5,5,0.01)
                dTn <- d(Tn)(x)
                ymax <- max(1/sqrt(2*pi), dTn)
                oldmar <- par("mar",no.readonly = TRUE)
                par(mfrow = c(1,2), mar = c(5.1,4.1,4.,2.1))
                plot(x, d(Tn)(x), ylim = c(0, ymax), type = "l", 
                     ylab = gettext("densities"), main ="", lwd = 4) 
                title(if(summands=="") gettextf("Sample size %i", k) else 
                     gettextf("Summands: %s", summands), line=0.6)
                lines(x, d(N)(x), col = "orange", lwd = 2)
                mtext(gettext("Distribution of Standardized Centered Sums"),
                      cex=2,adj=.5,outer=TRUE,padj=1.4)
                pn <- p(Tn)(x)
                pN <- p(N)(x)
                plot(x, pn, ylim = c(0, 1), type = "l", 
                     ylab = gettext("cdfs"), main="",
                     lwd = 4)
                title(paste("Sample size", k),line=0.6)     
                lines(x, pN, col = "orange", lwd = 2)
                kd <- round(max(abs(pn-pN)),4)
                text(1, 0.2, gettextf("Kolmogoroff-\nDistance:\n%1.4f",kd), 
                     adj = 0, cex = 0.8)
                legend("topleft", # x = -4.5, y = 1., 
                      legend = c(expression(italic(L)(frac(S[n]-E(S[n]), 
                                 sd(S[n])))), expression(italic(N)(0,1))), 
                       cex = .8, bty = "n", col = c("black", "orange"), 
                       lwd = c(4,2))
                par(mfrow = c(1,1), mar=oldmar)
       })



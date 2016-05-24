.lmomcohash <- new.env(hash=TRUE)


################# RICE DISTRIBUTION ######################
# START
"lmomrice2" <-
function(para, ...) {
    V   <- para$para[1]
    A   <- para$para[2]
    if(V == 0) {
      ray <- vec2par(c(0,A), type="ray")
      lmr <- lmomray(para=ray)
      lmr$source <- "lmomrice"
      return(lmr)
    }
    lmr <- theoLmoms.max.ostat(para=para,
                cdf=cdfrice, pdf=pdfrice,
                lower=0,
                upper=.Machine$double.max, ...)
    lmr$source <- "lmomrice"
    if(! are.lmom.valid(lmr)) {
      warning("The Rician parameters are producing invalid L-moments or L-moments outside of implementation of Rice distribution in lmomco")
      print(para)
      print(lmr)
      return()
    }
    return(lmr)
}

"Lhalf" <- function(x) {
   a <- x/2; I0 <- besselI(-a, n=0); I1 <- besselI(-a, n=1)
   return(exp(x/2)*((1-x)*I0 - x*I1))
}

# One should comment out the SNR > 24 test in pdfrice, cdfrice, quarice
# prior to running through the code herein. What you will see is the numerical
# problems for SNR > 24. The Normal distribution is switched to on the fly
# for Rice parameters having SNR > 24 to provide full support. As SNR gets big,
# the sqrt(pi/2)*Laguerre polynomial == SNR.

"RiceCooker" <- function(v, beg=NULL, end=NULL,
                         factor=5, step=NULL, nmom=4, tag=0) {
  if(is.null(step)) step <- v/100;
  if(is.null(beg))  beg <- v/100
  if(is.null(end))  end <- factor*v
  SNR  <- seq(beg,end, by=step);
  lcvs <- vector(mode="numeric")
  zetaLhalf <- zeta <- thevs <- tags <- lcvs
  t3 <- t4 <- lcvs
  j <- 0
  for(snr in SNR) {
    j <- j + 1
    mypara <- vec2par(c(v,v/snr), type="rice")
    lmr <- lmomrice2(mypara, nmom=nmom)
    if(is.null(lmr)) break
    lcv <- lmr$ratios[2]
    if(lcv <= 0) break
      tags[j] <- tag
     lcvs[j]  <- lcv
    thevs[j]  <- v
     zeta[j]  <- snr
       t3[j]  <- lmr$ratios[3]
       t4[j]  <- lmr$ratios[4]
    zetaLhalf[j] <- sqrt(pi/2)*Lhalf(-zeta[j]^2/2)
    cat(c(tags[j], lcvs[j], zeta[j], zetaLhalf[j], "\n"))
  }
  z <- list(tag=tags, LCV=lcvs, F.of.LCV = zeta,
            G.of.LCV=zetaLhalf, t3=t3, t4=t4)
  return(z)
}


try1 <- RiceCooker(1,   beg=0.035, end=1.2, step=0.02, nmom=4, tag=1)
try2 <- RiceCooker(10,  beg=1.1, end=6.1, step=0.05,  nmom=4, tag=2)
try3 <- RiceCooker(100, beg=5.1, end=25, step=0.1,  nmom=4, tag=3)

LCV      <- c(try1$LCV, try2$LCV, try3$LCV)
SNR      <- c(try1$F.of.LCV, try2$F.of.LCV, try3$F.of.LCV)
G.of.LCV <- c(try1$G.of.LCV, try2$G.of.LCV, try3$G.of.LCV)
T3       <- c(try1$t3, try2$t3, try3$t3)
T4       <- c(try1$t4, try2$t4, try3$t4)

# Now from the simulations, add in the Rayleigh end point
LCVray <- (.5*(sqrt(2)-1)*sqrt(pi))/(sqrt(pi/2))
raylmr <- lmomray(vec2par(c(0,1), type="ray"))
norlmr <- lmomnor(vec2par(c(0,1), type="nor"))
norT3 <- norlmr$TAU3
norT4 <- norlmr$TAU4
LCV      <- c(LCV, ray)
SNR      <- c(SNR, 0)
G.of.LCV <- c(G.of.LCV, sqrt(pi/2)*1)
T3       <- c(T3, raylmr$TAU3)
T4       <- c(T4, raylmr$TAU4)

xlim <- range(LCV)
ylim <- range(SNR)
plot(try1$LCV, try1$F.of.LCV,
     ylim=ylim, xlim=xlim,lwd=5,col=2,type="l")
lines(try2$LCV, try2$F.of.LCV, lwd=3,col=3)
lines(try3$LCV, try3$F.of.LCV, lwd=1,col=4)


ylim <- range(G.of.LCV)
plot(try1$LCV, try1$G.of.LCV,
     ylim=ylim, xlim=xlim,lwd=5,col=2,type="l")
lines(try2$LCV, try2$G.of.LCV, lwd=3,col=3)
lines(try3$LCV, try3$G.of.LCV, lwd=1,col=4)



RiceNomo <- data.frame(LCV=LCV, SNR=SNR, G=G.of.LCV, TAU3=T3, TAU4=T4)
row.names(RiceNomo) <- NULL
idx <- order(RiceNomo$LCV)
RiceNomo <- RiceNomo[idx,]
RiceNomo <- RiceNomo[RiceNomo$LCV <= ray,] # insure numerically <= Rayleigh
with(RiceNomo, plot(TAU3,TAU4))
RiceNomo <- RiceNomo[RiceNomo$TAU3 >= norT3,] # insure numerically >= Normal
with(RiceNomo, plot(TAU3,TAU4))
RiceNomo <- RiceNomo[RiceNomo$TAU4 <= norT4,] # insure numerically <= Normal
with(RiceNomo, plot(TAU3,TAU4)) # inspect plot!
RiceNomo <- RiceNomo[RiceNomo$G <= 24,] # now visually we have numeric problems
# still, with the TAU3,TAU4 diagram, looks like peak SNR is about 24
with(RiceNomo, plot(TAU3,TAU4)) # inspect plot!


LCVest <- seq(as.numeric(sprintf("%0.4f",min(RiceNomo$LCV))),
                                         ray, by=.0001)
LCVest <- c(LCVest,ray)
SNRest <- approx(RiceNomo$LCV, RiceNomo$SNR,  xout=LCVest, rule=2)$y
Gest   <- approx(RiceNomo$LCV, RiceNomo$G,    xout=LCVest, rule=2)$y
T3est  <- approx(RiceNomo$LCV, RiceNomo$TAU3, xout=LCVest, rule=2)$y
T4est  <- approx(RiceNomo$LCV, RiceNomo$TAU4, xout=LCVest, rule=2)$y
RiceNomoEst <- data.frame(LCV=LCVest, SNR=SNRest,
                          G=Gest, TAU3=T3est, TAU4=T4est)
with(RiceNomoEst, plot(TAU3,TAU4)) # inspect plot!
with(RiceNomoEst, plot(LCV,SNR))
with(RiceNomoEst, plot(LCV,G))

#T3seq <- seq(norT3, raylmr$TAU3,by=0.001)
#T4seq <- approx(T3,T4, xout=T3seq)$y
#T3seq <- c(T3seq, raylmr$TAU3)
#T4seq <- c(T4seq, raylmr$TAU4)
#.RiceT3T4 <- data.frame(TAU3=T3seq, TAU4=T4seq)
#assign("RiceT3T4", .RiceT3T4, .lmomcohash)

assign("RiceTable", RiceNomoEst, .lmomcohash)

assign("RiceTable.minLCV", min(RiceNomoEst$LCV), .lmomcohash)
assign("RiceTable.maxLCV", max(RiceNomoEst$LCV), .lmomcohash)
# END
################# RICE DISTRIBUTION ######################

save(.lmomcohash, file="sysdata.rda");


################# AEP DISTRIBUTION #######################
load(file.choose()) # go find the sysdata.rda

library(lmomco);

L1 <- L2 <- T3 <- T4 <- T5 <- vector(mode="numeric");
L1t <- L2t <- T3t <- T4t <- T5t <- L1;
Ur <- Ar <- Kr <- Hr <- L1;

i <- 0; step <- 0.05;
the.Kseq <- c(seq(-3,   0,   by=step),
              seq(step, 3.0, by=step));
the.Hseq <- c(seq(-3,   0,   by=step),
              seq(step, 3.0, by=step));
for(K in the.Kseq) { # outer loop on kappa parameter
   for(H in the.Hseq) { # inner loop on h parameter
      theK <- exp(K); theH <- exp(H);
      i <- i + 1
      cat("# K=",round(theK, digits=4),"  H=",round(theH, digits=4), "\n");
      PAR <- vec2par(c(0,1,theK,theH), type="aep4");

      lmr1 <- lmomaep4(PAR); # Delicado and Goria (2008)

      # numerical integration
      lmr2 <- theoLmoms(PAR, minF=1e-6, maxF=1-1e-6);

      if(! are.lmom.valid(lmr1)) next;
      if(! are.lmom.valid(lmr2)) next;

      T3[i]  <- lmr1$ratios[3];
      T4[i]  <- lmr1$ratios[4];
      T3t[i] <- lmr2$ratios[3];
      T4t[i] <- lmr2$ratios[4];
      T5t[i] <- lmr2$ratios[5];
      Kr[i]   <-    PAR$para[3];
      Hr[i]   <-    PAR$para[4];
  }
}

AEPD_kh2lmrTheo <- data.frame(K=Kr, H=Hr, T3=T3,  T4=T4,
                              T3t=T3t, T4t=T4t, T5t=T5t);

write.table(AEPD_kh2lmrTheo, file="AEPD_kh2lmrTheo.txt",
            quote=FALSE, row.names=FALSE);

save(AEPD_kh2lmrTheo, file="AEPD_kh2lmrTheo.RData");

assign("AEPkh2lmrTable", AEPD_kh2lmrTheo, .lmomcohash)
################# AEP DISTRIBUTION #######################

save(.lmomcohash, file="sysdata.rda");



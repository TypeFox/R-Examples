"lmrdia" <-
function() {
   step <- 0.005
   n <- 1
   aep4 <- matrix(nrow = 401, ncol = 2)
   for(t3 in seq(-1,1,step)) {
     a <- abs(t3)
     co <- c(0.7755464,  -3.3354852,  14.1955782, -29.9090294,
             37.2141451, -24.7411869,   6.7997646)
     T4.lowerbounds <- co[1]*a   + co[2]*a^2 + co[3]*a^3 +
                       co[4]*a^4 + co[5]*a^5 + co[6]*a^6 + co[7]*a^7
     aep4[n,1] <- t3
     aep4[n,2] <- T4.lowerbounds
     n <- n + 1
   }
   n <- 1
   lim <- matrix(nrow = 401, ncol = 2)
   gpa <- matrix(nrow = 401, ncol = 2)
   for(t3 in seq(-1,1,step)) {
     lim[n,1] <- t3
     lim[n,2] <- 0.25*(5*t3^2 - 1)
     gpa[n,1] <- t3
     gpa[n,2] <- (t3*(1+5*t3))/(5+t3)
     n <- n + 1
   }
   rgpa <- gpa; rgpa[,1] <- -rgpa[,1]

   n <- 1
   gev <- matrix(nrow = 582, ncol = 2)
   for(k in seq(-1,1,step)) {
     h <- -k
     gev[n,1] <- 2*(1-3^h)/(1-2^h) - 3
     gev[n,2] <- (5*(1-4^h)-10*(1-3^h)+6*(1-2^h))/(1-2^h)
     n <- n + 1
   }
   for(k in seq(1-step,10,0.05)) {
     h <- -k
     gev[n,1] <- 2*(1-3^h)/(1-2^h) - 3
     gev[n,2] <- (5*(1-4^h)-10*(1-3^h)+6*(1-2^h))/(1-2^h)
     n <- n + 1
   }
   gev <- gev[! is.nan(gev[,1]), ]
   wei <- gev; wei[,1] <- -wei[,1]

   n <- 1
   glo <- matrix(nrow = 401, ncol = 2)
   for(k in seq(-1,1,step)) {
     glo[n,1] <- -k
     glo[n,2] <- (1+5*k^2)/6
     n <- n + 1
   }

   n <- 1
   pIII <- matrix(nrow = 361, ncol = 2)
   for(t3 in seq(-.9,.9,step)) {
     pIII[n,1] <- t3
     pIII[n,2] <- 0.1224+0.30115*t3^2+0.95812*t3^4-0.57488*t3^6+0.19383*t3^8
     n <- n + 1
   }

   n <- 1
   ln <- matrix(nrow = 361, ncol = 2)
   for(t3 in seq(-.9,.9,step)) {
     ln[n,1] <- t3
     ln[n,2] <- 0.12282+0.77518*t3^2+0.12279*t3^4-0.13638*t3^6+0.11368*t3^8
     n <- n + 1
   }

   n <- 1
   t3s <- seq(-.5,1,0.005)
   gov <- matrix(nrow = length(t3s), ncol = 2)
   for(t3 in t3s) {
     gov[n,1] <- t3
     gov[n,2] <- -0.071430227 - 0.153032772*t3^1 + 1.049617013*t3^2 + 0.149620402*t3^3 + 0.021311828*t3^4 + 0.003912343*t3^5
     n <- n + 1
   }
   rgov <- gov; rgov[,1] <- -rgov[,1]

   exp <- matrix(nrow = 1, ncol = 2)
   exp[1,] <- c(1/3,1/6)
   gum <- matrix(nrow = 1, ncol = 2)
   gum[1,] <- c(log(9/8)/log(2),(16*log(2)-10*log(3))/log(2))
   nor <- matrix(nrow = 1, ncol = 2)
   nor[1,] <- c(0,30*pi^-1*atan(sqrt(2))-9)
   uni <- matrix(nrow = 1, ncol = 2)
   uni[1,] <- c(0,0)
   ray <- matrix(nrow = 1, ncol = 2)
   ray[1,] <- c(0.1139671, 0.1053695)
   cau <- matrix(nrow = 1, ncol = 2)
   cau[1,] <- c(0, 1)
   sla <- matrix(nrow = 1, ncol = 2)
   sla[1,] <- c(0, 0.3042045)
   z <- list(limits=lim, aep4=aep4, cau=cau, exp=exp,
             gev=gev, glo=glo, gpa=gpa, gum=gum, gno=ln, gov=gov,
             nor=nor, pe3=pIII, ray=ray, rgov=rgov, rgpa=rgpa,
             slash=sla, uniform=uni, wei=wei)
   return(z)
}


`modified.cvine.alg` <-
function(u1, u2, mu.x, phi.x, omega.x, mu.y, phi.y, omega.y, psi.x, psi.y,
                      rho.target, conv=0.001, margin.x, margin.y)
{
   N <- length(u1)
   rhoc <- rho.target
   lower <- -1
   upper <- 1
   lower.m <- -1
   upper.m <- 1

   #Burnin
   u1.in <- u1
   u1 <- h.inv(u1.in, u2, rhoc)
   u <- cbind(u1, u2)

   if (margin.x=="Poi")  x <- qpois(u[,1], mu.x)
   if (margin.x=="GP")   x <- pseudoinv.zigp(u[,1], mu.x, phi.x)
   if (margin.x=="ZIP")  x <- pseudoinv.zigp(u[,1], mu.x, 1, omega.x)
   if (margin.x=="ZIGP") x <- pseudoinv.zigp(u[,1], mu.x, phi.x, omega.x)
   if (margin.x=="NB")   x <- qnbinom(u[,1], mu=mu.x, size=psi.x)
   if (margin.y=="Poi")  y <- qpois(u[,2], mu.y)
   if (margin.y=="GP")   y <- pseudoinv.zigp(u[,2], mu.y, phi.y)
   if (margin.y=="ZIP")  y <- pseudoinv.zigp(u[,2], mu.y, 1, omega.y)
   if (margin.y=="ZIGP") y <- pseudoinv.zigp(u[,2], mu.y, phi.y, omega.y)
   if (margin.y=="NB")   y <- qnbinom(u[,2], mu=mu.y, size=psi.y)

   rho.m.alt <- cor(x,y)

   if (rho.m.alt > rho.target) { upper <- rhoc
                                 upper.m <- rho.m.alt }
   if (rho.m.alt < rho.target) { lower <- rhoc
                                 lower.m <- rho.m.alt }
   rhoc <- lower+(upper-lower)/2

   # Bisection
   rho.m.neu <- rho.m.alt
   konvergiert <- FALSE
   fehler <- FALSE
   notabbruch <- 0
   rausausfehlerhandling <- FALSE

   while(konvergiert==F & fehler==F) {
     while(abs(rho.m.neu - rho.target)>conv & fehler==F & rausausfehlerhandling==F) {
       notabbruch <- notabbruch + 1
       u1 <- h.inv(u1.in, u2, rhoc)
       u <- cbind(u1, u2)

       if (margin.x=="Poi")  x <- qpois(u[,1], mu.x)
       if (margin.x=="GP")   x <- pseudoinv.zigp(u[,1], mu.x, phi.x)
       if (margin.x=="ZIP")  x <- pseudoinv.zigp(u[,1], mu.x, 1, omega.x)
       if (margin.x=="ZIGP") x <- pseudoinv.zigp(u[,1], mu.x, phi.x, omega.x)
       if (margin.x=="NB")   x <- qnbinom(u[,1], mu=mu.x, size=psi.x)
       if (margin.y=="Poi")  y <- qpois(u[,2], mu.y)
       if (margin.y=="GP")   y <- pseudoinv.zigp(u[,2], mu.y, phi.y)
       if (margin.y=="ZIP")  y <- pseudoinv.zigp(u[,2], mu.y, 1, omega.y)
       if (margin.y=="ZIGP") y <- pseudoinv.zigp(u[,2], mu.y, phi.y, omega.y)
       if (margin.y=="NB")   y <- qnbinom(u[,2], mu=mu.y, size=psi.y)

       rho.m.neu <- cor(x,y)

       if (rho.m.neu > rho.target) { upper <- rhoc
                                     upper.m <- rho.m.neu }
       if (rho.m.neu < rho.target) { lower <- rhoc
                                     lower.m <- rho.m.neu }
       rhoc <- lower+(upper-lower)/2
       if (rho.target < min(lower.m, upper.m) | rho.target > max(lower.m, upper.m) | rho.m.neu==rho.m.alt | notabbruch >= 500) {
         fehler <- TRUE
         #print("fehler")
         if (notabbruch >= 500) {print(">500")}
         if (rho.m.neu==rho.m.alt) {
           #print("rho.m.neu==rho.m.alt")
           fehler <- FALSE
           konvergiert <- TRUE
           rausausfehlerhandling <- TRUE
        }
         #print(paste("rho.m.neu",rho.m.neu))
         #print(paste("rhoc",rhoc))
       }
       rho.m.alt <- rho.m.neu
     }
     if (abs(rho.m.neu - rho.target)<conv) { konvergiert <- TRUE }
   }

     return(list(u1=u1,rho.m.hat=rho.m.neu, rhoc=rhoc,fehler=fehler))
}


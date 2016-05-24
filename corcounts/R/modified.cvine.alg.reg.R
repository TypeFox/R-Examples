`modified.cvine.alg.reg` <-
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

   x <- double(N)
   y <- double(N)
   for (k in 1:N) {
     if (margin.x=="Poi")  x[k] <- qpois(u[k,1], mu.x[k])
     if (margin.x=="GP")   x[k] <- pseudoinv.zigp(u[k,1], mu.x[k], phi.x[k])
     if (margin.x=="ZIP")  x[k] <- pseudoinv.zigp(u[k,1], mu.x[k], 1, omega.x[k])
     if (margin.x=="ZIGP") x[k] <- pseudoinv.zigp(u[k,1], mu.x[k], phi.x[k], omega.x[k])
     if (margin.x=="NB")   x[k] <- qnbinom(u[k,1], mu=mu.x[k], size=psi.x[k])
     if (margin.y=="Poi")  y[k] <- qpois(u[k,2], mu.y[k])
     if (margin.y=="GP")   y[k] <- pseudoinv.zigp(u[k,2], mu.y[k], phi.y[k])
     if (margin.y=="ZIP")  y[k] <- pseudoinv.zigp(u[k,2], mu.y[k], 1, omega.y[k])
     if (margin.y=="ZIGP") y[k] <- pseudoinv.zigp(u[k,2], mu.y[k], phi.y[k], omega.y[k])
     if (margin.y=="NB")   y[k] <- qnbinom(u[k,2], mu=mu.y[k], size=psi.y[k])
   }

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

       for (k in 1:N) {
         if (margin.x=="Poi")  x[k] <- qpois(u[k,1], mu.x[k])
         if (margin.x=="GP")   x[k] <- pseudoinv.zigp(u[k,1], mu.x[k], phi.x[k])
         if (margin.x=="ZIP")  x[k] <- pseudoinv.zigp(u[k,1], mu.x[k], 1, omega.x[k])
         if (margin.x=="ZIGP") x[k] <- pseudoinv.zigp(u[k,1], mu.x[k], phi.x[k], omega.x[k])
         if (margin.x=="NB")   x[k] <- qnbinom(u[k,1], mu=mu.x[k], size=psi.x[k])
         if (margin.y=="Poi")  y[k] <- qpois(u[k,2], mu.y[k])
         if (margin.y=="GP")   y[k] <- pseudoinv.zigp(u[k,2], mu.y[k], phi.y[k])
         if (margin.y=="ZIP")  y[k] <- pseudoinv.zigp(u[k,2], mu.y[k], 1, omega.y[k])
         if (margin.y=="ZIGP") y[k] <- pseudoinv.zigp(u[k,2], mu.y[k], phi.y[k], omega.y[k])
         if (margin.y=="NB")   y[k] <- qnbinom(u[k,2], mu=mu.y[k], size=psi.y[k])
       }

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


library(spatsurv)

testIterators <- list(
                      m1 = mcmcLoop(N=15,burnin=3,thin=4),
                      m2 = mcmcLoop(N=15,burnin=0,thin=4,trim=FALSE),
                      m3 = mcmcLoop(N=15,burnin=1,thin=1,trim=FALSE)
                      )

for(tim in testIterators){
  print(tim)
  summary(tim)
  resetLoop(tim)
  summary(tim)
}

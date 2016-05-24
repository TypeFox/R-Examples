"randfn" <-
function(n, family, ...){
  args <- list(...)
  switch(family,
negbin= rnbinom(n, size=exp(args$size), mu=args$mu),        
poisson=rpois(n, lambda=args$mu), 
geometric=rgeom(n, prob= args$mu), 
binom=rbinom(n, size=args$size, prob=args$mu)
         
         )
}


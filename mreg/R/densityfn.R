"densityfn" <-
function(x, family, ...){
  
args <- list(...)

switch( family, 
"negbin"=  dnbinom(x,size=exp(args$size), mu=args$mu) ,
"negbin.ncar"= dnbinom(x, size=exp(args$size[1]), mu=args$mu)
       *exp( (length(x)>1)*(args$size[2]+x* args$size[3]  ))  /(1+exp( args$size[2]+ x*args$size[3])),
"poisson"= dpois(x, lambda=args$mu),
"geometric"=dgeom(x, prob= args$mu),
#"binom"=dbinom(x, size=args$size, prob=args$mu)      
)


}


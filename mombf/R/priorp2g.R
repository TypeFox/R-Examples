###
### priorp2g.R
###

priorp2g <- function(priorp,
                     q,
                     nu=1,
                     prior=c('iMom', 'normalMom','tMom')) {
  prior <- match.arg(prior)
  switch(EXPR=prior,
         normalMom = {
           e <- function(logg) {
             return((1-2*pmom(-abs(q), tau=exp(logg)) - priorp[i])^2)
           }

           ans <- double(length(priorp))
           for (i in 1:length(priorp)) {
             ans[i] <- exp(nlminb(start=0, objective=e)$par)
           }
           ans
         },
         tMom = {
           stop("prior=='tMom' is not currently implemented")
         },
         iMom = {
           p <- (1-priorp)/2
           qgamma(2*p,nu/2,1)*q^2
         })
}


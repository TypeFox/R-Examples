###
### mode2g.R
###

mode2g <- function(prior.mode,
                   prior=c('iMom', 'normalMom', 'tMom'),
                   nu=1,
                   dim=1) {
  prior <- match.arg(prior)
  switch(EXPR=prior,
         normalMom = {
           prior.mode / 2
         },
         iMom = {
           prior.mode * (nu+dim) / 2
         },
         tMom = {
           if (nu<3) {
             stop('tMom prior must have nu>2 degrees of freedom')
           }
           prior.mode * (nu-2+dim) / (2*nu)
         })
}


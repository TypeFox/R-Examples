###
### g2mode.R
###

g2mode <- function(g,
                   prior=c('iMom', 'normalMom', 'tMom'),
                   nu=1,
                   dim=1) {
  prior <- match.arg(prior)
  switch(EXPR=prior,
         normalMom = {
           2*g
         },
         tMom = {
           g*2*nu/(nu-2+dim)
         },
         iMom = {
           2*g/(nu+dim)
         })
}


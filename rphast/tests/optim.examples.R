library(rphast)

#' optim.rphast
# set up a tree model object using HKY85 and kappa=3
mod <- tm("((human:0.103064,(mouse:0.075131,rat:0.078848):0.279629):0.101528,cow:0.101528)", subst.mod="HKY85", backgd=rep(0.25, 4))
mod <- set.rate.matrix.tm(mod, params=3)
#'
# simulate an msa with 1000 bp under this model
align <- simulate.msa(mod, 1000, as.pointer=TRUE)
#'
# make a function which takes parameters and a model and returns likelihood
mylikelihood <- function(params, init.mod, align) {
  # first parameter is scale
  init.mod$tree <- rescale.tree(init.mod$tree, params[1])
  # second parameter is kappa
  init.mod <- set.rate.matrix.tm(init.mod, params=params[2])
  likelihood.msa(align, init.mod)
}
#'
result <- optim.rphast(mylikelihood, c(1.0, 3.0), lower=c(0, 0), init.mod=mod, align=align)
#'
# This should be about the same as:
compMod <- phyloFit(align, init.mod=mod, scale.only=TRUE)
# This is the estimated scale from phyloFit
branchlength.tree(compMod$tree)/branchlength.tree(mod$tree)
# This is the estimated kappa from phyloFit
get.rate.matrix.params.tm(compMod)
#'
# It is also the same as:
optim(c(1.0, 3.0), mylikelihood, method="L-BFGS-B", lower=c(0, 0), control=list(fnscale=-1.0),
      init.mod=mod, align=align)


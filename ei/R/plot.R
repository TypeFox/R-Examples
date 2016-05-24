#@x - ei.object to plot
#Plots in ei.

plot.ei <- function(x, ...){
  ei.object <- x
  lci.function.list <- list("tomogD" = .tomog, "tomog"=.tomogl)    # Fix for passing lci value
  function.list <- list("tomogCI"=.tomog80CI, "tomogCI95"=.tomog95CI,
                        "tomogE"=.tomogE, "tomogP" = .tomogP2,
                        "betab"=.betabd,
                        "betaw"=.betawd, "xt"=.xt, "xtc"=.xtc,
                        "xtfit"=.xtfit, "xtfitg"=.xtfitg,
                        "estsims"=.estsims, "boundXb"=.boundXb,
                        "boundXw"=.boundXw,
                        "truth"=.truthfn,"eiRxCtomog"=.bndplot,
                        "movieD"=.movieD, "movie"=.movie)
  arguments <- list(...)
  # Fix for passing lci value
  if ("lci" %in% names(arguments)){
    lci<-arguments$lci
    arguments$lci<-NULL
  } else {
    lci<-TRUE
  }

  results <- list()
  if (length(arguments)!=1) {row = ceiling(length(arguments)/2)
                             par(mfrow=c(row, 2))
                           }
  for (arg in arguments) {
    if (arg %in% names(function.list)){
      results[[arg]] <- function.list[[arg]] (ei.object=ei.object)
    } else if (arg %in% names(lci.function.list)) {
      results[[arg]] <- lci.function.list[[arg]] (ei.object=ei.object, lci=lci)    # Fix for passing lci value
    } else
      results[[arg]] <- NA
  }
}

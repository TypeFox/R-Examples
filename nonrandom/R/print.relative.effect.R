
print.relative.effect <- function(x,
                                  ...)
{

  object <- x
  
  cat("\n Treatment:",object$name.treat)
  cat("\n Outcome:",object$name.resp)
  cat("\n Covariates: ",object$name.sel, "\n\n")
  
  cat("\n Unadjusted treatment effect: ", round(object$unadj.treat,4),"\n", sep="")
  
  cat("\n Adjusted and relative effects: \n\n")
  
  rel.eff.tab <- matrix(c(as.numeric(object$adj.treat.cov),
                          as.numeric(object$rel.eff.treat)),
                        nrow=length(object$name.sel),
                        ncol=2,
                        dimnames=list(object$name.sel,
                          c("adj. treatment effect", "rel. effect")))

  if (dim(rel.eff.tab)[1] != 1)
  
    print(rel.eff.tab[order(object$rel.eff.treat, decreasing=TRUE),])

  else

    print(format(rel.eff.tab))
}

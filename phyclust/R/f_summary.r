### This file contains functions for print and summary.

### Print phyclust
print.phyclust <- function(x, digits = max(4, getOption("digits") - 3), ...){
  ret.phyclust <- x

  op.org <- options()
  options(digits = digits)
  init <- NULL

  if(!is.null(ret.phyclust$conv)){
    my.cat("Phyclust Results:\n",
           "code type: ", ret.phyclust$code.type,
           ", em method: ", ret.phyclust$em.method,
           ", boundary method: ", ret.phyclust$boundary.metho, ".\n",
           "init procedure: ", ret.phyclust$init.procedure,
           ", method: ", ret.phyclust$init.method, ".\n",
           "model substitution: ", ret.phyclust$substitution.model,
           ", distance: ", ret.phyclust$edist.model, ".\n",
           "label method: ", ret.phyclust$label.method, ".\n",
           "iter: ", ret.phyclust$conv$iter,
           " ", ret.phyclust$conv$inner.iter,
           " ", ret.phyclust$conv$cm.iter,
           ", convergence: ", ret.phyclust$conv$flag,
           ", check.param: ", ret.phyclust$conv$check.param, ".\n",
           "eps: ", ret.phyclust$conv$eps,
           ", error: ", ret.phyclust$conv$error, ".\n")
  }
      
  my.cat("N.X.org: ", ret.phyclust$N.X.org,
         ", N.X: ", ret.phyclust$N.X,
         ", L: ", ret.phyclust$L,
         ", K: ", ret.phyclust$K,
         ", p: ", ret.phyclust$p,
         ", N.seg.site: ", ret.phyclust$N.seg.site, ".\n",
         "logL: ", ret.phyclust$logL)
  if(!is.null(ret.phyclust$bic)) my.cat(", bic: ", ret.phyclust$bic)
  if(!is.null(ret.phyclust$aic)) my.cat(", aic: ", ret.phyclust$aic)
  if(!is.null(ret.phyclust$icl)) my.cat(", icl: ", ret.phyclust$icl)
  my.cat("\n")
  my.cat("identifier: ", ret.phyclust$QA$identifier, "\n")
  cat("  Eta:", ret.phyclust$Eta, "\n")
  if(!is.null(ret.phyclust$Q$pi)){
    my.cat("  pi:\n")
    my.print(ret.phyclust$Q$pi)
  }
  if(!is.null(ret.phyclust$Q$kappa)) cat("  kappa:", ret.phyclust$Q$kappa, "\n")
  cat("  Tt:", ret.phyclust$Q$Tt, "\n")
  if(!is.null(ret.phyclust$n.class)) cat("  n.class:", ret.phyclust$n.class, "\n")
  if(!is.null(ret.phyclust$SE)){
    my.cat("SE_model: ", ret.phyclust$SE$model,
           ", SE_constant: ", ret.phyclust$SE$constant, "\n")
    my.print(ret.phyclust$SE$f.err)
  }
  options(op.org)
  invisible()
} # End of print.phyclust().

summary.phyclust <- function(object, ...){
  ret.phyclust <- object

  print.phyclust(ret.phyclust)
  cat("Mu:\n")
  for(i in 1:ret.phyclust$K){
    cat("   ", as.character(.nucleotide$code[ret.phyclust$Mu[i,] + 1]), "\n")
  }
  cat("Class id:", ret.phyclust$class.id, "\n")
} # End of summary.phyclust().


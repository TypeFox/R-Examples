# Define summary() function for class "PerFit":
summary.PerFit <- function(object, #object = an object from 'PerFit' class
                           cutoff.obj=NULL, #cutoff.obj = an object from 'PerFit.cutoff' class
                           ModelFit="NonParametric", Nreps=1000, 
                           IP=object$IP, IRT.PModel=object$IRT.PModel, Ability=object$Ability, Ability.PModel=object$Ability.PModel,
                           mu=0, sigma=1, 
                           Blvl = 0.05, Breps = 1000, CIlvl = 0.95, 
                           UDlvl = NA, ...)
{
  x <- object
  # Sanity check - Class PerFit:
  Sanity.cls(x)  
  
  # Compute cutoff:
  if (is.null(cutoff.obj))
  {
    cutoff.res <- cutoff(x, ModelFit, Nreps, IP, IRT.PModel, Ability, Ability.PModel, mu, sigma, Blvl, Breps, CIlvl, UDlvl)
  } else
  {
    Sanity.clsPO(cutoff.obj)
    cutoff.res <- cutoff.obj
  }
  
  # Compute flagged:
  flagged.res <- flagged.resp(x, cutoff.res, scores = FALSE)$PFSscores[, 1]
  
  # Summarize results:
  cat(paste0("\nPFS = ", x$PFStatistic, "\n"))
  cat(paste0("Cutoff = ", cutoff.res$Cutoff, " (SE = ", cutoff.res$Cutoff.SE, ").\n"))
  cat(paste0("Tail = ", cutoff.res$Tail, ".\n"))
  cat(paste0("Proportion of flagged respondents = ", cutoff.res$Prop.flagged, ".\n"))
  cat("(N.B.: The cutoff varies each time cutoff() is run due to bootstrapping.)\n\n")
  # 
  cat(paste0("Identified respondents - ", length(flagged.res), " in total:\n"))
  cat("   ", flagged.res, "\n\n")
}

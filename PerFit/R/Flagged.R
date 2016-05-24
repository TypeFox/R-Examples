flagged.resp <- function(x, #x = an object from 'PerFit' class
                         cutoff.obj=NULL, #cutoff.obj = an object from 'PerFit.cutoff' class
                         scores=TRUE, ord=TRUE,
                         ModelFit="NonParametric", Nreps=1000, 
                         IP=x$IP, IRT.PModel=x$IRT.PModel, Ability=x$Ability, Ability.PModel=x$Ability.PModel, mu=0, sigma=1, 
                         Blvl = 0.05, Breps = 1000, CIlvl = 0.95, 
                         UDlvl=NA)
{
  # Sanity check - Class PerFit:
  Sanity.cls(x)  
  # 
  upp.PFS  <- c("Cstar", "C.Sato", "U3", "ZU3", "G", "Gnormed", "Gpoly", "Gnormed.poly", "U3poly", "D.KB")
  low.PFS  <- c("r.pbis", "NCI", "Ht", "A.KB", "E.KB", "lz", "lzstar", "lzpoly")
  # 
  if (is.null(cutoff.obj))
  {
    cutoff.res <- cutoff(x, ModelFit, Nreps, IP, IRT.PModel, Ability, Ability.PModel, mu, sigma, Blvl, Breps, CIlvl, UDlvl)
  } else
  {
    Sanity.clsPO(cutoff.obj)
    cutoff.res <- cutoff.obj
  }
  #   
  if (any(x$PFStatistic == upp.PFS)) 
  {
    flagged.subs <- which(x$PFscores[, 1] >= cutoff.res$Cutoff)
  }
  if (any(x$PFStatistic == low.PFS)) 
  {
    flagged.subs <- which(x$PFscores[, 1] <= cutoff.res$Cutoff)
  }
  Ps <- round(colMeans(x$Matrix, na.rm = TRUE), 3)
  # Not ordered by pvalue:
  if (ord == FALSE)
  {
    flagged.scores           <- x$Matrix[flagged.subs, ]
    colnames(flagged.scores) <- paste("It", 1:dim(x$Matrix)[2], sep = "")
  }
  # Ordered by pvalue:
  if (ord == TRUE)
  {
    matrix.ord               <- x$Matrix[, order(Ps, decreasing = TRUE)] # ordered from easy to difficult
    flagged.scores           <- matrix.ord[flagged.subs, ]
    colnames(flagged.scores) <- paste("It", order(Ps, decreasing = TRUE), sep = "")
    Ps                       <- sort(Ps, decreasing = TRUE)
  }
  flagged.scores           <- as.matrix(flagged.scores)
  rownames(flagged.scores) <- NULL
  res                      <- if (scores == FALSE)
  {
    list(PFSscores  = cbind(FlaggedID = flagged.subs, PFscores = x$PFscores[flagged.subs, 1]), 
         Cutoff.lst = cutoff.res, PFS = x[[2]])
  } else
  {
    list(Scores = cbind(FlaggedID = flagged.subs, flagged.scores, PFscores = x$PFscores[flagged.subs, 1]), 
         MeanItemValue = Ps, Cutoff.lst = cutoff.res, 
         PFS = x[[2]])
  }
  return(res)
}
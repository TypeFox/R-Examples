readLisrel <- function(x)
{
  
  # Read output:
  Out <- readLines(x)
  
  # Empty output structure (S3):
  Res <- list(
    fitIndices = list(), # list containing dataframes, 'global' and 'Group1'. 'Group2',...
    matrices = list( ), # List containing per matrix a list contaning 'est', 'par', 'se', 't' matrices
    variables = data.frame(), # data frame contaning information on each variable: manifest or latent and exo or endo.
    covariances = list(
      implied = NULL, # model Implied covariance matrix
      observed = NULL # observed covariance matrix
    ),
    output = Out
  )
  class(Res) <- "lisrel"
  class(Res$matrices) <- c("list","semMatrixModel")
  
  ### Find linenumbers of output structure:
  Struc <- list(
    parSpec = grep("Parameter Specifications",Out),
    est = grep("LISREL Estimates",Out),
    fit = grep("Goodness of Fit Statistics",Out),
    modInd = grep("Modification Indices and Expected Change",Out),
    std = grep("Standardized Solution",Out),
    stdComp = grep("Completely Standardized Solution",Out),
    stdcommet = grep("Common Metric Standardized Solution",Out),
    stdcommetComp = grep("Common Metric Completely Standardized Solution",Out)
  )
  
  Struc$std <- Struc$std[!Struc$std%in%c(Struc$stdComp,Struc$stdcommet,Struc$stdcommetComp)]
  Struc$stdComp <- Struc$stdComp[!Struc$stdComp%in%Struc$stdcommetComp]
  
  StrucUL <- unlist(Struc)
  
  ### Find linenumbers of matrices:
  Mats <- list(
    ObsCov = grep("^\\s*Covariance Matrix\\s*$",Out),
    ImpCov = grep("^\\s*Fitted Covariance Matrix\\s*$",Out),
    LX = grep("LAMBDA-X",Out),
    PH = grep("PHI",Out),
    TD = grep("THETA-DELTA",Out),
    GA = grep("GAMMA",Out),
    LY = grep("LAMBDA-Y",Out),
    PS = grep("PSI",Out),
    TE = grep("THETA-EPS",Out),
    BE = grep("BETA",Out),
    TX = grep("TAU-X",Out),
    TY = grep("TAU-Y",Out),
    AL = grep("ALPHA",Out),
    KA = grep("KAPPA",Out)
  )
  
  Mats$ObsCov <- Mats$ObsCov[!Mats$ObsCov%in%Mats$ImpCov]

  # Check for continued matrices:
  for (i in 1:length(Mats))
  {
    if (length(Mats[[i]])>0)
    {
      Inds <- lapply(Mats[[i]],matRange,Out=Out)
      Mats[[i]] <- Mats[[i]][!sapply(Mats[[i]],function(x)any(x>sapply(Inds,'[',1) & x<=sapply(Inds,'[',2)))]
    }
  }
  
  # Extract number of groups from parSpec matrices:
  Ng <- length(Struc$parSpec)
  
  ### EXTRACT MATRICES ###
  for (mat in c("LX","PH","TD","GA","TX","KA","LY","PS","TE","BE","TY","AL"))
  {
    Res$matrices[[mat]] <- list()
    for (g in Ng:1)
    {
      Res$matrices[[mat]][[g]] <- list()
      for (type in c("est","std","stdComp","stdcommet","stdcommetComp","parSpec"))
      {
        if (length(Struc[[type]])>0)
        {
          if (length(Struc[[type]])!=Ng) warning(paste("Number of",type,"matrices are not equal to number of parameter specification matrices. Resulting multigroup structure will NOT be reliable"))
          Res$matrices[[mat]][[g]][[type]] <- findMatrix(mat,type,Mats,Struc,Out,g)
          if (identical(Res$matrices[[mat]][[g]][[type]],"NextGroup"))
          {
            Res$matrices[[mat]][[g]][[type]] <- Res$matrices[[mat]][[g+1]][[type]]
            if (type=="est")
            {
              Res$matrices[[mat]][[g]][['se']] <- Res$matrices[[mat]][[g+1]][['se']]
              Res$matrices[[mat]][[g]][['t']] <- Res$matrices[[mat]][[g+1]][['t']]
            }
          } else if (type=="est" & is.list(Res$matrices[[mat]][[g]][[type]]))
          {
            Res$matrices[[mat]][[g]][['se']] <- Res$matrices[[mat]][[g]][[type]][['se']]
            Res$matrices[[mat]][[g]][['t']] <- Res$matrices[[mat]][[g]][[type]][['t']]
            Res$matrices[[mat]][[g]][[type]] <- Res$matrices[[mat]][[g]][[type]][['est']]
          }
        }
      }
    }
    names(Res$matrices[[mat]]) <- paste("group",seq_along(Res$matrices[[mat]]),sep="")
  }  

  Res$covariances$implied <- findCov("ImpCov",Mats,Out)
  if (length(Res$covariances$implied)>0)  names(Res$covariances$implied) <- paste("group",seq_along(Res$covariances$implied),sep="")
  Res$covariances$observed <- findCov("ObsCov",Mats,Out)
  if (length(Res$covariances$observed)>0) names(Res$covariances$observed) <- paste("group",seq_along(Res$covariances$observed),sep="")
  
#   
#   if (length(Struc$fit)==1)
#   {
#     Res$fitIndices$global <- findFit(Out,Struc$fit)
#   } else {
#     Res$fitIndices$global <- findFit(Out,Struc$fit[grep("Global",Out[Struc$fit])])
#     for (g in 1:(length(Struc$fit)-1))
#     {
#       Res$fitIndices[[paste("group",g,sep="")]] <- findFit(Out,Struc$fit[grep("Group",Out[Struc$fit])[g]])
#     }
#   }
#   
  # Return:
  return(Res)
}

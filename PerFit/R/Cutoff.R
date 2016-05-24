# Determine the cutoff for a PFS based on model-fitting item response vectors:
cutoff <- function(x, #x = an object from 'PerFit' class
                   ModelFit="NonParametric", Nreps=1000, 
                   IP=x$IP, IRT.PModel=x$IRT.PModel, Ability=x$Ability, Ability.PModel=x$Ability.PModel, mu=0, sigma=1, 
                   Blvl = 0.05, Breps = 1000, CIlvl = 0.95, 
                   UDlvl = NA)
{
  # Sanity check - Class PerFit:
  Sanity.cls(x)  
  # 
  N        <- dim(x$Matrix)[1]; I <- dim(x$Matrix)[2]; Ncat <- x$Ncat
  upp.PFS  <- c("Cstar", "C.Sato", "U3", "ZU3", "G", "Gnormed", "Gpoly", "Gnormed.poly", "U3poly", "D.KB")
  low.PFS  <- c("r.pbis", "NCI", "Ht", "A.KB", "E.KB", "lz", "lzstar", "lzpoly")
  dico.PFS <- c("Cstar", "C.Sato", "U3", "ZU3", "G", "Gnormed", "D.KB", "r.pbis", "NCI", "Ht", "A.KB", "E.KB", "lz", "lzstar")
  poly.PFS <- c("Gpoly", "Gnormed.poly", "U3poly", "lzpoly")
  perfectAllowed.PFS <- c("G", "Gnormed", "Gpoly", "Gnormed.poly", "U3poly", "NCI", "lz", "lzstar", "lzpoly")
  perfectNotAllowed.PFS <- c("Cstar", "C.Sato", "U3", "ZU3", "A.KB", "D.KB", "E.KB", "r.pbis", "Ht")
  # 
  pf.scores       <- x$PFscores[[1]]
  PFS.NA          <- is.na(pf.scores)
  pf.scores.noNAs <- pf.scores[!PFS.NA]
  
  if (ModelFit == "Parametric" & any(x$PFStatistic == dico.PFS))
  {
    # Generate model-fitting item score vectors based on a parametric model (1PL, 2PL, or 3PL).
    # 
    # Estimate item parameters if not provided:
    if (is.null(IRT.PModel)) 
    {
      IRT.PModel <- "2PL"
    } else
    {
      Sanity.IPm(IRT.PModel)
    }
    IP <- estIP(x$Matrix, IP, IRT.PModel)
    # Estimate ability parameters if not provided:
    if (is.null(Ability.PModel)) 
    {
      Ability.PModel <- "ML"
    } else
    {
      Sanity.Abm(Ability.PModel)
    }
    Ability <- estAb(x$Matrix, IP, Ability, Ability.PModel, mu, sigma)
    #
    Ability.gen <- sample(Ability[!is.na(Ability)], size = Nreps, replace = TRUE)
    #
    A   <- IP[, 1]; B <- IP[, 2]; C <- IP[, 3]
    P                   <- do.call(cbind, lapply(1:I, function (x) {C[x] + (1-C[x]) / (1 + exp(-A[x] * (Ability.gen - B[x])))}))
    matrix.modelfitting <- matrix(rbinom(length(P), 1, P), ncol = I)
    # Compute PFS values on the model-fitting item response vectors:
    if (any(x$PFStatistic == perfectNotAllowed.PFS))
    {
      nonperfectvects.mf <- (rowSums(matrix.modelfitting) %% I) != 0
      PFS.modelfitting <- eval(parse(text = x[[2]]))(matrix.modelfitting[nonperfectvects.mf,])$PFscores[[1]]
    } else
    {
      PFS.modelfitting <- eval(parse(text = x[[2]]))(matrix.modelfitting)$PFscores[[1]]
    }
  }
  
  if (ModelFit == "Parametric" & any(x$PFStatistic == poly.PFS))
  {
    # Generate model-fitting item score vectors based on a parametric model (PCM, GPCM, or GRM).
    # 
    # Estimate item parameters if not provided:
    if (is.null(IRT.PModel)) 
    {
      IRT.PModel <- "GRM"
    } else
    {
      Sanity.IPm.poly(IRT.PModel)
    }
    IP.res <- estIP.poly(x$Matrix, Ncat, IP, IRT.PModel)
    IP     <- IP.res[[1]]
    IP.ltm <- IP.res[[2]]
    # Estimate ability parameters if not provided:
    if (is.null(Ability.PModel)) 
    {
      Ability.PModel <- "EAP"
    } else
    {
      Sanity.Abm.poly(Ability.PModel)
    }
    Ability <- estAb.poly(x$Matrix, IP.ltm, Ability, Ability.PModel)
    #
    Ability.gen <- sample(Ability[!is.na(Ability)], size = Nreps, replace = TRUE)
    #
    P.CRF               <- estP.CRF(I, Ncat, IRT.PModel, IP, Ability.gen)
    P.CRF.ind           <- matrix(1:(Ncat * I), nrow = Ncat)
    matrix.modelfitting <- matrix(, nrow = Nreps, ncol = I)
    for (i in 1:I)
    {
      matrix.modelfitting[, i] <- rMultinom(P.CRF[, ((i - 1) * Ncat + 1) : (i * Ncat)], 1) - 1
    }
    # Compute PFS values on the model-fitting item response vectors:
    if (any(x$PFStatistic == perfectNotAllowed.PFS))
    {
      nonperfectvects.mf <- (rowSums(matrix.modelfitting) %% I) != 0
      PFS.modelfitting <- eval(parse(text = x[[2]]))(matrix.modelfitting[nonperfectvects.mf, ], Ncat)$PFscores[[1]]
    } else
    {
      PFS.modelfitting <- eval(parse(text = x[[2]]))(matrix.modelfitting, Ncat)$PFscores[[1]]
    }
  }
  
  # 
  if (ModelFit == "NonParametric" & any(x$PFStatistic == dico.PFS))
  {
    NC <- rowSums(x$Matrix, na.rm = TRUE)
    # 
    NC.gen <- sample(NC[!PFS.NA], size = Nreps, replace = TRUE)
    #
    uniqueNC            <- sort(unique(NC.gen))
    matrix.modelfitting <- matrix(, nrow = Nreps, ncol = I)
    # 
    for (i in 1:length(uniqueNC))
    {
      NC.i <- which(NC.gen == uniqueNC[i])
      pi.i <- colMeans(x$Matrix[NC == uniqueNC[i], , drop = FALSE], na.rm = TRUE)
      matrix.modelfitting[NC.i, ] <- rbinom(length(NC.i) * I, 1, rep(pi.i, each = length(NC.i)))
    }
    # Compute PFS values on the model-fitting item response vectors:
    if (any(x$PFStatistic == perfectNotAllowed.PFS))
    {
      nonperfectvects.mf <- (rowSums(matrix.modelfitting) %% I) != 0
      PFS.modelfitting <- eval(parse(text = x[[2]]))(matrix.modelfitting[nonperfectvects.mf,])$PFscores[[1]]
    } else
    {
      PFS.modelfitting <- eval(parse(text = x[[2]]))(matrix.modelfitting)$PFscores[[1]]
    }
  }
  
  # 
  if (ModelFit == "NonParametric" & any(x$PFStatistic == poly.PFS))
  {
    NC <- rowSums(x$Matrix, na.rm = TRUE)
    # 
    NC.gen <- sample(NC[!PFS.NA], size = Nreps, replace = TRUE)
    #
    uniqueNC            <- sort(unique(NC.gen))
    matrix.modelfitting <- matrix(, nrow = Nreps, ncol = I)
    # 
    for (i in 1:length(uniqueNC))
    {
      NC.i       <- which(NC.gen == uniqueNC[i])
      freq.abs.i <- apply(x$Matrix[NC == uniqueNC[i], , drop=FALSE], 2,
                          function(vec) {table(factor(vec, levels=0:(Ncat - 1)))})
      matrix.modelfitting[NC.i, ] <- 
        apply(freq.abs.i, 2, function(vect) 
        {
          which(rmultinom(length(NC.i),1,vect) == 1, arr.ind=TRUE)[,1] - 1
        })
    }
    # Compute PFS values on the model-fitting item response vectors:
    if (any(x$PFStatistic == perfectNotAllowed.PFS))
    {
      nonperfectvects.mf <- (rowSums(matrix.modelfitting) %% I) != 0
      PFS.modelfitting <- eval(parse(text = x[[2]]))(matrix.modelfitting[nonperfectvects.mf,], Ncat)$PFscores[[1]]
    } else
    {
      PFS.modelfitting <- eval(parse(text = x[[2]]))(matrix.modelfitting, Ncat)$PFscores[[1]]
    }
  }
  
  # Compute cutoff: 
  if (is.na(UDlvl)) # By means of bootstrap
  {
    if (any(x$PFStatistic == upp.PFS))
    {
      tail <- "upper"
      Blvl.use <- 1-Blvl
    }
    if (any(x$PFStatistic == low.PFS))
    {
      tail <- "lower"
      Blvl.use <- Blvl
    }   
    Bvec <- c()
    for (i in 1:Breps)
    {
      Bvec <- c(Bvec, quantile(sample(PFS.modelfitting, size=length(PFS.modelfitting), replace=TRUE), probs=Blvl.use))
    }
    cutoff.use <- round(median(Bvec), 4)
    cutoff.SE  <- round(sd(Bvec), 4)
    cutoff.CI  <- round(quantile(Bvec, probs=c((1-CIlvl)/2, (1+CIlvl)/2)), 4)
  } else
  { 
    if (any(x$PFStatistic == upp.PFS)) {tail <- "upper"}
    if (any(x$PFStatistic == low.PFS)) {tail <- "lower"}
    cutoff.use <- UDlvl
    cutoff.SE  <- NA
    cutoff.CI  <- NA
  }
  
  # Determine the proportion of flagged subjects:
  if (any(x$PFStatistic == upp.PFS))
  {
    prop.flagged <- round(sum(pf.scores.noNAs >= cutoff.use) / N, 4)
  }
  if (any(x$PFStatistic == low.PFS))
  {
    prop.flagged <- round(sum(pf.scores.noNAs <= cutoff.use) / N, 4)
  }
  #
  res        <- list(Cutoff=as.numeric(cutoff.use), Cutoff.SE=cutoff.SE, Prop.flagged=prop.flagged, Tail=tail, 
                     Cutoff.CI=cutoff.CI)
  class(res) <- "PerFit.cutoff"
  return(res)
}

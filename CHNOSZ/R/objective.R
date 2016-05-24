# CHNOSZ/objective.R
# objective functions for revisit(), findit()
# all objective functions in one file, added attributes  20121009 jmd 

get.objfun <- function(objective) {
  # get the function with the name given in 'objective'
  objfun <- get(objective)
  # perform a check on its usability
  if(!"optimum" %in% names(attributes(objfun)))
    stop(paste(objective, "is not a function with an attribute named 'optimum'"))
  # return the function
  return(objfun)
}

# input: a1 or loga1 is matrix with a column for each species and a row for each condition
# output: vector with length equal to number of rows of a1 or loga1
# attributes:
# optimum - the target for optimization (minimal or maximal)
## absolute - is the absolute value of the function used in the calculation? (TRUE or FALSE)

SD <- structure(
  function(a1) {
    # calculate standard deviation
    SD <- apply(a1, 1, sd)
    return(SD)
  },
  optimum="minimal"
)

CV <- structure(
  function(a1) {
    # get the coefficient of variation
    SD <- SD(a1)
    CV <- SD / rowMeans(a1)
    return(CV)
  },
  optimum="minimal"
)

shannon <- structure(
  function(a1) {
    # shannon entropy
    p <- a1/rowSums(a1)
    H <- rowSums(-p*log(p))
    return(H)
  },
  optimum="maximal"
)

DGmix <- structure(
  function(loga1) {
    # Gibbs energy/2.303RT of ideal mixing
    a1 <- 10^loga1
    DGmix <- rowSums(a1*loga1)
    return(DGmix)
  },
  optimum="minimal"
)

qqr <- structure(
  function(loga1) {
    # correlation coefficient of a q-q plot (for log-normality comparisons) 20100901
    # first, a function to calculate qqr
    qqrfun <- function(x) {
      # this is to catch errors from qqnorm
      qqn <- try(qqnorm(x, plot.it=FALSE), silent=TRUE)
      if(class(qqn)=="try-error") qqr <- NA
      else qqr <- cor(qqn[[1]], qqn[[2]])
    }
    # apply the function to the rows of loga1
    qqr <- apply(loga1, 1, qqrfun)
    return(qqr)
  },
  optimum="maximal"
)

logact <- structure(
  function(loga1, loga2) {
    # the logarithm of activity of the species identified in loga2 (species number)
    return(loga1[, loga2[1]])
  },
  optimum="maximal"
)

spearman <- structure(
  function(loga1, loga2) {
    # Spearman's rho (rank correlation coefficient) 20100912
    spearfun <- function(a, b) {
      # based on help(dSpearman) in package SuppDists
      if(length(a)!=length(b)) stop("a and b must be same length")
      if(any(is.na(a)) | any(is.na(b))) return(NA)
      ra <- rank(a)
      rb <- rank(b)
      dr <- rb - ra
      d <- sum(dr^2)
      r <- length(a)
      x <- 1-6*d/(r*(r^2-1))
      return(x)
    }
    rho <- apply(loga1, 1, spearfun, b=loga2)
    return(rho)
  },
  optimum="maximal"
)

pearson <- structure(
  function(loga1, loga2) {
    # pearson correlation coefficient
    H <- apply(loga1, 1, cor, y=loga2)
    return(H)
  },
  optimum="maximal"
)

RMSD <- structure(
  function(loga1, loga2) {
    rmsd <- function(a, b) {
      # to calculate the root mean square deviation
      d <- b - a
      sd <- d^2
      msd <- mean(sd)
      rmsd <- sqrt(msd)
      return(rmsd)
    }
    RMSD <- apply(loga1, 1, rmsd, b=loga2)
    return(RMSD)
  },
  optimum="minimal"
)

CVRMSD <- structure(
  function(loga1, loga2) {
    # coefficient of variation of the root mean squared deviation
    RMSD <- RMSD(loga1, loga2)
    CVRMSD <- RMSD/rowMeans(loga1)
    return(CVRMSD)
  },
  optimum="minimal"
)

DDGmix <- structure(
  function(loga1, loga2) {
    # difference in Gibbs energy/2.303RT of ideal mixing
    a2 <- 10^loga2
    DGmix2 <- sum(a2*loga2)
    DGmix1 <- DGmix(loga1)
    DDGmix <- DGmix2 - DGmix1
    return(DDGmix)
  },
  optimum="minimal"
)

DGinf <- structure(
  function(a1, a2) {
    dginf <- function(a1, a2) {
      # informatic Gibbs energy/2.303RT difference between assemblages
      p1 <- a1/sum(a1)
      p2 <- a2/sum(a2)
      sum(p2 * log10(p2/p1))
    }
    DGinf <- apply(a1, 1, dginf, a2=a2)
    return(DGinf)
  },
  optimum="minimal"
)

DGtr <- structure(
  function(loga1, loga2, Astar) {
    dgtr <- function(loga1, loga2, Astar) {
      # calculate the Gibbs energy/2.303RT of transformation
      # from an initial to a final assemblage at constant T, P and 
      # chemical potentials of basis species 20120917
      # loga1 - logarithms of activity of species in the initial assemblage
      # loga2 - logarithms of activity of species in the final assemblage
      # Astar - starved (of activity of the species of interest) values of chemical affinity/2.303RT
      # remove the logarithms
      a1 <- 10^loga1
      a2 <- 10^loga2
      # the molal Gibbs energy in the initial and final states
      # (derived from integrating A = Astar - RTln(a) from a1 to a2)
      G1 <- -a1*Astar + a1*loga1 - a1/log(10)
      G2 <- -a2*Astar + a2*loga2 - a2/log(10)
      # calculate the change in molal Gibbs energy for each species
      DG <- G2 - G1
      # return the sum
      return(sum(DG))
    }
    # we need to index both loga1 and Astar
    DGtr <- unlist(lapply(seq(nrow(loga1)), function(i) {
      dgtr(loga1[i, ], loga2, Astar[i, ])
    }))
    return(DGtr)
  },
  optimum="minimal"
)

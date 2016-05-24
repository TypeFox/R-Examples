#'Data simulation function
#'
#'Adapted from the supplementary material from Law \emph{et al}.
#'
#'@references Law, C. W., Chen, Y., Shi, W., & Smyth, G. K. (2014). voom: Precision
#'weights unlock linear model analysis tools for RNA-seq read counts. \emph{Genome
#'Biology}, 15(2), R29.
#'
#'
#'@examples
#'\dontrun{
#'set.seed(123)
#'data_sims <- data_sim_voomlike(maxGSsize=300)
#'}
#'@keywords internal
#'@importFrom utils data
#'@importFrom stats cor rchisq rgamma rnorm rpois model.matrix runif
#'@export
data_sim_voomlike <- function(seed=NULL, maxGSsize=400, minGSsize=30){


  ############################################################################
  # Get distribution function of abundance proportions
  # This distribution was generated from a real dataset
  utils::data("qAbundanceDist", envir = environment())
  qAbundanceDist <- get("qAbundanceDist")

  # Generate baseline proportions for desired number of genes -----
  ngenes <- 10000
  baselineprop <- qAbundanceDist( (1:ngenes)/(ngenes+1) )
  baselineprop <- baselineprop/sum(baselineprop)

  # Design ----
  n <- 18
  n1 <- n/2
  n2 <- n1
  nindiv <- 6
  ngroup <- 2
  ntime <- 3
  group <- factor(rep(1:2, each=n/ngroup))
  indiv <- factor(rep(1:nindiv, each=n/nindiv))
  time <- rep(1:ntime, nindiv)
  design <- stats::model.matrix(~ group + time)
  #design <- design[, 1, drop=FALSE]
  nlibs <- n

  # Library size ----
  # Use equal or unequal library sizes
  equal <- FALSE
  if(equal){
    expected.lib.size <- rep(11e6,nlibs)
  } else {
    expected.lib.size <- 20e6 * rep(c(1,0.1), n1)
  }

  # Set seed ----
  if(!is.null(seed)){set.seed(seed*54321)}
  u <- stats::runif(100)

  # Expected counts, group basis ----
  i <- sample(1:ngenes,200)
  i1 <- i[1:100]
  i2 <- i[101:200]
  fc <- 2
  baselineprop1 <- baselineprop2 <- baselineprop
  baselineprop1[i1] <- baselineprop1[i1]*fc
  baselineprop2[i2] <- baselineprop2[i2]*fc
  mu0.1 <- matrix(baselineprop1,ngenes,1) %*% matrix(expected.lib.size[1:n1],1,n1)
  mu0.2 <- matrix(baselineprop2,ngenes,1) %*% matrix(expected.lib.size[(n1+1):(n1+n2)],1,n2)
  mu0 <- cbind(mu0.1,mu0.2)
  status <- rep(0,ngenes)
  status[i1] <- -1
  status[i2] <- 1

  # Biological variation ----
  BCV0 <- 0.2+1/sqrt(mu0)

  # Use inverse chi-square or log-normal dispersion
  invChisq <- TRUE

  if(invChisq){
    df.BCV <- 40
    BCV <- BCV0*sqrt(df.BCV/stats::rchisq(ngenes ,df=df.BCV))
  } else {
    BCV <- BCV0*exp(stats::rnorm(ngenes,mean=0,sd=0.25)/2 )
  }
  if(NCOL(BCV)==1) BCV <- matrix(BCV, ngenes, nlibs)
  shape <- 1/BCV^2
  scale <- mu0/shape
  mu <- matrix(stats::rgamma(ngenes*nlibs, shape=shape,scale=scale), ngenes, nlibs)

  # Technical variation
  counts <- matrix(stats::rpois(ngenes*nlibs, lambda=mu), ngenes, nlibs)

  # Filter
  keep <- rowSums(counts)>=10
  nkeep <- sum(keep)
  counts2 <- counts[keep,]

  S_temp <- stats::cor(t(counts2))
  rownames(S_temp) <- as.character(1:nrow(S_temp))
  colnames(S_temp) <- as.character(1:ncol(S_temp))
  rownames(counts2) <- rownames(S_temp)
  GS <- list()
  for(i in 1:ncol(S_temp)){
    GS[[i]] <-which(abs(S_temp[,1])>0.8)
    #if(inherits(GS[[i]], "try-error")){browser()}
    S_temp <- S_temp[-GS[[i]], -GS[[i]]]
    if(is.null(dim(S_temp)) || nrow(S_temp)==0){
      break()
    }
  }
  gs_keep <- lapply(GS[sapply(GS, length)<maxGSsize & sapply(GS, length)>minGSsize], names)

  return(list("counts"=counts2, "design"=design, "gs_keep"=gs_keep, "indiv"=indiv))
}


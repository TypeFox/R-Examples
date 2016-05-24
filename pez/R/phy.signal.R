#' Calculate phylogenetic `signal'
#' 
#' @param data \code{\link{comparative.comm}} object
#' @param method what kind of signal to calculate, one of Pagel's
#' \eqn{$\lambda$}{lambda} (default), \eqn{$\delta$}{delta}, and
#' \eqn{$\kappa$}{kappa}, or Blomberg's K.
#' @details Phylogenetic `signal' is one of those concepts that is
#' said a lot in community ecology, but some do not full consider its
#' meaning. Think carefully before rushing to report a value whether:
#' (1) it makes sense to assess phylogenetic `signal' in your
#' datasets, and (2) what the phrase `phylogenetic signal' actually
#' means. This code makes use of \code{caper::\link{pgls}} to get
#' estimates of fit; alternatives that offer more flexibility exist
#' (see below).
#' @return Named numeric vector, where each element is a trait or community.
#' @author Will Pearse, Jeannine Cavender-Bares
#' @references Blomberg S.P., Garland T. & Ives A.R. Testing for phylogenetic signal in comparative data: behavioral traits are more labile. Evolution 57(4): 717--745.
#' @references R. P. Freckleton, P. H. Harvey, and M. Pagel. Phylogenetic analysis and comparative data: A test and review of evidence. American Naturalist, 160:712-726, 2002.
#' @references Mark Pagel (1999) Inferring the historical patterns of biological evolution. Nature 6756(401): 877--884.
#' @seealso fitContinuous fitDiscrete pgls phylosignal
#' @examples
#' data(laja)
#' data <- comparative.comm(invert.tree, river.sites, invert.traits)
#' phy.signal(data, "lambda")
#' @importFrom caper comparative.data pgls
#' @importFrom picante Kcalc
#' @importFrom stats setNames
#' @export
phy.signal <- function(data, method=c("lambda", "delta", "kappa", "blom.k")){
  #Assertions and argument handling
  if(!inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  method <- match.arg(method)
  if(is.null(data$data))
      stop("'data' must contain traits to calculate phylogenetic signal (of traits...!)")
  if(any(!sapply(data$data, is.numeric)))
      warning("'data' not all trait values are continuous\n\tinterpret signal of discrete data with caution; may also generate errors")

  traits <- numeric(ncol(data$data))
  for(i in seq(ncol(data$data))){
      #I know, this should be a case statement...
      if(method == "lambda"){
          model <- pgls(data$data[,i] ~ 1, data=data, lambda="ML")
          traits[i] <- model$param.CI$lambda$opt
      }
      if(method == "delta"){
          model <- pgls(data$data[,i] ~ 1, data=data, delta="ML")
          traits[i] <- model$param.CI$delta$opt
      }
      if(method == "kappa"){
        model <- pgls(data$data[,i] ~ 1, data=data, kappa="ML")
        traits[i] <- model$param.CI$kappa$opt
    }
      if(method == "blom.k"){
          #...better safe than sorry with names...
          traits[i] <- Kcalc(setNames(data$data[,i],rownames(data$data)), data$phy)
      }
  }
  names(traits) <- names(data$data)
  return(traits)
}

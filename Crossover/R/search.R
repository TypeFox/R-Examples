#' Search for a Cross-Over Design
#' 
#' Search for a Cross-Over Design
#' 
#' See the vignette of this package for further details.
#' 
#' @param s Number of sequences.
#' @param p Number of periods.
#' @param v Number of treatments.
#' @param model Model - one of the following: "Standard additive model" (2),
#' "Second-order carry-over effects" (3), "Full set of interactions" (3),
#' "Self-adjacency model" (3), "Placebo model" (2), "No carry-over into self
#' model" (2), "Treatment decay model" (2), "Proportionality model" (1), "No carry-over effects" (0).
#' The number in parentheses is the number of different efficiency factors that can
#' be specified.
#' @param eff.factor Weights for different efficiency factors. (Not used in the
#' moment.)
#' @param v.rep Integer vector specifying how often each treatment should be
#' assigned (sum must equal s*p).
#' @param balance.s Boolean specifying whether to allocate the treatments as
#' equally as possible to each sequence (can result in loss of efficiency).
#' @param balance.p Boolean specifying whether to allocate the treatments as
#' equally as possible to each period (can result in loss of efficiency).
#' @param verbose Level of verbosity, a number between 0 and 10. The default
#' \code{verbose=0} does not print any output, while \code{verbose=10} prints
#' any available notes.
#' @param model.param List of additional model specific parameters. In the
#' moment these are \code{ppp}, the proportionality parameter for the
#' proportionality model, and \code{placebos}, the number of placebo treatments
#' in the placebo model.
#' @param n \code{n=c(n1,n2)} with \var{n1} the number of hill climbing steps
#' per trial and \var{n2} the number of searches from random start matrices.
#' @param jumps To reduze the possibility of the hill-climbing algorithm to get
#' stuck in local extrema long jumps of distance \var{d} can be performed all
#' \var{k} steps. This can be specified as \code{long.jumps=c(d,k)}. If
#' \var{long.jumps} has only length 1 the default for \var{k} is 50.  If after
#' \var{k/2} hill-climbing steps the old design criterion is not enhanced (or
#' at least reached), the algorithm returns to the design from before the jump.
#' @param start.designs A single design or a list of start designs. If missing or to few start
#' designs are specified (with regard to parameter \var{n} which specifies a
#' number of 20 start designs as default) the start designs are generated
#' randomly with the sample function. Alternatively
#' \code{start.designs="catalog"} can be used to take start designs from the
#' catalog to which random designs are added till \var{n2} start designs are at
#' hand.
#' @param contrast Contrast matrix to be optimised. TODO: Example and better
#' explanation for contrast.
#' @param random.subject Should the subject effects be random (\code{random.subject=TRUE})
#' or fixed effects (\code{random.subject=FALSE}).
#' @param correlation Either a correlation matrix for the random subject effects or one 
#' of the following character strings: "equicorrelated", "autoregressive"
#' @param rho Parameter for the correlation if parameter \code{correlation} is a character string.
#' @return Returns the design as an integer matrix.
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @references John, J. A., Russell, K. G., & Whitaker, D. (2004). CrossOver:
#' an algorithm for the construction of efficient cross-over designs.
#' Statistics in medicine, 23(17), 2645-2658.
#' @keywords misc
#' @examples
#' 
#' \dontrun{
#' x <- searchCrossOverDesign(s=9, p=5, v=4, model=4)
#' 
#' jumps <- c(10000, 200) # Do a long jump (10000 changes) every 200 steps
#' n <- c(1000, 5)        # Do 5 trials with 1000 steps in each trial
#' result <- searchCrossOverDesign(s=9, p=5, v=4, model=4, jumps=jumps, n=n)
#' plot(result)
#' }
#' 
#' 
#' @export searchCrossOverDesign
searchCrossOverDesign <- function(s, p, v, model="Standard additive model", eff.factor=1,
                                  v.rep, balance.s=FALSE, balance.p=FALSE, verbose=0, model.param=list(), 
                                  n=c(5000, 20), jumps=c(5, 50), start.designs, random.subject=FALSE, contrast, correlation=NULL, rho=0) {
  #seed <<- .Random.seed #TODO Do not forget to remove this after testing! :)
  start.time <- proc.time()
  if (length(n)==1) {
    if (missing(start.designs)) { n <- c(n, 20) } else { n <- c(n, length(start.designs)) }
  }
  if (length(jumps)==1) jumps <- c(jumps, 50)
  if (jumps[2]==0) stop("The second component of 'jumps' must be a positive integer.")
  model <- getModelNr(model)
  
  H <- do.call( linkMatrix, c(list(model=model, v=v), model.param) )
  
  interchange <- TRUE
  if (missing(v.rep)) {
    v.rep <- rep((s*p) %/% v, v) + c(rep(1, (s*p) %% v), rep(0, v-((s*p) %% v)))
    interchange <- FALSE
  } else if (sum(v.rep)!=s*p) { # TODO Feature: Allow NA or sum(v.rep)<s*p
    stop("The sum of argument v.rep must equal s times p.")
  }
  if (balance.s || balance.p) interchange <- TRUE
  if (balance.s && balance.p) stop("Balancing sequences AND periods simultaneously is a heavy restriction and not supported (yet?).")  
  if (missing(contrast)) {
    Csub <- contrMat(n=rep(1, v), type="Tukey")
    class(Csub) <- "matrix" #TODO Package matrix can be improved here (IMO)!
    C <- appendZeroColumns(Csub, model, v)
  } else {
    if (is.matrix(contrast)) {
      C <- contrast
    } else {
      C <- contrMat2(type=contrast, v, model, eff.factor)
    }
  }
  if (!is.null(correlation) && !is.matrix(correlation)) correlation <- corMat(correlation, s=s, p=p, rho=rho)
  if (missing(start.designs)) { start.designs <- list() }  # In this list we save n[2] random start designs.
  if (isTRUE(start.designs %in% c("catalog","catalogue"))) { 
    st <- get(".summary_table", envir=Crossover.env)
    start.designs <- lapply(st[st$t==v & st$p==p & st$s==s,]$dataset, get, envir=Crossover.env)
    x <- getStartDesigns(s=s, p=p, v=v)
    if (length(x)!=0) start.designs <- c(start.designs, x)
  }
  if (!is.list(start.designs)) {
      start.designs <- list(start.designs)
  }
  i <- length(start.designs) + 1
  while (i <= n[2]) {    
    start.designs[[i]] <- randomDesign(s, p, v,  v.rep, balance.s=balance.s, balance.p=balance.p, model=model, C=C)
    i <- i + 1
  }
  if (length(start.designs)!=n[2]) { warning(paste("More start designs specified than number of search runs. Search runs are increased to ", length(start.designs), ".", sep="")) }

  CC <- t(C) %*% C

  if (model==8) { # Second-order carry-over effects
      r <- c(rep(s/v, v), rep((p-1)*s/v^2, v^2), rep((p-2)*s/v^3, v^3))
  } else if (model==9) {
      r <- c(rep(s*p/v, v))
  } else {
      r <- c(rep(s/v, v), rep((p-1)*s/v^2, v^2))
  }
  S2 <- sum(diag(ginv(t(H) %*% diag(r) %*% H) %*% CC))
  
  result <- .Call( "searchCOD", as.integer(s), as.integer(p), as.integer(v), 
                   start.designs, H, C, model, 
                   v.rep, balance.s, balance.p, verbose, 
                   as.integer(n), as.integer(jumps), S2, TRUE, random.subject, correlation, interchange, PACKAGE = "Crossover")
  
  design <- result$design
  
  if (!estimable(design, v=v, model=model, C=C)) {
      stop("Something went wrong. Specified contrasts are not estimable with this design.")
  }  
  time <- proc.time()-start.time
  class(time) <- NULL
  #varTrtPair <- paste(capture.output(print(general.carryover(t(design), model=model))), collapse = "\n")
  #return(list(design=design, varTrtPair=varTrtPair, eff=eff, search=, time=))
  return(new("CrossoverSearchResult", design=new("CrossoverDesign", design, v=v, model=model, description="Found by search algorithm"), startDesigns=start.designs, eff=result$eff,                   
             search=list(n=n, jumps=jumps), model=model, time=time, misc=list(designs=result$designs)))
}

getS1 <- function(design, v, model, C, randomS=FALSE, verbose=0) {
    if(missing(C)) {
        Csub <- contrMat(n=rep(1, v), type="Tukey")
        class(Csub) <- "matrix" #TODO Package matrix can be improved here (IMO)!
        C <- appendZeroColumns(Csub, model, v)
    }
    linkM <- linkMatrix(model, v)
    return(.Call( "getS12R", design, v, model, linkM, C, randomS))
}

appendZeroColumns <- function(Csub, model, v) {
  model <- getModelNr(model)
  if (model %in% c(2,8)) {
    C <- as.matrix(cbind(Csub, matrix(0, dim(Csub)[1], 2*v)))
  } else if (model %in% c(3,9)) {
    C <- Csub
  } else if (model == 7) {
      C <- as.matrix(cbind(Csub,matrix(0,dim(Csub)[1], v+v*v)))
  } else if (model %in% c(1,4,5,6) ) {
    C <- as.matrix(cbind(Csub, matrix(0,dim(Csub)[1], v)))
  }  
  return(C)
}

getTrtPair <- function(design, model=1) {
  gco <- general.carryover(t(design), model=model)
  return(triu(gco$Var.trt.pair))
}

getValues <- function(design, model=1, C, v, ppp=0.5, placebos=1) {
  model <- getModelNr(model)
  if (missing(C)) {
    Csub <- contrMat(n=rep(1, v), type="Tukey")
    class(Csub) <- "matrix" #TODO Package matrix can be improved here (IMO)!
    C <- appendZeroColumns(Csub, model, v)
    CC <- t(C) %*% C
  }
  rcDesign <- rcd(design, v, model=model)
  Ar <- infMatrix(rcDesign, v, model=model)
  H <- linkMatrix(model, v, ppp=ppp, placebos=placebos)
  return(diag(C %*% ginv(t(H) %*% Ar %*% H) %*% t(C)))  
}

getDesignMatrix <- function(design, v) {
    H <- linkMatrix(model="Standard additive model", v)
    rcDesign <- rcd(design, v=v, model=1)
    Xr <- rcdMatrix(rcDesign, v, model=1)
    return(Xr %*% H)
}
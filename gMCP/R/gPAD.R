#' EXPERIMENTAL: Evaluate conditional errors at interim for a pre-planned
#' graphical procedure
#' 
#' Computes partial conditional errors (PCE) for a pre-planned graphical
#' procedure given information fractions and first stage z-scores. -
#' Implementation of adaptive procedures is still in an early stage and may
#' change in the near future
#' 
#' For details see the given references.
#' 
#' @param graph A graph of class \code{\link{graphMCP}}.
#' @param z1 A numeric vector giving first stage z-scores.
#' @param v A numeric vector giving the proportions of pre-planned measurements
#' collected up to the interim analysis. Will be recycled of length different
#' than the number of elementary hypotheses.
#' @param alpha A numeric specifying the maximal allowed type one error rate.
#' @return An object of class \code{gPADInterim}, more specifically a list with
#' elements
#' \describe{
#' \item{\code{Aj}}{a matrix of PCEs for all elementary hypotheses in each
#' intersection hypothesis}
#' \item{\code{BJ}}{a numeric vector giving sum of PCEs per intersection
#' hypothesis}
#' \item{\code{preplanned}}{Pre planned test represented by an object of class}
#' \code{\link{graphMCP}}
#' }
#' @author Florian Klinglmueller \email{float@@lefant.net}
#' @seealso \code{\link{graphMCP}}, \code{\link{secondStageTest}}
#' @references Frank Bretz, Willi Maurer, Werner Brannath, Martin Posch: A
#' graphical approach to sequentially rejective multiple test procedures.
#' Statistics in Medicine 2009 vol. 28 issue 4 page 586-604.
#' \url{http://www.meduniwien.ac.at/fwf_adaptive/papers/bretz_2009_22.pdf}
#' 
#' Frank Bretz, Martin Posch, Ekkehard Glimm, Florian Klinglmueller, Willi
#' Maurer, Kornelius Rohmeyer (2011): Graphical approaches for multiple
#' comparison procedures using weighted Bonferroni, Simes or parametric tests.
#' Biometrical Journal 53 (6), pages 894-913, Wiley.
#' \url{http://onlinelibrary.wiley.com/doi/10.1002/bimj.201000239/full}
#' 
#' Posch M, Futschik A (2008): A Uniform Improvement of Bonferroni-Type Tests
#' by Sequential Tests JASA 103/481, 299-308
#' 
#' Posch M, Maurer W, Bretz F (2010): Type I error rate control in adaptive
#' designs for confirmatory clinical trials with treatment selection at interim
#' Pharm Stat 10/2, 96-104
#' @keywords htest graphs
#' @examples
#' 
#' 
#' ## Simple successive graph (Maurer et al. 2011)
#' ## two treatments two hierarchically ordered endpoints
#' a <- .025
#' G <- simpleSuccessiveI()
#' ## some z-scores:
#' 
#' p1=c(.1,.12,.21,.16)
#' z1 <- qnorm(1-p1)
#' p2=c(.04,1,.14,1)
#' z2 <- qnorm(1-p2)
#' v <- c(1/2,1/3,1/2,1/3)
#' 
#' intA <- doInterim(G,z1,v)
#' 
#' ## select only the first treatment 
#' fTest <- secondStageTest(intA,c(1,0,1,0))
#' 
#' 
#' 
#' @export doInterim
#' 
doInterim <- function(graph,z1,v,alpha=.025){
  g <- graph2matrix(graph)
  w <- getWeights(graph)
  ws <- generateWeights(g,w)
  n <- length(w)
  As <- t(apply(ws[,(1+n):(2*n)],1,partialCE,z1=z1,v=v,alpha=alpha))
  Bs <- rowSums(As)
  res <- new('gPADInterim',Aj=As,BJ=Bs,z1=z1,v=v,preplanned=graph,alpha=alpha)
  return(res)
}

#' EXPERIMENTAL: Construct a valid level alpha test for the second stage of an
#' adaptive design that is based on a pre-planned graphical MCP
#' 
#' Based on a pre-planned graphical multiple comparison procedure, construct a
#' valid multiple level alpha test that conserves the family wise error in the
#' strong sense regardless of any trial adaptations during an unblinded interim
#' analysis. - Implementation of adaptive procedures is still in an early stage
#' and may change in the near future
#' 
#' For details see the given references.
#' 
#' @param interim An object of class \code{\link{gPADInterim}}.
#' @param select A logical vector giving specifying which hypotheses are
#' carried forward to the second stage
#' @param matchCE Logical specifying whether second stage weights should be
#' computed proportional to corresponding PCEs
#' @param zWeights Either "reject","accept", or "strict" giving the rule what
#' should be done in cases where none of the selected hypotheses has positive
#' second stage weight.
#' @param G2 An object of class \code{\link{graphMCP}} laying down the rule to
#' compute second stage weights. Defaults to pre-planned graph.
#' @return A function of signature \code{function(z2)} with arguments
#' \code{z2} a numeric vector with second stage z-scores (Z-scores of
#' dropped hypotheses should be set no \code{NA})
#' that returns objects of class \code{\link{gMCPResult}}.
#'  
#' @author Florian Klinglmueller \email{float@@lefant.net}
#' @seealso \code{\link{graphMCP}}, \code{\link{doInterim}}
#' @references Frank Bretz, Willi Maurer, Werner Brannath, Martin Posch: A
#' graphical approach to sequentially rejective multiple test procedures.
#' Statistics in Medicine 2009 vol. 28 issue 4 page 586-604.
#' \url{http://www.meduniwien.ac.at/fwf_adaptive/papers/bretz_2009_22.pdf}
#' 
#' Bretz F., Posch M., Glimm E., Klinglmueller F., Maurer W., Rohmeyer K.
#' (2011): Graphical approaches for multiple endpoint problems using weighted
#' Bonferroni, Simes or parametric tests - to appear.
#' 
#' Posch M, Futschik A (2008): A Uniform Improvement of Bonferroni-Type Tests
#' by Sequential Tests JASA 103/481, 299-308
#' 
#' Posch M, Maurer W, Bretz F (2010): Type I error rate control in adaptive
#' designs for confirmatory clinical trials with treatment selection at interim
#' Pharm Stat 10/2, 96-104
#' @keywords htest graphs
#' @examples
#' 
#' 
#' ## Simple successive graph (Maurer et al. 2011)
#' ## two treatments two hierarchically ordered endpoints
#' a <- .025
#' G <- simpleSuccessiveI()
#' ## some z-scores:
#' 
#' p1=c(.1,.12,.21,.16)
#' z1 <- qnorm(1-p1)
#' p2=c(.04,1,.14,1)
#' z2 <- qnorm(1-p2)
#' v <- c(1/2,1/3,1/2,1/3)
#' 
#' intA <- doInterim(G,z1,v)
#' 
#' ## select only the first treatment 
#' fTest <- secondStageTest(intA,c(1,0,1,0))
#' 
#' 
#' 
#' @export secondStageTest
#' 
secondStageTest <- function(interim,select,matchCE=TRUE,zWeights="reject",G2=interim@preplanned){
  n <- nhyp(interim@preplanned)
  w2s <- t(sapply(1:(2^n-1),function(J) adaptWeights(to.binom(J,n),select,G2,zWeights)))
  Cs <- w2s*interim@BJ
  if(matchCE){
    Cs <- t(apply(cbind(interim@BJ,w2s),1,function(Bw){
      matchCE(Bw[-1],Bw[1],interim@z1,interim@v,interim@alpha)
    }))
  }
  return(function(z) {
    decideTest(z,Cs)
  })
}
                      
nhyp <- function(graph){
  return(nrow(graph2matrix(graph)))
}

            
validPartialCEs <- function(object) {
  ## if(all(rowSums(object@Aj)==BJ)){
  ##   return(TRUE)
  ## } else {
  ##   stop("Invalid interim results PCEs do not match corresponding sums")
  ## }
  return(TRUE)
}
                        

partialCE <- function(w,z1,v,alpha){
  ## conditional error for an elementary hypothesis with weight at level alpha and first stage z-score and proportion v for the first stage
  ## also works for vectors
  ## returns the A(i,J) or if called with a vector the vector A(J,J)
  1-pnorm((qnorm(1-w*alpha)-(sqrt(v)*z1))/sqrt(1-v))
}

matchCE <- function(w2,B,z1,v,alpha,enhanced=T){
  ## find a suitable alpha level that matches the sum of PCEs for the selected hypotheses and adapted weights to that of the pre-planned procedure
  if(all(w2==0)){
    return(w2)
  }
  ## enhanced for B>1 we can reject the intersection at interim
  if(B>1){
    return(rep(1,length(w2)))
  }
  d <- function(alpha,w2,z1,v,B){
    sum(partialCE(w2,z1,v,alpha))-B
  }
  ## catch zero's
  r <- uniroot(d,c(0,1),w2=w2,z1=z1,v=v,B=B)$root
  partialCE(w2,z1,v,r)
}

adaptWeights <- function(J,select,G2,dw='reject'){
  ## adapt the weights this is basically a wrapper to mtp.weights that handles dropped hypotheses
  w <- getWeights(G2)
  g <- graph2matrix(G2)
  ## in case we only include selected hypotheses
  if(all((J-select)>0)){
    return(mtp.weights(J,g,w))
  }
  ## canonical rule number 1
  if(dw=='reject'){
    return(mtp.weights(J * select,g,w))
  }
  ## canonical rule number 2 with fallback to 1 in case all weights are zero
  if(dw=='accept'){
    if(all(w <- mtp.weights(J,g,w)*select)==0){
      return(mtp.weights(J * select,g,w))
    } else {
      return(w)
    }
  }
  ## strict rule may produce all zero weights
  if(dw=='strict'){
    return(mtp.weights(J,g,w)*select)
  }
  ## 
  stop('Invalid rule to determine second stage weights')
}

decideTest <- function(z,bounds){
  p <- 1-pnorm(z)
  dm <- t(sapply(1:nrow(bounds),function(n) {
    ## check whether z is larger than boundary
    m <- ncol(bounds)
    b <- bounds[n,]
    J <- to.binom(n,m)
    d <- rep(NA,length(b))
    d[which(J==1)] <- (b[which(J==1)]>=p[which(J==1)])
    return(d)
  }))
  d <- apply(dm,2,function(h) {
    ## closed testing
    all(apply(dm[!is.na(h),],1,any,na.rm=T))
  })
  d
}

to.binom <- function(int,n=floor(log2(int))+1){
  ## 6 times faster than the old function (Thankyou!)
  if(n+2<=floor(log2(int))){
    stop('Vector length to small to hold binary number')
  }
  ((int)%/% 2^((n:1)-1))%%2
}



parse.intersection <- function(binom){
  paste("H(",paste(which(binom==1),collapse=','),")",sep="")
}

to.intersection <- function(int){
  maxn <- floor(log2(max(int)))+1
  if(length(int)>1){
    unlist(lapply(lapply(int,to.binom,n=maxn),parse.intersection))
  } else {
    parse.intersection(to.binom(int,n=maxn))
  }
}
         
## test.to.binom <- function(v){
##   sum(2^(which(v)-1))
## }

## all((1:1000-sapply(lapply(1:1000,to.binom),test.to.binom))==0)

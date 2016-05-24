pps.sampling <-
function(z, n, id=1:N, method="sampford", return.PI=FALSE){
    ## The function provides a sample which is proportional to the size
    ## of a quantity z, say.
    ## z = vector of size N. quantity which determines the proportionality
    ## n = sample size
    ## id = identification number, default is 1:N
    ## method = sampling method. Options are "sampford", "tille", "midzuno" or "madow"
    ## return.PI = if TRUE the pairwise inclusion probabilities for all individuals
    ## in the population is returned
    if(missing(z)) stop("Wrong input: vector ", sQuote("z")," of quantities needs to be given.")
    if(any(z==0)) stop("Wrong input: quantities in vector ", sQuote("z")," not allowed to be zero.")
    N <- length(z)
    if(missing(n) || n<=0) stop("Wrong input: sample size ", sQuote("n")," has to be positive integer.")
    else if(!n < N) stop("Wrong input: Sample size ", sQuote("n")," has to be smaller than length of ", sQuote("z"),", which is population size ", sQuote("N"),".")
    if(!(method=="sampford" || method=="tille" || method=="midzuno" || method=="madow" )){
      stop("Wrong input: Only types ", sQuote("sampford"),", ", sQuote("tille"),", ", sQuote("midzuno")," or ", sQuote("madow")," allowed")
    }
    if(method=="sampford" & (n/N) > 0.3) stop("Wrong input: For using method = ", sQuote("sampford")," the relation ", sQuote("n/N")," needs to be smaller than 0.3. Use other methods.")
    if(method=="sampford" && N>200) stop("Long run: method = ", sQuote("sampford")," has for ", sQuote("N > 200")," a long run. Better use method ", sQuote("midzuno")," or ", sQuote("madow"),".")
    if(method=="tille" && N>500) stop("Long run: method = ", sQuote("tille")," has for ", sQuote("N > 500")," a long run. Better use method ", sQuote("midzuno")," or ", sQuote("madow"),".")
    if(method=="sampford")
      {
        PI.full <- sampfordpi(z,n)
        index <-  sampford(z,n)
        pps.sample <-  id[index]
        PI <- PI.full[index,index]
      }
    if (method=="tille")
      {
        pik <-  inclusionprobabilities(z,n)
        PI.full <-  UPtillepi2(pik)
        index <-  (1:N)[UPtille(pik)==1]
        pps.sample <- id[index]
        PI <- PI.full[index,index]
     }
    if (method=="midzuno")
      {
        pik <-  inclusionprobabilities(z,n)
        PI.full <-  UPmidzunopi2(pik)
        index <-  (1:N)[UPmidzuno(pik)==1]
        pps.sample <- id[index]
        PI <- PI.full[index,index]
      }
   if (method=="madow")
      {
        pik <-  inclusionprobabilities(z,n)
        PI.full <-  UPsystematicpi2(pik)
        index <-  (1:N)[UPsystematic(pik)==1]
        pps.sample <- id[index]
        PI <- PI.full[index,index]
        warning("Systematic Sample with zeros in ", sQuote("PI"),": For calculating estimates use approximate methods.")
     }                               
    pik <- diag(PI.full)
    if ( return.PI == FALSE)  PI.full <-  NULL
## return argument  
    ret <- list()
    ret$call <- list(z=z,n=n,id=id,method=method)
    ret$sample <- pps.sample 
    ret$pik <- pik
    ret$PI <- PI 
    ret$PI.full <- PI.full
    structure(ret,class="pps.sampling")
}

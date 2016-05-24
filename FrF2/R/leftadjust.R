leftadjust <- function(k, gen, early=NULL, show=10){
  ## function that leftadjusts as long as the current factors remain generators
  ## it is however possible to find a better changing design with modified generator columns
  ##  ?? how can these be chosen ??
  ## for practical purposes, this will presumably be sufficient most of the time
  
  hilf <- c(k,abs(gen),show)
  if (!is.null(early)) hilf <- c(hilf,early)
  if (!is.numeric(hilf))
      stop ("All inputs to leftadjust must be numeric.")
  if (!all(hilf == floor(hilf) & hilf > 0))
      stop ("All inputs to leftadjust must contain integer numbers.")
  if (!k >= 3) stop ("leftadjust requires k>=3.")
  g <- length(gen)
  minus <- sign(gen)
  gen <- abs(gen)
  if (any(gen %in% 2^(0:(k-1)))) 
        stop ("gen must not contain column numbers of base factors.")
  if (any(!gen %in% 3:(2^k-1))) 
        stop ("Column numbers in gen must be smaller than ", 2^k,".")
  
  if (!is.null(early))
   if(early > k + g) stop ("early must not be larger than the total number of factors.")
  
  perm <- permutations(k)
  hilf <- digitsBase(gen,ndigits=k)
  ergeb <- matrix(0,nrow=nrow(perm),ncol=length(gen))
  if (!is.null(early)) {
         maxpos <- rep(NA,nrow(perm)) 
         k.early <- rep(NA,nrow(perm)) 
      }
      else {
         maxpos <- NULL
         k.early <- NULL
      }
  for (i in 1:factorial(k)){
    ergeb[i,] <- sort(as.integer(reord(hilf,perm[i,])))
    if (!is.null(early)) {
        for (j in ceiling(log2(early+1)):k){
           if (early==j) zwischen <- 2^(j-1) else
           zwischen <- max(ergeb[i,early-j],2^(j-1))
           if (is.na(maxpos[i]) | zwischen < maxpos[i]){
                  maxpos[i] <- zwischen
                  k.early[i] <- j
               }
           }
        }
  }
  if (is.null(early)) reorder <- ord(ergeb)
  else reorder <- ord(cbind(k.early,maxpos,ergeb))
  
  list(orig=gen*minus,basics=c(nruns=2^k, nfactors=k+g, early=early),
         perms=perm[reorder[1:show],,drop=FALSE],
         maxpos=maxpos[reorder[1:show]],k.early=k.early[reorder[1:show]],
         gen=ergeb[reorder[1:show],,drop=FALSE])
}
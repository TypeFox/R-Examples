splitpick <- function(k, gen, k.WP, nfac.WP, show=10){
  ## function that picks design with right number of factors within first 
  ## so many columns as long as the current factors remain generators
  ## it may be possible to find an appropriate design within 
  ##    the same base design with modified generator columns
  ##  ?? how can these be chosen ??
  ##  ?? is this worth the extra effort ??

  ## k.WP is the dimension for WP factors
  ## 2^k.WP - 1 first columns are used exclusively for whole plots
  ## and should not contain SP factors

  ## k.WP <= nfac.WP is needed
  ##                  if not true in the subject matter situation,
  ##                  extra whole plot generation factors are required
  
  ## potential future improvement:
  ## enhance output by alias structure, if requested
  ## alias for whole plot separate as well
  
  hilf <- c(k,abs(gen),k.WP,nfac.WP,show)
  if (!is.numeric(hilf))
      stop ("All inputs to splitpick must be numeric.")
  if (!all(hilf == floor(hilf) & hilf > 0))
      stop ("All inputs to splitpick must contain integer numbers.")
  minus <- which(gen<0)
  gen <- abs(gen)
  if (!k >= 3) stop ("splitpick requires k>=3.")
  if (!k.WP < k) stop ("splitpick requires k.WP < k.")
  if (!nfac.WP < 2^k.WP) stop ("nfac.WP >= 2^k.WP is not permitted.")
  if (!nfac.WP >= k.WP) stop ("nfac.WP < k.WP is not permitted. 
       You must increase nfac.WP by ", k.WP - nfac.WP, 
       " in order to support generation of ", 2^k.WP, " whole plots.")
  if (any(gen %in% 2^(0:(k-1)))) 
        stop ("gen must not contain column numbers of base factors.")
  if (any(!gen %in% 3:(2^k-1))) 
        stop ("Column numbers in gen must be smaller than ", 2^k,".")
  g <- length(gen)
  if (!nfac.WP < k+g) stop ("nfac.WP must be smaller than the total number of factors!")
  
  perm <- permutations(k)
  hilf <- digitsBase(gen,ndigits=k) ## always k rows
  ergeb <- matrix(0,nrow=nrow(perm),ncol=length(gen))
  nfacs.WP <- rep(NA, nrow=nrow(perm))
  for (i in 1:factorial(k)){
    ergeb[i,] <- sort(as.integer(reord(hilf,perm[i,])))
    nfacs.WP[i] <- sum(ergeb[i,]<2^k.WP)
  }
  pick <- which(nfacs.WP==nfac.WP-k.WP)
  if (length(pick)<show) show <- length(pick) 
  if (show==0) stop("no adequate split-plot design found\n")
  else {
    gens <- ergeb[pick[1:show],,drop=FALSE]
    if (nfac.WP==k.WP) res.WP <- rep(Inf,show)
    else {
         if (nfac.WP > 2^(k.WP-1)) res.WP <- rep(3,show)
         else {
         hilf <- gens[,1:(nfac.WP-k.WP),drop=FALSE]
         res.WP <- apply(hilf,1,function(obj) min(sapply(words.all(k.WP,obj)[[2]],length)))
         }
    } 
    reorder <- sort(res.WP,index.return=TRUE,decreasing=TRUE)$ix
    ## if > 7, res.WP=Inf

    gen[minus] <- -gen[minus]
    ## for documentation purposes only
    ## does not affect the resulting design whether or not there is a minus there
    
    list(orig=gen, basics=c(nruns=2^k, nWPs=2^k.WP, nfac.WP=nfac.WP, nfac.SP=k+g-nfac.WP),
         perms=perm[pick[1:show][reorder],,drop=FALSE],
         res.WP=res.WP[reorder], gen=gens[reorder,,drop=FALSE])
  }
}
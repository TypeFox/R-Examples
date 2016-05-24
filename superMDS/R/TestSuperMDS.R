TestSuperMDS <-
function(trout, xte=NULL, dtetr=NULL){
  if(is.matrix(trout$z)) S <- ncol(trout$z)
  if(!is.matrix(trout$z)) S <- 1
  if(!is.null(xte) && !is.null(dtetr)) stop("Please give xte or dtetr, but not both.")
  if(is.null(xte) && is.null(dtetr)) stop("Must give either xte or dtetr.")
  if(is.null(dtetr)){
    if(is.null(trout$x)) stop("If dtetr is null, x and xte must both be given.")
    dtetr <- distmat(rbind(trout$x, xte))[(nrow(trout$x)+1):(nrow(trout$x)+nrow(xte)),1:nrow(trout$x)]
  }
  crit1s <- crit2s <- rep(NA, nrow(dtetr))
  znew1s <- znew2s <- matrix(NA, nrow(dtetr), ncol=ncol(trout$z))
  for(i in 1:nrow(dtetr)){ 
    teout <- TestSuperMDSSingleObs(dtetr[i,], trout$alpha, trout$z, trout$y, S)
    crit1s[i] <- teout$crit1
    crit2s[i] <- teout$crit2
    znew1s[i,] <- teout$znew1
    znew2s[i,] <- teout$znew2
  }
  if(is.null(trout$cutpoint)) trout$cutpoint <- 0
  ytehat <- rep(2, nrow(dtetr))
  ytehat[crit1s-crit2s < trout$cutpoint] <- 1
  zte <- znew2s
  zte[crit1s-crit2s < trout$cutpoint,] <- znew1s[crit1s-crit2s < trout$cutpoint,]
  return(list(crit1s=crit1s, crit2s=crit2s, znew1s=znew1s, znew2s=znew2s, ytehat=ytehat, zte=zte))
}


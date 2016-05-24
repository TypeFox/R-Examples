diagnostics <-
function(x, s2=1, ar.LjungB=c(1,0.025),
  arch.LjungB=c(1,0.025), normality.JarqueB=NULL,
  verbose=FALSE)
{
  diagnostics.chk <- TRUE
  if(s2 == 1){ zhat <- x }else{ zhat <- x/sqrt(s2) }

  ##verbose:
  ##to do: diagnostic tables etc.

  ##serial correlation:
  if(!is.null(ar.LjungB)){
    ar.LjungBox <- Box.test(zhat, lag=ar.LjungB[1], type="L")
    if(ar.LjungBox$p.value <= ar.LjungB[2]){
      diagnostics.chk <- FALSE
    }
  } #end if(!is.null(..))

  ##arch:
  if(diagnostics.chk && !is.null(arch.LjungB)){
    zhat2 <- zhat^2
    arch.LjungBox <- Box.test(zhat2, lag=arch.LjungB[1], type="L")
    if(arch.LjungBox$p.value <= arch.LjungB[2]){
      diagnostics.chk <- FALSE
    }
  } #end arch

  ##normality:
  if(diagnostics.chk && !is.null(normality.JarqueB)){
    n <- length(zhat)
    avgzhat <- mean(zhat)
    zhat.avgzhat <- zhat-avgzhat
    zhat.avgzhat2 <- zhat.avgzhat^2
    K <- n*sum(zhat.avgzhat^4)/(sum(zhat.avgzhat2)^2)
    S <- (sum(zhat.avgzhat^3)/n)/(sum(zhat.avgzhat2)/n)^(3/2)
    JB <- (n/6)*(S^2 + 0.25*((K-3)^2))
    pval <- pchisq(JB, df = 2, lower.tail=FALSE)
    if(pval <= normality.JarqueB){
      diagnostics.chk <- FALSE
    }
  } #end normality

  #return result:
  if(!verbose){ return(diagnostics.chk) }
}

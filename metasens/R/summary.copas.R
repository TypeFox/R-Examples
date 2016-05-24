summary.copas <- function(object, level=0.95, sign.rsb=object$sign.rsb, ...){
  
  meta:::chkclass(object, "copas")
  
  seTE <- object$seTE
  TE.random <- object$TE.random
  seTE.random <- object$seTE.random
  gamma0.slope <- object$gamma0.slope
  gamma1.slope <- object$gamma1.slope
  TE.slope <- object$TE.slope
  seTE.slope <- object$seTE.slope
  publprob <- object$publprob
  pval.rsb <- object$pval.rsb
  N.unpubl <- object$N.unpubl
  ##
  ci.random <- ci(TE.random, seTE.random, level)
  ##
  if (is.null(sign.rsb))
    sign.rsb <- 0.1
  else
    meta:::chklevel(sign.rsb)
  
  
  ci.lab <- paste(round(100*level, 1), "%-CI", sep="")
  
  
  ord <- rev(order(publprob)) 
  pom <- publprob[ord]
  ##
  TE.slope <- TE.slope[ord]
  seTE.slope <- seTE.slope[ord]
  pval.rsb <- pval.rsb[ord]
  N.unpubl <- N.unpubl[ord]
  ##
  ci.slope <- ci(TE.slope, seTE.slope, level)
  
  
  ##
  ## Copas estimate adjusted for selection bias (added by sc, 24.09.2007):
  ##
  tres <- data.frame(seq=seq(along=pval.rsb),
                     cumsum=cumsum(pval.rsb <= sign.rsb),
                     diff=seq(along=pval.rsb) - cumsum(pval.rsb <= sign.rsb))
  pval.rsb.sign.all <- all(tres$diff==0)
  pval.rsb.sign <- ifelse(sum(tres$diff==0)>0, TRUE, FALSE)
  ##
  if (pval.rsb.sign.all){
    TE.adj <- NA
    seTE.adj <- NA
    pval.rsb.adj <- NA
    N.unpubl.adj <- NA
  }
  else{
    if(pval.rsb.sign){
      sel.adj <- tres$seq[tres$diff>0][1]
      TE.adj <- TE.slope[sel.adj]
      seTE.adj <- seTE.slope[sel.adj]
      pval.rsb.adj <- pval.rsb[sel.adj]
      N.unpubl.adj <- N.unpubl[sel.adj]
    }
    else{
      TE.adj <- TE.slope[1]
      seTE.adj <- seTE.slope[1]
      pval.rsb.adj <- pval.rsb[1]
      N.unpubl.adj <- N.unpubl[1]
    }
  }
  ##
  adjust <- ci(TE.adj, seTE.adj, level)
  
  res <- list(slope=ci.slope,
              publprob=pom,
              pval.rsb=pval.rsb,
              N.unpubl=N.unpubl,
              adjust=adjust,
              sign.rsb=sign.rsb,
              pval.rsb.adj=pval.rsb.adj,
              N.unpubl.adj=N.unpubl.adj,
              random=ci.random,
              sm=object$sm,
              ci.lab=ci.lab
              )
  
  class(res) <- c("summary.copas")
  
  res$complab <- object$complab
  res$outclab <- object$outclab
  res$title   <- object$title
  
  res$backtransf <- object$backtransf
  
  res$version <- utils::packageDescription("metasens")$Version
  
  res
}

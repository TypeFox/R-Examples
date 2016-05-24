
meboot.pdata.frame <- function(x, reps=999, trim=0.10, reachbnd=TRUE,
  expand.sd=TRUE, force.clt=TRUE, scl.adjustment = FALSE, sym = FALSE, elaps=FALSE,
  colsubj, coldata, coltimes, ...)
{
  ref1 <- x[,colsubj]
  ref2 <- unique(ref1)
  
  xens <- NULL
  for(i in seq(along=ref2)){
    ir <- which(ref1==ref2[i])
    xs <- x[ir,coldata]  
    bxs <- meboot(xs, reps=reps, trim=0.10, reachbnd=reachbnd,
      expand.sd=expand.sd, force.clt=force.clt, scl.adjustment=scl.adjustment, sym=sym, 
      elaps=elaps, ...)
    xens <- rbind(xens, bxs$ensemble) 
  }

  xens <- data.frame(xens)
  dimnames(xens)[[2]] <- paste("Pseries", seq.int(reps))
  if(!missing(coltimes))
    dimnames(xens)[[1]] <- paste("Series", seq.int(reps), x[,coltimes])

  xens
}

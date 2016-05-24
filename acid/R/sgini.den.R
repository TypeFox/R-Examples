sgini.den <-
function(incs,dens,nu=2, pm0 =NA, lower = NULL, upper=NULL){
  if (length(incs) != length(dens)) 
    print("incs and dens must have the same length!")
  incs <- as.numeric(sort(incs)) 
  if (is.null(lower)) lower <- min(incs)
  if (is.null(upper)) upper <- max(incs)
  if (min(incs) < lower) print("WARNING: min(incs)<lower int limit!")
  if (!is.na(pm0) & incs[1] == 0) incs <- incs[-1]
  if (!is.na(pm0) & dens[1] == pm0) dens <- dens[-1]
  
  n.inc <- length(incs)
  incs.diff <- diff(incs) 
  av.incs <- (incs[-1] + incs[-n.inc])/2
  av.dens <- (dens[-1] + dens[-n.inc])/2
  fi.inner <- av.dens * incs.diff
  fi.outer <- c(dens[1]/2 * (incs[1] - lower), dens[n.inc]/2 * 
                  (upper - incs[n.inc]), pm0)
  if (!is.na(pm0)){
    #   print("pm0")
    probs <- c(fi.outer[3],fi.outer[1],incs.diff*av.dens,fi.outer[2])
    if(round(sum(probs),digits=1)!=1) print(paste("sum(probs)=",round(sum(probs),digits=3)," - probs don't sum up to 1!"))
    probs <- probs/sum(probs) 
    Fx.den    <- frac.ranks(c(0,incs,upper),probs)
    m.inner.int <- t(av.incs) %*% (fi.inner)
    m.outer.int <- incs[1] * fi.outer[1] + incs[n.inc] * fi.outer[2]
    mux         <- m.inner.int + m.outer.int
    cov.den <- weighted.mean(c(0,incs,upper)/mux*(1-Fx.den)^(nu-1),probs)-weighted.mean(c(0,incs,upper)/mux,probs)*weighted.mean((1-Fx.den)^(nu-1),probs)
    Gini    <- -nu *cov.den
  }else{
    probs <- c(fi.outer[1],incs.diff*av.dens,fi.outer[2])
    if(round(sum(probs),digits=1)!=1) print(paste("sum(probs)=",round(sum(probs),digits=3)," - probs don't sum up to 1!"))
    probs <- probs/sum(probs)
    Fx.den    <- frac.ranks(c(incs,upper),probs)
    m.inner.int <- t(av.incs) %*% (fi.inner)
    m.outer.int <- incs[1] * fi.outer[1] + incs[n.inc] * fi.outer[2]
    mux         <- m.inner.int + m.outer.int
    cov.den <- weighted.mean(c(incs,upper)/mux*(1-Fx.den)^(nu-1),probs)-weighted.mean(c(incs,upper)/mux,probs)*weighted.mean((1-Fx.den)^(nu-1),probs)
    Gini    <- -nu *cov.den
  } 
  list(Gini = Gini, pm0 = pm0, lower = lower, upper = upper)
}

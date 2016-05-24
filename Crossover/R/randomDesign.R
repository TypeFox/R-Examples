randomDesign <- function(s, p, v,  v.rep, balance.s=FALSE, balance.p=FALSE, model, C) {
  if (missing(v.rep)) {
    v.rep <- rep((s*p) %/% v, v) + c(rep(1, (s*p) %% v), rep(0, v-((s*p) %% v)))
  }
  design <- randomDesignWithoutCheck(s, p, v,  v.rep, balance.s, balance.p, model)
  i <- 0
  # We should disable the really rare warnings estimable_R could throw.
  while (!estimable_R(design, v, model, C)) {   
    i <- i + 1
    if (i>1000) stop("Could not find design that allows estimation of contrasts after 1000 tries.")
    design <- randomDesignWithoutCheck(s, p, v,  v.rep, balance.s, balance.p, model)
  } 
  return(design)
}

randomDesignWithoutCheck <- function(s, p, v,  v.rep, balance.s=FALSE, balance.p=FALSE, model) {
  if (balance.s) {
    design <- matrix(unlist(tapply(rep(1:v, v.rep), as.factor(rep(1:s,p)), sample)), p, s)
  } else if (balance.p) {
    design <- matrix(unlist(tapply(rep(1:v, v.rep), as.factor(rep(1:p,s)), sample)), p, s, byrow=TRUE)
  } else {
    design <- matrix(sample(rep(1:v, v.rep)), p, s)
  }
  return(design)
}
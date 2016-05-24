ei.reg <- function(formula, data, ...) {  
  D <- model.frame(formula, data = data)
  G <- D[[2]]
  T <- D[[1]]
  idx.r <- apply(G,1,sum)
  idx.c <- apply(T,1,sum)

    countG <- countT <- propG <- propT <- FALSE

  if(all(as.integer(T)==T) && all(T)>=0){
    countT <- TRUE
  }
  else{
    propT <- TRUE
  }
  
  if(all(as.integer(G)==G) & all(G)>=0){
    countG <- TRUE
  }
  else{propG <- TRUE}
  
  if(countT & countG){
    if(all(idx.r == idx.c)){
      G <- G/idx.r
      T <- T/idx.c
      countT <- countG <- FALSE
      propT <- propG <- TRUE
    }
    else{
      stop("row and column count totals unequal in some precincts -
please respecify data")
    }
  }
  
  if(countT & propG){
    if(!all(0 <= G && G <= 1)){
      stop("row proportions are not within [0,1] - please respecify
data")
    }
    else{
      T <- T/idx.c
      propT <- TRUE
      countT <- FALSE
    }
  }

  if(propT & countG){
    if(!all(0 <= T && T <= 1)){
      stop("column proportions are not within [0,1] - please respecify
data")
    }
    else{
      G <- G/idx.r
      propG <- TRUE
      countG <- FALSE
    }
  }

  if(propT & propG){
    idx.r <- apply(G,1,sum)
    idx.c <- apply(T,1,sum)

    if(!all(round(idx.r, digits=3)==1)){
      stop("row marginals are proportions that do not sum to 1 - please
respecify data")
    }

    if(!all(round(idx.c, digits=3)==1)){
      stop("column marginals are proportions that do not sum to 1 -
please respecify data")
    }
  }
  
  out <- list()
  se <- list()
  cov.out <- list()
  for (i in 1:(ncol(T))) {
    lm.out <- lm(T[,i] ~ G - 1, ...)
    out[[colnames(T)[i]]] <- lm.out$coef
    se[[colnames(T)[i]]] <- summary(lm.out)$coef[,2]
    cov.out[[colnames(T)[i]]] <- summary(lm.out)$sigma *
      summary(lm.out)$cov.unscaled
    colnames(cov.out[[colnames(T)[i]]]) <- colnames(G)
    rownames(cov.out[[colnames(T)[i]]]) <- colnames(G)
  }
  
  tab <- cbind(out[[1]], out[[2]])
  se.tab <- cbind(se[[1]], se[[2]])
  if (length(out)>2){
    for(i in 3:length(out)){
      tab <- cbind(tab,out[[i]])
      se.tab <- cbind(se.tab,se[[i]])
    }
  }
  colnames(tab) <- colnames(se.tab) <- colnames(T) 
  rownames(tab) <- rownames(se.tab) <- colnames(G)
  
  out <- list(call=match.call(), coefficients=tab, se=se.tab,
              cov.matrices=cov.out)
  
  class(out) <- "eiReg"
  out
}


MTMM <-
function(x, traits, methods) {
  MAT <- cor(x, use="pair")
  startindx <- c(0,seq(methods, traits*methods, by=methods)) + 1
  endindx <- startindx - 1
  SameTrait.vals <- c()
  OUTER.vals <- c()
  for(i in 1:traits) {  
    SameTrait.vals <- c(SameTrait.vals, MAT[startindx[i]:endindx[i+1],startindx[i]:endindx[i+1]][upper.tri(MAT[startindx[i]:endindx[i+1],startindx[i]:endindx[i+1]])])
    OUTER.vals <- c(OUTER.vals, MAT[-startindx[i]:-endindx[i+1],startindx[i]:endindx[i+1]])
  }   
  MethIndx <- seq(1, methods*traits, methods)
  MethodL <- list()
  for(i in 1:methods) {     
    MethodL[[i]] <- data.frame(x[,MethIndx])
    MethIndx <- MethIndx + 1
  }
  SameMethod.vals <- unlist(lapply(MethodL, function(x) cor(x, use="pair")[upper.tri(cor(x, use="pair"))]))
  SameTrait <- mean(SameTrait.vals)
  SameMethod <- mean(SameMethod.vals)
  DiffDiff <- ((sum(OUTER.vals)/2) - sum(SameMethod.vals)) / ((length(OUTER.vals)/2) - length(SameMethod.vals))
  out <- cbind(SameTrait, SameMethod, DiffDiff)
  rownames(out) <- "Results"
  return(out)
}

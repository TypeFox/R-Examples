naive <- function(x, block.vars, id.vars, vcov, n.tr, l2, l1names, valid, validvar, validlb, validub, verbose, ismahal, dist){
  vcd <- as.matrix(x[, block.vars])
  if(nrow(vcd) == 1){
    result <- data.frame(matrix(c(x[, id.vars], rep(NA, n.tr)), ncol = n.tr + 1))
  }
  else{
    vcovi <- solve(vcov)
    if(is.numeric(validvar) && length(validvar)==1){
      p = (length(unique(l1names)) %/% n.tr) + as.integer((length(unique(l1names)) %% n.tr) > 0)
    }
    else{
      p = nrow(x)
    }
    
    out <- .C("naive",
              data = as.double(vcd),
              vec = as.double(dist[row(dist) < col(dist)]),
              nrow = as.integer(nrow(x)),
              ncol = as.integer(ncol(vcd)),
              vcovi = as.double(vcovi),
              ntr = as.integer(n.tr),
              l2 = as.integer(l2),
              l1names = as.integer(as.factor(l1names)),
              valid = as.integer(valid),
              validvar = as.double(validvar),
              validlb = as.double(validlb),
              validub = as.double(validub),
              verbose = as.integer(verbose),
              pairdist = numeric(p),
              ismahal = as.integer(ismahal),
              result = integer(p * n.tr), 
              p = as.integer(p))
    
    result <- data.frame(matrix(out$result, ncol=(n.tr), byrow = TRUE), out$pairdist)
  }
  return(result)
}

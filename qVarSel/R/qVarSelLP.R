library(lpSolveAPI)
library(Rcpp)

locL2 <- function(i, j, k, m, p){
  # assume 111, 121, 131,   1p1, 211, 221, ....2p1, 311, ...., mp1
  # assume 112, 122, 132,   1n2, 212, 222, ....2n2, 312, ...., mp2
  # eccetera
  # m = numero unità statistiche
  lb1 = (k-1)*m*p
  lb2 = (i - 1)*p
  lb3 = j
  loc = lb1 + lb2 + lb3
  return(loc)
}


# This is the covering version of the VarSelModel
qVarSelLP<- function(d, q, binary = FALSE, write = FALSE){
  m = dim(d)[1]
  p = dim(d)[2]
  n = dim(d)[3]
  ncol = m*p*n + n
  vnames = rep(" ",ncol)
  for (i in 1:m){
    for (j in 1:p){
      for (k in 1:n){
        z = locL2(i, j, k, m, p)
        vnames[z] = paste("w", i, j, k, sep = ",")
      }
    }
  }
  for ( k in 1:n){
    vnames[ (m*p*n) + k] = paste("Z", k, sep = ",")
  }
  lps.model <- make.lp(0, ncol)
  obj <- rep(0, ncol)
  for (i in 1:m){
    for (j in 1:p){
      for (k in 1:n){
        z = locL2(i, j, k, m, p)
        obj[z] = d[i, j, k]
      }
    }
  }
  set.objfn(lps.model, obj)	
  for (i in 1:m){
    for (k in 1:n){
      for (j in 1:p){
        if (j == 1){
          a = 1
          ind = locL2(i, j, k, m, p)
        } else {
          a = c(a, 1)
          ind = c(ind, locL2(i, j, k, m, p))
        }}
      a = c(a, -1)
      ind = c(ind, m*p*n + k )
      add.constraint(lps.model, a, type = "=",
                     0, indices = ind) 
    }
  }
  for (i in 1:m){
    for (j in 1:p){
      for (k in 1:n){
        for (t in 1:n){
          if (t == 1){
            a = 1
            ind = locL2(i, j, t, m, p)
          } else {
            a = c(a, 1)
            ind = c(ind, locL2(i, j, t, m, p))
          } # chiude if
        } # chiude for t
        a[k] = -(q - 1)
        add.constraint(lps.model, a, type = ">=", 0, indices = ind)
      }
    }
  }
  
  for (k in 1:n){
    if (k == 1){
      a = 1
      ind = m*p*n + k
    } else {
      a = c(a, 1)
      ind = c(ind, m*p*n + k)
    }}
  add.constraint(lps.model, a, type = "=", q, indices = ind)
  set.bounds(lps.model, lower = rep(0.0, ncol), upper = rep(1.0, ncol), columns = 1:ncol)
  if (binary == T) {
     set.type(lps.model, (m*p*n + 1):ncol, type = "binary") }
  n1 = m*n
  cn1 = paste("C", 1:n1, sep = "")
  n2 = m*p*n 
  cn2 = paste("D", 1:n2, sep = "")
  dimnames(lps.model) = list( a = c(cn1, cn2, "F1"), b = vnames )
  if ( write == T)
            write.lp(lps.model, filename = "qdistsel.lp", type = "lp", 
            use.names = c("TRUE", "TRUE"))
  lp.control(lps.model, bb.depthlimit = 0)
  sol = solve(lps.model)
  obj.sol = get.objective(lps.model)
  z.sol = get.variables(lps.model)[(m*p*n + 1):ncol]
  sl = list(status = sol, obj = obj.sol, x = z.sol)
  return(sl)
}

qVarSelH <- function(d, q, maxit = 100){
  n = dim(d)[1]
  p = dim(d)[2]
  m = dim(d)[3]
  sl = VarSelH(d, n, p, m, q, maxit)
  return(sl)
}

PrtDist <- function(a, p){
  numunits = dim(a)[1]
  nummedians = dim(p)[1]
  numvariables = dim(a)[2]
  nvarprot = dim(p)[2]
  if (numvariables != nvarprot)
    return("Number of columns incorrect")
  dist <- array(0, dim = c(numunits, nummedians, numvariables))
  for (i in 1:numunits){
    for (j in 1:nummedians){      
      dist[i,j, ] = (a[i, ] - p[j, ])^2.0
    }}
  return(dist)
}

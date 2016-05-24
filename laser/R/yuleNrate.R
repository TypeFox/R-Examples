# 24.4.2013 Klaus Schliep
yule2rate <- function(x, verbose = FALSE, ints = NULL, file = "out_yule2rate.txt") #3 parameters; 2 speciation rates
  {
    checkbasal(x)
    if (!is.numeric(x)) stop("object x not of class 'numeric'")
    x <- rev(sort(x))
    if (is.null(ints))
      stvec <- x[3:(length(x)-2)]
    else
      stvec <- seq(x[3], x[length(x)-2], length.out = ints)
    #stvec is a vector of possible shift times, determined either by
    # observed branching times or by specified intervals.
    
   
    #verbose TRUE prints LH, r1, r2 for each shift rate.
    if (verbose == TRUE)
      cat("i",  "LH",  "r1",  "r2",  "st1\n", file = file, sep = "\t")
    
    res = numeric(5)
    res[1] = -Inf
        
    for (i in 1:(length(stvec)))
    {
      v1 <- yuleint2(x, x[1], stvec[i])
      v2 <- yuleint2(x, stvec[i], 0)
      LHtemp <- v1[2] + v2[2]
      if(LHtemp>res[1]){
        res[1] <- LHtemp
        res[2] <- v1[1]
        res[3] <- v2[1]
        res[4] <- stvec[i]
      }
      if (verbose == TRUE)
        res[1] <- LHtemp
        res[2] <- v1[1]
        res[3] <- v2[1]
        res[4] <- stvec[i]   
        cat(i, res, "\n", file = file, append = TRUE, sep = "\t")
    }
    res[5] <- (-2*res[1]) + 6
    names(res) = c("LH",  "r1",  "r2","st1", "AIC") 
    return(res)
  }


yule3rate <- function(x, ints = NULL, verbose = FALSE, file = "out_yule3rate.txt")
  {
    if (!is.numeric(x)) stop("object x not of class 'numeric'")
    x <- rev(sort(x))
    N <- length(x)+1
    
    if (is.null(ints))
    {
      tv <- x[2:(N-2)]
      tv <- unique(tv)
    }
    else
    {
      inc <- (x[2] - x[length(x)])/ints
      tv <- seq(x[2], length.out = ints, by = -inc)
    }  
    if (verbose == TRUE)
    cat("i",  "LH",  "r1",  "r2", "r3", "st1", "st2\n", file = file, sep = "\t")
  
    res = numeric(7)
    res[1]=-Inf
    
    l <- length(tv)   
    pos = integer(l)
    for(i in 1:l) pos[i] = sum(x>tv[i])  
    
    tv = c(x[1], tv, 0)
    
    ltv = length(tv)
    LH = matrix(NA, ltv, ltv) 
    V =  matrix(NA, ltv, ltv)
    
    for(i in 1:(ltv-1)){        
      for(j in (i+1):ltv){
        v <- yuleint2(x, tv[i], tv[j])
        LH[i,j] = v[2]  
        V[i,j] = v[1]
      }
    } 
    
    
    for(i in 2:(ltv-2)){
      l1 <- LH[1,i]
      for(j in (i+1):(ltv-1)){ 
        l2 <- LH[i,j]
        l3 <- LH[j,ltv]
        LHtemp=l1 + l2 + l3
      
        if (is.finite(LHtemp) == TRUE && LHtemp > res[1]){
          res[1] <- LHtemp
          res[2] <- V[1,i]
          res[3] <- V[i,j]
          res[4] <- V[j,ltv]
          res[5] <- tv[i]
          res[6] <- tv[j]
        }
        
        if (verbose == TRUE){
          res[1] <- LHtemp
          res[2] <- V[1,i]
          res[3] <- V[i,j]
          res[4] <- V[j,ltv]
          res[5] <- tv[i+1]
          res[6] <- tv[j+1]

          cat(i, res, "\n", file = file, append = TRUE, sep = "\t")
        }
      }
      
    }
    res[7] <- (-2*res[1]) + 10
    names(res) = c("LH",  "r1",  "r2", "r3", "st1", "st2", "AIC")   
    res
  }


yule4rate <- function(x, ints = NULL) #7 parameters: 4 speciation rates, 3 shift times
  {
    #uses observed shift times only; computational time very high otherwise#
    if (!is.numeric(x)) stop("object x not of class 'numeric'")
    x <- rev(sort(x))
    N <- length(x)+1
    
    if (is.null(ints))
    {
      tv <- x[2:(N-2)]
      tv <- unique(tv)
    }
    else
    {
      inc <- (x[2] - x[length(x)])/ints
      tv <- seq(x[2], length.out = ints, by = -inc)
    }  

    res = numeric(9)
    res[1]=-Inf
    
    l <- length(tv)   
    pos = integer(l)
    for(i in 1:l) pos[i] = sum(x>tv[i])  
    
    tv = c(x[1], tv, 0)
    
    ltv = length(tv)
    LH = matrix(NA, ltv, ltv) 
    V =  matrix(NA, ltv, ltv)
    
    for(i in 1:(ltv-1)){        
      for(j in (i+1):ltv){
        v <- yuleint2(x, tv[i], tv[j])
        LH[i,j] = v[2]  
        V[i,j] = v[1]
      }
    }  

    
    for(i in 2:(ltv-3)){
      l1 <- LH[1,i]
#      v1 <- V[1,i]
      for(j in (i+1):(ltv-2)){
        l2 <- LH[i,j]
        tmp2 <- l1 + l2
#        v2 <- V[i,j]
        for(k in (j+1):(ltv-1)){
        l3 <- LH[j,k]
        l4 <- LH[k,ltv]
#        v3 <- V[j,k] 
#        v4 <- V[k,ltv]          
#        LHtemp=l1 + l2 + l3 + l4
        LHtemp <- tmp2 + l3 + l4
        if (is.finite(LHtemp) == TRUE && LHtemp > res[1]){
          res[1] <- LHtemp
          res[2] <- V[1,i]
          res[3] <- V[i,j]
          res[4] <- V[j,k]
          res[5] <- V[k,ltv]  
          res[6] <- tv[i]
          res[7] <- tv[j]
          res[8] <- tv[k]
        }

        }
      }      
    }
    res[9] <- (-2*res[1]) + 14
    names(res) = c("LH",  "r1",  "r2", "r3", "r4", "st1", "st2", "st3", "AIC")   
    res
  }    
    


yule5rate <- function (x, ints = NULL) 
{
  if (!is.numeric(x)) 
    stop("object x not of class 'numeric'")
  x <- rev(sort(x))
  N <- length(x) + 1
  if (is.null(ints)) {
    tv <- x[2:(N - 2)]
    tv <- unique(tv)
  }
  else {
    inc <- (x[2] - x[length(x)])/ints
    tv <- seq(x[2], length.out = ints, by = -inc)
  }
  res = numeric(9)
  res[1] = -Inf
  l <- length(tv)
  pos = integer(l)
  for (i in 1:l) pos[i] = sum(x > tv[i])
  tv = c(x[1], tv, 0)
  ltv = length(tv)
  LH = matrix(NA, ltv, ltv)
  V = matrix(NA, ltv, ltv)
  for (i in 1:(ltv - 1)) {
    for (j in (i + 1):ltv) {
      v <- yuleint2(x, tv[i], tv[j])
      LH[i, j] = v[2]
      V[i, j] = v[1]
    }
  }
  for (i in 2:(ltv - 4)) {
    l1 <- LH[1, i]
    #       v1 <- V[1, i]
    for (j in (i + 1):(ltv - 3)) {
      l2 <- LH[i, j]
      tmp2 <- l1 + l2 
      #           v2 <- V[i, j]
      for (k in (j + 1):(ltv - 2)) {
        l3 <- LH[j, k]
        tmp3 <- tmp2 + l3 
        #               v3 <- V[j, k]
        for (m in (k + 1):(ltv - 1)) {
          l4 <- LH[k, m]
          #                 v4 <- V[k, m]
          l5 <- LH[m, ltv]
          #                 v5 <- V[m, ltv]
          #                  LHtemp = l1 + l2 + l3 + l4 + l5
          LHtemp = tmp3 + l4 + l5
          if (is.finite(LHtemp) && LHtemp > res[1]) {
            res[1] <- LHtemp
            res[2] <- V[1, i]
            res[3] <- V[i, j]
            res[4] <- V[j, k]
            res[5] <- V[k, m]
            res[6] <- V[m, ltv]
            res[7] <- tv[i]
            res[8] <- tv[j]
            res[9] <- tv[k]
            res[10] <- tv[m]
          }
        }
      }
    }
  }
  res[11] <- (-2 * res[1]) + 18
  names(res) = c("LH", "r1", "r2", "r3", "r4", "r5", "st1", "st2", "st3", "st4", "AIC")
  return(res)
}

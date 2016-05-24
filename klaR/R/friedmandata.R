friedman.data <- function(setting = 1, p = 6, samplesize = 40, asmatrix = FALSE)
#  setting    :  1, ..., 6
#  p          :  6, 10, 20, 40
#  samplesize :  40, 100
{
  stopifnot(is.element(setting, 1:6), is.element(p, c(6, 10, 20, 40)), samplesize >= 6)
  friedman.mean <- function(setting, group, p)
  #  setting :  1,...,6    
  #  group   :  1,2,3      
  #  p       :  6,10,20,40 
  {
    eigen34 <- function(i, p)
    { 
      return((9 * ((i-1) / (p-1)) + 1)^2) 
    }

    if(group == 1) result <- rep(0, p)
    else switch(setting,
        "1" = {# "equal spherical" 
                if (group == 2) result <- c(3, rep(0, p-1))
                else result <- c(0, 3, rep(0, p-2))
              },
        "2" = {# "unequal spherical" 
                if (group == 2) result <- c(3, rep(0, p-1))
                else result <- c(0, 4, rep(0, p-2))
              },
        "3" = {# "equal highly ellipsoidal, low-var subspace" 
                ev <- apply(matrix(1:p, ncol=1), 1, eigen34, p=p)
                result <- 2.5 * sqrt(ev/p) * ((p - (1:p)) / (p/2 - 1))
                if (group == 3) result <- result * rep(c(-1,1), p/2)
              },
        "4" = {# "equal highly ellipsoidal, high-var subspace" 
                ev <- apply(matrix(1:p,ncol=1), 1, eigen34, p=p)
                result <- 2.5 * sqrt(ev/p) * (((1:p) - 1) / (p/2 - 1))
                if (group == 3) result <- result * rep(c(-1,1), p/2)
              },
        "5" = {# "unequal highly ellipsoidal, low-var subspace" 
                result <- rep(0, p)  # (equal means (!)) 
              },
        "6" = {# "equal highly ellipsoidal, high-var subspace" 
                result <- rep(14/sqrt(p), p)
                if (group == 3) result <- result * rep(c(-1, 1), p/2)
              })
    return(result)
  }

  friedman.cov <- function(setting, group, p)
  #  setting :  1,...,6    
  #  group   :  1,2,3      
  #  p       :  6,10,20,40 
  {
    eigen34 <- function(i,p)
    { 
      return((9*((i-1)/(p-1))+1)^2) 
    }
    eigen56 <- function(i,p,group)
    { 
      if (group == 1) result <- eigen34(i,p)
      else if (group == 2) result <- (9 * ((p-i) / (p-1)) + 1)^2
      else result <- (9 * ((i - (p-1) / 2) / (p-1)) + 1)^2
      return(result) 
    }
    if (setting == 1) result <- rep(1, p)                             # "equal spherical" 
    else if (setting == 2) result <- rep(group^2, p)                  # "unequal spherical" 
    else {
                                                                       # !! Aufbau anders als bei Friedman !! 
       # "equal highly ellipsoidal" (TRUE) or  "unequal highly ellipsoidal" (FALSE):
      if (setting <= 4) result <- apply(matrix(1:p, ncol = 1), 1, eigen34, p = p)  
      else result <- apply(matrix(1:p, ncol = 1), 1, eigen56, p = p, group = group)
    }
    return(result)
  }
  
  sizes <- tabulate(sample(1:3, samplesize - 6, replace = TRUE), nbins = 3) + 2
  data <- NULL
  for (i in 1:3) {
    me <- friedman.mean(setting = setting, group = i, p = p)
    co <- friedman.cov(setting = setting, group = i, p = p)
    rfried <- function(param) {rnorm(sizes[i], param[1], sqrt(param[2]))}
    data <- rbind(data, apply(rbind(me, co), 2, rfried))
  }
  grp <- rep(1:3, sizes)
  colnames(data) <- paste("x", as.character(1:p), sep = "")
  if (asmatrix) result <- cbind("class" = grp, data)
  else result <- cbind.data.frame("class" = factor(grp), data)
  return(result)
}

qqrank <-
function(X, INDEX, alternative = c("two.sided", "less", "greater"), absrank=TRUE
                 , N=NA, b=NA, plotpenalty = TRUE, allowed.error = 0.005) { 
  
  alternative = match.arg(alternative)
  
  Xpop = mean(X)
  
  R = data.frame('n' = tapply(X, INDEX, length),
                 'x' = tapply(X, INDEX, mean))
  
  if (sum(ifelse(X==0|X==1, 0, 1))==0) {
    Vpop = Xpop * (1 - Xpop)
    R = cbind(R,'e' = tapply(X, INDEX, sum)
              ,'v' = R$x * (1 - R$x))
    R = cbind(R,'t' = (R$e - (R$n * Xpop))/sqrt((R$n * Xpop)*(1 - Xpop))
              ,'pval' = as.numeric(NA))
    for (i in 1:length(R[,1])) {    
      R$pval[i] = 1 - binom.test(R$e[i], R$n[i], Xpop, alternative = alternative)$p.value
    }
  } else {
    Npop = as.numeric(length(X))
    Vpop = var(X)
    R = cbind(R,'v' = tapply(X, INDEX, var))
    R = cbind(R,'t' = (R$x - Xpop)/ sqrt((R$v/R$n) + (Vpop/Npop))
              ,'pval' = as.numeric(NA))
    
    if (alternative == "less") {
      R$pval = 1 - (pt(R$t, R$n-1))
      R$t = ifelse(R$t>0, 0, R$t)
    } else if (alternative == "greater") {
      R$pval = 1 - (pt(R$t, R$n - 1, lower.tail = FALSE))
      R$t = ifelse(R$t<0 ,0, R$t)
    } else {
      R$pval = 1 - (2 * pt(-abs(R$t), R$n-1))
    }  
    
  }
  
  if (sum(ifelse(X==0|X==1, 0, 1))==0) {
    test.used = "binomial"
  } else {
    test.used = "welch unequal variance"
  }
  
  
  R = R[order(R$n),]
  
  if (is.na(b)==FALSE) {
    glpenalty = bglpenalty(R$n, N, b, "growth", plotpenalty, allowed.error)
  } else {
    glpenalty = rep(1, length(R$n))
  }
  
  if (absrank==TRUE) {
    R$t = abs(R$t)
  }
  
  R = cbind(R,
            'glpenalty' = glpenalty)
  
  R = cbind(R,'r' = R$t * R$pval * R$glpenalty)
  R = R[order(-R$r, -R$n), ]
  
  
  
  structure(list(rankmatrix = data.frame('n' = R$n
                                         ,'mean' = R$x
                                         ,'sd' = sqrt(R$v)
                                         ,'qqrank' = R$r)
                 ,test.used = test.used
                 ,pop.mean = Xpop
                 ,pop.sd = sqrt(Vpop))
            ,class = 'qqrank')
}

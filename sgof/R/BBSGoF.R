BBSGoF<-function (u, alpha = 0.05, gamma = 0.05, kmin = 2, kmax = min(length(u)%/%10, 100), tol = 10, adjusted.pvalues = FALSE, blocks = NA) 
{

    if (missing(blocks) & adjusted.pvalues == T) {
        stop("blocks argument is required to compute the adjusted p-values")  }
   
    if (missing(u)) { stop("data argument is required") }
    
    if (kmin > kmax) {
        stop("kmax should be larger than kmin") }

    n = length(u)

    if (kmax == n) {
        stop("kmax should be lower than n")
    }

    bbsgof <- function(u, alpha = 0.05, gamma = 0.05, kmin = 2, 
        kmax = min(length(u)%/%10, 100), tol = 10, adjusted.pvalues = FALSE, 
        blocks = NA) {
        k <- blocks
        BBSGoF.ap <- function(u, k) {
            n = length(u)
            prob = vector(length = n)
            ro = vector(length = n)
            beta1 = vector(length = n)
            g1 = vector(length = n)
            sb1 = vector(length = n)
            low = vector(length = n)
            tlow = vector(length = n)
            effects = vector(length = n)
            vb1 = vector(length = n)
            AA = list()
            BB = list()
            nnI = list()
            for (i in 1:n) {
                v = as.numeric(u <= u[i])
                p = mean(v)
                A = vector(length = k)
                n1 = n%/%k
                A <- sapply(1:k, function(m) A <- sum(v[((m - 
                  1) * n1 + 1):(m * n1)]))
                A[k] = sum(v[((k - 1) * n1 + 1):n])
                AA[[i]] <- A
                B = c(rep(n1, k - 1), length(v[((k - 1) * n1 + 
                  1):n]))
                BB[[i]] <- B
                l1 = vector(length = (length(A)))
                L1 <- function(pe, rho) {
                  for (j in 1:length(A)) {
                    if (A[j] == 0) {
                      l1[j] = 0
                    }
                    else {
                      l1[j] <- sum(log(pe + (-rho/(rho - 1)) * 
                        (0:(A[j] - 1))))
                    }
                  }
                  return(sum(l1))
                }
                l2 = vector(length = (length(B - A)))
                L2 <- function(pe, rho) {
                  for (j in 1:length(B - A)) {
                    if ((B[j] - A[j]) == 0) {
                      l2[j] = 0
                    }
                    else {
                      l2[j] <- sum(log(1 - pe + (-rho/(rho - 
                        1)) * (0:(B[j] - A[j] - 1))))
                    }
                  }
                  return(sum(l2))
                }
                l3 = vector(length = (length(B)))
                L3 <- function(pe, rho) {
                  for (j in 1:length(A)) {
                    if (B[j] == 0) {
                      l1[j] = 0
                    }
                    else {
                      l3[j] <- sum(log(1 + (-rho/(rho - 1)) * 
                        (0:(B[j] - 1))))
                    }
                  }
                  return(sum(l3))
                }
                L <- function(pe, rho) {
                  L <- L1(pe, rho) + L2(pe, rho) - L3(pe, rho)
                  return(L)
                }
                Max <- function(x) {
                  -L(x[1], x[2])
                }
                opt <- nlminb(c(0.01, 0.01), Max, lower = rep(0.001, 
                  2), upper = rep(0.999, 2), control = list(rel.tol = 1e-06))
                prob[i] = opt$par[1]
                ro[i] = opt$par[2]
                l111 = vector(length = (length(A)))
                L111 <- function(pe, rho) {
                  for (j in 1:length(A)) {
                    if (A[j] == 0) {
                      l111[j] = 0
                    }
                    else {
                      l111[j] <- sum((-1)/((pe + (-rho/(rho - 
                        1)) * (0:(A[j] - 1)))^2))
                    }
                  }
                  return(sum(l111))
                }
                l112 = vector(length = (length(B - A)))
                L112 <- function(pe, rho) {
                  for (j in 1:length(B - A)) {
                    if ((B[j] - A[j]) == 0) {
                      l112[j] = 0
                    }
                    else {
                      l112[j] <- sum((-1)/((1 - pe + (-rho/(rho - 
                        1)) * (0:(B[j] - A[j] - 1)))^2))
                    }
                  }
                  return(sum(l112))
                }
                der11 <- function(pe, rho) {
                  der11 <- L111(pe, rho) + L112(pe, rho)
                  return(der11)
                }
                l121 = vector(length = (length(A)))
                L121 <- function(pe, rho) {
                  for (j in 1:length(A)) {
                    if (A[j] == 0) {
                      l121[j] = 0
                    }
                    else {
                      l121[j] <- sum((-(0:(A[j] - 1)))/((pe + 
                        (-rho/(rho - 1)) * (0:(A[j] - 1)))^2))
                    }
                  }
                  return(sum(l121))
                }
                l122 = vector(length = (length(B - A)))
                L122 <- function(pe, rho) {
                  for (j in 1:length(B - A)) {
                    if ((B[j] - A[j]) == 0) {
                      l122[j] = 0
                    }
                    else {
                      l122[j] <- sum((0:(B[j] - A[j] - 1))/((1 - 
                        pe + (-rho/(rho - 1)) * (0:(B[j] - A[j] - 
                        1)))^2))
                    }
                  }
                  return(sum(l122))
                }
                der12 <- function(pe, rho) {
                  der12 <- L121(pe, rho) + L122(pe, rho)
                  return(der12)
                }
                l221 = vector(length = (length(A)))
                L221 <- function(pe, rho) {
                  for (j in 1:length(A)) {
                    if (A[j] == 0) {
                      l221[j] = 0
                    }
                    else {
                      l221[j] <- sum((-((0:(A[j] - 1))^2))/((pe + 
                        (-rho/(rho - 1)) * (0:(A[j] - 1)))^2))
                    }
                  }
                  return(sum(l221))
                }
                l222 = vector(length = (length(B - A)))
                L222 <- function(pe, rho) {
                  for (j in 1:length(B - A)) {
                    if ((B[j] - A[j]) == 0) {
                      l222[j] = 0
                    }
                    else {
                      l222[j] <- sum((-((0:(B[j] - A[j] - 1))^2))/((1 - 
                        pe + (-rho/(rho - 1)) * (0:(B[j] - A[j] - 
                        1)))^2))
                    }
                  }
                  return(sum(l222))
                }
                l223 = vector(length = (length(B)))
                L223 <- function(pe, rho) {
                  for (j in 1:length(B)) {
                    if ((B[j]) == 0) {
                      l223[j] = 0
                    }
                    else {
                      l223[j] <- sum(((0:(B[j] - 1))^2)/((1 + 
                        (-rho/(rho - 1)) * (0:(B[j] - 1)))^2))
                    }
                  }
                  return(sum(l223))
                }
                der22 <- function(pe, rho) {
                  der22 <- L221(pe, rho) + L222(pe, rho) + L223(pe, 
                    rho)
                  return(der22)
                }
                nI <- matrix(-1 * c(der11(prob[i], ro[i]), der12(prob[i], 
                  ro[i]), der12(prob[i], ro[i]), der22(prob[i], 
                  ro[i])), nrow = 2, ncol = 2)
                nnI[[i]] <- nI
                beta1[i] <- log(prob[i]/(1 - prob[i]))
                g1[i] = (exp(2 * beta1[i]))/((1 + exp(beta1[i]))^4)
                vb1[i] <- (solve(matrix(nnI[[i]], 2, 2))[1, 1])/g1[i]
                if (vb1[i] < 0) 
                  sb1[i] <- 0
                else sb1[i] = sqrt(vb1[i])
            }
            low <- beta1 - sb1 * qnorm(1 - u)
            tlow <- (exp(low)/(1 + exp(low))) - u
            effects1 = floor(pmax(n * tlow, 0, na.rm = T))
            effects2 <- numeric(n)
            for (i in 1:n) {
                effects2[i] <- max(which(as.integer(n * ecdf(su)(su)) <= 
                  effects1[i]), 0)
            }
            effects <- pmin(effects1, effects2)
            umax2 <- which.max(effects)
            au2 = rep(1, n)
            out <- outer(effects, as.integer(n * ecdf(u)(u)), 
                ">=")
            ind <- which(sapply(1:n, function(i) ind <- length(effects[out[, 
                i]])) != 0)
            au2[ind] <- sapply(ind, function(i) au2[i] <- min(u[which(as.integer(n * 
                ecdf(u)(u[i])) <= effects)]))
            return(sort(au2))
        }

        v = as.numeric(u <= gamma)
        n = length(v)
        p = mean(v)
        jmin = kmin - 1
        jmax = kmax - 1
        ro = vector(length = jmax)
        prob = vector(length = jmax)
        AA <- list()
        BB <- list()
        nnI = list()
        for (j in jmin:jmax) {
            k = j + 1
            A = vector(length = k)
            n1 = n%/%k
            A <- sapply(1:k, function(i) A <- sum(v[((i - 1) * 
                n1 + 1):(i * n1)]))
            A[k] = sum(v[((k - 1) * n1 + 1):n])
            AA[[j]] = A
            B = c(rep(n1, k - 1), length(v[((k - 1) * n1 + 1):n]))
            BB[[j]] = B
            l1 = vector(length = (length(A)))
            L1 <- function(pe, rho) {
                for (i in 1:length(A)) {
                  if (A[i] == 0) {
                    l1[i] = 0
                  }
                  else {
                    l1[i] <- sum(log(pe + (-rho/(rho - 1)) * 
                      (0:(A[i] - 1))))
                  }
                }
                return(sum(l1))
            }
            l2 = vector(length = (length(B - A)))
            L2 <- function(pe, rho) {
                for (i in 1:length(B - A)) {
                  if ((B[i] - A[i]) == 0) {
                    l2[i] = 0
                  }
                  else {
                    l2[i] <- sum(log(1 - pe + (-rho/(rho - 1)) * 
                      (0:(B[i] - A[i] - 1))))
                  }
                }
                return(sum(l2))
            }
            l3 = vector(length = (length(B)))
            L3 <- function(pe, rho) {
                for (i in 1:length(A)) {
                  if (B[i] == 0) {
                    l1[i] = 0
                  }
                  else {
                    l3[i] <- sum(log(1 + (-rho/(rho - 1)) * (0:(B[i] - 
                      1))))
                  }
                }
                return(sum(l3))
            }
            L <- function(pe, rho) {
                L <- L1(pe, rho) + L2(pe, rho) - L3(pe, rho)
                return(L)
            }
            Max <- function(x) {
                -L(x[1], x[2])
            }
            opt <- nlminb(c(0.01, 0.01), Max, lower = rep(0.001, 
                2), upper = rep(0.999, 2), control = list(rel.tol = 1e-06))
            prob[j] = opt$par[1]
            ro[j] = opt$par[2]
            l111 = vector(length = (length(A)))
            L111 <- function(pe, rho) {
                for (i in 1:length(A)) {
                  if (A[i] == 0) {
                    l111[i] = 0
                  }
                  else {
                    l111[i] <- sum((-1)/((pe + (-rho/(rho - 1)) * 
                      (0:(A[i] - 1)))^2))
                  }
                }
                return(sum(l111))
            }
            l112 = vector(length = (length(B - A)))
            L112 <- function(pe, rho) {
                for (i in 1:length(B - A)) {
                  if ((B[i] - A[i]) == 0) {
                    l112[i] = 0
                  }
                  else {
                    l112[i] <- sum((-1)/((1 - pe + (-rho/(rho - 
                      1)) * (0:(B[i] - A[i] - 1)))^2))
                  }
                }
                return(sum(l112))
            }
            der11 <- function(pe, rho) {
                der11 <- L111(pe, rho) + L112(pe, rho)
                return(der11)
            }
            l121 = vector(length = (length(A)))
            L121 <- function(pe, rho) {
                for (i in 1:length(A)) {
                  if (A[i] == 0) {
                    l121[i] = 0
                  }
                  else {
                    l121[i] <- sum((-(0:(A[i] - 1)))/((pe + (-rho/(rho - 
                      1)) * (0:(A[i] - 1)))^2))
                  }
                }
                return(sum(l121))
            }
            l122 = vector(length = (length(B - A)))
            L122 <- function(pe, rho) {
                for (i in 1:length(B - A)) {
                  if ((B[i] - A[i]) == 0) {
                    l122[i] = 0
                  }
                  else {
                    l122[i] <- sum((0:(B[i] - A[i] - 1))/((1 - 
                      pe + (-rho/(rho - 1)) * (0:(B[i] - A[i] - 
                      1)))^2))
                  }
                }
                return(sum(l122))
            }
            der12 <- function(pe, rho) {
                der12 <- L121(pe, rho) + L122(pe, rho)
                return(der12)
            }
            l221 = vector(length = (length(A)))
            L221 <- function(pe, rho) {
                for (i in 1:length(A)) {
                  if (A[i] == 0) {
                    l221[i] = 0
                  }
                  else {
                    l221[i] <- sum((-((0:(A[i] - 1))^2))/((pe + 
                      (-rho/(rho - 1)) * (0:(A[i] - 1)))^2))
                  }
                }
                return(sum(l221))
            }
            l222 = vector(length = (length(B - A)))
            L222 <- function(pe, rho) {
                for (i in 1:length(B - A)) {
                  if ((B[i] - A[i]) == 0) {
                    l222[i] = 0
                  }
                  else {
                    l222[i] <- sum((-((0:(B[i] - A[i] - 1))^2))/((1 - 
                      pe + (-rho/(rho - 1)) * (0:(B[i] - A[i] - 
                      1)))^2))
                  }
                }
                return(sum(l222))
            }
            l223 = vector(length = (length(B)))
            L223 <- function(pe, rho) {
                for (i in 1:length(B)) {
                  if ((B[i]) == 0) {
                    l223[i] = 0
                  }
                  else {
                    l223[i] <- sum(((0:(B[i] - 1))^2)/((1 + (-rho/(rho - 
                      1)) * (0:(B[i] - 1)))^2))
                  }
                }
                return(sum(l223))
            }
            der22 <- function(pe, rho) {
                der22 <- L221(pe, rho) + L222(pe, rho) + L223(pe, 
                  rho)
                return(der22)
            }
            MM <- function(pe, rho) {
                matrix(c(der11(pe, rho), der12(pe, rho), der12(pe, 
                  rho), der22(pe, rho)), 2, 2)
            }
            nI <- matrix(-1 * c(der11(prob[j], ro[j]), der12(prob[j], 
                ro[j]), der12(prob[j], ro[j]), der22(prob[j], 
                ro[j])), nrow = 2, ncol = 2)
            nnI[[j]] <- nI
        }
        beta1 <- sapply(jmin:jmax, function(j) beta1 <- log(prob[j]/(1 - 
            prob[j])))


        beta2 <- sapply(jmin:jmax, function(j) beta2 <- log(ro[j]/(1 - 
            ro[j])))
        varp <- sapply(jmin:jmax, function(i) varp <- solve(matrix(nnI[[i]], 
            2, 2))[1, 1])
        varo <- sapply(jmin:jmax, function(i) varo <- solve(matrix(nnI[[i]], 
            2, 2))[2, 2])
        g1 = (exp(2 * beta1))/((1 + exp(beta1))^4)
        g2 = (exp(2 * beta2))/((1 + exp(beta2))^4)
        vb1 = varp/g1
        vb2 = varo/g2
        sb1 = sqrt(vb1[vb1 >= 0])
        sb2 = sqrt(vb2[vb2 >= 0])
         
        low<-beta1 -  sb1 * qnorm(1 - alpha)

        tlow.original <- exp(low)/(1 + exp(low)) - gamma
        
        ii=unique(c( which(vb2 <= 0), which(vb1 <= 0)))


if((kmax-kmin+1)==length(ii)){stop("All the blocks in the grid have been removed because they provided negative variances: change kmin and/or kmax")}



        if (length(ii) != 0) 
            tlow = exp(low[-ii])/(1 + exp(low[-ii])) - gamma
        else tlow = exp(low)/(1 + exp(low)) - gamma
        effects1 = pmax(n * tlow, 0, na.rm = T)
        effects <- numeric(length(effects1))
        for (i in 1:length(effects1)) {
            effects[i] <- min(effects1[i], sum(as.integer(n * 
                ecdf(u)(u)) <= effects1[i]))
        }
        Mini = min(tlow, na.rm = T)
        mineffects = floor(max(n * Mini, 0))
        mineffects <- min(mineffects, sum(as.integer(n * ecdf(u)(u)) <= 
            mineffects))
        
         kN=  jmin + which(tlow.original == Mini) 
        SGoF = floor(max(min(n * (mean(v) - gamma) - n * sqrt(mean(v) * 
            (1 - mean(v))/n) * qnorm(1 - alpha) + 1, sum(as.integer(n * 
            ecdf(u)(u)) <= n * (mean(v) - gamma) - n * sqrt(mean(v) * 
            (1 - mean(v))/n) * qnorm(1 - alpha) + 1)), 0))
        su <- sort(u)
        jj <- which(u == 1)
        if (length(jj) != 0) 
            pi0 <- 1
        else pi0 <- min((-1/n) * sum(log(1 - u)),1)
        if (mineffects == 0) {
            FDR_BB <- 0
        }
        else {
            FDR_BB <- round((pi0 * su[mineffects])/(ecdf(u)(su[mineffects])), 
                4)
        }


        S <- sapply(1:jmax, function(j) S <- sum((AA[[j]] - 
            prob[j] * BB[[j]])^2)/(prob[j] * (1 - prob[j])))
        Zvalue <- sapply(1:jmax, function(j) Zvalue <- (S[j] - 
            sum(BB[[j]]))/sqrt(2 * sum(BB[[j]] * (BB[[j]] - 1))))
        pvalue <- sapply(1:jmax, function(j) pvalue <- 1 - 
            pnorm(Zvalue[j]))
        Tarone.pvalue.auto = round(pvalue[kN - 1], 4)
       

       a = (1 - ro[kN - 1]) * prob[kN - 1]/ro[kN - 1]
        b = (1 - ro[kN - 1]) * (1 - prob[kN - 1])/ro[kN - 1]
       
 n.blocks = (jmin:jmax) + 1
        
if (length(ii) == 0) 
            n.blocks = n.blocks
        else n.blocks = n.blocks[-ii]

        if (length(ii) == 0) 
            deleted.blocks = NA
        else deleted.blocks = (ii + jmin)

        if (length(ii) == 0)

            cor = ro[jmin:jmax]
        else cor = ro[jmin:jmax][-ii]

        if (length(ii) == 0) 
      
 Tarone.pvalues = pvalue[jmin:jmax]
  else Tarone.pvalues = pvalue[jmin:jmax][-ii]


        if (adjusted.pvalues == TRUE) {
            return(c(list(Rejections = mineffects, FDR = min(FDR_BB,1), 
                Adjusted.pvalues = BBSGoF.ap(u, blocks), effects = floor(effects), 
                SGoF = SGoF, automatic.blocks = kN, deleted.blocks = deleted.blocks, 
                n.blocks = n.blocks, p = prob[kN - 1], cor = cor, 
                Tarone.pvalues = Tarone.pvalues, Tarone.pvalue.auto = Tarone.pvalue.auto, 
                beta.parameters = c(round(a, 4), round(b, 4)), 
                betabinomial.parameters = c(round(prob[kN - 1], 
                  4), round(ro[kN - 1], 4)), sd.betabinomial.parameters = c(round(sqrt(varp[kN - 
                  jmin]), 4), round(sqrt(varo[kN - jmin]), 4)))))
        }
        else {
            return(c(list(Rejections = mineffects, FDR = min(FDR_BB,1), 
                effects = floor(effects), SGoF = SGoF, automatic.blocks = kN, 
                deleted.blocks = deleted.blocks, n.blocks = n.blocks, 
                p = prob[kN - 1], cor = cor, Tarone.pvalues = Tarone.pvalues, 
                Tarone.pvalue.auto = Tarone.pvalue.auto, beta.parameters = c(round(a, 
                  4), round(b, 4)), betabinomial.parameters = c(round(prob[kN - 
                  1], 4), round(ro[kN -1], 4)), sd.betabinomial.parameters = c(round(sqrt(varp[kN-jmin]), 4), round(sqrt(varo[kN-jmin]), 4)))))}
    }



u<-as.vector(u)
res<-bbsgof(u,alpha,gamma,kmin,kmax,tol,adjusted.pvalues,blocks)
res$data<-sort(u)
res$adjusted.pvalues<-adjusted.pvalues
res$blocks<-blocks

res$n<-length(u)
res$alpha<-alpha
res$gamma<-gamma
res$kmin<-kmin
res$kmax<-kmax

res$tol<-tol
res$call<-match.call()
class(res)<-"BBSGoF"
return(res)
}

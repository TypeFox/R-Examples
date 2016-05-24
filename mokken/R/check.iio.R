"check.iio" <- 
function (X, method="MIIO", minvi = default.minvi, minsize = default.minsize, alpha = .05, item.selection=TRUE, verbose=FALSE)
{
   # CHECK DATA
   X <- check.data(X)

   # READ 
   J <- ncol(X)
   N <- nrow(X)
   ncat <- max(X) + 1 # CHECK VOOR IT EN MS-CPM (WAARSCHIJNLIJK GOED) VOOR MIIO WAARSCHIJNLIJK m <- max(X)

   # DETERMINE METHOD
   if (substr(method,1,1)=="I" || substr(method,1,1)=="i") method="IT" else
   if (substr(method,1,2)=="MS"|| substr(method,1,2)=="ms" ||substr(method,1,1)=="C" ||substr(method,1,1)=="c") method="MSCPM" else method <- "MIIO"
   
   # DETERMINE MINVI and MINSIZE
   default.minvi <- ifelse(method=="MIIO",(ncat-1)*.03,.03)
   default.minsize <- ifelse(N > 500, floor(N/10), floor(N/5))
   default.minsize <- ifelse(N <= 250, floor(N/3), default.minsize)
   default.minsize <- ifelse(N < 150, 50, default.minsize)
   # minvi <- default.minvi
   # minsize <- default.minsize
    
   # STOP IF THERE ARE NOT ENOUGH ITEMS OR RESPONDENTS
   if (N < minsize) stop("Sample size less than Minsize")
   if (J < 3) stop("Less than 3 items. Restscore cannot be computed")

   # INITIAL VALUES
   results <- list()
   vi.matrix <- matrix(0,J,J)
   item.order <- rev(order(colMeans(X)))
   if(is.null(dimnames(X)[[2]])) dimnames(X)[[2]] <- paste("V", 1:ncol(X),sep="")
   X <- X[, item.order]
   I.labels <- dimnames(X)[[2]]
   dimnames(vi.matrix) <- list(I.labels,I.labels)
   Hi <- coefH(X, FALSE)$Hi
   g <- 0
   gg <- 0
   h <- 0
   i <- 1
   j <- 2
   k <- 0
   item.mean <- colMeans(X)
   
   # METHOD IT
   if (method=="IT"){
     ii <- jj <- kk <- 0
     uv <- matrix(0,nrow=choose(ncat,2),ncol=2)
     for (ii in 0:(ncat-2)) for (jj in (ii+1):(ncat-1)) {kk <- kk + 1; uv[kk,] <- c(ii,jj)}
     for (i in 1:(J - 1)) {
        R <- as.matrix(X) %*% (matrix(1, J, J) - diag(J)) - X[,i] # R_(ij) for j != i
        R[, i] <- 0
        for (j in (i + 1):J) {
           k <- k + 1
           rvm <- nrow(uv) + 1 #D#
           violation.matrix <- matrix(0, nrow = rvm, ncol = 8)
           dimnames(violation.matrix) <- list(c(t(paste("P(X",i,"=",uv[,2],",X",j,"=",uv[,1],") P(X",i,"=",uv[,1],",X",j,"=",uv[,2],")" ,sep="")),"Total"), #D#
                                             c("#ac", "#vi", "#vi/#ac", "maxvi", "sum", "sum/#ac","X^2 max", "#X^2 sig"))
           results[[k]] <- list()
           results[[k]][[1]] <- list()
           results[[k]][[1]][1] <- I.labels[i]
           results[[k]][[1]][2] <- I.labels[j]
           sorted.R <- sort(R[, j])
           group <- max(which(sorted.R == sorted.R[minsize]))
           repeat {
              if (N - max(group) < minsize) break
              group <- c(group, max(which(sorted.R == sorted.R[minsize + max(group)])))
           } 
           group <- group[-length(group)]
           summary.matrix <- matrix(nrow = length(group) + 1, ncol = 5 + 2*rvm)
           dimnames(summary.matrix)[[2]] <- c("Group", "Lo","Hi", "N", "n", paste("E(X", i, ")", sep = ""), paste("E(X",j, ")", sep = ""), paste("P(X",i,"=",uv[,1],",X",j,"=",uv[,2],")",sep=""), paste("P(X",i,"=",uv[,2],",X",j,"=",uv[,1],")" ,sep=""))
           summary.matrix[, 1] <- 1:nrow(summary.matrix)
           summary.matrix[, 4] <- c(group, N) - c(0, group)
           group <- c(sorted.R[group], max(sorted.R))
           L <- length(group)
           summary.matrix[, 3] <- group
           summary.matrix[, 2] <- c(min(sorted.R), group[-L] + 1)
           member <- apply(1 - outer(R[, j], group, "<="), 1, sum) + 1
           for (g in 1:L) {
              summary.matrix[g, 6] <- mean(X[member == g, i])
              summary.matrix[g, 7] <- mean(X[member == g, j])
              u <- 0
              for (u in 1:nrow(uv)){ 
                summary.matrix[g,7+u]           <- length(X[member==g & X[,i]==uv[u,1] & X[,j]==uv[u,2],])/J
                summary.matrix[g,7+nrow(uv)+ u] <- length(X[member==g & X[,i]==uv[u,2] & X[,j]==uv[u,1],])/J 
              }  
           }# end g-loop
           summary.matrix[,5] <- apply(summary.matrix[,8:(7+2*nrow(uv))],1,sum)
           results[[k]][[2]] <- summary.matrix
           ac <- violation.matrix[1:(rvm - 1), 1] <- L - 1
           for (g in 1:nrow(uv)) {
              n.uv <- summary.matrix[, 7 + g]
              n.vu <- summary.matrix[, 7 + nrow(uv) + g]
                Ng <- summary.matrix[,4]
                d <- (n.uv - n.vu)/Ng
                d[d <= minvi] <- 0
                vi <- length(d[d > minvi/2])
                sum.vi <- sum(d)
                chi.statistic <- rep(0, L)
                chi.pvalue <- rep(1, L)
                if (any(d > 0)) {
                    for (gg in 1:L) {
                    if (d[gg] > 0) {
                        e <- (n.uv[gg] + n.vu[gg])/2  
                        chi.statistic[gg] <- (n.uv[gg] - e)^2/e + (n.vu[gg] - e)^2/e
                        chi.pvalue[gg] <- 1- pchisq(chi.statistic[gg],1)
                    }#endif
                    }#end gg-loop
                }#endif
                violation.matrix[g, 2:8] <- c(vi,vi/ac, max(d), sum.vi, sum.vi/ac, max(chi.statistic),sum(chi.pvalue < alpha))
              vi.matrix[i,j] <- vi.matrix[i,j] + sum(chi.pvalue < alpha)
           }#end g-loop
           violation.matrix[rvm,c(1,2,5,8)] <- apply(violation.matrix[1:(rvm-1),c(1,2,5,8)],2,sum)
           violation.matrix[rvm,3] <- violation.matrix[rvm,2]/violation.matrix[rvm,1]
           violation.matrix[rvm,6] <- violation.matrix[rvm,5]/violation.matrix[rvm,1]
           violation.matrix[rvm,c(4,7)] <- apply(violation.matrix[1:(rvm-1),c(4,7)],2,max)
           results[[k]][[3]] <- violation.matrix
        }#end j-loop
     }#end i-loop
   }#end if (method IT)


   # METHOD MS-CPM
   if (method=="MSCPM"){
      for (i in 1:(J - 1)) {
         R <- as.matrix(X) %*% (matrix(1, J, J) - diag(J)) - X[,i]
         R[, i] <- 0
         for (j in (i + 1):J) {
            k <- k + 1
            rvm <- (ncat - 1) + 1
            violation.matrix <- matrix(0, nrow = rvm, ncol = 8)
            dimnames(violation.matrix) <- list(c(t(paste("P(X",i, ">=", 1:(ncat - 1), ")  P(X",j, ">=", 1:(ncat - 1), ")", sep = "")), "Total"), 
                                               c("#ac", "#vi", "#vi/#ac", "maxvi", "sum", "sum/#ac","zmax", "#zsig"))
            results[[k]] <- list()
            results[[k]][[1]] <- list()
            results[[k]][[1]][1] <- I.labels[i]
            results[[k]][[1]][2] <- I.labels[j]
            sorted.R <- sort(R[, j])
            group <- max(which(sorted.R == sorted.R[minsize]))
            repeat {
               if (N - max(group) < minsize) break
               group <- c(group, max(which(sorted.R == sorted.R[minsize + max(group)])))
            } 
            group <- group[-length(group)]
            summary.matrix <- matrix(nrow = length(group) + 1, ncol = 6 + 2 * (ncat - 1))
            dimnames(summary.matrix)[[2]] <- c("Group", "Lo","Hi", "N", paste("E(X", i, ")", sep = ""), paste("E(X",j, ")", sep = ""), paste("P(X", i, ">=", 1:(ncat -1), ")", sep = ""), paste("P(X", j, ">=", 1:(ncat -1), ")", sep = ""))
            summary.matrix[, 1] <- 1:nrow(summary.matrix)
            Ng <- summary.matrix[, 4] <- c(group, N) - c(0, group)
            group <- c(sorted.R[group], max(sorted.R))
            L <- length(group)
            summary.matrix[, 3] <- group
            summary.matrix[, 2] <- c(min(sorted.R), group[-L] + 1)
            member <- apply(1 - outer(R[, j], group, "<="), 1, sum) + 1
             
            for (g in 1:L) {
               summary.matrix[g, 5] <- mean(X[member == g, i])
               summary.matrix[g, 6] <- mean(X[member == g, j])
               freqi <- tabulate(X[member == g, i] + 1, ncat)
               freqj <- tabulate(X[member == g, j] + 1, ncat)
               cum.freqi <- rev(cumsum(rev(freqi))/Ng[g])
               cum.freqj <- rev(cumsum(rev(freqj))/Ng[g])
               summary.matrix[g, 7:(5 + ncat)] <- cum.freqi[2:ncat]
               summary.matrix[g, (6 + ncat):(4 + 2 * ncat)] <- cum.freqj[2:ncat]
            }# end g-loop
            results[[k]][[2]] <- summary.matrix

            ac <- violation.matrix[1:(rvm - 1), 1] <- L - 1
            for (g in 1:(ncat - 1)) {
               p.1 <- summary.matrix[, 6 + g]
               p.2 <- summary.matrix[, 5 + ncat + g]
               d <- p.2 - p.1
               d[d <= minvi] <- 0
               vi <- length(d[d > minvi/2])
               sum.vi <- sum(d)
               z.statistic <- rep(0, L)
               if (any(d > 0)) {
                 for (gg in 1:L) {
                   if (d[gg] > 0) {
                     Xgg <- X[member == gg, ]
                     f.01 <- length(which(Xgg[, i]>=g& Xgg[, j] < (g+ncat)))
                     f.10 <- length(which(Xgg[, i]< g& Xgg[, j]>= (g+ncat)))
                     f.k <- min(f.01, f.10)
                     f.n <- f.01 + f.10
                     f.b <- ((2 * f.k + 1 - f.n)^2 - 10 * f.n)/(12 * f.n)
                     z.statistic[gg] <- abs(sqrt(2 * f.k + 2 + f.b) - sqrt(2 * f.n - 2 * f.k + f.b))
                   }#endif
                 }#end gg-loop
               }#endif
               violation.matrix[g, 2:8] <- c(vi,vi/ac, max(d), sum.vi, sum.vi/ac, max(z.statistic),length(z.statistic[abs(z.statistic) > qnorm(1-alpha)]))
               vi.matrix[i,j] <- vi.matrix[i,j] + length(z.statistic[abs(z.statistic) > qnorm(1-alpha)])
            }#end g-loop
            violation.matrix[rvm,c(1,2,5,8)] <- apply(violation.matrix[1:(rvm-1),c(1,2,5,8)],2,sum)
            violation.matrix[rvm,3] <- violation.matrix[rvm,2]/violation.matrix[rvm,1]
            violation.matrix[rvm,6] <- violation.matrix[rvm,5]/violation.matrix[rvm,1]
            violation.matrix[rvm,c(4,7)] <- apply(violation.matrix[1:(rvm-1),c(4,7)],2,max)
            results[[k]][[3]] <- violation.matrix
         }#end j-loop
      }#end i-loop
   }#end if (method MS-CPM)

   # METHOD MIIO
   dich=1
   if (method=="MIIO"){
      for (i in 1:(J-1)){
         R <- as.matrix(X) %*% (matrix(1,J,J) - diag(J)) - X[,i]
         R[,i] <- 0
         for (j in (i+1):J){
            k <- k + 1
            rvm <- 2
            violation.matrix <- matrix(0, nrow = 1, ncol = 8)
            dimnames(violation.matrix) <- list(c(t(paste("E(X",i,")  E(X",j,")", sep = ""))), 
                                               c("#ac", "#vi", "#vi/#ac", "maxvi", "sum", "sum/#ac",ifelse(ncat == dich,"zmax", "tmax"), ifelse(ncat == dich, "#zsig", "#tsig")))
            results[[k]] <- list()
            results[[k]][[1]] <- list()
            results[[k]][[1]][1] <- I.labels[i]
            results[[k]][[1]][2] <- I.labels[j]

            sorted.R <- sort(R[,j])
            group <- max(which(sorted.R==sorted.R[minsize]))
            repeat{
              if(N - max(group) < minsize)break
              group <- c(group,max(which(sorted.R==sorted.R[minsize+max(group)])))
            }
            group <- group[-length(group)]  
            summary.matrix <- matrix(nrow = length(group) + 1, ncol = 8)
            dimnames(summary.matrix)[[2]] <- c("Group", "Lo","Hi", "N", paste("E(X", i, ")", sep = ""), paste("E(X",j, ")", sep = ""),paste("SD(X", i, ")", sep = ""), paste("SD(X",j, ")", sep = ""))
            summary.matrix[, 1] <- 1:nrow(summary.matrix)
            summary.matrix[, 4] <- c(group, N) - c(0, group)
            group <- c(sorted.R[group],max(sorted.R))
            L <- length(group)
            summary.matrix[, 3] <- group
            summary.matrix[, 2] <- c(min(sorted.R), group[-L] + 1)
            member <- apply(1 - outer(R[,j], group, "<="),1,sum) + 1

            for (g in 1:L){
               summary.matrix[g, 5] <- mean(X[member == g, i])
               summary.matrix[g, 6] <- mean(X[member == g, j])
               summary.matrix[g, 7] <- sd(X[member == g, i])
               summary.matrix[g, 8] <- sd(X[member == g, j])
            }# end g-loop
            results[[k]][[2]] <- summary.matrix

            ac <- violation.matrix[1:(rvm - 1), 1] <- L - 1
            E.1 <- summary.matrix[, 5]
            E.2 <- summary.matrix[, 6]
            d <- E.2 - E.1
            d[d <= minvi] <- 0
            vi <- length(d[d > minvi/2])
            sum.vi <- sum(d)
            t.statistic <- rep(0, L)
            t.pvalue <- rep(1, L)
            if (any(d > 0)) {
              for (gg in 1:L) {
                if (d[gg] > 0) {
                   if (ncat > dich){
                      tt <- t.test(X[member==gg,i],X[member==gg,j],alternative="less")
                      t.statistic[gg] <- abs(tt$statistic)
                      t.pvalue[gg] <- tt$p.value
                   } else {
                      Xgg <- X[member == gg, ]
                      f.01 <- length(which(Xgg[, i]>=1& Xgg[, j] < 1)) 
                      f.10 <- length(which(Xgg[, i]< 1& Xgg[, j] >= 1))
                      f.k <- min(f.01, f.10)
                      f.n <- f.01 + f.10
                      f.b <- ((2 * f.k + 1 - f.n)^2 - 10 * f.n)/(12 * f.n)
                      t.statistic[gg] <- abs(sqrt(2 * f.k + 2 + f.b) - sqrt(2 * f.n - 2 * f.k + f.b)) # even though it is a z-statistic
                      t.pvalue[gg] <- 1-pnorm(t.statistic[gg])
                   } 
                }#endif
              }#end gg-loop
            }#endif
            violation.matrix[1, 2:8] <- c(vi,vi/ac, max(d), sum.vi, sum.vi/ac, max(t.statistic),sum(t.pvalue < alpha))
            vi.matrix[i,j] <- vi.matrix[i,j] + sum(t.pvalue < alpha)
           results[[k]][[3]] <- violation.matrix
         }# end j-loop
      }# end i-loop
   }# end if(method="MIIO") 


   vi.matrix <- vi.matrix + t(vi.matrix)
   vi.matrix.a <- vi.matrix
   VI <- matrix(apply(sign(vi.matrix),1,sum)) 
   items.removed <- NULL
   if (item.selection){
     repeat{
       nvi <- apply(sign(vi.matrix.a),2,sum)
       if (sum(nvi)==0) break
       maxvi <- which(nvi ==max(nvi))
       if(length(maxvi) > 1){ 
          H <- rep(0,length(maxvi))
          for (i in 1:length(maxvi)) H[i] <- coefH(X[,-c(items.removed,maxvi[i])],FALSE)$H + rnorm(1,0,1e-5)
          maxvi <- maxvi[which(H==max(H))[1]]
       }# end if 
       vi.matrix.a[maxvi,] <- vi.matrix.a[,maxvi] <- 0
       items.removed <- c(items.removed,maxvi)
       VI. <- matrix(apply(sign(vi.matrix.a),1,sum)) 
       VI.[items.removed,] <- NA
       VI <- cbind(VI,VI.)
     }
   }# end if(item.selection)  
   dimnames(VI) <- list(I.labels,paste("step",1:ncol(VI)))
   if(verbose) print(VI)
   
   # Computation of HT

   if (length(items.removed) > 0) X <- X[,-items.removed]
   HT <- coefHT(X)
   
   iio.list <- list(results = results, violations = VI, items.removed = items.removed, Hi = Hi, HT = HT, method = method, item.mean = item.mean, m = ncat-1)
   class(iio.list) <- "iio.class"
   return(iio.list)
}

 
coefHT <- function(Y){
   eq.var <- apply(Y, 1, sd)
   Y <- Y[eq.var > 0, ]
   N <- nrow(Y)
   S    <- try(var(t(Y)) , silent = TRUE)
   Smax <- try(var(apply(t(Y), 2, sort)), silent = TRUE)
  
   if (class(S)!="try-error" & class(Smax)!="try-error"){
      diag(S) <- diag(Smax) <- 0
      HT <- sum(S)/sum(Smax)
   } else {  
      Y    <- Y[sample(1:N,min(N,5000)),]
      S    <- try(var(t(Y)) , silent = TRUE)
      Smax <- try(var(apply(t(Y), 2, sort)), silent = TRUE)
      if (class(S)!="try-error" & class(Smax)!="try-error"){
          diag(S) <- diag(Smax) <- 0
          HT <- sum(S)/sum(Smax)
      } else {
           HT <- NA      
      }
   }    
   return(HT)
}

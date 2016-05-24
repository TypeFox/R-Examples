"check.restscore" <-
function(X, minvi = .03, minsize = default.minsize){

  X <- check.data(X)
  J <- ncol(X)
  N <- nrow(X)
  m <- max(X) + 1
  I.labels <- dimnames(X)[[2]]
  if(length(I.labels)==0) I.labels <- paste("C",1:ncol(X))
  default.minsize <- ifelse(N > 500, floor(N/10), floor(N/5))
  default.minsize <- ifelse(N <= 250, floor(N/3), default.minsize)
  default.minsize <- ifelse(N <  150, 50, default.minsize)

  # Is the sample size large enough?
  if (N < minsize) stop("Sample size less than Minsize")

  # Are there enough items?
  if (J < 3) stop("Less than 3 items. Restscore cannot be computed")

  # Initial computation
  results <- list()
  g <- 0; gg <- 0; h <- 0; i <- 0; j <- 0; k <- 0
  # Non-intersection checks per item
  for (i in 1:(J-1)){
    R <- as.matrix(X) %*% (matrix(1,J,J) - diag(J)) - X[,i]
    R[,i] <- 0
    for (j in (i+1):J){
      k <- k + 1
      rvm <- (m-1)*(m-1)+1
      violation.matrix <- matrix(0,nrow=rvm,ncol=8)
      dimnames(violation.matrix) <-  list(
       c(t(outer(paste("P(X",i,">=",1:(m-1),")",sep=""),paste("P(X",j,">=",1:(m-1),")",sep=""),paste)),"Total"),
       c("#ac","#vi","#vi/#ac","maxvi","sum","sum/#ac","zmax","#zsig"))
      results[[k]] <- list()
      results[[k]][[1]] <- list()
      results[[k]][[1]]$FIRST.ITEM <- I.labels[i]
      results[[k]][[1]]$SECOND.ITEM <- I.labels[j]
      sorted.R <- sort(R[,j])
      group <- max(which(sorted.R==sorted.R[minsize]))
      repeat{
        if(N - max(group) < minsize)break
        group <- c(group,max(which(sorted.R==sorted.R[minsize+max(group)])))
      }
      group <- group[-length(group)]

      summary.matrix <- matrix(nrow = length(group)+1,ncol = 6 + 2*(m-1))
      dimnames(summary.matrix)[[2]] <- c("Group", "Lo", "Hi", "N", paste("E(X",i,")",sep=""), paste("E(X",j,")",sep=""), paste("P(X",i,">=",1:(m-1),")",sep=""), paste("P(X",j,">=",1:(m-1),")",sep=""))
      summary.matrix[,1] <- 1:nrow(summary.matrix)
      summary.matrix[,4] <- c(group,N) - c(0,group)
      group <- c(sorted.R[group],max(sorted.R))
      L <- length(group)
      summary.matrix[,3] <- group
      summary.matrix[,2] <- c(min(sorted.R),group[-L]+1)
      member <- apply(1 - outer(R[,j], group, "<="),1,sum) + 1



      for (g in 1:L){
        Ng <- summary.matrix[g,4]
        summary.matrix[g,5] <- mean(X[member==g,i])
        summary.matrix[g,6] <- mean(X[member==g,j])
        freqi <- tabulate(X[member==g,i]+1,m)
        freqj <- tabulate(X[member==g,j]+1,m)
        cum.freqi <- rev(cumsum(rev(freqi))/Ng)
        cum.freqj <- rev(cumsum(rev(freqj))/Ng)
        summary.matrix[g,7:(5+m)] <- cum.freqi[2:m]
        summary.matrix[g,(6+m):(4+2*m)] <- cum.freqj[2:m]
      }
      results[[k]]$SUMMARY.MATRIX <- summary.matrix

      ac <- violation.matrix[1:(rvm-1),1] <- L-1

      z.scores <- matrix(NA,(m-1)^2,L)
      dimnames(z.scores) <- list(t(outer(paste("P(X",i,">=",1:(m-1),")",sep=""),paste("P(X",j,">=",1:(m-1),")",sep=""),paste)),
                                  paste("Rest-score", apply(summary.matrix[,2:3],1,paste, collapse="-")))

#      for (g in 1:(m-1)) for (h in 1:(m-1)){
#         violation.matrix[(g-1)*(m-1)+h,2:8] <- compute.violations(summary.matrix[,6+g],summary.matrix[,5+m+h],L-1,minvi)
#      }  

      hh <- 0
      for (g in 1:(m-1)) for (h in 1:(m-1)){
        hh <- hh + 1
        p.1 <- summary.matrix[,6+g]
        p.2 <- summary.matrix[,5+m+h]
        n <- summary.matrix[,4]
        if(sum(n*p.1) >  sum(n*p.2)) d <- p.2-p.1
        if(sum(n*p.1) <= sum(n*p.2)) d <- p.1-p.2
        d[d <= minvi] <- 0
        vi <- length(d[d > minvi/2])
        sum.vi <- sum(d)
        z <- rep(0,L)
        if(any(d > 0)){
           for (gg in 1:L){
              if(d[gg] > 0){
                 Xgg <-  X[member==gg,]
                 f.01 <- length(Xgg[Xgg[,i] >= g & Xgg[,j] <  h,])/J
                 f.10 <- length(Xgg[Xgg[,i] <  g & Xgg[,j] >= h,])/J
                 f.k <- min(f.01,f.10); f.n <- f.01 + f.10; f.b <- ((2 * f.k + 1 - f.n)^2 - 10 * f.n) /(12 * f.n)
                 z[gg] <- abs(sqrt(2*f.k + 2 + f.b) - sqrt(2*f.n-2*f.k+f.b))          
              }  
           }
        }        
        violation.matrix[(g-1)*(m-1)+h,2:8] <- c(vi,vi/(L-1),max(d),sum.vi,sum.vi/ac,max(z),length(z[abs(z) > qnorm(.95)]))
        z.scores[hh,] <- z
      }  

#     violation.matrix[(rvm-1),2:8] <- compute.violations(,summary.matrix[,6],L-1,minvi)
      if(sum(summary.matrix[,5]) > sum(summary.matrix[,6])) d <- summary.matrix[,6]-summary.matrix[,5]
      if(sum(summary.matrix[,5]) <= sum(summary.matrix[,6])) d <- summary.matrix[,5]-summary.matrix[,6]

      d[d <= minvi] <- 0
      vi <- length(d[d > minvi/2])
      sum.vi <- sum(d)
#      # BEGIN NEW AANPASSEN ALS var(X[member==gg,i]) == var(X[member==gg,j]) == 0  geen t-test mogelijk.
#       T.statistic <- rep(0,L)  # T statistic
#       t.res <- list(NULL)
#       p <- rep(1,L)             # p-value
#       if(any(d > 0)){
#         for (gg in 1:L){
#            if(d[gg] > 0){
#              t.res[[gg]] <- t.test(X[member==gg,i],X[member==gg,j],alternative="less")
#              T.statistic[gg] <- abs(t.res[[gg]]$statistic)
#              p[gg] <- t.res[[gg]]$p.value
#            } else t.res[[gg]] <- 0
#         }
#       }
#       violation.matrix[(rvm-1),2:8] <- c(vi,vi/ac,max(d),sum.vi,sum.vi/ac,max(T.statistic),length(p[p < .05]))
#      # END NEW

      if(rvm > 2){
         violation.matrix[rvm,c(1,2,5,8)] <- apply(violation.matrix[1:(rvm-1),c(1,2,5,8)],2,sum)
         violation.matrix[rvm,c(4,7)] <- apply(violation.matrix[1:(rvm-1),c(4,7)],2,max)
      }else{
         violation.matrix[rvm,c(1,2,4,5,7,8)] <- violation.matrix[1,c(1,2,4,5,7,8)]
      }
      violation.matrix[rvm,3] <- violation.matrix[rvm,2]/violation.matrix[rvm,1]
      violation.matrix[rvm,6] <- violation.matrix[rvm,5]/violation.matrix[rvm,1]
      results[[k]]$VIOLATION.MATRIX <- violation.matrix
#      # BEGIN TMP
      results[[k]]$Z.SCORES <- z.scores
#      # END TMP
   }
 }
 Hi <- coefH(X,FALSE)$Hi
 restscore.list <- list(results=results,I.labels=I.labels,Hi=Hi,m=m)
 class(restscore.list) <- "restscore.class"
 return(restscore.list)
}

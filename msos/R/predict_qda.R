predict_qda <-
function(qd,newx) {
	 newx <- c(newx)
     disc <- NULL
     K <- length(qd$c)
     for(k in 1:K) {
          dk <- -t(newx-qd$Mean[k,])%*%
                  solve(qd$Sigma[k,,],newx-qd$Mean[k,])+qd$c[k]
          disc <- c(disc,dk)
     }
     disc
}

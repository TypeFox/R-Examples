"plot.restscore.class" <-
function(x, item.pairs = all.pairs, plot.ci = TRUE, color.ci = c("orange","yellow"), alpha.ci = .05, ask = TRUE, ...){

cn2n <- function(y) c(rev(diff(rev(y))),y[length(y)])

up.lo.bound.ISRF <- function(n,alpha=.05){
   n[n < 1e-10] <- 1e-10
   n <- matrix(n)
   p <- length(n)
   A1 <- upper.tri(matrix(,p,p),diag=TRUE)*1
   A2 <- cbind(matrix(-1,p-1,1),diag(p-1))
   g.1 <- A1 %*% n
   g <- exp(A2 %*% log(g.1))
   G <- as.numeric(exp(A2 %*% log(g.1))) * (A2 %*% (A1/as.numeric(g.1)))
   ase <- sqrt(diag((G*as.numeric(n)) %*% t(G)))
   z <- qnorm(1 - alpha/2)
   matrix(c(g - z * ase,g + z * ase),p-1,2)
}

  def.par <- par(no.readonly = TRUE)
  J <- length(x$Hi)
  max.item.pairs <- J*(J-1)/2
  all.pairs <- 1:max.item.pairs
  results <- x$results
  m <- x$m
  if (ask==TRUE) par("ask"=TRUE) else par("ask"=FALSE)
  i <- 0; j <- 0

  for (j in item.pairs){
     plot.matrix <- results[[j]][[2]]

     x.labels <- paste(plot.matrix[,2],"-",plot.matrix[,3],sep="")
     if(plot.ci){ 
        cn1 <- plot.matrix[,4] * cbind(1,plot.matrix[,6 + 1:(m-1)])
        n1  <- t(apply(cn1,1,cn2n))
        cn2 <- plot.matrix[,4] * cbind(1,plot.matrix[,6 + (m-1) + 1:(m-1)])
        n2  <- t(apply(cn2,1,cn2n))

        up.lo1 <- apply(n1,1,up.lo.bound.ISRF, alpha.ci)
        up.lo2 <- apply(n2,1,up.lo.bound.ISRF, alpha.ci)
        lo1 <- up.lo1[1:(m-1),]
        lo2 <- up.lo2[1:(m-1),]
        up1 <- up.lo1[m:(2*m-2),]
        up2 <- up.lo2[m:(2*m-2),]
      }   
      mi1 <- t(plot.matrix[,6 + 1:(m-1)])
      mi2 <- t(plot.matrix[,6 + (m-1) + 1:(m-1)])

      plot(plot.matrix[,1],mi1[1,],
      ylim=c(0,1),
      xaxt = 'n',
      xlab = "Rest score group",
      ylab = "Item step response functions",
      type = "n")
      title(paste(results[[j]][[1]][1],"(solid)",results[[j]][[1]][2],"(dashed)"))
      axis(1, at=1:nrow(plot.matrix),labels=x.labels)

      if(m==2){
        if(plot.ci){
           polygon(c((1:length(up1))[!is.na(up1)],rev((1:length(lo1))[!is.na(lo1)])),c(up1[!is.na(up1)],rev(lo1[!is.na(lo1)])),col=color.ci[1], border=NA)
           polygon(c((1:length(up2))[!is.na(up2)],rev((1:length(lo2))[!is.na(lo2)])),c(up2[!is.na(up2)],rev(lo2[!is.na(lo2)])),col=color.ci[2], border=NA)
           lines(plot.matrix[,1],up1, lwd=1, lty=1)
           lines(plot.matrix[,1],up2, lwd=1, lty=2)
           lines(plot.matrix[,1],lo1, lwd=1, lty=1)
           lines(plot.matrix[,1],lo2, lwd=1, lty=2)
        }   
        lines(plot.matrix[,1],mi1, lwd=4, lty=1)
        lines(plot.matrix[,1],mi2, lwd=4, lty=2)
      } 
      if(m>2){
        if(plot.ci) for(i in 1:(m-1)){
           polygon(c((1:length(up1[i,]))[!is.na(up1[i,])],rev((1:length(lo1[i,]))[!is.na(lo1[i,])])),c(up1[i,!is.na(up1[i,])],rev(lo1[i,!is.na(lo1[i,])])),col=color.ci[1], border=NA)
           polygon(c((1:length(up2[i,]))[!is.na(up2[i,])],rev((1:length(lo2[i,]))[!is.na(lo2[i,])])),c(up2[i,!is.na(up2[i,])],rev(lo2[i,!is.na(lo2[i,])])),col=color.ci[2], border=NA)
        }
        if(plot.ci) for(i in 1:(m-1)){
           lines(plot.matrix[,1],up1[i,], lwd=1, lty=1)
           lines(plot.matrix[,1],up2[i,], lwd=1, lty=2)
           lines(plot.matrix[,1],lo1[i,], lwd=1, lty=1)
           lines(plot.matrix[,1],lo2[i,], lwd=1, lty=2)
        }
        for(i in 1:(m-1)) lines(plot.matrix[,1],mi1[i,], lwd=3, lty=1)
        for(i in 1:(m-1)) lines(plot.matrix[,1],mi2[i,], lwd=3, lty=2)
      }
    }  
 invisible()
 par(def.par)
}

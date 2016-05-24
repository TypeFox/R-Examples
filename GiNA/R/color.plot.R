color.plot <- function(x, cex.axis=.7, cex=5){
  
  transp <-  function (col, alpha = 0.5)  {
    res <- apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255, 
                                                  c[3]/255, alpha))
    return(res)
  }
  
  
  nombres <- strsplit(rownames(x),"\\.")
  to.extract <- lapply(nombres, length)
  real.names <- numeric()
  for(h in 1:length(nombres)){
    real.names[h] <- nombres[[h]][to.extract[[h]]]
  }
  
  yylab <- unique(unlist(lapply(nombres, function(x){x[1]})))
  
  real.names <- as.numeric(real.names)
  w <- which(real.names == 1)
  w2 <- c((w-1)[-1], length(real.names))#rev(length(real.names)-(w-1))
  ww <- cbind(w,w2)
  w3 <- list(NA)
  for(f in 1:dim(ww)[1]){
    rt <- seq(ww[f,1],ww[f,2],by=1); 
    w3[[f]] <- rep(f,length(rt))
  }
  ys <- unlist(w3)
  ## which are the initial plants
  x$xx <- real.names
  x$yy <- ys
  
  avera <- list(NA)
  for(u in 1:length(w)){
    prov <- x[which(x[,"yy"] == u),]
    avera[[u]] <- apply(prov, 2, mean)
  }
  
  avera2 <- lapply(avera, function(x){x[15] <- 0; return(x)})
  avera3 <- data.frame(matrix(unlist(avera2), nrow=length(avera), byrow=TRUE))
  names(avera3) <- names(x)
  
  x2 <- rbind(x,avera3)
  with(x2, plot(xx,yy,col=transp(grDevices::rgb(x2[,c("Red","Green","Blue")]), .8), pch=20, main="", cex=cex, yaxt="n", xaxt="n", xlab="Fruits analyzed", ylab="Genotypes analyzed"))
  
  axis(1, seq(0,500,1),  c("Aver",seq(1,500,1)), cex.axis=cex.axis)
  axis(3, seq(0,500,1),  c("Aver",seq(1,500,1)), cex.axis=cex.axis)
  axis(2, seq(1,length(yylab)), yylab, cex.axis=cex.axis)
  axis(4, seq(1,length(yylab)), yylab, cex.axis=cex.axis)
  points(x=avera3$xx, y=avera3$yy, pch=20, cex=cex+2, col=transp(grDevices::rgb(avera3[,c("Red","Green","Blue")]), 1))
}
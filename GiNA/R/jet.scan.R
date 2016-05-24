jet.scan <- function(mydata, var=1){
  
  nombres <- strsplit(rownames(mydata),"\\.")
  to.extract <- lapply(nombres, length)
  real.names <- numeric()
  for(h in 1:length(nombres)){
    real.names[h] <- nombres[[h]][to.extract[[h]]]
  }
  
  yylab <- unique(unlist(lapply(nombres, function(mydata){mydata[1]})))
  
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
  mydata$xx <- real.names
  mydata$yy <- ys
  
  list22 <- list(NA)
  for(g in 1:max(mydata$yy)){
    list22[[g]] <- mydata[which(mydata$yy == g),]
  }
  names(list22) <- yylab
  
  ## -----------------------------------------------------------
  list.jet <- lapply(list22,function(x){y <- x[,var];return(y)})
  list.jet <- lapply(list.jet, sort, decreasing =T)
  m <- max(unlist(lapply(list.jet, length))) # maximum number of fruits
  mai <- max(unlist(lapply(list.jet, max))) # maximum number of fruits
  # fill the gaps for the plants that do not have the same number of fruits
  list.jet2 <- lapply(list.jet, function(x,m){if(length(x) < m){y <- vector(length=m, mode="double"); y[1:length(x)] <- x}else{y <- x}; return(y)}, m)
  z <- matrix(unlist(list.jet2), ncol=length(list.jet2))
  classif <- seq(mai/100, mai, by=mai/100)
  plo.mat <- matrix(0,nrow=length(classif), ncol=length(list.jet))
  #################################################
  # find to which classification they below
  for(i in 1:dim(z)[1]){
    for(j in 1:dim(z)[2]){
      difs <- abs(classif- z[i,j])
      v <- which(difs == min(difs))
      plo.mat[v,j] <- plo.mat[v,j]+1
    }
  }
  plo.mat[1,] <- 0 # missing values were coded as ero and most be removed
  #a <- apply(z,2,function(x){classif-x}) # rows classifs, col fruits
  #apply(abs(a),1, function(x){which(x == max(x))})
  ##################################################
  # function to generate jet colors
  jet.colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red"))
  palette <- jet.colors(25)
  graphics::image(t(plo.mat), col=palette, xaxt = "n", yaxt = "n", main=names(list22[[1]])[var])
  if(length(list.jet) > 1){
    axis(side=3,at=seq(0,1,by=(1/(length(list.jet2)-1))),names(list22), cex.axis=0.5) #above
    axis(side=2,at=seq(.01,1, by=.01),round(classif, digits=2), las=2, cex.axis=0.5) #left
  }else{
    axis(side=3,at=0,names(list22), cex.axis=0.5) # above
    axis(side=2,at=seq(.01,1, by=.01),round(classif, digits=2), las=2, cex.axis=0.5) #left
  }
  z2 <- data.frame(z); names(z2) <- names(list22)
  return(z2)
}
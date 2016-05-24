
TP.prep <- function(markers=NULL, vp.names=NULL, tp.size=seq(50,200,25), method="sim-dissim", npop=300, nelite=1, mutprob=.5, niterations=100){
  
  if(!is.character(vp.names)){
    cat("'vp.names' argument needs to be a vector with names of the plants in character format\n")
    stop()
  }
  if(is.null(markers)){
    cat("'markers' argument is neccesary to create the TP\n")
    stop()
  }
  #if(length(vp.names) > min(tp.size)){
  #  cat("The size of the TP size has to be greater than the VP size\n")
  #  stop()
  #}
  # markers have to have rownames for individuals
  # also colnames are the marker names
  #PCA
  if(method == "sim" | method == "sim-dissim"){
  mark6.centered <- scale(markers,center=T,scale=F)
  cat(paste("Performing principal components with", dim(markers)[2],"markers \n\n"))
  mark.svd <- svd(mark6.centered)
  layout(matrix(1:2,1,2))
  ######## plot1
  po <- mark.svd$d^2/sum(mark.svd$d^2)
  plot(po,ylab="Percent variation explained", xlab="Principal components", col="firebrick", pch=20, las=2, bty="n")#breaks=dim(po)[2]
  grid()
  legend("topright", paste("Var explained 3 PCs = ",round(sum(po[1:3]),3)), bty="n", cex=.5)
  ix <- which(mark.svd$d<1e-9)
  #Project markers onto the PCs
  mark6.PCbasis <- mark6.centered %*% mark.svd$v[,-ix]
  #library(RColorBrewer)
  #col.scheme <- brewer.pal(6,"Set1")
  ####
  
  VP <- which(rownames(markers) %in% vp.names)
  #VP <- sample(1:dim(mark6.PCbasis)[1], 25);
  names(VP) <- rownames(mark6.PCbasis)[VP]
  #or
  #VP <- which(mark6.PCbasis[,1] < -98 & mark6.PCbasis[,2] < 0)
  #VP <- VP[1:25]
  
  ALLDIST <- dist(mark6.PCbasis[,1:2])
  ALLDIST2 <- as.matrix(ALLDIST)
  ALLDIST2[1:5,1:5]
  
  # potential plants to predict the VP based on euclidian distance
  
  ## for each tp size  find the best
  remove.list <- list(NA)
  
  for(k in 1:length(tp.size)){
    
    plot(mark6.PCbasis[,1:2],col="cadetblue", pch=20, cex=1.3, xlab = "PC1", ylab="PC2", bty="n", las=2, main=paste("TP size =",tp.size[k]))
    
    remove <- 1
    cat(paste("Looking for the best",tp.size[k],"individuals to form TP\n"))
    
    starto <- 1
    is.even <- function(x){x %% 2 == 0}
    while(length(remove) < tp.size[k]){
      potential <- apply(as.data.frame(VP),1, function(x,y,z){
        
        if(method=="sim"){
          closer <- (sort(y[x,], decreasing = FALSE))[-1]; best <- closer[1:z];
        }
        else if(method=="sim-dissim"){
          closer1 <- (sort(y[x,], decreasing = FALSE))[-1]; 
          closer2 <- (sort(y[x,], decreasing = TRUE))[-1]; 
          best <- c( closer1[1:z], closer2[1:z]);
        }
        return(names(best))}, y=ALLDIST2, z=starto)
      
      TP <- unique(as.vector(unlist(potential))) # can be duplicated plants predicting
      VPP <- names(VP) # get only the names of the VP
      # remove from TP plants already in the VP
      remove <- setdiff(TP, VPP)
      starto <- starto + 1
    }
    
    actual <- length(remove) - tp.size[k]
    if(actual > 0){
      #remove[(length(remove) - actual + 1):length(remove)]
      xxx <- tp.size[k]
      remove <- remove[1:xxx]
    }
    ###
    ## show selected TP and the VP
    for(i in 1:length(remove)){
      yo <- remove[i]
      points(x=mark6.PCbasis[yo,1],y=mark6.PCbasis[yo,2], col="red", pch=20)
    }
    for(i in 1:length(VPP)){
      yo <- VPP[i]
      points(x=mark6.PCbasis[yo,1],y=mark6.PCbasis[yo,2], col="blue", pch=20)
    }
    legend("topleft", col=c("red","blue"), pch=20, legend = c("Train pop selected","Val pop provided"), bty="n", cex=.7)
    #cat(paste(" Total", length(remove),"individuals selected for TP\n"))
    remove.list[[k]] <- remove
  }
  tp.list <- remove.list
  }
  else if(method == "random"){
    cat(paste("Picking randomly plants from population\n"))
    rest <- (1:dim(markers)[1])[-which(rownames(markers) %in% vp.names)]
    tp.list <- apply(data.frame(tp.size), 1, function(x,y,z){z[sample(y, x)]}, y=rest, z=rownames(markers))
    if(length(tp.size) == 1){
      tp.list <- list(tp.list[,1])
    }
  }
  else if(method == "PEV"){
    mark6.centered <- scale(markers,center=T,scale=F)
    cat(paste("Performing principal components with", dim(markers)[2],"markers \n\n"))
    mark.svd <- svd(mark6.centered, nu=40,nv=40)
    layout(matrix(1:2,1,2))
    ######## plot1
    po <- mark.svd$d^2/sum(mark.svd$d^2)
    plot(po,ylab="Percent variation explained", xlab="Principal components", col="firebrick", pch=20, las=2, bty="n")#breaks=dim(po)[2]
    grid()
    legend("topright", paste("Var explained 3 PCs = ",round(sum(po[1:3]),3)), bty="n", cex=.5)
    ix <- which(mark.svd$d<1e-9)
    #Project markers onto the PCs
    mark6.PCbasis <- mark6.centered %*% mark.svd$v[,-ix]
    candidates<-setdiff(rownames(markers),vp.names)
    LambdaTrait<-1/ncol(markers)
    tp.list <-list(NA)
    for(h in 1:length(tp.size)){
     cat(paste("Performing TP selection based in PEV for pop size=",tp.size[h],sep=""))
     tp.list[[h]] <- PEV(PCAs=mark6.PCbasis,candidates=candidates,Test=vp.names,ntoselect=tp.size[h], npop=npop, nelite=nelite, mutprob=mutprob, niterations=niterations, lambda=LambdaTrait)
    }
  }
  return(tp.list)
}
`weighted_richclub_tm` <-
function(net, NR=1000, seed=NULL, projection.method="Newman", nbins=30){
  # Ensure that the network conforms to the tnet standard
  if(is.null(attributes(net)$tnet))                 net <- as.tnet(net, type="binary two-mode tnet")
  if(attributes(net)$tnet!="binary two-mode tnet")  stop("Network not loaded properly")

  if(!is.null(seed))
    set.seed(as.integer(seed))
  #Internal function: the non-normalised coefficient
  `phi` <- function(net){
    output <- cbind(x=xlevels,num=NaN,den=NaN,y=NaN)
    net <- cbind(net, rc.i=prominence[net[,1],"r"], rc.j=prominence[net[,2],"r"])
    net <- cbind(net, rc=pmin.int(net[,"rc.i"],net[,"rc.j"]))
    Er <- sapply(output[,"x"], function(a) which(net[,"rc"]>a))
    output[,"num"] <- unlist(lapply(Er, function(a) sum(net[a,"w"])))
    net <- net[order(-net[,"w"]),]
    output[,"den"] <- unlist(lapply(Er, function(a) sum(net[1:length(a),"w"])))
    output <- output[,"num"]/output[,"den"]
    return(output)
  }
  dimnames(net)[[2]] <- c("i","p")
  net.2mode <- net;
  net.2mode.list <- split(net.2mode[,2], net.2mode[,1])
  net.1mode <- projecting_tm(net.2mode, method=projection.method)
  
  #Defining prominence
  prominence <- unique(net.2mode[,1])
  prominence <- cbind(node=prominence, r=unlist(lapply(net.2mode.list, length)))
  #Creating log bins
  tmp1 <- prominence[,2]
  tmp1 <- tmp1[order(tmp1)]
  tmp1 <- tmp1[1:(length(tmp1)-1)]
  xlevels <- vector()
  xlevels[1] <- tmp1[1]-0.00001
  xlevels[nbins] <- tmp1[length(tmp1)]-0.00001
  tmp2 <- (log(xlevels[nbins])-log(xlevels[1]))/(nbins-1)
  for(i in 2:(nbins-1))
    xlevels[i] <- exp(log(xlevels[i-1])+tmp2)
  
  #The non-normalised coefficient
  ophi <- data.frame(x=xlevels, y=phi(net.1mode))
  #Calculating phi_NULL
  rphi <- matrix(data=0, nrow=nrow(ophi), ncol=NR)
  for(i in 1:NR) {
    if(i/10 == round(i/10) )
      cat(paste("Random network ", i, "/", NR, " @ ", date(), "\n", sep=""))
    rdm.2mode <- rg_reshuffling_tm(net.2mode)
    rdm.1mode <- projecting_tm(rdm.2mode)
    rphi[,i] <- phi(rdm.1mode)
  }
  #Defining the normalised coefficient
  rho <- data.frame(x=ophi[,"x"], y=0, l99=0, l95=0, h95=0, h99=0)
  rho[,"y"] <- ophi[,"y"]/rowMeans(rphi)
  rho <- rho[!is.na(rho[,"y"]),]
  
  #Finding the significance levels (see my Ph.D. thesis; only if NR>100)
  if(NR>100) {
    for(i in 1:nrow(rho)) {
      rphi[i,] <- rphi[i,order(rphi[i,])]
      rho[i,"l99"] <- rphi[i,round(NR/100*00.5)]/mean(rphi[i,])
      rho[i,"l95"] <- rphi[i,round(NR/100*02.5)]/mean(rphi[i,])
      rho[i,"h95"] <- rphi[i,round(NR/100*97.5)]/mean(rphi[i,])
      rho[i,"h99"] <- rphi[i,round(NR/100*99.5)]/mean(rphi[i,])
    }
  }
  return(rho)
}
`weighted_richclub_w` <-
function(net,rich="k", reshuffle="weights", NR=1000, nbins=30, seed=NULL, directed=NULL){
  # Ensure that the network conforms to the tnet standard
  if (is.null(attributes(net)$tnet))                      net <- as.tnet(net, type = "weighted one-mode tnet")
  if (attributes(net)$tnet != "weighted one-mode tnet")   stop("Network not loaded properly")

  # Set seed
  if(!is.null(seed))
    set.seed(as.integer(seed))
  # Check if network is directed
  if(is.null(directed)) {
    tmp <- symmetrise_w(net, method = "MAX")
    directed <- (nrow(tmp) != nrow(net) | sum(tmp[,"w"]) != sum(net[,"w"]))
  }
  # Internal function: the non-normalised coefficient
  `phi` <- function(net){
    output <- cbind(x=xlevels,num=NaN,den=NaN,y=NaN)
    net <- cbind(net, rc.i=prominence[net[,1],"r"], rc.j=prominence[net[,2],"r"])
    net <- cbind(net, rc=pmin.int(net[,"rc.i"],net[,"rc.j"]))
    Er <- lapply(output[,"x"], function(a) which(net[,"rc"]>a))
    output[,"num"] <- unlist(lapply(Er, function(a) sum(net[a,"w"])))
    net <- net[order(-net[,"w"]),]
    output[,"den"] <- unlist(lapply(Er, function(a) sum(net[1:length(a),"w"])))
    output <- output[,"num"]/output[,"den"]
    return(output)
  }
  # Defining prominence
  prominence <- degree_w(net)
  if(rich=="k") {
    prominence <- cbind(prominence, r=prominence[,"degree"])
  } else if(rich=="s") { 
    prominence <- cbind(prominence, r=prominence[,"output"])
  } else {
    stop("The rich-parameter is not properly specified; only k and s.")
  }
  # Creating log bins
  tmp1 <- prominence[prominence[,"degree"]>0,"r"]
  tmp1 <- tmp1[order(tmp1)]
  tmp1 <- tmp1[1:(length(tmp1)-1)]
  xlevels <- vector()
  xlevels[1] <- tmp1[1]-0.00001
  xlevels[nbins] <- tmp1[length(tmp1)]-0.00001
  tmp2 <- (log(xlevels[nbins])-log(xlevels[1]))/(nbins-1)
  for(i in 2:(nbins-1))
    xlevels[i] <- exp(log(xlevels[i-1])+tmp2)
  # The non-normalised coefficient
  ophi <- data.frame(x=xlevels, y=phi(net))
  # Calculating phi_NULL
  rphi <- matrix(data=0, nrow=nrow(ophi), ncol=NR)
  for(i in 1:NR) {
    if(i/10 == round(i/10) )
      cat(paste("Random network ", i, "/", NR, " @ ", date(), "\n", sep=""))
    rnet <- rg_reshuffling_w(net, option=reshuffle)
    rphi[,i] <- phi(rnet)
  }
  # Defining the normalised coefficient
  rho <- data.frame(x=ophi[,"x"], y=0, l99=0, l95=0, h95=0, h99=0)
  rho[,"y"] <- ophi[,"y"]/rowMeans(rphi)
  rho <- rho[!is.na(rho[,"y"]),]
  
  # Finding the significance levels (see my Ph.D. thesis; only if NR>100)
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


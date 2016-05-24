`infit.outfit` <-
function(Pxji, n.x, resp, x.i){
  
  x.i <- x.i + 1       # handdles that 0 category is needed

  n.i  <- ncol(x.i) 
  j.weight <- n.x*resp
  residual <- x.i
  
  Eji <- array(0,dim=dim(Pxji)[2:3])
  Wji <- array(0,dim=dim(Pxji)[2:3])
  Cji <- array(0,dim=dim(Pxji)[2:3])

  item.in.out <- array(dim=c(dim(Pxji)[3],4))
  person.in.out <- array(dim=c(dim(Pxji)[2],4))

  colnames(item.in.out) <- c("infit", "in.Z", "outfit", "out.Z")

  for(i in 1:n.i){
    steps    <- sum(! is.na(Pxji[,1,i]))   
    stepper  <- 1:(steps+1)

    Pxj <- rbind( 1 - colSums(matrix(Pxji[1:steps,,i],nrow=steps)), Pxji[1:steps,,i])

    kPxj    <- matrix(stepper*Pxj,nrow=nrow(Pxj))
    Eji[,i] <- colSums(kPxj)
    Wji[,i] <- colSums(((outer(stepper,Eji[,i],"-")^2 )*Pxj),na.rm=T)
    Cji[,i] <- colSums(((outer(stepper,Eji[,i],"-")^4 )*Pxj),na.rm=T)
  }
  
  WiSumj <- colSums(j.weight*Wji)
  N.col  <- colSums(j.weight)
  
  qi <- sqrt(colSums(j.weight*(Cji - (Wji^2)))/(WiSumj^2))

  q2i <- sqrt(colSums(j.weight*(Cji/(Wji^2)))/(N.col^2) - (1/N.col))

  residual <- x.i - Eji
  
  # Small but obvious differeces with Winsteps #
  item.in.out[,1] <- colSums((j.weight*residual^2),na.rm=TRUE)/WiSumj
  item.in.out[,2] <- (item.in.out[,1]^(1/3) - 1)*(3/qi) + (qi/3)
  item.in.out[,3] <- colSums(j.weight*(residual^2)/Wji,na.rm=TRUE)/N.col
  item.in.out[,4] <- (item.in.out[,3]^(1/3) - 1)*(3/q2i) + (q2i/3)

  WjSumi <- rowSums(j.weight*Wji)
  N.row  <- rowSums(j.weight)
  
  qj <- sqrt(rowSums(j.weight*(Cji - (Wji^2)))/(WjSumi^2))

  q2j <- sqrt(rowSums(j.weight*(Cji/(Wji^2)))/(N.row^2) - (1/N.row))

  person.in.out[,1] <- rowSums((j.weight*residual^2),na.rm=TRUE)/WjSumi
  person.in.out[,2] <- (person.in.out[,1]^(1/3) - 1)*(3/qj) + (qj/3)
  person.in.out[,3] <- rowSums(j.weight*(residual^2)/Wji,na.rm=TRUE)/N.row
  person.in.out[,4] <- (person.in.out[,3]^(1/3) - 1)*(3/q2j) + (q2j/3)

list(item.in.out,person.in.out)
}


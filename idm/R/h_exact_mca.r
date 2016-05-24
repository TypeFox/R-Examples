h_exact_mca <- function(data1, data2,nchunk, disk = TRUE) {
  ## data1 the starting matrix; data2 the upcoming data
  ## is.update: states if there are one or more updates
  ## nchunk: number of chunks to split data2

  ## This is equivalent to MCA on the indicator matrix 
  ## mjca(data,lambda="indicator")
 
  if(disk==TRUE){
    suppressWarnings(dir.create("./PCstory"))
  }else{
    allCoords=list()
    allCoordsU=list()
    allctr=list()
    allctrU=list()
    allcor=list()
    allcorU=list()
  }
  
  data1 <- data.frame(lapply(data.frame(data1), factor))
  data2 <- data.frame(lapply(data.frame(data2), factor))
  
  mods1 = apply(data1, 2, unique)
  if(is.null(dim(mods1))== FALSE){
    mods1=split(t(mods1),colnames(mods1))
  }    
  
  n.mods1 = sapply(mods1, length)
  
  bigMat = rbind(data1, data2)
  fn <- rep(names(bigMat), unlist(lapply(bigMat, nlevels))) 
  ln <- unlist(lapply(bigMat, levels))
  labs <- paste(fn, ln, sep = ":")
  
  n1 = nrow(data1)
  n = nrow(bigMat)
  Q = ncol(bigMat)
  bigMatFac=apply(bigMat,2,factor)
  bigMatTab=apply(bigMatFac,2,table)
  
  if(is.null(dim(bigMatTab))== F){
    bigMatTab=split(t(bigMatTab),colnames(bigMatTab))
  }
  
  c=bigMatTab[[1]]
  for(j in 2:Q){
    c=c(c,bigMatTab[[j]])
  }
  
  c = c/(n*Q)
  r = rep(1/n,n)
  tZ1 = transform_z(data1, is.weight = T,is.exact=T, r = r, c = c)
  
  sZ1 = tZ1$SZ
  dims=ncol(tZ1$SZ) - Q
  
  if ((length(nchunk) > 1 ) & (sum(nchunk) != nrow(data2))) {
    stop("\n The overall size of given 'nchunk' blocks does not match the size of 'data2'")
  }
  
  eg1 = do_es(sZ1)
  r1 = r[1:n1]
  PC1 = t(t(eg1$v / sqrt(c)) * as.vector((eg1$d)))
  PCu1 = t(t(eg1$u / sqrt(r1)) * as.vector((eg1$d)))
  PC1.ctr=eg1$v^2
  PC1.cor=(PC1^2)/apply(PC1^2,1,sum)
  PCu1.ctr = t(((1/n)*t(PCu1^2))*as.vector(1/(eg1$d)^2) )
  PCu1.cor = (PCu1^2)/apply(PCu1^2,1,sum) 
  
  if(disk==TRUE){
    fnameA=paste("./PCstory/PCstart",1,".txt",sep="")
    fnameB=paste("./PCstory/PCEnd",1,".txt",sep="")
    fnameC=paste("./PCstory/PCstartUnit",1,".txt",sep="")
    fnameD=paste("./PCstory/PCendUnit",1,".txt",sep="")
    fnameE=paste("./PCstory/PCctr",1,".txt",sep="")
    fnameF=paste("./PCstory/PCctrUnit",1,".txt",sep="")
    fnameG=paste("./PCstory/PCcor",1,".txt",sep="")
    fnameH=paste("./PCstory/PCcorUnit",1,".txt",sep="")
    #   write.table(file=fnameA, matrix(0,dim(PC1[,1:dims])[1],dims))
    write.table(file=fnameA, PC1[,1:dims])
    write.table(file=fnameB, PC1[,1:dims])
    #     write.table(file=fnameC, matrix(0,dim(PCu1[,1:dims])[1],dims))
    write.table(file=fnameC, PCu1[,1:dims])
    write.table(file=fnameD, PCu1[,1:dims])
    write.table(file=fnameE, PC1.ctr[,1:dims])
    write.table(file=fnameF, PCu1.ctr[,1:dims])
    write.table(file=fnameG, PC1.cor[,1:dims])
    write.table(file=fnameH, PCu1.cor[,1:dims])
  }
  
  if(disk==FALSE){
    allCoords[[1]] = PC1[,c(1:dims)]
    allCoordsU[[1]] = PCu1[,c(1:dims)]
    allctr[[1]] = PC1.ctr[,c(1:dims)]
    allctrU[[1]] = PCu1.ctr[,c(1:dims)]
    allcor[[1]] = PC1.cor[,c(1:dims)]
    allcorU[[1]] = PCu1.cor[,c(1:dims)]
  }
  
  out.split = mat_split(data2, (nchunk))
  mat.story = out.split$splitMat
  
  #if block sizes are given, switch back to number
  if (length(nchunk) > 1) {
    nchunk = length(nchunk)
  } 
  
  for (q in 1:length(mat.story)) {
    mat.chu = data.matrix(mat.story[[q]])
    n.chu = nrow(mat.chu)
    mods.up = apply(mat.chu, 2, unique)
    if(is.null(dim(mods.up))== F){
      mods.up=split(t(mods.up),colnames(mods.up))
    }    
    
    n.mods.up = sapply(mods.up, length)   
    catDiff = n.mods1 - n.mods.up
    
    if (any(catDiff) != 0) {
      fake.row = fake_row_make(mods1, mods.up, n.mods1, n.mods.up)
      mat.chu = rbind(mat.chu, fake.row)
      tZ2 = transform_z(mat.chu, is.weight = T,is.exact=T, r = r, c = c)
      sZ2 = tZ2$SZ[1:n.chu, ]
    } else {
      tZ2 = transform_z(mat.chu, is.weight = T, is.exact=T,r = r, c = c)
      sZ2 = tZ2$SZ
    }
    
    n2 = n.chu
    eg2 = do_es(sZ2)
    eg12 = add_es(eg1, eg2,method="esm")
    
    if (q > 1) {
      #######################
      PCu1 = PCuall
      PC1 = PCall
      #######################
    }
    
    SCall = eg12$v / sqrt(c)
    PCall = t(t(eg12$v / sqrt(c)) * as.vector((eg12$d)))
    r2 = r[1:(n1+n2)]
    SRall = eg12$u / sqrt(r2)
    PCuall = t(t(eg12$u / sqrt(r2)) * as.vector((eg12$d)))
    
    PCall = sign_match(PC1, PCall)
    PCuall = sign_match(PCu1, PCuall)
    
    PCall.ctr = eg12$v^2
    PCall.cor = (PCall^2)/apply(PCall^2,1,sum)
    PCuall.ctr = t(((1/n)*t(PCuall^2))*as.vector(1/(eg12$d)^2) )
    PCuall.cor = (PCuall^2)/apply(PCuall^2,1,sum)
    
    if(disk==FALSE){      
      allCoords[[q+1]]=PCall[,c(1:dims)]
      allCoordsU[[q+1]]=PCuall[,c(1:dims)]
      allctr[[q+1]] = PCall.ctr[,c(1:dims)]
      allctrU[[q+1]] = PCuall.ctr[,c(1:dims)]
      allcor[[q+1]] = PCall.cor[,c(1:dims)]
      allcorU[[q+1]] = PCuall.cor[,c(1:dims)]
    } 
    
    eg1 = eg12
    n1 = n1 + n2
    r1 = r2
    
    if(disk==TRUE){      
      fnameA=paste("./PCstory/PCstart",q,".txt",sep="")
      fnameB=paste("./PCstory/PCEnd",q+1,".txt",sep="")
      fnameC=paste("./PCstory/PCstartUnit",q,".txt",sep="")
      fnameD=paste("./PCstory/PCendUnit",q+1,".txt",sep="")
      fnameE=paste("./PCstory/PCctr",q+1,".txt",sep="")
      fnameF=paste("./PCstory/PCctrUnit",q+1,".txt",sep="")
      fnameG=paste("./PCstory/PCcor",q+1,".txt",sep="")
      fnameH=paste("./PCstory/PCcorUnit",q+1,".txt",sep="")
      write.table(file=fnameA, PC1[,1:dims])
      write.table(file=fnameB, PCall[,1:dims])
      write.table(file=fnameC, PCu1[,1:dims])
      write.table(file=fnameD, PCuall[,1:dims])
      write.table(file=fnameE, PCall.ctr[,1:dims])
      write.table(file=fnameF, PCuall.ctr[,1:dims])
      write.table(file=fnameG, PCall.cor[,1:dims])
      write.table(file=fnameH, PCuall.cor[,1:dims])
    }
    
  }

  #calculates adjusted inertia
  J = length(eg12$d)
  inertia0 = eg12$d[1:(J-Q)]^4
  alldim <- sum(sqrt(inertia0) >= 1/Q)
  inertia.adj  <- ((Q/(Q-1))^2 * (sqrt(inertia0)[1:alldim] - 1/Q)^2)
  inertia.t    <- (Q/(Q-1)) * (sum(inertia0) - ((J - Q) / Q^2))
    
  out = list()
#  out$colpcoordStart = PC1[,c(1:dims)]
  out$colcoord = SCall[,c(1:dims)]
  out$colpcoord = PCall[,c(1:dims)]
 # out$rowpcoordStart = PCu1[,c(1:dims)]
  out$rowcoord = SRall[,c(1:dims)]
  out$rowpcoord = PCuall[,c(1:dims)]
  out$levelnames = labs 
  out$colctr = PCall.ctr[,c(1:dims)]
  out$colcor = PCall.cor[,c(1:dims)]
  out$rowctr = PCuall.ctr[,c(1:dims)]
  out$rowcor = PCuall.cor[,c(1:dims)]
  out$sv = eg12$d[c(1:(J-Q))]
  out$inertia_e = inertia.adj / inertia.t
  
  if(disk==FALSE){
    out$allrowcoords=allCoordsU
    out$allcolcoords=allCoords
    out$allcolctr=allctr
    out$allrowctr=allctrU
    out$allcolcor=allcor
    out$allrowcor=allcorU
  }
  
  out$rowmass = r
  out$colmass = c
  # out$sZ1 = sZ1
  out$nchunk = nchunk
  out$disk = disk
  out$ff = 0
  out
}

r_live_mca <-  function(data1, data2, nchunk, current_rank,ff = 0,disk = TRUE) {
  
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
  #calculate frequency tables for each variable
  Q = ncol(data1)
  if (missing("current_rank")) {
    #full rank
    current_rank = sum(n.mods1) - Q
  }
  
  # labs=names(unlist(lapply(data1,levels)))
  fn <- rep(names(data1), unlist(lapply(data1, nlevels))) 
  ln <- unlist(lapply(data1, levels))
  labs <- paste(fn, ln, sep = ":")
  
  tZ1 = transform_z(data1, is.weight = F)#, r = r, c = c)
  sZ1 = tZ1$SZ
  Z = tZ1$dZ
  c = tZ1$c
  r = tZ1$r
  
  A <- sZ1
  n1 = nrow(A)
  # center first chunk and calculate eigenspace
  # m should be > k
  m <- dim(A)[1]
  orgn <- colMeans(A)#apply(A,2,sum) / m
  #  Ac <- A - t(as.matrix(orgn) %*% as.matrix(t(rep(1,m))))
  #  Ac = A - rep(orgn, rep.int(nrow(A), ncol(A)))
  Ac = A 
  eg <- fast.svd(Ac,0)
  #  eg <- do_eig(Ac)
  eg$m <- m
  eg$orgn <- orgn
  eg$u = eg$u[,1:current_rank]
  eg$d = eg$d[1:current_rank]
  eg$v = eg$v[,1:current_rank]
  SRall <- eg$u / sqrt(r)
  #column standard coordinates
  SCall <- eg$v / sqrt(c)
  #column principal coordinates
  PC1 <- (eg$v / sqrt(c)) %*% diag(as.vector((eg$d[1:current_rank])))
  #signPC1 = sign(PC1)
  # row principal coordinates
  PCu1 = (eg$u / sqrt(r)) %*% diag(as.vector((eg$d[1:current_rank])))
  dims = current_rank
  
  if ((length(nchunk) > 1 ) & (sum(nchunk) != nrow(data2))) {
    stop("\n The overall size of given 'nchunk' blocks does not match the size of 'data2'")
  }
  PC1.ctr=eg$v^2
  PC1.cor=(PC1^2)/apply(PC1^2,1,sum)
  PCu1.ctr = t(((1/n1)*(PCu1^2))*as.vector(1/(eg$d[1:current_rank])^2) )
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
    # write.table(file=fnameA, matrix(0,dim(PC1[,1:dims])[1],dims))
    write.table(file=fnameA, PC1[,1:dims])
    write.table(file=fnameB, PC1[,1:dims])
    #  write.table(file=fnameC, matrix(0,dim(PCu1[,1:dims])[1],dims))
    write.table(file=fnameC, PCu1[,1:dims])
    write.table(file=fnameD, PCu1[,1:dims])
    write.table(file=fnameE, PC1.ctr[,1:dims])
    write.table(file=fnameF, PCu1.ctr[,1:dims])
    write.table(file=fnameG, PC1.cor[,1:dims])
    write.table(file=fnameH, PCu1.cor[,1:dims])
  }
  
  if(disk==FALSE){
    allCoords[[1]]=PC1[,c(1:dims)]
    allCoordsU[[1]]=PCu1[,c(1:dims)]
    allctr[[1]] = PC1.ctr[,c(1:dims)]
    allctrU[[1]] = PCu1.ctr[,c(1:dims)]
    allcor[[1]] = PC1.cor[,c(1:dims)]
    allcorU[[1]] = PCu1.cor[,c(1:dims)]
  }
  #this is in order to keep the rank of the following chunks stable
  current_rank =  dims  ## This is J-Q
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
      tZ2 = transform_z(mat.chu, is.weight = T, is.exact = F, c=c)#, r = r, c = c)
      sZ2 = tZ2$SZ[1:n.chu, ]
      #    Z = rbind(Z,tZ2$dZ[1:n.chu,])
      c2=tZ2$c
      #  r2=tZ2$r[1:n.chu]    
    } else {
      tZ2 = transform_z(mat.chu, is.weight = T, is.exact = F, c=c)#, r = r, c = c)
      sZ2 = tZ2$SZ[1:n.chu, ]
      #      Z = rbind(Z,tZ2$dZ)
      c2=tZ2$c
      #   r2=tZ2$r[1:n.chu]
    }
    ###### END OF THE NEW CODE ########################
    n2 = n.chu
    n12= n1 + n2
    c12 = (c*n1 + c2*n2)/n12
    r12 = rep(1/n12,n12)
    #   sZ2 = sZ2 - rep(eg$orgn, rep.int(nrow(sZ2), ncol(sZ2)))
    
    eg = add_es(eg,sZ2,current_rank,ff=ff,method="isvd")
    #eg,eg2,current_rank,orgn,ff = 0,method=c("esm","isvd")
    ### coordinate computation
    ## column computation (modalities)
    if (q > 1) {
      #######################
      PCu1 = PCuall
      PC1 = PCall
      #######################
    }
    
    #update column standard coordinates
    SCall <- (eg$v / sqrt(c12))
    #update column principal coordinates
    PCall <- (eg$v / sqrt(c12)) %*% diag(as.vector((eg$d[1:current_rank])))
    #update row standard coordinates
    SRall <- (eg$u / sqrt(r12))
    #update row principal coordinates
    PCuall <- (eg$u / sqrt(r12)) %*% diag(as.vector((eg$d[1:current_rank])))
    #  PCuall <- Z %*% as.matrix(SCall) * Q^(-1)
    
    n1 = n12
    c = c12
    r = r12
    
    PCall = sign_match(PC1, PCall[,1:current_rank])
    PCuall = sign_match(PCu1, PCuall[,1:current_rank])
    
    PCall.ctr=eg$v^2
    PCall.cor=(PCall^2)/apply(PCall^2,1,sum)
    PCuall.ctr =  suppressWarnings(t(((1/n1)*t(PCuall^2))*as.vector(1/(eg$v)^2)))
    PCuall.cor = (PCuall^2)/apply(PCuall^2,1,sum) 
    
    if(disk==FALSE){      
      allCoords[[q+1]]=PCall[,c(1:dims)]
      allCoordsU[[q+1]]=PCuall[,c(1:dims)]
      allctr[[q+1]] = PCall.ctr[,c(1:dims)]
      allctrU[[q+1]] = PCuall.ctr[,c(1:dims)]
      allcor[[q+1]] = PCall.cor[,c(1:dims)]
      allcorU[[q+1]] = PCuall.cor[,c(1:dims)]
    }
    
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
  J = length(eg$d) ###this should be calculated somehat else
  #get (almost) same eigenvalues with exact
  dd = eg$d/sqrt(nchunk+1)
  inertia0 = dd^4
  alldim <- sum(sqrt(inertia0) >= 1/Q)
  inertia.adj <- ((Q/(Q-1))^2 * (sqrt(inertia0)[1:alldim] - 1/Q)^2)
  inertia.t <- (Q/(Q-1)) * (sum(inertia0) - ((J - Q) / Q^2))
  out = list()
  # out$colpcoordStart = PC1[,c(1:dims)]
  out$colpcoord = PCall[,c(1:dims)]/sqrt(nchunk+1)
  #  out$rowpcoordStart = PCu1[,c(1:dims)]
  out$rowcoord = SRall[,c(1:dims)]#/sqrt(nchunk+1)
  out$colcoord = SCall[,c(1:dims)]
  out$rowpcoord = PCuall[,c(1:dims)]/sqrt(nchunk+1)
  out$levelnames=labs 
  out$colctr=PCall.ctr[,c(1:dims)]
  out$colcor=PCall.cor[,c(1:dims)]
  out$rowctr=PCuall.ctr[,c(1:dims)]
  out$rowcor=PCuall.cor[,c(1:dims)]
  out$sv=eg$d
  #out$inertia_e=(eg$d[c(1:alldim)])^2/inertia_t
  out$inertia_e <- inertia.adj / inertia.t
  
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
  out$ff = ff
  out$nchunk = nchunk
  out$disk = disk
  out
  
}   

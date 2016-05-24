add_svd <- function(eg,B,m,current_rank,orgn,ff = 0) {
  out = list()
  #new data block
  B = data.matrix(B) 
  #columns (fixed)
  c = dim(B)[2]
  #num of new rows
  n = dim(B)[1]
  if (missing("current_rank")) {
    #full rank
    print(current_rank)
    current_rank = c
  }
  #for convenience
  ff = 1 -ff
  #get low-rank matrices
  Pk = eg$u[,1:current_rank]
  Sk = eg$d[1:current_rank]
  Qk = eg$v[,1:current_rank]
  if (missing("orgn")) {
    #without orgn data is assumed as zero-mean
    Bc = B
    t = 0
  } else {
    #mean update
    orgnb = colMeans(B) 
    #center data
    Bc = B - rep(orgnb, rep.int(nrow(B), ncol(B)))
    #  print(head(Bc))
    #  Bc <- B - t(as.matrix(orgnb) %*% as.matrix(t(rep(1,n))))
    #    print(head(Bc))
    #account for the variance of the mean
    Bc <- rbind(Bc,t(sqrt((n*m)/(n+m))*as.matrix((orgnb-orgn))))
    orgnc <- (ff*m*orgn + n*orgnb)/(n+ff*m)
    #indicates the extra row/col needed for the mean
    t = 1
    out$orgn <- orgnc
  }
  #QR-decomposition of (I-Qk'Qk)B
  qrstr = qr((diag(c) - Qk%*%t(Qk))%*%t(Bc))
  Qt = qr.Q(qrstr)
  L = qr.R(qrstr)
  
  #Form middle matrix K
  #fix this when k=1
  K = rbind(cbind(ff*diag(Sk),matrix(0,current_rank,c)),cbind(Bc%*%Qk,t(L)))
  #get svd of K
  eg = fast.svd(K)
  #keep current_rank singular values and vectors
  Uk = eg$u[,1:current_rank]
  Sk = eg$d[1:current_rank]
  Vk = eg$v[,1:current_rank]
  #calculate U and V
  Vk = cbind(Qk,Qt)%*%Vk
  Uk = rbind(cbind(Pk,matrix(0,m,n+t)),cbind(matrix(0,n,current_rank+t),diag(n)))%*%Uk
  #update number of columns m processed so far
  m <- m + n
  #but the correct is 
  #  m <- ff*m + n
  
  #after thousands of updates you may need this to preserve orthogonality
  #URQ = qr(Uk)
  #UQ = qr.Q(URQ)
  #UR = qr.R(URQ)
  #VRQ= qr(Vk)
  #VQ = qr.Q(VRQ)
  #VR = qr.R(VRQ)
  #egk <- fast.svd(UR %*% diag(Sk) %*% t(VR))
  #Up = UQ %*% egk$u
  #Vp = VQ %*% egk$v
  #Sk = egk$d
  
  #output
  out$u <- Uk
  out$d <- Sk
  out$v <- Vk
  out$m <- m 
  
  
  out
  
}


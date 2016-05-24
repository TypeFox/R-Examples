## The commands perform the test data of the robust of the
## reconstructions from ridges in wavelet transform
## The command used to generate run of testing 128 run for each db in
## c(25,20,10,5,0) chbat signal with 5 noct, 10 voice, 5 compr
## tt <- RunRec(chbat,5,10,5,128,c(25,20,10,5,0))

## plot.ts(HOWAREYOU)
## cgtHOWAREYOU <- cgt(HOWAREYOU,70,0.01,100)
## clHOWAREYOU <- crc(Mod(cgtHOWAREYOU),nbclimb=1000)
## cfHOWAREYOU <- cfamily(clHOWAREYOU,ptile=0.001)
## image(cfHOWAREYOU$ordered > 0) 
## HOWord <- cfHOWAREYOU$ordered > 0 

robustrec <- function(x, nvoice, freqstep, scale, db, ptile) {
  ## HACK!
  InverseDB <- function(db) 10^(db/10)
  
  lng <- length(x)
  nx <- x
  wn <- rnorm(lng, 0, sqrt(var(nx)/InverseDB(db))) 
  nx <- nx + wn
  
  cgtnx <- cgt(nx, nvoice, freqstep, scale, plot = FALSE)
  crcnx <- crc(Mod(cgtnx), nbclimb = 1000)
  cfnx <- cfamily(crcnx,, ptile = ptile)
  nxordered <- (cfnx$ordered > 0)
  
  nxordered
}

## db is a vector of DB
## number of experiments (nrun) for each db (ndb) and each radius (nr)
## vdb <- c(20,15,10,5,0,-5)
## pvec <- c(0.001,0.001,0.001,0.01,0.01,0.02)
## pr <- c(0,1)
## tt <- RunRec(HOWord,HOWAREYOU,70,0.01,100,1,vdb,pvec,pr)
## MM <- apply(tt,c(1,2),mean)

RunRec <- function(HOWord, x, nvoice, freqstep, scale, nrun, vdb, pvec, pr) {
  ndb <- length(vdb)
  nr <- length(pr)
  Run <- array(0,c(ndb,nr,nrun))
  
  for(ii in 1:nrun) {
    cat("The number of run:",ii,"\n")
    
    for(jj in 1:ndb) {
      db <- vdb[jj]
      ptile <- pvec[jj]
      tt <- robustrec(x,nvoice,freqstep,scale,db,ptile)
      for(kk in 1:nr) {
        rr <- pr[kk]
        Run[jj,kk,ii] <- RidgeDist(tt,HOWord,rr)
      }
    }
  }
  Run
}

## wn <- rnorm(2048,0,sqrt(var(CLICKS)/InverseDB(1)))
## 10*log10(var(CLICKS)/var(wn))
## plot.ts(HOWAREYOU)
## cgtHOWAREYOU <- cgt(HOWAREYOU,70,0.01,100)
## clHOWAREYOU <- crc(Mod(cgtHOWAREYOU),nbclimb=1000)
## cfHOWAREYOU <- cfamily(clHOWAREYOU,ptile=0.001)
## image(cfHOWAREYOU$ordered > 0) 
## HOWord <- cfHOWAREYOU$ordered > 0 

band <- function(x, r) {
  xx <- x
  lng <- length(xx)
  yy <- logical(lng)
  x.rts <- as.ts(xx)
  for(k in 1:r) {
    y.rts <- lag(x.rts,-k) 
    zz <- as.logical(y.rts)
    yy[(1+(k)):lng] <- zz[1:(lng-(k))]
    xx <- yy | xx
    y.rts <- lag(x.rts,k)
    zz <- as.logical(y.rts)
    yy[1:(lng-(k))] <- zz[(1+(k)):lng]
    yy[(lng-(k)+1):lng] <- FALSE
    xx <- yy | xx
  }
  xx
}

Sausage <- function(A, r) {
  if( abs(r) > 0 ) {
    r <- abs(r)
    B <- band(c(A),r)
    dim(B) <- c(dim(A)[1],dim(A)[2])
    C <- band(c(t(B)),r)
    ## C <- band(c(t(A)),r)
    dim(C) <- c(dim(t(A))[1],dim(t(A))[2])
    B <- B | t(C)
  }
  else {
    B <- A
  }
  B
}

## Measure the distance of the two set M1 and M2: [0,1]
## (#(M1 - (M1&Sausage(M2))) +  #(M2 - (M2&Sausage(M1))))/(#M1+#M2)

RidgeDist <- function(m1, m2, r) {
  ms1 <- Sausage(m1,r)
  ms2 <- Sausage(m2,r)
  m3 <- m1 & ms2
  m4 <- m1 & (!m3)
  m5 <- m2 & ms1
  m6 <- m2 & (!m5)
  d <- (sum(m4) + sum(m6))/(sum(m1)+sum(m2))
  d
}

## confidence level:

confident <- function(x, low = 0.05, high = 0.95) {
  ndb <- dim(x)[1]
  nr <- dim(x)[2]
  nrun <- dim(x)[3]
  for(jj in 1:ndb) {
    for(kk in 1:nr) {
      tmp <- x[jj,kk,]
      ttmp <- quantile(c(tmp), c(low,high))
      cat("<dB r> = ", jj, kk, "<low high>=", ttmp, "\n")
    }
  }
}

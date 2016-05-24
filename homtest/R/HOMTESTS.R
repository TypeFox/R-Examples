


# ------------------------------------------------------------------------- #

HW.tests <- function(x,cod,Nsim=500) {

  # x = vettore contenente i dati di tutte le stazioni da sottoporre al test
  # cod = vettore tipo factor associato ad x con gli identificativi delle stazioni
  # Nsim = numero di generazioni

  if (length(x)!=length(cod)) {stop('x and cod must have the same length')}

  fac <- factor(cod)
  ni <- tapply(x,fac,length)
  k <- nlevels(fac)

  Lm <- sapply(split(x,fac),Lmoments)
  rLm <- regionalLmoments(x,fac)

  ti <- as.numeric(Lm[3,])
  t3i <- as.numeric(Lm[4,])
  t4i <- as.numeric(Lm[5,])
  lambda1Reg <- as.numeric(rLm[1])
  lambda2Reg <- as.numeric(rLm[2])
  tauReg <- as.numeric(rLm[3])
  tau3Reg <- as.numeric(rLm[4])
  tau4Reg <- as.numeric(rLm[5])

  V1 <- (sum(ni*(ti-tauReg)^2)/sum(ni))^0.5
  #V2 <- (sum(ni * ((ti - tauReg)^2 + (t3i - tau3Reg)^2))/sum(ni))^0.5
  V2 <- sum(ni * ((ti - tauReg)^2 + (t3i - tau3Reg)^2)^0.5)/sum(ni)

  parkappa <- par.kappa(lambda1Reg,lambda2Reg,tau3Reg,tau4Reg)

  xi <- parkappa$xi
  alfa <- parkappa$alfa
  kappa <- parkappa$k
  hacca <- parkappa$h

  V1s <- rep(NA,Nsim)
  V2s <- rep(NA,Nsim)
  for (i in 1:Nsim) {
    ti.sim <- rep(NA,k)
    t3i.sim <- rep(NA,k)
    for (j in 1:k) {
      campione <- rand.kappa(ni[j],xi,alfa,kappa,hacca)
      campione.ad <- campione/mean(campione)
      lmom <- Lmoments(campione.ad)
      ti.sim[j] <- lmom[3]
      t3i.sim[j] <- lmom[4]
    }
    tauReg.sim <- sum(ni*ti.sim)/sum(ni)
    tau3Reg.sim <- sum(ni*t3i.sim)/sum(ni)
    V1s[i] <- (sum(ni*(ti.sim-tauReg.sim)^2)/sum(ni))^0.5
    #V2s[i] <- (sum(ni * ((ti.sim - tauReg.sim)^2 + (t3i.sim - tau3Reg.sim)^2))/sum(ni))^0.5
    V2s[i] <- sum(ni * ((ti.sim - tauReg.sim)^2 + (t3i.sim - tau3Reg.sim)^2)^0.5)/sum(ni)
  }

  muV1 <- mean(V1s)
  stdV1 <- sd(V1s)
  muV2 <- mean(V2s)
  stdV2 <- sd(V2s)

  H1 <- (V1 - muV1)/stdV1
  H2 <- (V2 - muV2)/stdV2

  output <- c(H1,H2)
  names(output) <- c("H1","H2")

  return(output)
}


# ------------------------------------------------------------------------- #

ADbootstrap.test <- function (x,cod,Nsim=500,index=2) {

  # x = vettore contenente i dati di tutte le stazioni da sottoporre al test
  # fac = vettore tipo factor associato ad x con gli identificativi delle stazioni
  # Nsim = numero di generazioni
  # index = 1: grandezza indice è la media
  #         2: grandezza indice è la mediana (default)

  if (length(x)!=length(cod)) {stop('x and cod must have the same length')}

  fac <- factor(cod)
  ni <- tapply(x,fac,length)
  k <- nlevels(fac)
  N <- sum(ni)

  if(index==1) {
    indexflood <- function(x) {m <- mean(x); return(m)}
  }
  else if(index==2) {
    indexflood <- function(x) {m <- median(x); return(m)}
  }

  # adimensionalizzo
  med <- tapply(x,fac,indexflood)
  regione.adim <- x/unsplit(med,fac)

  # Applico Anderson-Darling
  A2kN <- ksampleA2(regione.adim,fac)

  # bootstrap di Nsim regioni e adimensionalizzazione
  A2kNs <- rep(NA,Nsim)
  for (i in 1:Nsim) {
    regione.simul <- nonparboot(regione.adim)
    med.simul <- tapply(regione.simul,fac,indexflood)
    regione.simul.adim <- regione.simul/unsplit(med.simul,fac)
    A2kNs[i] <- ksampleA2(regione.simul.adim,fac)
  }

  #quantili <- quantile(A2kNs,c(0.90,0.95,0.99))
  #output <- c(A2kN,quantili)
  ecdfA2kNs <- ecdf(A2kNs)
  probabilita <- ecdfA2kNs(A2kN)
  output <- c(A2kN, probabilita)
  names(output) <- c("A2kN","P")

  return(output)
}


# ------------------------------------------------------------------------- #

DK.test <- function (x,cod) {

  # Francesco Laio (2005)
  # INPUT:
  # x = vettore contenente i dati di tutte le stazioni da sottoporre al test
  # cod = vettore associato ad x con gli identificativi delle stazioni

  if (length(x)!=length(cod)) {stop('x and cod must have the same length')}

  fac <- factor(cod)
  Y <- x
  k <- nlevels(fac)
  nn <- tapply(Y,fac,length)
  N <- sum(nn)
  NN <- cumsum(nn)

  X <- matrix(0,nrow=N,ncol=k)
  X[1:nn[1],1] <- 1
  for (j in 2:k) {
    X[(NN[j-1]+1):NN[j],j] <- 1
  }

  Z <- cbind(Y,X)

  U <- Z[(sort(Z[,1],index.return=TRUE))$ix,]
  U[,1] <- (1:N)/N
  AA <- matrix(NA, nrow=max(nn), ncol=k)
  for (j in 1:k) {
    AA[1:nn[j],j] <- U[U[,j+1]==1,1]
  }

  BB=cos(2*pi*AA)
  CC=colSums(BB, na.rm=TRUE)*(2/nn)^0.5
  Ak=sum(CC^2)

  #quantili <- qchisq(c(0.90,0.95,0.99), k-1)
  #output <- c(Ak,quantili)
  probabilita <- pchisq(Ak, df=k-1, lower.tail=TRUE)
  output <- c(Ak, probabilita)
  names(output) <- c("Ak","P")

  return(output)
}


# ------------------------------------------------------------------------- #

ksampleA2 <- function (x,cod) {

  # Francesco Laio (2005)
  # INPUT
  # x = vettore contenente i dati di tutte le stazioni da sottoporre al test
  # cod = vettore associato ad x con gli identificativi delle stazioni

  if (length(x)!=length(cod)) {stop('x and cod must have the same length')}

  fac <- factor(cod)
  Y <- x
  N <- length(Y)
  k <- nlevels(fac)
  nn <- tapply(Y,fac,length)
  NN <- cumsum(nn)

  X <- matrix(0,nrow=N,ncol=k)
  X[1:nn[1],1] <- 1                            #   Y[1:NN[1]]
  for (j in 2:k) {
    X[(NN[j-1]+1):NN[j],j] <- 1                #   Y[(NN[j-1]+1):NN[j]]
  }

  Z <- cbind(Y,X)

  U <- Z[(sort(Z[,1],index.return=TRUE))$ix,]

  M <- as.matrix(cumsum(as.data.frame(U[1:(N-1),2:(k+1)])))

  J <- c(1:(N-1))

  numeratore2 <- (N*M - J %*% t(nn))^2
  denominatore <- t(matrix(J*(N-J), ncol=(N-1), nrow=k, byrow=TRUE))

  sommaj <- colSums(numeratore2/denominatore)/nn

  A2k <- sum(sommaj)/N

  output <- A2k

  return(output)
}


# ------------------------------------------------------------------------- #

nonparboot <- function (z,n=length(z)) {

  # dall'analogo algoritmo di Francesco Laio
  # INPUT:
  # z = original sample
  # n = derived sample length

  m <- length(z)
  R <- round(runif(n, min=0, max=1)*m + 0.5)
  z2 <- sort(z)
  A <- z2[R]
  return(A)
}



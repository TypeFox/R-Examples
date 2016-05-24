mySimulations <- function(
  nSim,
  model,
  n.ind,
  mu.range=c(0,0),
  add.eff1.range=c(0,0),
  dom.eff1.range=c(0,0),
  add.eff2.range=c(0,0),
  dom.eff2.range=c(0,0),
  beta21.range=c(0,0),
  beta1h.range=c(0,0),
  beta2h.range=c(0,0),
  sig2.1.range=c(1,1),
  sig2.2.range=c(1,1),
  sig2.h.range=c(1,1),
  eq.spacing=FALSE,
  cross.type="f2",
  thr,
  peak.dist=2,
  normalize=FALSE)
{
  cor12 <- matrix(NA,nSim,1, 
    dimnames=list(c(1:nSim),c("cor.Y1.Y2")))
  R2s <- matrix(NA,nSim,2,
    dimnames=list(c(1:nSim),c("R2.Y1~Q","R2.Y2~Q")))
  BICs <- matrix(NA,nSim,4,
    dimnames=list(c(1:nSim),c("BIC.1","BIC.2","BIC.3","BIC.4")))
  AICs <- matrix(NA,nSim,4,
    dimnames=list(c(1:nSim),c("AIC.1","AIC.2","AIC.3","AIC.4")))
  z.scores <- matrix(NA,nSim,6,
    dimnames=list(c(1:nSim),c("z.12","z.13","z.14","z.23","z.24","z.34")))
  sig2s <- matrix(NA,nSim,6,
    dimnames=list(c(1:nSim),c("s.12.12","s.13.13","s.14.14","s.23.23","s.24.24","s.34.34")))
  pval.par.cmst.joint.BIC <- matrix(NA,nSim,4,
    dimnames=list(c(1:nSim),c("pval.1","pval.2","pval.3","pval.4")))
  pval.par.cmst.iu.BIC <- matrix(NA,nSim,4,
    dimnames=list(c(1:nSim),c("pval.1","pval.2","pval.3","pval.4")))
  pval.non.par.cmst.iu.BIC <- matrix(NA,nSim,4,
    dimnames=list(c(1:nSim),c("pval.1","pval.2","pval.3","pval.4")))
  pval.par.cmst.joint.AIC <- matrix(NA,nSim,4,
    dimnames=list(c(1:nSim),c("pval.1","pval.2","pval.3","pval.4")))
  pval.par.cmst.iu.AIC <- matrix(NA,nSim,4,
    dimnames=list(c(1:nSim),c("pval.1","pval.2","pval.3","pval.4")))
  pval.non.par.cmst.iu.AIC <- matrix(NA,nSim,4,
    dimnames=list(c(1:nSim),c("pval.1","pval.2","pval.3","pval.4")))
  pval.cit <- matrix(NA,nSim,2,
    dimnames=list(c(1:nSim),c("pval.1","pval.2")))

  if(model=="A"){
    k <- 1
    while(k <= nSim){
      mu <- runif(1,mu.range[1],mu.range[2])
      beta21 <- runif(1,beta21.range[1],beta21.range[2])
      add.eff1 <- runif(1,add.eff1.range[1],add.eff1.range[2])
      dom.eff1 <- runif(1,dom.eff1.range[1],dom.eff1.range[2])
      sig2.1 <- runif(1,sig2.1.range[1],sig2.1.range[2])
      sig2.2 <- runif(1,sig2.2.range[1],sig2.2.range[2])
      Cross <- sim.cross.1(n.ind, mu, beta21, add.eff1, dom.eff1, sig2.1, 
        sig2.2, eq.spacing, cross.type, normalize)
      Cross <- calc.genoprob(Cross, step=2)
      cq <- get.common.qtls(Cross, "y1", "y2", thr, peak.dist)
      if(!is.na(cq[1])){
        print(k)
        cor12[k,] <- cor(Cross$pheno[,1],Cross$pheno[,2])
        aux <- try(cmstTests(Cross, "y1", "y2", Q.chr=cq[1,2], Q.pos=cq[1,3], 
          , , , , cross.type),silent=TRUE)
        if(class(aux) != "try-error"){
          R2s[k,] <- aux$R2
          BICs[k,] <- aux$BIC.stats[1:4]
          AICs[k,] <- aux$AIC.stats[1:4]
          z.scores[k,] <- aux$BIC.stats[5:10]
          sig2s[k,] <- aux$Sig.stats[1:6]
          pval.par.cmst.joint.BIC[k,] <- aux$pvals.par.cmst.joint.BIC
          pval.par.cmst.iu.BIC[k,] <- aux$pvals.par.cmst.iu.BIC
          pval.non.par.cmst.iu.BIC[k,] <- aux$pvals.non.par.cmst.iu.BIC
          pval.par.cmst.joint.AIC[k,] <- aux$pvals.par.cmst.joint.AIC
          pval.par.cmst.iu.AIC[k,] <- aux$pvals.par.cmst.iu.AIC
          pval.non.par.cmst.iu.AIC[k,] <- aux$pvals.non.par.cmst.iu.AIC
        }
        cit.mar <- find.marker(Cross,cq[1,2],cq[1,3])
        LL <- pull.geno(Cross)[,cit.mar]
        GG <- Cross$pheno[,1]
        TT <- Cross$pheno[,2]
        aux2 <- try(citTests(LL, GG, TT),silent=TRUE)
        if(class(aux2) != "try-error")
          pval.cit[k,] <- aux2
        k <- k + 1
      }
    }
  }
  if(model=="B"){
    k <- 1
    while(k <= nSim){
      mu <- runif(1,mu.range[1],mu.range[2])
      beta21 <- runif(1,beta21.range[1],beta21.range[2])
      beta1h <- runif(1,beta1h.range[1],beta1h.range[2])
      beta2h <- runif(1,beta2h.range[1],beta2h.range[2])
      add.eff1 <- runif(1,add.eff1.range[1],add.eff1.range[2])
      dom.eff1 <- runif(1,dom.eff1.range[1],dom.eff1.range[2])
      sig2.1 <- runif(1,sig2.1.range[1],sig2.1.range[2])
      sig2.2 <- runif(1,sig2.2.range[1],sig2.2.range[2])
      sig2.h <- runif(1,sig2.h.range[1],sig2.h.range[2])
      Cross <- sim.cross.2(n.ind, mu, beta21, beta1h, beta2h, add.eff1, 
        dom.eff1, sig2.1, sig2.2, sig2.h, eq.spacing, cross.type, normalize)
      Cross <- calc.genoprob(Cross, step=2)
      cq <- get.common.qtls(Cross, "y1", "y2", thr, peak.dist)
      if(!is.na(cq[1])){
        print(k)
        cor12[k,] <- cor(Cross$pheno[,1],Cross$pheno[,2])
        aux <- try(cmstTests(Cross, "y1", "y2", Q.chr=cq[1,2], Q.pos=cq[1,3], 
          , , , , cross.type),silent=TRUE)
        if(class(aux) != "try-error"){
          R2s[k,] <- aux$R2
          BICs[k,] <- aux$BIC.stats[1:4]
          AICs[k,] <- aux$AIC.stats[1:4]
          z.scores[k,] <- aux$BIC.stats[5:10]
          sig2s[k,] <- aux$Sig.stats[1:6]
          pval.par.cmst.joint.BIC[k,] <- aux$pvals.par.cmst.joint.BIC
          pval.par.cmst.iu.BIC[k,] <- aux$pvals.par.cmst.iu.BIC
          pval.non.par.cmst.iu.BIC[k,] <- aux$pvals.non.par.cmst.iu.BIC
          pval.par.cmst.joint.AIC[k,] <- aux$pvals.par.cmst.joint.AIC
          pval.par.cmst.iu.AIC[k,] <- aux$pvals.par.cmst.iu.AIC
          pval.non.par.cmst.iu.AIC[k,] <- aux$pvals.non.par.cmst.iu.AIC
        }
        cit.mar <- find.marker(Cross,cq[1,2],cq[1,3])
        LL <- pull.geno(Cross)[,cit.mar]
        GG <- Cross$pheno[,1]
        TT <- Cross$pheno[,2]
        aux2 <- try(citTests(LL, GG, TT),silent=TRUE)
        if(class(aux2) != "try-error")
          pval.cit[k,] <- aux2
        k <- k + 1
      }
    }
  }
  if(model=="C"){
    k <- 1
    while(k <= nSim){
      mu <- runif(1,mu.range[1],mu.range[2])
      beta21 <- runif(1,beta21.range[1],beta21.range[2])
      add.eff1 <- runif(1,add.eff1.range[1],add.eff1.range[2])
      dom.eff1 <- runif(1,dom.eff1.range[1],dom.eff1.range[2])
      add.eff2 <- runif(1,add.eff2.range[1],add.eff2.range[2])
      dom.eff2 <- runif(1,dom.eff2.range[1],dom.eff2.range[2])
      sig2.1 <- runif(1,sig2.1.range[1],sig2.1.range[2])
      sig2.2 <- runif(1,sig2.2.range[1],sig2.2.range[2])
      Cross <- sim.cross.3(n.ind, mu, beta21, add.eff1, dom.eff1, add.eff2, 
        dom.eff2, sig2.1, sig2.2, eq.spacing, cross.type, normalize)
      Cross <- calc.genoprob(Cross, step=2)
      cq <- get.common.qtls(Cross, "y1", "y2", thr, peak.dist)
      if(!is.na(cq[1])){
        print(k)
        cor12[k,] <- cor(Cross$pheno[,1],Cross$pheno[,2])
        aux <- try(cmstTests(Cross, "y1", "y2", Q.chr=cq[1,2], Q.pos=cq[1,3], 
          , , , , cross.type),silent=TRUE)
        if(class(aux) != "try-error"){
          R2s[k,] <- aux$R2
          BICs[k,] <- aux$BIC.stats[1:4]
          AICs[k,] <- aux$AIC.stats[1:4]
          z.scores[k,] <- aux$BIC.stats[5:10]
          sig2s[k,] <- aux$Sig.stats[1:6]
          pval.par.cmst.joint.BIC[k,] <- aux$pvals.par.cmst.joint.BIC
          pval.par.cmst.iu.BIC[k,] <- aux$pvals.par.cmst.iu.BIC
          pval.non.par.cmst.iu.BIC[k,] <- aux$pvals.non.par.cmst.iu.BIC
          pval.par.cmst.joint.AIC[k,] <- aux$pvals.par.cmst.joint.AIC
          pval.par.cmst.iu.AIC[k,] <- aux$pvals.par.cmst.iu.AIC
          pval.non.par.cmst.iu.AIC[k,] <- aux$pvals.non.par.cmst.iu.AIC
        }
        cit.mar <- find.marker(Cross,cq[1,2],cq[1,3])
        LL <- pull.geno(Cross)[,cit.mar]
        GG <- Cross$pheno[,1]
        TT <- Cross$pheno[,2]
        aux2 <- try(citTests(LL, GG, TT),silent=TRUE)
        if(class(aux2) != "try-error")
          pval.cit[k,] <- aux2
        k <- k + 1
      }
    }
  }
  if(model=="D"){
    k <- 1
    while(k <= nSim){
      mu <- runif(1,mu.range[1],mu.range[2])
      add.eff1 <- runif(1,add.eff1.range[1],add.eff1.range[2])
      dom.eff1 <- runif(1,dom.eff1.range[1],dom.eff1.range[2])
      add.eff2 <- runif(1,add.eff2.range[1],add.eff2.range[2])
      dom.eff2 <- runif(1,dom.eff2.range[1],dom.eff2.range[2])
      sig2.1 <- runif(1,sig2.1.range[1],sig2.1.range[2])
      sig2.2 <- runif(1,sig2.2.range[1],sig2.2.range[2])
      Cross <- sim.cross.4(n.ind, mu, add.eff1, dom.eff1, add.eff2, 
        dom.eff2, sig2.1, sig2.2, eq.spacing, cross.type, normalize)
      Cross <- calc.genoprob(Cross, step=2)
      cq <- get.common.qtls(Cross, "y1", "y2", thr, peak.dist)
      if(!is.na(cq[1])){
        print(k)
        cor12[k,] <- cor(Cross$pheno[,1],Cross$pheno[,2])
        aux <- try(cmstTests(Cross, "y1", "y2", Q.chr=cq[1,2], Q.pos=cq[1,3], 
          , , , , cross.type),silent=TRUE)
        if(class(aux) != "try-error"){
          R2s[k,] <- aux$R2
          BICs[k,] <- aux$BIC.stats[1:4]
          AICs[k,] <- aux$AIC.stats[1:4]
          z.scores[k,] <- aux$BIC.stats[5:10]
          sig2s[k,] <- aux$Sig.stats[1:6]
          pval.par.cmst.joint.BIC[k,] <- aux$pvals.par.cmst.joint.BIC
          pval.par.cmst.iu.BIC[k,] <- aux$pvals.par.cmst.iu.BIC
          pval.non.par.cmst.iu.BIC[k,] <- aux$pvals.non.par.cmst.iu.BIC
          pval.par.cmst.joint.AIC[k,] <- aux$pvals.par.cmst.joint.AIC
          pval.par.cmst.iu.AIC[k,] <- aux$pvals.par.cmst.iu.AIC
          pval.non.par.cmst.iu.AIC[k,] <- aux$pvals.non.par.cmst.iu.AIC
        }
        cit.mar <- find.marker(Cross,cq[1,2],cq[1,3])
        LL <- pull.geno(Cross)[,cit.mar]
        GG <- Cross$pheno[,1]
        TT <- Cross$pheno[,2]
        aux2 <- try(citTests(LL, GG, TT),silent=TRUE)
        if(class(aux2) != "try-error")
          pval.cit[k,] <- aux2
        k <- k + 1
      }
    }
  }
  if(model=="E"){
    k <- 1
    while(k <= nSim){
      mu <- runif(1,mu.range[1],mu.range[2])
      add.eff1 <- runif(1,add.eff1.range[1],add.eff1.range[2])
      dom.eff1 <- runif(1,dom.eff1.range[1],dom.eff1.range[2])
      add.eff2 <- runif(1,add.eff2.range[1],add.eff2.range[2])
      dom.eff2 <- runif(1,dom.eff2.range[1],dom.eff2.range[2])
      beta1h <- runif(1,beta1h.range[1],beta1h.range[2])
      beta2h <- runif(1,beta2h.range[1],beta2h.range[2])
      sig2.1 <- runif(1,sig2.1.range[1],sig2.1.range[2])
      sig2.2 <- runif(1,sig2.2.range[1],sig2.2.range[2])
      sig2.h <- runif(1,sig2.h.range[1],sig2.h.range[2])
      Cross <- sim.cross.5(n.ind, mu, add.eff1, dom.eff1, add.eff2, dom.eff2, 
        beta1h, beta2h, sig2.1, sig2.2, sig2.h, eq.spacing, cross.type,
        normalize)
      Cross <- calc.genoprob(Cross, step=2)
      cq <- get.common.qtls(Cross, "y1", "y2", thr, peak.dist)
      if(!is.na(cq[1])){
        print(k)
        cor12[k,] <- cor(Cross$pheno[,1],Cross$pheno[,2])
        aux <- try(cmstTests(Cross, "y1", "y2", Q.chr=cq[1,2], Q.pos=cq[1,3], 
          , , , , cross.type),silent=TRUE)
        if(class(aux) != "try-error"){
          R2s[k,] <- aux$R2
          BICs[k,] <- aux$BIC.stats[1:4]
          AICs[k,] <- aux$AIC.stats[1:4]
          z.scores[k,] <- aux$BIC.stats[5:10]
          sig2s[k,] <- aux$Sig.stats[1:6]
          pval.par.cmst.joint.BIC[k,] <- aux$pvals.par.cmst.joint.BIC
          pval.par.cmst.iu.BIC[k,] <- aux$pvals.par.cmst.iu.BIC
          pval.non.par.cmst.iu.BIC[k,] <- aux$pvals.non.par.cmst.iu.BIC
          pval.par.cmst.joint.AIC[k,] <- aux$pvals.par.cmst.joint.AIC
          pval.par.cmst.iu.AIC[k,] <- aux$pvals.par.cmst.iu.AIC
          pval.non.par.cmst.iu.AIC[k,] <- aux$pvals.non.par.cmst.iu.AIC
        }
        cit.mar <- find.marker(Cross,cq[1,2],cq[1,3])
        LL <- pull.geno(Cross)[,cit.mar]
        GG <- Cross$pheno[,1]
        TT <- Cross$pheno[,2]
        aux2 <- try(citTests(LL, GG, TT),silent=TRUE)
        if(class(aux2) != "try-error")
          pval.cit[k,] <- aux2
        k <- k + 1
      }
    }
  }
  list(cor12=cor12,
       R2s=R2s,
       BICs=BICs,
       AICs=AICs,
       z.scores=z.scores,
       sig2s=sig2s,
       pval.par.cmst.joint.BIC=pval.par.cmst.joint.BIC,
       pval.par.cmst.iu.BIC=pval.par.cmst.iu.BIC,
       pval.non.par.cmst.iu.BIC=pval.non.par.cmst.iu.BIC,
       pval.par.cmst.joint.AIC=pval.par.cmst.joint.AIC,
       pval.par.cmst.iu.AIC=pval.par.cmst.iu.AIC,
       pval.non.par.cmst.iu.AIC=pval.non.par.cmst.iu.AIC,
       pval.cit=pval.cit)
}
#########################################################################
SimCrossCausal <- function(n.ind, len, n.mar, beta, add.eff, dom.eff, 
                           sig2.1 = 1, sig2.2 = 1, eq.spacing = FALSE, 
                           cross.type = c("bc", "f2"), normalize = FALSE) {
  n.traits <- length(beta)
  beta <- matrix(rep(beta, each = n.ind), n.ind, n.traits)
  Map <- sim.map(len, n.mar, eq.spacing = eq.spacing, include.x = FALSE)
  Cross <- sim.cross(map = Map, n.ind = n.ind, type = cross.type)
  mygeno <- pull.geno(Cross)
  q <- mygeno[, "D1M51"]
  
  cross.type <- match.arg(cross.type)
  if (cross.type == "bc") {
    add.q <- q - 1.5
    y1 <- add.q * add.eff + rnorm(n.ind, 0, sqrt(sig2.1))
  }
  if (cross.type == "f2") {
    add.q <- q - 2
    dom.q <- (1 + add.q) * (1 - add.q) - 0.5
    y1 <- add.q * add.eff + dom.q * dom.eff + rnorm(n.ind, 0, sqrt(sig2.1))
  }
  y <- beta * y1 + matrix(rnorm(n.ind * n.traits, 0, sqrt(sig2.2)), n.ind, n.traits)
  y <- data.frame(y1, y)
  names(y) <- paste("y", 1 : (n.traits + 1), sep = "")
  if (normalize) {
    apply(y, 2, normal.trans)
  }
  Cross$pheno <- y
  Cross
}
################################################################################
SimCross1 <- function(n.ind, mu, beta21, add.eff1, dom.eff1, 
                      sig2.1 = 1, sig2.2 = 1, eq.spacing = FALSE, 
                      cross.type = "f2", normalize = FALSE) {
  Map <- sim.map(len = rep(100,3), n.mar = 101, eq.spacing = eq.spacing, 
                 include.x = FALSE)
  Cross <- sim.cross(map = Map, n.ind = n.ind, type = cross.type)
  mygeno <- pull.geno(Cross)
  q <- mygeno[, "D1M51"]
  if (cross.type == "bc") {
    add.q <- q - 1.5
    y1 <- mu + add.q * add.eff1 + rnorm(n.ind, 0, sqrt(sig2.1))
  }
  if (cross.type == "f2") {
    add.q <- q - 2
    dom.q <- (1 + add.q) * (1 - add.q) - 0.5
    y1 <- mu + add.q * add.eff1 + dom.q * dom.eff1 + 
          rnorm(n.ind, 0, sqrt(sig2.1))
  }
  y2 <- mu + beta21 * y1 + rnorm(n.ind, 0, sqrt(sig2.2))
  if (normalize) {
    y1 <- normal.trans(y1)
    y2 <- normal.trans(y2)
  }
  phenos <- data.frame(y1, y2)
  Cross$pheno <- phenos
  Cross
}
################################################################################
SimCross2 <- function(n.ind, mu, beta21, beta1h, beta2h, add.eff1, dom.eff1, 
                      sig2.1 = 1, sig2.2 = 1, sig2.h = 1, eq.spacing = FALSE, 
                      cross.type = "f2", normalize = FALSE) {
  Map <- sim.map(len = rep(100,3), n.mar = 101, eq.spacing = eq.spacing, 
                 include.x = FALSE)
  Cross <- sim.cross(map = Map, n.ind = n.ind, type = cross.type)
  mygeno <- pull.geno(Cross)
  q <- mygeno[, "D1M80"]
  h <- mu + rnorm(n.ind, 0, sqrt(sig2.h))
  if (cross.type == "bc") {
    add.q <- q - 1.5
    y1 <- mu + add.q * add.eff1 + beta1h * h + rnorm(n.ind, 0, sqrt(sig2.1))
  }
  if (cross.type == "f2") {
    add.q <- q - 2
    dom.q <- (1 + add.q) * (1 - add.q) - 0.5
    y1 <- mu + add.q * add.eff1 + dom.q * dom.eff1 + beta1h * h + 
          rnorm(n.ind, 0, sqrt(sig2.1))
  }
  y2 <- mu + beta21 * y1 + beta2h * h + rnorm(n.ind, 0, sqrt(sig2.2))
  if (normalize) {
    y1 <- normal.trans(y1)
    y2 <- normal.trans(y2)
  }
  phenos <- data.frame(y1, y2)
  Cross$pheno <- phenos
  Cross
}
##############################################################################
SimCross3 <- function(n.ind, mu, beta21, add.eff1, dom.eff1, add.eff2,
                      dom.eff2, sig2.1 = 1, sig2.2 = 1, eq.spacing = FALSE, 
                      cross.type = "f2", normalize = FALSE) {
  Map <- sim.map(len = rep(100, 3), n.mar = 101, eq.spacing = eq.spacing, 
                 include.x = FALSE)
  Cross <- sim.cross(map = Map, n.ind = n.ind, type = cross.type)
  mygeno <- pull.geno(Cross)
  q <- mygeno[, "D1M80"]
  if (cross.type == "bc") {
    add.q <- q - 1.5
    y1 <- mu + add.q * add.eff1 + rnorm(n.ind, 0, sqrt(sig2.1))
    y2 <- mu + add.q * add.eff2 + beta21 * y1 + rnorm(n.ind, 0, sqrt(sig2.2))
  }
  if (cross.type == "f2") {
    add.q <- q - 2
    dom.q <- (1 + add.q) * (1 - add.q) - 0.5
    y1 <- mu + add.q * add.eff1 + dom.q * dom.eff1 + 
          rnorm(n.ind, 0, sqrt(sig2.1))
    y2 <- mu + add.q * add.eff2 + dom.q * dom.eff2 + beta21 * y1 + 
          rnorm(n.ind, 0, sqrt(sig2.2))
  }
  if (normalize) {
    y1 <- normal.trans(y1)
    y2 <- normal.trans(y2)
  }
  phenos <- data.frame(y1, y2)
  Cross$pheno <- phenos
  Cross
}
##############################################################################
SimCross4 <- function(n.ind, mu, add.eff1, dom.eff1, add.eff2, dom.eff2, 
                      sig2.1 = 1, sig2.2 = 1, eq.spacing = FALSE, 
                      cross.type = "f2", normalize = FALSE) {
  Map <- sim.map(len = rep(100, 3), n.mar = 101, eq.spacing = eq.spacing, 
                 include.x = FALSE)
  Cross <- sim.cross(map = Map, n.ind = n.ind, type = cross.type)
  mygeno <- pull.geno(Cross)
  q <- mygeno[, "D1M80"]
  if (cross.type == "bc") {
    add.q <- q - 1.5
    y1 <- mu + add.q * add.eff1 + rnorm(n.ind, 0, sqrt(sig2.1))
    y2 <- mu + add.q * add.eff2 + rnorm(n.ind, 0, sqrt(sig2.2))
  }
  if (cross.type == "f2") {
    add.q <- q - 2
    dom.q <- (1 + add.q) * (1 - add.q) - 0.5
    y1 <- mu + add.q * add.eff1 + dom.q * dom.eff1 + 
          rnorm(n.ind, 0, sqrt(sig2.1))
    y2 <- mu + add.q * add.eff2 + dom.q * dom.eff2 + 
          rnorm(n.ind, 0, sqrt(sig2.2))
  }
  if (normalize) {
    y1 <- normal.trans(y1)
    y2 <- normal.trans(y2)
  }
  phenos <- data.frame(y1, y2)
  Cross$pheno <- phenos
  Cross
}
##############################################################################
SimCross5 <- function(n.ind, mu, add.eff1, dom.eff1, add.eff2, dom.eff2, 
                      beta1h, beta2h, sig2.1 = 1, sig2.2 = 1, sig2.h = 1, 
                      eq.spacing = FALSE, cross.type = "f2", 
                      normalize = FALSE) {
  Map <- sim.map(len = rep(100, 3), n.mar = 101, eq.spacing = eq.spacing, 
                 include.x = FALSE)
  Cross <- sim.cross(map = Map, n.ind = n.ind, type = cross.type)
  mygeno <- pull.geno(Cross)
  q <- mygeno[, "D1M80"]
  h <- mu + rnorm(n.ind, 0, sqrt(sig2.h))
  if (cross.type == "bc") {
    add.q <- q - 1.5
    y1 <- mu + add.q * add.eff1 + h * beta1h + rnorm(n.ind, 0, sqrt(sig2.1))
    y2 <- mu + add.q * add.eff2 + h * beta2h + rnorm(n.ind, 0, sqrt(sig2.2))
  }
  if (cross.type == "f2") {
    add.q <- q - 2
    dom.q <- (1 + add.q) * (1 - add.q) - 0.5
    y1 <- mu + add.q * add.eff1 + dom.q * dom.eff1 + h * beta1h + 
          rnorm(n.ind, 0, sqrt(sig2.1))
    y2 <- mu + add.q * add.eff2 + dom.q * dom.eff2 + h * beta2h + 
          rnorm(n.ind, 0, sqrt(sig2.2))
  }
  if (normalize) {
    y1 <- normal.trans(y1)
    y2 <- normal.trans(y2)
  }
  phenos <- data.frame(y1, y2)
  Cross$pheno <- phenos
  Cross
}
##############################################################################
SimCross6 <- function(n.ind, mu, add.eff, dom.eff, beta1h, beta2h, 
                      sig2.1 = 1, sig2.2 = 1, sig2.h = 1, eq.spacing = FALSE, 
                      cross.type = "f2", normalize = FALSE) {
  Map <- sim.map(len = rep(100, 3), n.mar = 101, eq.spacing = eq.spacing, 
                 include.x = FALSE)
  Cross <- sim.cross(map = Map, n.ind = n.ind, type = cross.type)
  mygeno <- pull.geno(Cross)
  q <- mygeno[, "D1M80"]
  if (cross.type == "bc") {
    add.q <- q - 1.5
    h <- mu + add.q * add.eff + rnorm(n.ind, 0, sqrt(sig2.h))
  }
  if (cross.type == "f2") {
    add.q <- q - 2
    dom.q <- (1 + add.q) * (1 - add.q) - 0.5
    h <- mu + add.q * add.eff + dom.q * dom.eff + 
         rnorm(n.ind, 0, sqrt(sig2.h))
  }
  y1 <- mu + h * beta1h + rnorm(n.ind, 0, sqrt(sig2.1))
  y2 <- mu + h * beta2h + rnorm(n.ind, 0, sqrt(sig2.2))
  if (normalize) {
    y1 <- normal.trans(y1)
    y2 <- normal.trans(y2)
  }
  phenos <- data.frame(y1, y2)
  Cross$pheno <- phenos
  Cross
}

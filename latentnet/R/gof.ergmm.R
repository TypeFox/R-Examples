gof.ergmm <- function (object, ..., nsim=100,
                      GOF=~idegree+odegree+distance, 
		      verbose=FALSE) {

  formula <- object[["model"]][["formula"]]

  trms <- ergm.getterms(formula)
  if(length(trms)>2){
    nw <- eval(trms[[2]], sys.parent())
  }else{
    stop("A network object on the RHS of the formula argument must be given")
  }

  nsim <- max(nsim, dim(object[["Z"]])[3])

  all.gof.vars <- all.vars(GOF)

# match variables

  for(i in seq(along=all.gof.vars)){
   all.gof.vars[i] <- match.arg(all.gof.vars[i],
    c('distance', 'espartners', 'dspartners', 'odegree', 'idegree', 
      'degree','triadcensus','model'
     )
                               )
  }
  GOF <- as.formula(paste("~",paste(all.gof.vars,collapse="+")))

  if(!is.network(nw)){
    stop("A network object on the RHS of the formula argument must be given")
  }

  pval.model<-pval.triadcensus<-pval.dist<-pval.deg<-pval.espart<-pval.espart<-NULL
#
  obs.model<-pobs.model<-sim.model<-psim.model<-pval.model<-bds.model<-NULL
  obs.triadcensus<-pobs.triadcensus<-sim.triadcensus<-psim.triadcensus<-pval.triadcensus<-bds.triadcensus<-NULL
  obs.dist<-pobs.dist<-sim.dist<-psim.dist<-pval.dist<-bds.dist<-NULL
  obs.deg<-pobs.deg<-sim.deg<-psim.deg<-pval.deg<-bds.deg<-NULL
  obs.espart<-pobs.espart<-sim.espart<-psim.espart<-pval.espart<-bds.espart<-NULL
  obs.dspart<-pobs.dspart<-sim.dspart<-psim.dspart<-pval.dspart<-bds.dspart<-NULL

  obs.ideg<-pobs.ideg<-sim.ideg<-psim.ideg<-pval.ideg<-bds.ideg<-pval.ideg<-NULL
  obs.odeg<-pobs.odeg<-sim.odeg<-psim.odeg<-pval.odeg<-bds.odeg<-pval.odeg<-NULL

  n <- network.size(nw)

  # Calculate network statistics for the observed graph
  # Set up the output arrays of sim variables

  if ('model' %in% all.gof.vars) {
   obs.model <- summary(formula)
   sim.model <- array(0,dim=c(nsim,length(obs.model)))
   dimnames(sim.model) <- list(paste(c(1:nsim)),names(obs.model))
  }

  if ('distance' %in% all.gof.vars) {
   obs.dist <- ergm.geodistdist(nw)
   obs.dist[obs.dist==Inf] <- n
   sim.dist <-array(0,dim=c(nsim,n))
   dimnames(sim.dist)  <- list(paste(c(1:nsim)),paste(1:n))
  }

  if ('odegree' %in% all.gof.vars) {
    mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
    obs.odeg <- summary(as.formula(paste('nw ~ odegree(',mesp,')',sep="")),drop=FALSE)
   sim.odeg <- array(0,dim=c(nsim,n))
   dimnames(sim.odeg)   <- list(paste(c(1:nsim)),paste(0:(n-1)))
   names(obs.odeg) <- dimnames(sim.odeg)[[2]]
  }

  if ('idegree' %in% all.gof.vars) {
    mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
    obs.ideg <- summary(as.formula(paste('nw ~ idegree(',mesp,')',sep="")),drop=FALSE)
   sim.ideg <- array(0,dim=c(nsim,n))
   dimnames(sim.ideg)   <- list(paste(c(1:nsim)),paste(0:(n-1)))
   names(obs.ideg) <- dimnames(sim.ideg)[[2]]
  }

  if ('degree' %in% all.gof.vars) {
     mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
     obs.deg <- summary(as.formula(paste('nw ~ degree(',mesp,')',sep="")),drop=FALSE)
   sim.deg <- array(0,dim=c(nsim,n))
   dimnames(sim.deg)   <- list(paste(c(1:nsim)),paste(0:(n-1)))
   names(obs.deg) <- dimnames(sim.deg)[[2]]
  }
 
  if ('espartners' %in% all.gof.vars) {
    mesp <- paste("c(",paste(0:(network.size(nw)-2),collapse=","),")",sep="")
    obs.espart <- summary(as.formula(paste('nw ~ esp(',mesp,')',sep="")), drop=FALSE)
   sim.espart <- array(0,dim=c(nsim,n-1))
   dimnames(sim.espart) <- list(paste(c(1:nsim)),paste(0:(n-2)))
  }
 
  if ('dspartners' %in% all.gof.vars) {
    mesp <- paste("c(",paste(0:(network.size(nw)-2),collapse=","),")",sep="")
    obs.dspart <- summary(as.formula(paste('nw ~ dsp(',mesp,')',sep="")), drop=FALSE)
   sim.dspart <- array(0,dim=c(nsim,n-1))
   dimnames(sim.dspart) <- list(paste(c(1:nsim)),paste(0:(n-2)))
  }

  if ('triadcensus' %in% all.gof.vars) {
   if(is.directed(nw)){
    triadcensus <- 0:15
    namestriadcensus <- c("003","012", "102", "021D", "021U", "021C",
      "111D", "111U", "030T",
      "030C", "201", "120D", "120U", "120C", "210", "300")
    triadcensus.formula <- "~ triadcensus(0:15)"
   }else{
    triadcensus <- 0:3
    namestriadcensus <- c("0","1","2", "3")
    triadcensus.formula <- "~ triadcensus(0:3)"
   }
   obs.triadcensus <- summary(as.formula(paste('nw',triadcensus.formula,sep="")), drop=FALSE)
   sim.triadcensus <- array(0,dim=c(nsim,length(triadcensus)))
   dimnames(sim.triadcensus) <- list(paste(c(1:nsim)), namestriadcensus)
   names(obs.triadcensus) <- namestriadcensus
  }
 
  # Simulate an exponential family random graph model

  SimNetworkSeriesObj <- simulate(object,nsim=nsim)

  if(verbose){cat("\nCollating simulations\n")}

  for (i in 1:nsim)
  { 
    if(verbose){
     cat("\nCalculating statistics for simulation",i,"\n")
    }

    if ('model' %in% all.gof.vars) {
     sim.model[i,] <- summary(update(formula,SimNetworkSeriesObj[["networks"]][[i]] ~ .))
    }

    if ('distance' %in% all.gof.vars) {
     sim.dist[i,] <- ergm.geodistdist(SimNetworkSeriesObj[["networks"]][[i]])
    }
    if ('idegree' %in% all.gof.vars) {
     mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
     gi <- SimNetworkSeriesObj[["networks"]][[i]]
     sim.ideg[i,] <- summary(as.formula(paste('gi ~ idegree(',mesp,')',sep="")),drop=FALSE)
#    temp <- table(degreedist(SimNetworkSeriesObj[["networks"]][[i]], print=verbose)[1,])
#    sim.ideg[i,] <- c(temp, rep(0, n-length(temp)))
    }
    if ('odegree' %in% all.gof.vars) {
     mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
     gi <- SimNetworkSeriesObj[["networks"]][[i]]
     sim.odeg[i,] <- summary(as.formula(paste('gi ~ odegree(',mesp,')',sep="")),drop=FALSE)
    }
    if ('degree' %in% all.gof.vars) {
     gi <- SimNetworkSeriesObj[["networks"]][[i]]
     if(is.bipartite(gi)){
      temp <- degreedist(gi, print=FALSE)[["event"]]
      sim.deg[i,] <- c(temp,rep(0,n-length(temp)))
     }else{
      mesp <- paste("c(",paste(0:(n-1),collapse=","),")",sep="")
      sim.deg[i,] <- summary(as.formula(paste('gi ~ degree(',mesp,')',sep="")),drop=FALSE)
     }
    }
    if ('espartners' %in% all.gof.vars) {
     gi <- SimNetworkSeriesObj[["networks"]][[i]]
     mesp <- paste("c(",paste(0:(network.size(gi)-2),collapse=","),")",sep="")
     sim.espart[i,] <- summary(as.formula(paste('gi ~ esp(',mesp,')',sep="")), drop=FALSE)
    }
    if ('dspartners' %in% all.gof.vars) {
     gi <- SimNetworkSeriesObj[["networks"]][[i]]
     mesp <- paste("c(",paste(0:(network.size(gi)-2),collapse=","),")",sep="")
     sim.dspart[i,] <- summary(as.formula(paste('gi ~ dsp(',mesp,')',sep="")), drop=FALSE)
    }
    if ('triadcensus' %in% all.gof.vars) {
     gi <- SimNetworkSeriesObj[["networks"]][[i]]
     sim.triadcensus[i,] <- summary(as.formula(paste('gi',triadcensus.formula,sep="")), drop=FALSE)
    }
  }

  # calculate p-values

 if ('model' %in% all.gof.vars) {
  pval.model <- apply(sim.model <= obs.model[col(sim.model)],2,mean)
  pval.model.top <- apply(sim.model >= obs.model[col(sim.model)],2,mean)
  pval.model <- cbind(obs.model,apply(sim.model, 2,min), apply(sim.model, 2,mean),
                apply(sim.model, 2,max), pmin(1,2*pmin(pval.model,pval.model.top)))
  dimnames(pval.model)[[2]] <- c("obs","min","mean","max","MC p-value")
  pobs.model <- pval.model.top
  psim.model <- apply(sim.model,2,rank)/nrow(sim.model)
  bds.model <- apply(psim.model,2,quantile,probs=c(0.025,0.975))
 }

 if ('distance' %in% all.gof.vars) {
  pval.dist <- apply(sim.dist <= obs.dist[col(sim.dist)],2,mean)
  pval.dist.top <- apply(sim.dist >= obs.dist[col(sim.dist)],2,mean)
  pval.dist <- cbind(obs.dist,apply(sim.dist, 2,min), apply(sim.dist, 2,mean),
                apply(sim.dist, 2,max), pmin(1,2*pmin(pval.dist,pval.dist.top)))
  dimnames(pval.dist)[[2]] <- c("obs","min","mean","max","MC p-value")
  pobs.dist <- obs.dist/sum(obs.dist)
  psim.dist <- sweep(sim.dist,1,apply(sim.dist,1,sum),"/")
  psim.dist[is.na(psim.dist)] <- 1
  bds.dist <- apply(psim.dist,2,quantile,probs=c(0.025,0.975))
 }

# cat("\nGoodness-of-fit for minimum geodesic distance\n\n")
# print(pval.dist)

 if ('idegree' %in% all.gof.vars) {
  pval.ideg <- apply(sim.ideg <= obs.ideg[col(sim.ideg)],2,mean)
  pval.ideg.top <- apply(sim.ideg >= obs.ideg[col(sim.ideg)],2,mean)
  pval.ideg <- cbind(obs.ideg,apply(sim.ideg, 2,min), apply(sim.ideg, 2,mean),
                apply(sim.ideg, 2,max), pmin(1,2*pmin(pval.ideg,pval.ideg.top)))
  dimnames(pval.ideg)[[2]] <- c("obs","min","mean","max","MC p-value")
  pobs.ideg <- obs.ideg/sum(obs.ideg)
  psim.ideg <- sweep(sim.ideg,1,apply(sim.ideg,1,sum),"/")
  psim.ideg[is.na(psim.ideg)] <- 1
  bds.ideg <- apply(psim.ideg,2,quantile,probs=c(0.025,0.975))
 }

 if ('odegree' %in% all.gof.vars) {
  pval.odeg <- apply(sim.odeg <= obs.odeg[col(sim.odeg)],2,mean)
  pval.odeg.top <- apply(sim.odeg >= obs.odeg[col(sim.odeg)],2,mean)
  pval.odeg <- cbind(obs.odeg,apply(sim.odeg, 2,min), apply(sim.odeg, 2,mean),
                apply(sim.odeg, 2,max), pmin(1,2*pmin(pval.odeg,pval.odeg.top)))
  dimnames(pval.odeg)[[2]] <- c("obs","min","mean","max","MC p-value")
  pobs.odeg <- obs.odeg/sum(obs.odeg)
  psim.odeg <- sweep(sim.odeg,1,apply(sim.odeg,1,sum),"/")
  psim.odeg[is.na(psim.odeg)] <- 1
  bds.odeg <- apply(psim.odeg,2,quantile,probs=c(0.025,0.975))
 }

 if ('degree' %in% all.gof.vars) {
  pval.deg <- apply(sim.deg <= obs.deg[col(sim.deg)],2,mean)
  pval.deg.top <- apply(sim.deg >= obs.deg[col(sim.deg)],2,mean)
  pval.deg <- cbind(obs.deg,apply(sim.deg, 2,min), apply(sim.deg, 2,mean),
                apply(sim.deg, 2,max), pmin(1,2*pmin(pval.deg,pval.deg.top)))
  dimnames(pval.deg)[[2]] <- c("obs","min","mean","max","MC p-value")
  pobs.deg <- obs.deg/sum(obs.deg)
  psim.deg <- sweep(sim.deg,1,apply(sim.deg,1,sum),"/")
  psim.deg[is.na(psim.deg)] <- 1
  bds.deg <- apply(psim.deg,2,quantile,probs=c(0.025,0.975))
 }

# cat("\nGoodness-of-fit for degree\n\n")
# print(pval.deg)

 if ('espartners' %in% all.gof.vars) {
  pval.espart <- apply(sim.espart <= obs.espart[col(sim.espart)],2,mean)
  pval.espart.top <- apply(sim.espart >= obs.espart[col(sim.espart)],2,mean)
  pval.espart <- cbind(obs.espart,apply(sim.espart, 2,min), apply(sim.espart, 2,mean),
                apply(sim.espart, 2,max), pmin(1,2*pmin(pval.espart,pval.espart.top)))
  dimnames(pval.espart)[[2]] <- c("obs","min","mean","max","MC p-value")
  pobs.espart <- obs.espart/sum(obs.espart)
  psim.espart <- sweep(sim.espart,1,apply(sim.espart,1,sum),"/")
  psim.espart[is.na(psim.espart)] <- 1
  bds.espart <- apply(psim.espart,2,quantile,probs=c(0.025,0.975))
 }

# cat("\nGoodness-of-fit for edgewise shared partner\n\n")
# print(pval.espart)

 if ('dspartners' %in% all.gof.vars) {
  pval.dspart <- apply(sim.dspart <= obs.dspart[col(sim.dspart)],2,mean)
  pval.dspart.top <- apply(sim.dspart >= obs.dspart[col(sim.dspart)],2,mean)
  pval.dspart <- cbind(obs.dspart,apply(sim.dspart, 2,min), apply(sim.dspart, 2,mean),
                apply(sim.dspart, 2,max), pmin(1,2*pmin(pval.dspart,pval.dspart.top)))
  dimnames(pval.dspart)[[2]] <- c("obs","min","mean","max","MC p-value")
  pobs.dspart <- obs.dspart/sum(obs.dspart)
  psim.dspart <- sweep(sim.dspart,1,apply(sim.dspart,1,sum),"/")
  psim.dspart[is.na(psim.dspart)] <- 1
  bds.dspart <- apply(psim.dspart,2,quantile,probs=c(0.025,0.975))
 }

# cat("\nGoodness-of-fit for dyadwise shared partner\n\n")
# print(pval.dspart)

 if ('triadcensus' %in% all.gof.vars) {
  pval.triadcensus <- apply(sim.triadcensus <= obs.triadcensus[col(sim.triadcensus)],2,mean)
  pval.triadcensus.top <- apply(sim.triadcensus >= obs.triadcensus[col(sim.triadcensus)],2,mean)
  pval.triadcensus <- cbind(obs.triadcensus,apply(sim.triadcensus, 2,min), apply(sim.triadcensus, 2,mean),
                apply(sim.triadcensus, 2,max), pmin(1,2*pmin(pval.triadcensus,pval.triadcensus.top)))
  dimnames(pval.triadcensus)[[2]] <- c("obs","min","mean","max","MC p-value")
  pobs.triadcensus <- obs.triadcensus/sum(obs.triadcensus)
  psim.triadcensus <- sweep(sim.triadcensus,1,apply(sim.triadcensus,1,sum),"/")
  psim.triadcensus[is.na(psim.triadcensus)] <- 1
  bds.triadcensus <- apply(psim.triadcensus,2,quantile,probs=c(0.025,0.975))
 }

# cat("\nGoodness-of-fit for edgewise shared partner\n\n")
# print(pval.espart)
# Return

  returnlist <- list(n,
   pval.model, pval.triadcensus, pval.dist, pval.ideg, pval.odeg, pval.deg, pval.espart, pval.dspart,
   obs.model, pobs.model, sim.model, psim.model, pval.model, bds.model,
   obs.triadcensus, pobs.triadcensus, sim.triadcensus, psim.triadcensus, pval.triadcensus, bds.triadcensus,
   obs.dist, pobs.dist, sim.dist, psim.dist, pval.dist, bds.dist,
   obs.ideg, pobs.ideg, sim.ideg, psim.ideg, pval.ideg, bds.ideg,
   obs.odeg, pobs.odeg, sim.odeg, psim.odeg, pval.odeg, bds.odeg,
   obs.deg, pobs.deg, sim.deg, psim.deg, pval.deg, bds.deg,
   obs.espart, pobs.espart, sim.espart, psim.espart, pval.espart, bds.espart,
   obs.dspart, pobs.dspart, sim.dspart, psim.dspart, pval.dspart, bds.dspart,
   GOF
                   )

  names(returnlist) <- c(
  "network.size",
  "summary.model",
  "summary.triadcensus",
  "summary.dist",
  "summary.ideg",
  "summary.odeg",
  "summary.deg",
  "summary.espart",
  "summary.dspart",
  "obs.model", "pobs.model", "sim.model", "psim.model", "pval.model", "bds.model",
  "obs.triadcensus", "pobs.triadcensus", "sim.triadcensus", "psim.triadcensus", "pval.triadcensus", "bds.triadcensus",
  "obs.dist", "pobs.dist", "sim.dist", "psim.dist", "pval.dist", "bds.dist",
  "obs.ideg", "pobs.ideg", "sim.ideg", "psim.ideg", "pval.ideg", "bds.ideg",
  "obs.odeg", "pobs.odeg", "sim.odeg", "psim.odeg", "pval.odeg", "bds.odeg",
  "obs.deg", "pobs.deg", "sim.deg", "psim.deg", "pval.deg", "bds.deg",
  "obs.espart", "pobs.espart", "sim.espart", "psim.espart", "pval.espart", "bds.espart",
  "obs.dspart", "pobs.dspart", "sim.dspart", "psim.dspart", "pval.dspart", "bds.dspart",
  "GOF"
                        )
  class(returnlist) <- "gofobject"
  returnlist
  }

equilibrium <-
function(start, model, data, tolerance=1e-5, max.iter=100, coal=0, alpha=0, 
margin=NULL, fixed=NULL, gamma=0, boot=0, MC=0,self.var="self", prox.var="prox", position=NULL, votes=NULL,quadratic=TRUE){
    ccc <- match.call()
    continue <- TRUE
    iter <- 0
    tmp <- NULL
    k <- length(unique(data$alt))
    if(missing(start))
     start <- sample(data[,self.var], k)
    tmp <- start

    basic <- equilibrium.internal(start=start, model=model, data=data, tolerance=tolerance, max.iter=max.iter, coal=coal, alpha=alpha,margin=margin, fixed=fixed, gamma=gamma, 
        self.var=self.var, prox.var=prox.var,quadratic=quadratic)
    bootstrap <- NULL
    montecarlo <- NULL
    
    est <- NULL
    mP <- NULL
    chid <- unique(data$chid)
    nvoters <- length(chid)


    if(boot>0){
        cat("\nBootstrapping \n\n")
        n.sim <- boot
        est <- matrix(,n.sim, k)
        mP <- matrix(,n.sim, k)
        pb <- txtProgressBar(min=1, max=n.sim, style=3)
        for(sim in 1:n.sim){
            
         newchid <- sample(chid, nvoters, replace=TRUE)
         iidd <- match(data$chid, newchid)
         iidd <- which(!is.na(iidd))
         new.election <- data[iidd,]
         ans <- equilibrium.internal(start=start, model=model, data=new.election, tolerance=tolerance, max.iter=max.iter, coal=coal, alpha=alpha, margin=margin, fixed=fixed, gamma=gamma,self.var=self.var, prox.var=prox.var,quadratic=quadratic)
      #   cat("o")
         est[sim,] <- ans$est
         mP[sim,] <- ans$mP  
         setTxtProgressBar(pb, sim)
            #     print(ans)
        }
        close(pb)
        
        est.mean <- apply(est,2,mean)
        est.sd <- apply(est,2,sd)
        names(est.mean) <- names(ans$est)
        names(est.sd) <- names(ans$est)

        mP.mean <- apply(mP,2,mean)
        mP.sd <- apply(mP,2,sd)
        names(mP.mean) <- names(ans$mP)
        names(mP.sd) <- names(ans$mP)
   
        bootstrap <-  list(est.mean=est.mean, est.sd=est.sd, mP.mean=mP.mean, mP.sd=mP.sd,replications=boot)  
        
        cat("\n")
    }  
     
   

    if(MC>0){
        cat("\nDoing Monte Carlo \n\n")
        n.sim <- MC
        est <- matrix(,n.sim, k)
        mP <- matrix(,n.sim, k)
      
        IC <- summary(model)$CoefTable
        pb <- txtProgressBar(min=1, max=n.sim, style=3)
      
        for(sim in 1:n.sim){
        
         newpar <- apply(IC, 1, function(x) rnorm(1, x[1],x[2]))
         idx <- match(names(newpar), names(model$coefficients))
         mod1 <- model
         mod1$coefficients[idx] <- newpar
         ans <- equilibrium.internal(start=start, model=mod1, data=data, tolerance=tolerance, max.iter=max.iter, coal=coal, alpha=alpha, margin=margin, fixed=fixed, gamma=gamma,
             self.var=self.var, prox.var=prox.var,quadratic=quadratic)
#         cat("o")
         est[sim,] <- ans$est
         mP[sim,] <- ans$mP  
         setTxtProgressBar(pb, sim) 
        }
         close(pb)
    
    est.mean <- apply(est,2,mean)
    est.sd <- apply(est,2,sd)
    names(est.mean) <- names(ans$est)
    names(est.sd) <- names(ans$est)

    mP.mean <- apply(mP,2,mean)
    mP.sd <- apply(mP,2,sd)
    names(mP.mean) <- names(ans$mP)
    names(mP.sd) <- names(ans$mP)
     montecarlo <-  list(est.mean=est.mean, est.sd=est.sd, mP.mean=mP.mean, mP.sd=mP.sd,replications=MC)  
        cat("\n")
    }
    
    obj <- list(basic=basic, MC=montecarlo,boot=bootstrap, position=position, votes=votes, call=ccc)
    class(obj) <- "nash.eq"
    obj
}


# ottimizzare spostando margin, gamma, alpha ecc fuori da nash.hip e prepararli
# dentro equilibrium.internal

equilibrium.internal <- function(start, model, data, tolerance=1e-5, max.iter=100, coal=0, alpha=0, margin, fixed=NULL, gamma=0,
    self.var, prox.var,quadratic=TRUE){
    continue <- TRUE
    iter <- 0
    tmp <- start
    names(tmp) <- levels(attr(data,"index")$alt)
    while(continue){
        tmp1 <- nash.hip(tmp, model, data, coal=coal, alpha=alpha, margin=margin, fixed=fixed, gamma=gamma,
            self.var=self.var, prox.var=prox.var,quadratic=quadratic)
        iter <- iter + 1
        #    if(iter==2) continue <- FALSE
        if( (sum(abs(tmp-tmp1$est))<tolerance) | (iter>max.iter) )
        continue <- FALSE
        tmp <- tmp1$est  
     #   cat(".") # print(tmp)
    }
    #  cat("\n")
    tmp1
}


# calculate nash equilibrium
# tutti max voti
# start : posizione del partito
# election: data set agginrato per la distanza elettore-partito
nash.hip <- function(start, model, data, coal=0, alpha=0, margin=NULL, fixed=NULL, gamma=0, self.var, prox.var, quadratic=TRUE){
    if(is.null(coal))
     cat("\nno coalitions")
    
    
    # aggiornamento posizione dei partiti
    #cat("\n")
    #    print(start)
    true.order <- attr(attr(data, "reshapeLong")$varying,"times")
    nvoters <- length(unique(data$chid))
    aa <-  levels(attr(data,"index")$alt) 
    # print(aa)
    newpos <- rep(start[aa], nvoters)
    new.election <- data
    if(quadratic){
     new.election[prox.var] <- -(data[self.var]-newpos)^2
    } else {
     new.election[prox.var] <- -abs(data[self.var]-newpos)
    }
    #print(head(new.election))
    P <- predict(model, newdata=new.election)
    K <- NCOL(P)
  #  print(head(P))
  #  print(str(data))
  #  print(str(new.election))
    
    #print(true.order)
    
    P <- P[, true.order]
    mP <- colMeans(P)
    #  print(mP)
    #    print(colnames(P))
    #print(names(mP))
    
    
    npar <- colnames(P)
    
    
    CC <- rep(0, length(npar) )
    names(CC) <- npar
    if(is.list(coal)){
        C1 <- unlist(coal)
        idx <- match( names(C1), names(CC))
        CC[idx] <- C1
    } else {
        for(i in 1:length(CC))
        CC[i] <- coal[1]
    }
    
    AA <- rep(0, length(npar) ) # default value for alpha=0
    names(AA) <- npar
    
    if(is.list(alpha)){
        A1 <- unlist(alpha)
        idx <- match( names(A1), names(AA))
        AA[idx] <- A1    
    } else {
        for(i in 1:length(AA))
        AA[i] <- alpha[1]
    }
    
    AA[ which(CC==0) ] <- 0  # lonely parties
    
    MM <- vector( mode="list", length(npar))
    names(MM) <- npar
    
    if(is.list(margin)){
        idx <- match( names(margin), names(MM))
        MM[idx] <- margin    
    }
    
    
    # qua sotto formula per max. solo voti del partito
    v2 <- matrix(, nvoters, K )  # Pi*(1-Pi)
    v1 <- matrix(, nvoters, K )  # Pi*(1-Pi)*X
    
    #   for(i in 1:K){
    #    v2[,i] <-  P[,i]*(1-P[,i])
    #    v1[,i] <-  P[,i]*(1-P[,i])*(data$auto[K*(1:nvoters)-(K-1)])  ###
    #}
    
    # formula per coalizioni
    
    for(i in 1:K){
        mm <- 0
        cln <- CC[i] # number of coalition the party belongs to
        idx <- which(CC == cln) # index of parties in the same coalition
        cl <- CC[idx]
        
        if(!is.null(MM[[i]])){
            
            mm <- rowSums(matrix(P[,idx],nvoters,length(idx)))   
            idd <- match(MM[[i]], npar)
            tt <- rowSums(matrix(P[,idd],nvoters,length(idd)))   
            mm <- mm*tt
        }
        
        if(length(idx)>1){
            idx <- idx[-which(names(idx)==names(CC)[i])]
            idx <- as.numeric(idx)
            v2[,i] <- P[,i]*(1-P[,i] - AA[i]*rowSums(matrix(P[, idx],nvoters,length(idx)) ) + mm)
            v1[,i] <- (P[,i]*(1-P[,i] - AA[i]*rowSums(matrix(P[, idx],nvoters,length(idx))) + mm))*(data[K*(1:nvoters)-(K-1),self.var])
            
        } else {
            v2[,i] <-  P[,i]*(1-P[,i] + mm)
            v1[,i] <-  P[,i]*(1-P[,i] + mm)*(data[K*(1:nvoters)-(K-1),self.var])  ###
        }
        
    }
    
    
    
    v1m <- colMeans(v1)
    v2m <- colMeans(v2)
    
    est <- v1m/v2m # nuova posizione partiti
    names(est) <- names(mP)
    
    BB <- numeric(length(npar))
    names(BB) <- npar
    
    if(is.list(gamma)){
        idx <- match( names(gamma), names(BB))
        BB[idx] <- unlist(gamma)    
    } else {
        for(i in 1:length(BB))
        BB[i] <- gamma[1]
    }
    
    
    if(!is.null(fixed)){
        for(i in names(fixed)){
            est[i] <- fixed[[i]] *(1-BB[i]) + est[i]*BB[i]                        
        }
    }
    
    return(list(est=est, mP=mP))
}




EM <- function(y, X=NULL, ETA=NULL, R=NULL, init=NULL, iters=50, REML=TRUE, draw=TRUE, silent=FALSE, forced=NULL){
  nanana <- names(ETA)
  n.var.com <- length(ETA)
  #ETA <- lapply(ETA, function(x){if(length(x) == 2){y <- list(Z=as(x[[1]],Class="sparseMatrix"), K=as(x[[2]],Class="sparseMatrix"))}else{y <- list(Z=as(x[[1]],Class="sparseMatrix"))}; return(y)})
  ### make sure fixed effects are there
  if(is.null(X) & is.null(ETA)){ # nothing in the model
    tn = length(y); xm <- matrix(1,tn,1) 
    output <- lm(y~xm-1) # intercept model
  }else{ # user provided some structures
    ############################ FIX LIST LEVEL STRUCTURE
    # if was provided as a one level list
    if(is.list(ETA)){
      if(is.list(ETA[[1]])){
        ETA=ETA
      }else{
        ETA=list(ETA)
      }
    }else{
      stop;
      cat("\nThe random effects need to be provided in a list format, please see examples")
    }
    ########################### DONE...
    Zs <- lapply(ETA, function(x){x[[1]]})
    Gs <- lapply(ETA, function(x){x[[2]]})
    Zsp <- as(do.call("cbind", Zs),Class="sparseMatrix") # column bind Z=[Z1 Z2 Z3]
    Ksp <- as(do.call("adiag1", Gs),Class="sparseMatrix") # G as diagonal
    ############################
    if(is.null(X) & !is.null(ETA)){ # only random effects present
      ETA0 <- list(X=matrix(rep(1,times=length(y)), ncol=1),K=diag(dim(matrix(rep(1,times=length(y))))[2]))
      ETA <- c(list(ETA0),ETA)
    }
    if(!is.null(X) & !is.null(ETA)){ # both present, extract xm from X, check double list
      if(is.list(X)){
        if(is.list(X[[1]])){
          ETA0 <- list(X=X[[1]][[1]],K=diag(dim(X[[1]][[1]])[2]))
          ETA <- c(list(ETA0),ETA)
        }else{
          ETA0 <- list(X=X[[1]],K=diag(dim(X[[1]])[2]))
          ETA <- c(list(ETA0),ETA)
        }
      }else{
        X=X
        K <- diag(dim(X)[2])
        ETA0 <- list(X=X,K=K)
        ETA <- c(list(ETA0),ETA)
      }
    }
    ## R matrix
    if(is.null(R)){R <- diag(length(y))} # dim(x[[1]])[2]
    # markers do not allow missing data
    fixed <- which(unlist(lapply(ETA, function(x){names(x)[1]})) == "X")
    if(length(fixed) == 0){
      X <- matrix(rep(1,times=length(y)), ncol=1)
      K <- diag(dim(X)[2])
      ETA0 <- list(X=X,K=K)
      ETA <- c(list(ETA0),ETA)
    }
    # add an identity matrix to all effects in ETA that did not provide a var-cov matrix
    ETA <- lapply(ETA, function(x){if(length(x) == 1){x[[2]] <- as(diag(dim(x[[1]])[2]), Class = "sparseMatrix")}else{x <- x}; return(x)})
    # convert everything to sparse matrices
    ##
    eta.or <- ETA
    eta.or <- lapply(eta.or, function(x){lapply(x, as.matrix)}) # put back everything as matrices again
    ##
    ETA2 <- ETA; y2 <- y ; good <- which(!is.na(y)) # make a copy and identify no missing data
    ETA <- lapply(ETA, function(x,good){x[[1]] <- x[[1]][good,]; x[[2]]<- x[[2]]; return(x)}, good=good)
    y <- y[good]
    R <- R[good,good]
    ################
    nran <- length(which(unlist(lapply(ETA, function(x){names(x)[1]})) == "Z")) # number of random effects
    var.com <- vector(mode="list", length=nran+1) # empty list for values of variance components plus error
    if(is.null(init)){var.y=var(y, na.rm=TRUE); var.com <- lapply(var.com, function(x){x=var.y/nran})}else{var.com <- as.list(init)}
    var.com.sto <- var.com # list to store the values of variance components for each iteration
    if(is.null(init)){
      var.com[[length(var.com)]]=var.y
    }
    ETA <- lapply(ETA, function(x){lapply(x, as.matrix)}) # put back everything as matrices again
    # ETA[[1]] is incidence matrix and ETA[[2]] is var-cov matrix
    
    ####################
    ## initialize the progress bar
    if(!silent){
      count <- 0
      tot <- iters
      pb <- txtProgressBar(style = 3)
      setTxtProgressBar(pb, 0)
    }
    #####################
    if(is.null(init)){varE = var(y, na.rm=TRUE)}else{varE=init[length(init)]}#initial value
    conv=0
    wi=0
    change=1
    #ready=2
    if(is.null(forced)){ # &&&&&&&&&&&&&&&&& IF NOT FORCED &&&&&&&&&&&&&&&&&&&&&
      
      while (conv==0) { ## START!!!!! ===============================
        wi=wi+1
        ###################
        if(!silent){
          count <- count + 1
        }
        ###################
        axs = lapply(var.com, function(x){varE/x}) #varE/varGCA1
        CC = numeric()
        ################### ROWS OF MME
        for(j in 1:length(ETA)){ # for each variance component
          prov <- numeric()
          for(k in 1:length(ETA)){ #multiply it for all other variance components
            if(j == k & names(ETA[[j]])[1] != "X"){ ## diagonal element of CC and not fixed
              ranked <- rankMatrix(as.matrix(ETA[[j]][[2]]))[1]
              #if(ranked == dim(as.matrix(ETA[[j]][[2]]))[2]){
                solvedK <- solve(as.matrix(ETA[[j]][[2]]))
              #}else{
              #  stop(cat("\nYour covariance matrices are not full rank, please try the AI algorithm\n")
              #       ,call. = FALSE)
                #solvedK <- ginv(as.matrix(ETA[[j]][[2]]))
              #}
              res <- crossprod(ETA[[j]][[1]], ETA[[k]][[1]]) + (as.vector(axs[[j-1]]) * solvedK) # var(e)/var(x) K-
              ## we used j-1 because the 2nd element of ETA is the first var.component and the var(e) is never used in here
            }else{
              res <- crossprod(ETA[[j]][[1]], ETA[[k]][[1]])
            }
            prov <- cbind(prov,res) # C-BINDING
          }
          CC <- rbind(CC,prov) # R-BINDING
        } # which((C == CC) == FALSE)
        ##############################
        l <- lapply(ETA, function(x,y){t(x[[1]])%*%y}, y=y)
        l2 <- as.matrix(unlist(l))
        #rownames(l2) <- c(unlist(lapply(ETA, function(x){colnames(x[[1]])})))
        ## inverse C, which is the coefficient matrix
        CInv<-solve(CC) 
        thetaHat<-CInv%*%l2 
        #Mstep 
        nn <- lapply(ETA, function(x){dim(x[[1]])[2]}) # number of coefficients to estimate
        nn <- lapply(nn, function(x){if(is.null(x)){x=1}else{x=x};return(x)}) # correct for vectors such as intercepts with only one column present as sinlge vectors sometime
        pairs.a = list(NA)
        for(h in 1:length(nn)){
          pairs.a[[h]] <- ((sum(unlist(nn[1:h])) - (unlist(nn[h])-1) ) : sum(unlist(nn[1:h])))
        }
        ## ADJUST VARIANCE COMPONENTS FOR THIS ITERATION
        # adjust error variance
        now <- 0
        for(f in 1:length(pairs.a)){now <- now + crossprod(as.matrix(ETA[[f]][[1]])%*%thetaHat[pairs.a[[f]],],y)}
        varE = ( crossprod(y) - now ) / (length(y)-nn[[1]])
        # adjust the rest
        indexK <- which(unlist(lapply(ETA, function(x){names(x)[1]})) == "Z")
        ##
        track <- var.com ## keep track of old variance component before update
        ##
        for(k in 1:(length(var.com)-1)){ # adjust var comps except ERROR VARIANCE
          rrr <- indexK[k]
          
          #ranked <- rankMatrix(ETA[[rrr]][[2]])
          #if(ranked == dim(ETA[[rrr]][[2]])[2]){
            Kinv <- solve(ETA[[rrr]][[2]])
          #}else{
          #  Kinv <- ginv(ETA[[rrr]][[2]])
          #}
          
          var.com[[k]] = ( ( t(thetaHat[pairs.a[[rrr]],]) %*% Kinv  %*% thetaHat[pairs.a[[rrr]],] ) + matrixcalc::matrix.trace(Kinv%*%CInv[pairs.a[[rrr]],pairs.a[[rrr]]]*as.numeric(varE)))/nrow(ETA[[rrr]][[2]]) 
        }
        var.com[[length(var.com)]] = varE # last variance component
        #######################################
        ## store the results for each iteration
        for(k in 1:length(var.com.sto)){var.com.sto[[k]] <- c(var.com.sto[[k]],var.com[[k]])}
        ###
        if(is.null(nanana)){
          varosss <- c(paste("u.",1:n.var.com, sep=""))
        }else{
          varosss <- nanana
        }
        
        lege2 <- list()
        for(k in 1:length(var.com)){
          if(k == length(var.com)){
            lege2[[k]] <- paste("Var(Error):")
          }else{
            lege2[[k]] <- paste("Var(",varosss[k],"):",sep="")
          }
        }
        ##
        if(draw){
          ylim <- max(unlist(var.com.sto), na.rm=TRUE)
          my.palette <- brewer.pal(7,"Accent")
          
          layout(matrix(1,1,1))
          plot(var.com.sto[[1]],ylim=c(0,ylim),type="l",col=my.palette[1],lwd=3, xlim=c(0,iters), xaxt="n", las=2, main="Expectation-Maximization algorithm results", ylab="Value of the variance component", xlab="Iteration to be processed to reach convergence", cex.axis=.8) 
          axis(1, las=1, at=0:10000, labels=0:10000, cex.axis=.8)
          for(t in 1:length(var.com)){lines(var.com.sto[[t]],col=my.palette[t],lwd=3)} 
          ww <- length(var.com.sto[[1]])
          lege <- list()
          #lege2 <- list()
          for(k in 1:length(var.com)){
            if(k == length(var.com)){
              lege[[k]] <- paste("Var(e):",round(var.com.sto[[k]][ww],4), sep="")
              #lege2[[k]] <- paste("Var(e):")
            }else{
              lege[[k]] <- paste("Var(",varosss[k],"):",round(var.com.sto[[k]][ww],4), sep="")
              #lege2[[k]] <- paste("Var(u",k,"):",sep="")
            }
          }
          legend("topright",bty="n", col=my.palette, lty=1, legend=unlist(lege), cex=.75)
        }
        ################################
        if(!silent){
          setTxtProgressBar(pb, (count/tot))### keep filling the progress bar
        }
        ################################
        if(wi > 1){
          #change=abs(sum(unlist(var.com) - unlist(track)))
          uuu1 <- cbind(unlist(var.com), unlist(track))
          #uuu2 <- apply(uuu1,1,function(x){if(abs(x[1] - x[2]) < 0.0001){y <- TRUE}else{y <- FALSE}; return(y)})
          uuu2 <- apply(uuu1,1,function(x){if(abs(x[1] - x[2]) < x[1]/100){y <- TRUE}else{y <- FALSE}; return(y)})
          change <- length(which(uuu2))
        }
        ####
        if (change == length(var.com)  | wi == iters ){ ## CONVERGENCE change < 0.000000001
          conv=1
          if(!silent){
            setTxtProgressBar(pb, (tot/tot))### keep filling the progress bar
          }
          if(wi == iters){ # if we didn't reach convergence
            if(!silent){
              ylim <- max(unlist(var.com.sto), na.rm=TRUE)
              my.palette <- brewer.pal(7,"Accent")
              layout(matrix(1,1,1))
              plot(var.com.sto[[1]],ylim=c(0,ylim),type="l",col=my.palette[1],lwd=3, xlim=c(0,iters), xaxt="n", las=2, main="Expectation-Maximization algorithm results", ylab="Value of the variance component", xlab="Iteration to be processed to reach convergence", cex.axis=.8, cex.main=.8) 
              axis(1, las=1, at=0:10000, labels=0:10000, cex.axis=.8)
              for(t in 1:length(var.com)){lines(var.com.sto[[t]],col=my.palette[t],lwd=3)}
              ww <- length(var.com.sto[[1]])
              lege <- list()
              #lege2 <- list()
              for(k in 1:length(var.com)){
                if(k == length(var.com)){
                  lege[[k]] <- paste("Var(e):",round(var.com.sto[[k]][ww],4), sep="")
                  #lege2[[k]] <- paste("Var(e):")
                }else{
                  lege[[k]] <- paste("Var(",varosss[k],"):",round(var.com.sto[[k]][ww],4), sep="")
                  #lege2[[k]] <- paste("Var(u",k,"):",sep="")
                }
              }
              legend("topleft",bty="n", col=my.palette, lty=1, legend=unlist(lege), cex=.75)
              cat("\nMaximum number of iterations reached with no convergence using the EM algorithm, look at the variance components change over iterations (plot) and be cautious using the variance components estimated if they don't look steady")
            } 
          }########### END of if we didn't reach convergence
        }
        ################################
      } ### ======================= END while loop
      
    }else{ ## &&&&&&&&&&&&&&&&& IF FORCED &&&&&&&&&&&&&&&&&&&&&
      setTxtProgressBar(pb, (tot/tot))
      
      cat("\nVariance components forced\n")
      var.com <- as.vector(forced)
      lege2 <- list()
      for(k in 1:length(var.com)){
        if(k == length(var.com)){
          lege2[[k]] <- paste("Var(e):")
        }else{
          lege2[[k]] <- paste("Var(u",k,"):",sep="")
        }
      }
      ##
      if(is.null(nanana)){
        varosss <- c(paste("u.",1:n.var.com, sep=""))
      }else{
        varosss <- nanana
      }
    } ## &&&&&&&&&&&&&&&&& END OF FORCING &&&&&&&&&&&&&&&&&&&&&
    ################################
    ################################
    ################################
    ################################
    ################################
    ################################
    ## axs are lambdas var(e)/var(x)
    # calculate phenotypic variance var(y)= var(Zu+e) = ZGZ+R
    xm <- as.matrix(ETA[[1]][[1]])
    # RANDOM VARIABLES
    zvar <- which(unlist(lapply(ETA, function(x){names(x)[1]})) == "Z")
    
    V <- matrix(0, nrow=dim(xm)[1], ncol=dim(xm)[1])
    V2 <- R * (as.numeric(varE))
    for (k in zvar) {
      #V + ZZK2 or V + ZZK1 ... , where each K is the covariance matrix for each random effect
      V2 <- V2 + tcrossprod(ETA[[k]][[1]], ETA[[k]][[1]] %*% (ETA[[k]][[2]]*as.numeric(var.com[k-1])) )
      V <- V + tcrossprod(ETA[[k]][[1]], ETA[[k]][[1]] %*% (ETA[[k]][[2]]) ) + (R * (as.numeric(varE)/as.numeric(var.com[k-1])))
    }# then V is just the sum of all Z'KZ including the error
    # therefore V == ZGZ + IR
    # Inverse of total variance
    Vinv <- solve(as(V, Class = "sparseMatrix"))
    Vinv2 <- solve(as(V2,Class = "sparseMatrix"))
    
    X <- ETA[[fixed]][[1]]
    xvx=t(X)%*%Vinv%*%X # X'V-X
    xvxi=solve(xvx) # (X'V-X)-
    pm=Vinv-Vinv%*%X%*%xvxi%*%t(X)%*%Vinv # P = Vinv - [Vinv X (X'V-X)- X Vinv]
    if(REML==TRUE){
      #logL=-0.5*(log(det(V))+log(det(xvx))+t(y)%*%pm%*%y) # log likelihood, problem
      logL=as.numeric(-0.5*((determinant(V, logarithm = TRUE)$modulus[[1]])+determinant(xvx, logarithm = TRUE)$modulus[[1]]+t(y)%*%(pm%*%y)))
    }else{
      #logL=-0.5*(log(det(V)) + t(y)%*%pm%*%y ) # log likelihood, problem
      logL=as.numeric(-0.5*((determinant(V, logarithm = TRUE)$modulus[[1]]) + t(y)%*%(pm%*%y )))
    }
    #logL=NULL
    # Estimates of coefficients of the linear mixed model
    ## Fixed effects
    # B=(X'X V- XVy)-
    # FIXED VARIABLES
    xvar <- which(unlist(lapply(ETA, function(x){names(x)[1]})) == "X")
    beta <- solve(crossprod(ETA[[xvar]][[1]], Vinv %*% ETA[[xvar]][[1]])) %*% crossprod(ETA[[xvar]][[1]], Vinv %*% y)
    #dasa <- lapply(ETA, function(x){dim(x[[1]])[2]})
    #l <- lapply(ETA, function(x,y){t(x[[1]])%*%y}, y=y)
    #l2 <- as.matrix(unlist(l))
    #rownames(l2) <- c(unlist(lapply(ETA, function(x){colnames(x[[1]])})))
    ## inverse C, which is the coefficient matrix
    #theta<-CInv%*%l2
    #beta <- as.matrix(theta[c(1:dasa[[1]]),1]) 
    ## Random effects
    u <- list()
    for (k in zvar) {
      #u = KZ'V- (y- XB)
      u[[k]] <- ETA[[k]][[2]] %*% t(ETA[[k]][[1]]) %*% Vinv %*% (y - ETA[[xvar]][[1]] %*% beta)
    }
    u <- u[zvar]
    #a0 <- unlist(dasa)
    #a1 <- a0[1]
    #if(length(a0) > 1){
    #for(r in 2:length(a0)){
    # a1[r] <- sum(a0[1:r])
    #}}
    #aa <- a1+1; aa <- aa[-c(length(aa))]
    #a2 <- c(1,aa)
    #u <- list()
    #for (k in zvar) {
    #  b1 <- a1[k]
    #  b2 <- a2[k]
    #  u[[k]] <- as.matrix(theta[c(b2:b1),1] )
    #}
    #u <- u[zvar]
    # variance for the variance or 2nd moment sigma2^2 = sigma4
    Var.u.hat <- list()
    PEV.u.hat <- list()
    #for(i in zvar){
    # Var.u.hat[[i]] <- diag(CInv[pairs.a[[i]],pairs.a[[i]]])
    #  PEV.u.hat[[i]] <- diag(CInv[pairs.a[[i]],pairs.a[[i]]]) * as.numeric(varE) # https://aipl.arsusda.gov/publish/jds/1988/71_27.pdf
    #}
    for(i in zvar){ # we need "pm" matrix
      Var.u.hat[[i]] <- (as.numeric(var.com[i-1])^2) *  ( crossprod(ETA[[i]][[1]]%*%ETA[[i]][[2]], pm)  %*%  (ETA[[i]][[1]]%*%ETA[[i]][[2]])   ) # sigma^4 ZKP ZK
      PEV.u.hat[[i]] <- as.numeric(var.com[i-1]) * ETA[[i]][[2]] - Var.u.hat[[i]]  # standard errors (SE) for each individual
    }
    
    Var.u.hat <- Var.u.hat[zvar]
    PEV.u.hat <- PEV.u.hat[zvar]
    Var.beta.hat <- solve(crossprod(xm, Vinv2 %*% xm)) 
    #Var.beta.hat=NULL
    ###############
    AIC = as.vector((-2 * logL ) + ( 2 * dim(X)[2]))
    BIC = as.vector((-2 * logL ) + ( log(length(y)) * dim(X)[2]))
    
    fitted.y <- eta.or[[1]][[1]] %*% beta
    fitted.u <- 0
    for(h in 2:length(eta.or)){
      fitted.y <- fitted.y + (eta.or[[h]][[1]] %*% u[[h-1]])
      if(h > 1){
        fitted.u <- fitted.u + (eta.or[[h]][[1]] %*% u[[h-1]])
      }
    }
    fitted.y.good <- fitted.y[good]
    residuals2 <- y - (ETA[[1]][[1]] %*% beta) # regular y - Xb
    residuals3 <- y - fitted.y[good] # conditional residuals
    out1 <- unlist(var.com); names(out1) <- unlist(lege2) # variance components
    ## add names to fixed and random effects
    rownames(beta) <- colnames(ETA[[1]][[1]])
    for(i in 2:length(ETA)){
      rownames(u[[i-1]]) <- colnames(ETA[[i]][[1]])
    }
    #########################################
    if(!is.null(nanana)){
      names(u) <- nanana
      names(Var.u.hat) <- nanana
      names(PEV.u.hat) <- nanana
    }
    logL <- as.vector(logL)
    ##################################
    VCC <- (as.matrix(out1,ncol=1))
    rownames(VCC) <- c(paste("Var(",varosss,")",sep=""), "Var(Error)")
    output <- list(var.comp=VCC, V.inv = Vinv, u.hat=u, Var.u.hat=Var.u.hat, 
                   PEV.u.hat=PEV.u.hat, beta.hat=beta, Var.beta.hat=Var.beta.hat,
                   LL=logL, AIC=AIC, BIC=BIC, X=ETA[[xvar]][[1]], 
                   fitted.y=fitted.y, fitted.u=fitted.u, residuals=residuals2, 
                   cond.residuals=residuals3, fitted.y.good=fitted.y.good, Z=Zsp, K=Ksp)
    #cat("\n\nVariance components:\n\n")
    #print(out1)
  }
  return(output)
}

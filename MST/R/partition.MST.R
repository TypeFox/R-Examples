partition.MST <-
function(dat, test=NULL, name="1", method=c("marginal", "gamma.frailty", "exp.frailty"),     		# method= Choose between marginal, gamma.frailty, & exp.frailty
                          col.time, col.status, col.id, col.split.var, col.ctg=NULL,
                          minsplit=20, min.nevents=5, max.depth=10, mtry=length(col.split.var),
                          cont.split=c("distinct","percentiles"), delta=0.05, nCutPoints=50,details=FALSE){
  method<-match.arg(method,c("marginal", "gamma.frailty", "exp.frailty"))
  cont.split<-match.arg(cont.split,c("distinct", "percentiles"))

  suppressWarnings(warning("coxph")); suppressWarnings(warning("coxpenal.fit"))
  call <- match.call(); out <- match.call(expand.dots = FALSE)
  out$info <- out$name.l <- out$name.r <- out$left <- out$right <- out$... <- NULL
  n <- nrow(dat);  vnames <- colnames(dat); var <- NA; cut <- NA;
  max.score <- -1e50; score.test <- NA; 
  depth <- nchar(name) # CONTROL THE MAX TREE DEPTH
  time <- dat[, col.time]; status <- dat[, col.status]; id <- dat[, col.id]
  n.event <- sum(status==1)  # COMPUTE NUMBER OF EVENT TIMES
  if (details) print(paste("Tree ", name, ": the sample size: ", n, " and num of events: ", n.event))
  if (details) print(depth)
  if (sum(!is.element(col.ctg, col.split.var)) >0) warning("col.ctg should be a subset of col.split.var.")	
  if (depth <= max.depth && n >= minsplit && n.event > min.nevents) {
    if(method=="exp.frailty"){
      # FIT THE NULL MODEL
      X.i <- aggregate(x=time, by=list(id), FUN=sum)$x
      Delta.i <- aggregate(x=status, by=list(id), FUN=sum)$x
      dat1 <- data.frame(X.i=X.i, Delta.i=Delta.i)
      x <- optim(par=c(1, 1), fn=loglik0, gr=gr0, method = "SANN", control = list(maxit=200), hessian=FALSE, dat=dat1)    	############## GLOBAL OPTIMIZATION (SIMULATED ANNEALING)
      x <- optim(par=x$par, fn=loglik0, gr=gr0, method = "BFGS", hessian = TRUE, control = list(maxit=50), dat=dat1)  		############## 
      # x <- optim(par=x$par, fn=loglik0, method = "BFGS", hessian = TRUE, dat=dat1)  							############## NUMERIC DERIVATIVES
      if (!is.null(x$message)) print(paste("At the ", i, "-th run, there are ", x$message, " in fitting the null model.")) 
      beta0 <- x$par[1]; v0 <- x$par[2]; I0.inv <- solve(x$hessian)
      Ai <- 1/v0 + Delta.i; Bi <- 1/v0 + exp(beta0)*X.i
    }
    
    # SEARCH FOR BEST SPLIT OF THE TRAINING DAT
    for(i in sort(sample(col.split.var, size=mtry, replace=FALSE))) {
      z <- dat[,i]; v.name <- vnames[i]; temp <- sort(unique(na.omit(z)));
      if(length(temp) > 1) {
        if (is.element(i,col.ctg)){ zcut <- power.set(temp)    				############################ CLASS VARIABLE
        } else {
          if(cont.split=="distinct" | length(temp) <= 50){zcut <- temp[-length(temp)]
          } else if(cont.split=="percentiles" & (delta>0.2 | delta<0.01 | nCutPoints < 5)){stop("Choice of percentile cutpoints too small")
          } else if(cont.split=="percentiles"){zcut <- quantile(temp, probs = seq(delta,1-delta,length=nCutPoints))} ## TAKE PERCENTILEs OF CUTOFF POINT
        }
        for (j in zcut) {
          score <- NA
          if (is.element(i,col.ctg)) {grp <- sign(is.element(z, j)); cut1 <- paste(j, collapse=" ")      ##########################
          } else  {grp <- sign(z <= j); cut1 <- as.character(j)}
          if (method=="marginal"){ score <- splitting.stat.MST1(time, status, id, z=grp, min.nevents)
          } else if (method=="gamma.frailty"){ score <- splitting.stat.MST2(time, status, id, z=grp, min.nevents, method="wald.test")
          } else if (method=="exp.frailty"){
            n.R1 <- sum(grp==0&status==1); n.L1 <- sum(grp==1&status==1);
            if (min(n.R1, n.L1)>= min.nevents){
              # COMPUTE THE SCORE TEST STATISTIC
              mi <- aggregate(x=time*grp, by=list(id), FUN=sum)$x
              zdi <- aggregate(x=status*grp, by=list(id), FUN=sum)$x
              U <- sum(zdi) - exp(beta0)*sum(mi*Ai/Bi)
              I11 <- exp(beta0)*sum(mi*Ai*(Bi - exp(beta0)*mi)/Bi^2)
              I1theta <- c(exp(beta0)*sum(mi*Ai*(Bi-exp(beta0)*X.i)/Bi^2), exp(beta0)/v0^2*sum(mi*(Ai-Bi)/Bi^2))
              score <- U^2/(I11 - t(I1theta)%*%I0.inv%*%I1theta)
              score <- max(0, score)
            }
          }
          if (identical(score, numeric(0))) score <- NA		# TO DEAL WITH THE PROBLEM THAT THE WALD TEST COULD RETURN numeric(0). 
          if (!is.na(score) && score >= max.score) {max.score <- score; var <- i; vname <- v.name; cut <- cut1; best.cut<-j; grp.best <- grp}
          if (details) {print(cbind(var=i, v.name=v.name, cut=j, score=score, max.score=max.score));} # print(is.null(score)); print(score)}
        }
      }
    }
  }
  
  # THE TEST SAMPLE
  if (!(is.null(test))) {
    n.test <- nrow(test);
    if (!(is.na(var)) && max.score!=0) {
      time.test <- test[, col.time]; status.test <- test[, col.status]; id.test <- test[, col.id]	
      # COMPUTE THE SCORE STAT BASED ON THE TEST SAMPLE
      if (is.element(var,col.ctg)){ grp.test <- sign(is.element(test[,var], best.cut))                       ############################
      } else {grp.test <- sign(test[,var] <= best.cut)}
      if (method=="marginal") {score.test <- splitting.stat.MST1(time.test, status.test, id.test, z=grp.test, min.nevents)
      } else if (method=="gamma.frailty") {score.test <- splitting.stat.MST2(time.test, status.test, id.test, z=grp.test, min.nevents, method="wald.test")
      } else if (method=="exp.frailty") {
        nL1.test <- sum(grp.test==1&status.test==1); nR1.test <- sum(grp.test==0&status.test==1);
        if (min(nL1.test, nR1.test) >= min.nevents) {
          # fit the null model with the test sample
          X.i <- aggregate(x=time.test, by=list(id.test), FUN=sum)$x
          Delta.i <- aggregate(x=status.test, by=list(id.test), FUN=sum)$x
          dat2 <- data.frame(X.i=X.i, Delta.i=Delta.i)
          x1 <- optim(c(2,1.5), fn=loglik0, gr=gr0, method = "Nelder-Mead", hessian = TRUE, dat=dat2)
          beta0 <- x1$par[1]; v0 <- x1$par[2]
          I0.inv <- solve(x1$hessian)
          Ai <- 1/v0 + Delta.i; Bi <- 1/v0 + exp(beta0)*X.i
          # compute the score test
          mi <- aggregate(x=time.test*grp.test, by=list(id.test), FUN=sum)$x
          zdi <- aggregate(x=status.test*grp.test, by=list(id.test), FUN=sum)$x
          U <- sum(zdi) - exp(beta0)*sum(mi*Ai/Bi)
          I11 <- exp(beta0)*sum(mi*Ai*(Bi - exp(beta0)*mi)/Bi^2)
          I1theta <- c(exp(beta0)*sum(mi*Ai*(Bi-exp(beta0)*X.i)/Bi^2), exp(beta0)/v0^2*sum(mi*(Ai-Bi)/Bi^2))
          score.test <- U^2/(I11 - t(I1theta)%*%I0.inv%*%I1theta)
          score.test <- ifelse(score.test<0, NA, score.test)
        }
      }
    }
  } else {n.test <- n; score.test <- max.score}
  
  if (!is.na(score.test) && !is.na(var)){
    out$name.l <- paste(name, 1, sep=""); out$name.r <- paste(name, 2, sep="")
    if(!is.null(test)) {out$left.test <- test[grp.test==1,  ]; out$right.test <- test[grp.test==0,]}
    else {out$left.test <- out$right.test <- NULL}
    out$left  <- dat[grp.best==1,];  out$right <- dat[grp.best==0,] 
  } else {var <- NA; vname <- NA; cut <- NA;  max.score <- score.test <- NA}
  out$info <- data.frame(node=name, size = n, var = var, vname=vname, cut= cut, 
                         score=ifelse(max.score==-1e10, NA, max.score), size.test=n.test, score.test)
  out 
}

rs.diff <- function (formula = formula(data), data = parent.frame(), ratetable = relsurv::slopop, 
                     na.action,precision=1) 
  
  #formula: for example Surv(time,cens)~sex
  #data: the observed data set
  #ratetable: the population mortality tables
{
  
  call <- match.call()
  rform <- rformulate(formula, data, ratetable, na.action)  		#get the data ready 
  data <- rform$data								#the data set
  p <- rform$m								#number of covariates
  if (p > 0) 									#if covariates 
    data$Xs <- strata(rform$X[, ,drop=FALSE ])	#make groups according to covariates
  else data$Xs <- rep(1, nrow(data))        					#if no covariates, just put 1
  # Xs is a  vector of factors determining the groups we wish to compare
  strats <- rform$strata.keep # added for strata
  str_num <- length(levels(strats)) # number of strata

  out <- NULL
  out$n <- table(data$Xs)	#table of groups
  out$time <- out$n.risk <- out$n.event <- out$n.censor <- out$surv <- out$std.err <- out$groups <- NULL
  
  #TIMES ARE EQUAL FOR ALL GROUPS  								
  if(!precision)tis <- sort(unique(rform$Y))					#unique times
  else{
  extra <- as.numeric(seq(1,max(rform$Y),by=precision))
  tis <- sort(union(extra,rform$Y))					#1-day long intervals used - to take into the account the continuity of the pop. part
  }
  # start working
  kgroups <- length(out$n)  						#number of groups
  
  if (kgroups == 1) stop("There is only one group in your data. You should choose another variable.")
  
  w.risk <- w.event <- dnisisq <- array(NA,dim=c(length(tis),length(out$n),str_num))			#MATRIX - COLUMNS ARE GROUPS, ROWS ARE TIMES,levels are strata
  #numOfSmallGrps <- 0
  numOfFewEvents <- 0
  for (s in 1:str_num){ # added for strata
    for (kt in 1:kgroups) {						#for each group
      inx <- which(data$Xs == names(out$n)[kt] & strats == levels(strats)[s])				#individuals within this group
      #if (length(inx)<10)numOfSmallGrps <- numOfSmallGrps + 1
      temp <- exp.prep(rform$R[inx,,drop=FALSE],rform$Y[inx],ratetable,rform$status[inx],times=tis,fast=TRUE)	#calculate the values for each interval of time
      out$time <- c(out$time, tis)						#add times
      out$n.risk <- c(out$n.risk, temp$yi)					#add number at risk for each time
      out$n.event <- c(out$n.event, temp$dni)					#add number of events for each time
      if (sum(temp$dni) < 10) numOfFewEvents <- numOfFewEvents + 1
      out$n.censor <- c(out$n.censor,  c(-diff(temp$yi),temp$yi[length(temp$yi)]) - temp$dni) 	#add number of censored for each time
      
      w.risk[,kt,s] <- temp$yisi						#Y_h^w
      w.event[,kt,s] <- temp$dnisi - temp$yidlisi				#dN_eh^w
      dnisisq[,kt,s] <- temp$dnisisq						#dN/S_p^2
      out$groups <- c(out$groups, length(tis))				#number of times in this group
    }
  }
  #if (numOfSmallGrps > 0) warning(numOfSmallGrps, " out of ", kgroups*str_num, " groups is/are smaller than 10.")
  if (numOfFewEvents > 0) warning("In ", numOfFewEvents, " out of ", kgroups*str_num, " groups there are less than 10 events.")
    
  w.risk.total <- apply(w.risk,c(1,3),sum)					#sum over all individuals at each time point ## Y_{.,s}^w
  w.event.total <- apply(w.event,c(1,3),sum)       #sum over all individuals at each time point ## dN_{E,.,s}^w
  
  zs <- rep(0,kgroups) # added for strata
  for (s in 1:str_num){
    # znotraj danega stratuma
    inx_str <- which(w.risk.total[,s] > 0)
    zhst <- w.event[inx_str,,s,drop=FALSE] - w.risk[inx_str,,s,drop=FALSE]/w.risk.total[inx_str,s]*w.event.total[inx_str,s]   #value under the integral of zh
    # integriramo po casu - sestejemo po casih dogodkov
    
    zhs <- apply(zhst,2,sum)    # the vector of test statistics
    zs <- zs + zhs 
  }
  # cat("vektor testnih statistik je = \n")
  # print(zs)
    
  #covariance matrix:
  covmats <- matrix(0,nrow=kgroups,ncol=kgroups)
  d <- diag(kgroups)  					#identity matrix of groups size (for the kronecker deltas)
  for (s in 1:str_num){
    underint <- 0
    inx_str <- which(w.risk.total[,s] > 0)
    for(kt in 1:kgroups){						#matrix calculation through the groups
      ys <- matrix(d[kt,],nrow=length(inx_str),ncol=kgroups,byrow=T) - w.risk[inx_str,,s]/w.risk.total[inx_str,s]	#preparing the matrix for the first two terms	
      #yslist <- apply(apply(ys,1,list),unlist)			#a list, each row of ys (each time point) represents one item
      yslist <- as.list(data.frame(t(ys)))				#a list, each row of ys (each time point) represents one item
      yprod <- lapply(yslist,function(x)outer(x,x))			#a list of matrices with y products through all the time points, 
      yproda <- array(unlist(yprod),dim=c(kgroups,kgroups,length(inx_str)))#y terms transformed to an array
      dnisisqa <- array(rep(dnisisq[,kt,s],each=kgroups^2),dim=c(kgroups,kgroups,length(inx_str)))	#dnisisq terms transformed into an array of equal size
      underint <- underint +  yproda * dnisisqa			#the terms under the integral
    }
    covmat <- apply(underint,1:2,sum)						#summing down the array
    covmats <- covmats + covmat
  }
  # cat("kovariancna matrika je = \n")
  # print(covmats)
  # del za testiranje
  
  zs <-    zs[-kgroups]							# the last one is deleted
  zs <- matrix(zs,nrow=1)
  # print(covmats)
  covmats <- covmats[-kgroups,-kgroups,drop=F]
  # print(covmats)
  test.stat <- zs %*% solve(covmats) %*% t(zs)
  p.value <- 1-pchisq(test.stat,df=kgroups-1)
  
  
  names(out$groups) <- names(out$n)
  if (p == 0) out$groups <- NULL						#if no covariates
  out$n <- as.vector(out$n)
  out$call <- call
  #class(out) <- c("survdiff", "rs.surv")
  #cat(zh)
  out$zh <- zs  
  out$covmat <- covmats
  out$test.stat <- test.stat
  out$p.value <- p.value
  out$df <- kgroups-1
  class(out) <- "rsdiff"
  out
}
 print.rsdiff <- function(x,...){
 invisible(cat("Value of test statistic:", x$test.stat, "\n"))
 invisible(cat("Degrees of freedom:", x$df, "\n"))
 invisible(cat("P value:", x$p.value, "\n"))
 }

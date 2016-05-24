
recode <- function(sym_list, ref_sym){
  the.list <- setdiff(sort(unique(sym_list)),ref_sym)
  the.mat <- matrix(0,nrow=length(sym_list),ncol=length(the.list))
  for(i in 1:length(the.list)) the.mat[,i] <- as.numeric(sym_list==the.list[i])
  return(the.mat)
}

##  this is just to estimate permuation based variance
affun_temp <- function(the.form=NULL,the.data=NULL,nsample_perm = 10000,prev=.015,allperm=TRUE){
  
  d=model.frame(the.form,the.data)
    
  y <- d[,1]
  w <- rep(1,length(y))
  if(!is.na(prev)){ 
    
    thenas <- apply(d[,2:ncol(d)],1,function(x){sum(is.na(x))>0})
    w[y==0] <- (sum(y==1 & !thenas)/sum(y==0 & !thenas))*((1-prev)/prev)
    w[y==1] <-  1
    
  }
  m=the.form
  
  m$coefficients <- m$coefficients[!is.na(m$coefficients)]
  mean_beta <- m$coefficients
  cov_beta <- summary(m)$cov.unscaled
  
  ##  remove columns that are 0
  d <- d[,apply(d,2,function(x){sum(abs(x))})!=0]
  ## remove weight column
  if(colnames(d)[ncol(d)]=="(weights)") d <- d[,1:(ncol(d)-1)]
  
  
  startcol <- grep("(c|o)[0-9]_",colnames(d),perl=TRUE)[1]
  
  vars=length(unique(substring(colnames(d)[startcol:ncol(d)],1,3)))
  n=nrow(d)
  obs.cases=sum(d[,1])
  continuous_indices <-  (startcol-1) + grep("c[0-9]+_*",unique(substring(colnames(d)[startcol:ncol(d)],1,3)),perl=TRUE)
  lc <- length(continuous_indices)
  ordinal_indices <-  (startcol-1) + grep("o[0-9]+_*",unique(substring(colnames(d)[startcol:ncol(d)],1,3)),perl=TRUE)
  lo <- length(ordinal_indices)
  combined_indices <- c(ordinal_indices, continuous_indices)
  df_ordinal <- numeric(length(ordinal_indices))
  if(length(df_ordinal) > 0){
    for(j in 1:length(df_ordinal)) df_ordinal[j] <- sum(grep(paste("o",j,"_",sep=""),substring(colnames(d),1,nchar(paste("o",j,"_",sep=""))),perl=TRUE)>0)
  }
  combined_df <- c(df_ordinal)
  combined_df <- combined_df[combined_df > 0] ###  remove 0's
  
  if(allperm) indices=perm(vars,vars,startcol:(vars+startcol-1))   ### create sampling frame of every possible permtuation
  ### in the case of not using all the permutations - just sample some of them
  if(!allperm){
    
    indices <- matrix(0,nrow=nsample_perm,ncol=vars)
    for(i in 1:nsample_perm){
      indices[i,] <- sample(startcol:(vars+startcol-1),vars,replace=FALSE)
    }
    
  }
  
  
  if(sum(combined_df) > 0){
    the.dim <- ncol(indices)+sum(combined_df-1)
    indices_w <- matrix(0,nrow(indices),the.dim)
    
    for(i in 1:nrow(indices)){
      o_vec <- indices[i,]
      new_vec <- numeric(the.dim)
      break_p <- (1:length(o_vec))[o_vec %in% combined_indices]
      nbreak_p <- setdiff((1:length(o_vec)),break_p)
      
      # order the break_p interms of order of old vec
      the_order <- order(o_vec[break_p])  #what is the position of the breakpoint for x1,the breakpoint for x2 ....
      combined_df_ordered <- combined_df[rank(o_vec[break_p])]  #which df comes first, which second and which last
      break_p <- break_p[order(o_vec[break_p])]
      count_col <- 1
      addfactor <- 0
      
      for(k in 1:length(break_p)){
        
         if(the_order[k]>1) start_i <- sum(combined_df_ordered[1:(the_order[k]-1)]-1)+break_p[k]
        if(the_order[k]==1) start_i <- break_p[k]
        end_i <- sum(combined_df_ordered[1:the_order[k]]-1)+break_p[k]
        new_vec[start_i:end_i] <- (o_vec[break_p[k]]+addfactor):(o_vec[break_p[k]]+sum(combined_df[1:k] -1))
        addfactor <- sum(combined_df[1:k]-1)
      }
      new_vec[new_vec == 0] <- o_vec[nbreak_p]
      
      if((100*i/nrow(indices))%%10 == 0){flush.console()
                                         #print(paste('permutations ', 100*i/nrow(indices),'% done',sep=''))
      }
      indices_w[i,] <- new_vec
    }
    
    indices <- indices_w
  }
  
  
  d1=d
  ##  an error if indices matrix loses its dimensionality
  if(is.null(nrow(indices))) indices <- as.matrix(indices,nrow=length(indices),ncol=1)
  pred.cases.m=matrix(NA,nrow=nrow(indices),ncol=ncol(indices))
  prev.cases.m=matrix(NA,nrow=nrow(indices),ncol=ncol(indices))
  
  ### set minimum values for all the variables
  for(i in 1:nrow(indices)){
    minvalues <- numeric(ncol(d))
    minvalues[1] <- NA
    addfactor <- 0
    
    if(length(df_ordinal)>0) addfactor <- sum(df_ordinal-1)
    
    for(k in 1:ncol(indices)){
      d[,indices[i,k]]=rep(minvalues[indices[i,k]],nrow(d))
      
      eta <- as.matrix(cbind(1,d[,2:ncol(d)]),nrow=nrow(d))%*%as.vector(m$coefficients)
      eta <- exp(eta)
      pred.cases.m[i,k]<- sum(w*eta/(1+eta))
    }
    d=d1
    if((100*i/nrow(indices))%%10 == 0){flush.console()
                                       #print(paste('calculations ', 100*i/nrow(indices),'% done',sep=''))
    }
    
  }
   pred.cases.m=cbind(rep(obs.cases,nrow(indices)),pred.cases.m)
  for(i in 1:nrow(indices)){
    for(k in 1:ncol(indices)){ prev.cases.m[i,indices[i,k]-(startcol-1)]=pred.cases.m[i,k]-pred.cases.m[i,k+1]
    }
  }
  prev.cases.m.new <- as.data.frame(prev.cases.m)
  
  if(length(combined_indices)>0){
    if(length(nbreak_p)>0) prev.cases.m.new <- as.data.frame(cbind(prev.cases.m[,1:length(nbreak_p)],matrix(0,nrow(indices),ncol=length(break_p))))
    if(length(nbreak_p)==0) prev.cases.m.new <- as.data.frame(matrix(0,nrow(indices),ncol=length(break_p)))
    
    ord_colnames <- c()
    if(lo>=1) ord_colnames <- paste('o',1:lo,sep="")
    cont_colnames <- c()
    if(lc>=1) cont_colnames <- paste('c',1:lc,sep="")
    colnames(prev.cases.m.new)[1:length(break_p)] <- c(ord_colnames,cont_colnames)
    
    addfactor <-1
    for(k in 1:length(break_p)){
      if(combined_df[k] == 1) prev.cases.m.new[,length(nbreak_p)+k] <- prev.cases.m[,(length(nbreak_p) + addfactor):(length(nbreak_p) + addfactor)]
      if(combined_df[k] > 1) prev.cases.m.new[,length(nbreak_p)+k] <- apply(prev.cases.m[,(length(nbreak_p) + addfactor):(length(nbreak_p) + sum(combined_df[1:k]))],1,sum)
      addfactor <- sum(combined_df[1:k])+1
    }
  }
  prev.cases.m.new <- cbind(prev.cases.m.new,total=apply(prev.cases.m.new,1,sum))
  colnames(prev.cases.m.new)[ncol(prev.cases.m.new)] <- "total"  
  return(apply(prev.cases.m.new,2,sd)/obs.cases)
}


affun <- function(the.form=NULL,the.data=NULL,nsample_perm = 10000,prev=.015,approx_error=0.001,allperm=TRUE){
  
  d=model.frame(the.form,the.data)
  
  
   
  y <- d[,1]
  w <- rep(1,length(y))
  if(!is.na(prev)){ 
    
    thenas <- apply(d[,2:ncol(d)],1,function(x){sum(is.na(x))>0})
    w[y==0] <- (sum(y==1 & !thenas)/sum(y==0 & !thenas))*((1-prev)/prev)
    w[y==1] <-  1
    
  }
  m=the.form
  
  m$coefficients <- m$coefficients[!is.na(m$coefficients)]
  mean_beta <- m$coefficients
  cov_beta <- summary(m)$cov.unscaled
  
  ##  remove columns that are 0
  d <- d[,apply(d,2,function(x){sum(abs(x))})!=0]
  ## remove weight column
  if(colnames(d)[ncol(d)]=="(weights)") d <- d[,1:(ncol(d)-1)]
  
  
  startcol <- grep("(c|o)[0-9]_",colnames(d),perl=TRUE)[1]
  
  vars=length(unique(substring(colnames(d)[startcol:ncol(d)],1,3)))
  n=nrow(d)
  obs.cases=sum(d[,1])
  continuous_indices <-  (startcol-1) + grep("c[0-9]+_*",unique(substring(colnames(d)[startcol:ncol(d)],1,3)),perl=TRUE)
  lc <- length(continuous_indices)
  ordinal_indices <-  (startcol-1) + grep("o[0-9]+_*",unique(substring(colnames(d)[startcol:ncol(d)],1,3)),perl=TRUE)
  lo <- length(ordinal_indices)
  combined_indices <- c(ordinal_indices, continuous_indices)
  df_ordinal <- numeric(length(ordinal_indices))
  if(length(df_ordinal) > 0){
    for(j in 1:length(df_ordinal)) df_ordinal[j] <- sum(grep(paste("o",j,"_",sep=""),substring(colnames(d),1,nchar(paste("o",j,"_",sep=""))),perl=TRUE)>0)
  }
  combined_df <- c(df_ordinal)
  combined_df <- combined_df[combined_df > 0] ###  remove 0's
  
  if(allperm) indices=perm(vars,vars,startcol:(vars+startcol-1))   ### create sampling frame of every possible permtuation
  ### in the case of not using all the permutations - just sample some of them
  if(!allperm){
    
    if(!is.na(approx_error)){
      est_sd <- affun_temp(the.form,the.data=the.data,nsample_perm = 100,  prev=prev, allperm=allperm)
      nsample_perm <- max(ceiling(est_sd^2*1.96^2/approx_error^2))+1 
      print(paste(nsample_perm, " permutations will be performed",sep=""))
    }
    
    indices <- matrix(0,nrow=nsample_perm,ncol=vars)
    for(i in 1:nsample_perm){
      indices[i,] <- sample(startcol:(vars+startcol-1),vars,replace=FALSE)
    }
    
  }
  
  
  if(sum(combined_df) > 0){
    the.dim <- ncol(indices)+sum(combined_df-1)
    indices_w <- matrix(0,nrow(indices),the.dim)
    
    for(i in 1:nrow(indices)){
      o_vec <- indices[i,]
      new_vec <- numeric(the.dim)
      break_p <- (1:length(o_vec))[o_vec %in% combined_indices]
      nbreak_p <- setdiff((1:length(o_vec)),break_p)
      
       the_order <- order(o_vec[break_p])  #what is the position of the breakpoint for x1,the breakpoint for x2 ....
      combined_df_ordered <- combined_df[rank(o_vec[break_p])]  #which df comes first, which second and which last
      break_p <- break_p[order(o_vec[break_p])]
      count_col <- 1
      addfactor <- 0
      
      for(k in 1:length(break_p)){
        
          if(the_order[k]>1) start_i <- sum(combined_df_ordered[1:(the_order[k]-1)]-1)+break_p[k]
        if(the_order[k]==1) start_i <- break_p[k]
        end_i <- sum(combined_df_ordered[1:the_order[k]]-1)+break_p[k]
        new_vec[start_i:end_i] <- (o_vec[break_p[k]]+addfactor):(o_vec[break_p[k]]+sum(combined_df[1:k] -1))
        addfactor <- sum(combined_df[1:k]-1)
      }
      new_vec[new_vec == 0] <- o_vec[nbreak_p]
      
      if((100*i/nrow(indices))%%10 == 0){flush.console()
                                         print(paste('permutations ', 100*i/nrow(indices),'% done',sep=''))
      }
      indices_w[i,] <- new_vec
    }
    
    indices <- indices_w
  }
  
  
  ####
  d1=d
  ##  an error if indices matrix loses its dimensionality
  if(is.null(nrow(indices))) indices <- as.matrix(indices,nrow=length(indices),ncol=1)
  pred.cases.m=matrix(NA,nrow=nrow(indices),ncol=ncol(indices))
  prev.cases.m=matrix(NA,nrow=nrow(indices),ncol=ncol(indices))
  
  ### set minimum values for all the variables
  for(i in 1:nrow(indices)){
    minvalues <- numeric(ncol(d))
    minvalues[1] <- NA
    addfactor <- 0
    
    if(length(df_ordinal)>0) addfactor <- sum(df_ordinal-1)
    
    for(k in 1:ncol(indices)){
      d[,indices[i,k]]=rep(minvalues[indices[i,k]],nrow(d))
      
      eta <- as.matrix(cbind(1,d[,2:ncol(d)]),nrow=nrow(d))%*%as.vector(m$coefficients)
      eta <- exp(eta)
      pred.cases.m[i,k]<- sum(w*eta/(1+eta))
    }
    d=d1
    if((100*i/nrow(indices))%%10 == 0){flush.console()
                                       print(paste('calculations ', 100*i/nrow(indices),'% done',sep=''))
    }
    
  }
  pred.cases.m=cbind(rep(obs.cases,nrow(indices)),pred.cases.m)
  for(i in 1:nrow(indices)){
    for(k in 1:ncol(indices)){ prev.cases.m[i,indices[i,k]-(startcol-1)]=pred.cases.m[i,k]-pred.cases.m[i,k+1]
    }
  }
  prev.cases.m.new <- as.data.frame(prev.cases.m)
   
  if(length(combined_indices)>0){
    if(length(nbreak_p)>0) prev.cases.m.new <- as.data.frame(cbind(prev.cases.m[,1:length(nbreak_p)],matrix(0,nrow(indices),ncol=length(break_p))))
    if(length(nbreak_p)==0) prev.cases.m.new <- as.data.frame(matrix(0,nrow(indices),ncol=length(break_p)))
    
    ord_colnames <- c()
    if(lo>=1) ord_colnames <- paste('o',1:lo,sep="")
    cont_colnames <- c()
    if(lc>=1) cont_colnames <- paste('c',1:lc,sep="")
    colnames(prev.cases.m.new)[1:length(break_p)] <- c(ord_colnames,cont_colnames)
    
    addfactor <-1
    for(k in 1:length(break_p)){
      if(combined_df[k] == 1) prev.cases.m.new[,length(nbreak_p)+k] <- prev.cases.m[,(length(nbreak_p) + addfactor):(length(nbreak_p) + addfactor)]
      if(combined_df[k] > 1) prev.cases.m.new[,length(nbreak_p)+k] <- apply(prev.cases.m[,(length(nbreak_p) + addfactor):(length(nbreak_p) + sum(combined_df[1:k]))],1,sum)
      addfactor <- sum(combined_df[1:k])+1
    }
  }
  prev.cases.m.new <- cbind(prev.cases.m.new,total=apply(prev.cases.m.new,1,sum))
  colnames(prev.cases.m.new)[ncol(prev.cases.m.new)] <- "total" 
  #if(!is.na(approx_error)) print(paste("estimated approximation error is within", 1.96*apply(prev.cases.m.new,2,sd)/(obs.cases*sqrt(nsample_perm))))
  prev.cases=apply(prev.cases.m.new,2,mean)
  return(prev.cases/obs.cases)
}

affun_ci <- function(the.form=NULL,the.data=NULL,nsample_perm = 10000, approx_error=NA, nsample_var = 100, prev=0.015, conf_level=0.99, sep_est=FALSE, correction_factor=TRUE,allperm=TRUE, quantile_est=TRUE){
  if(!is.na(approx_error)) sep_est=TRUE
  d=model.frame(the.form,the.data)
  
  chunksize = nsample_perm/nsample_var
  y <- d[,1]
  w <- rep(1,length(y))
  if(!is.na(prev)){ 
    
    thenas <- apply(d[,2:ncol(d)],1,function(x){sum(is.na(x))>0})
    w[y==0] <- (sum(y==1 & !thenas)/sum(y==0 & !thenas))*((1-prev)/prev)
    w[y==1] <-  1  
  }
  m=the.form
  
  m$coefficients <- m$coefficients[!is.na(m$coefficients)]
  mean_beta <- m$coefficients
  cov_beta <- summary(m)$cov.unscaled
  ## remove weight column
  if(colnames(d)[ncol(d)]=="(weights)") d <- d[,1:(ncol(d)-1)]
  
  if(!is.na(prev)){ 
    thenas <- apply(d[,2:ncol(d)],1,function(x){sum(is.na(x))>0})
    newd <- cbind(1,d[!thenas,2:ncol(d)])
    newd <- as.matrix(newd)
    #cov_beta <- solve(t(newd)%*%solve(diag(w[!thenas]))%*%diag(m$weights)%*%newd)
    newm <- glm(d[,1]~as.matrix(d[,2:ncol(d)]),family='binomial')
    cov_beta <-  summary(newm)$cov.unscaled
  }
  ##  remove columns that are 0
  d <- d[,apply(d,2,function(x){sum(abs(x))})!=0]
  
  
  startcol <- grep("(c|o)[0-9]_",colnames(d),perl=TRUE)[1]
  
  vars=length(unique(substring(colnames(d)[startcol:ncol(d)],1,3)))
  n=nrow(d)
  obs.cases=sum(d[,1])
  continuous_indices <-  (startcol-1) + grep("c[0-9]+_*",unique(substring(colnames(d)[startcol:ncol(d)],1,3)),perl=TRUE)
  lc <- length(continuous_indices)
  ordinal_indices <-  (startcol-1) + grep("o[0-9]+_*",unique(substring(colnames(d)[startcol:ncol(d)],1,3)),perl=TRUE)
  lo <- length(ordinal_indices)
  combined_indices <- c(ordinal_indices, continuous_indices)
  df_ordinal <- numeric(length(ordinal_indices))
  if(length(df_ordinal) > 0){
    for(j in 1:length(df_ordinal)) df_ordinal[j] <- sum(grep(paste("o",j,"_",sep=""),substring(colnames(d),1,nchar(paste("o",j,"_",sep=""))),perl=TRUE)>0)
  }
  combined_df <- c(df_ordinal)
  combined_df <- combined_df[combined_df > 0] ###  remove 0's
  
  
  if(allperm){
    PARF <- AF_exact(the.form=the.form,the.data=the.data,prev=prev)
    the.ests <- matrix(0,nrow=nsample_var,ncol=vars+1) 
    beta_to_use <- MASS::mvrnorm(n = nsample_var, mu=mean_beta , Sigma=cov_beta, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    for(l in 1:nsample_var){
      m$coefficients <- beta_to_use[l,]
      d <- d[,apply(d,2,function(x){sum(abs(x))})!=0]
      if(colnames(d)[ncol(d)]=="(weights)") d <- d[,1:(ncol(d)-1)]
      
      
      startcol <- grep("(c|o)[0-9]_",colnames(d),perl=TRUE)[1]
      
      vars=length(unique(substring(colnames(d)[startcol:ncol(d)],1,3)))
      n=nrow(d)
      obs.cases=sum(d[,1])
      continuous_indices <-  (startcol-1) + grep("c[0-9]+_*",unique(substring(colnames(d)[startcol:ncol(d)],1,3)),perl=TRUE)
      lc <- length(continuous_indices)
      ordinal_indices <-  (startcol-1) + grep("o[0-9]+_*",unique(substring(colnames(d)[startcol:ncol(d)],1,3)),perl=TRUE)
      lo <- length(ordinal_indices)
      combined_indices <- c(ordinal_indices, continuous_indices)
      df_ordinal <- numeric(length(ordinal_indices))
      if(length(df_ordinal) > 0){
        for(j in 1:length(df_ordinal)) df_ordinal[j] <- sum(grep(paste("o",j,"_",sep=""),substring(colnames(d),1,nchar(paste("o",j,"_",sep=""))),perl=TRUE)>0)
      }
      combined_df <- c(df_ordinal)
      combined_df <- combined_df[combined_df > 0] ###  remove 0's
      AFvec <- length(vars)
      minvalues <- numeric(ncol(d))
      minvalues[1] <- NA
      addfactor <- 0
      
      if(length(df_ordinal)>0) addfactor <- sum(df_ordinal-1)
      AFvec <- numeric(vars)
      for(j in 1:vars){
        theframe <- create_frame(vars-1)
        SAF <- as.numeric(nrow(theframe)) ## this stores the sequential estimates
        the.w <- theframe[,1]
        tomin <- matrix(as.logical(theframe[,2:ncol(theframe)]),nrow=nrow(theframe),byrow=FALSE)
         tempdf <- combined_df[-j]
        tomin1 <- matrix(0,nrow(tomin),ncol=0)
        tomin2 <- matrix(0,nrow(tomin),ncol=0)
        count <- 1
          the.col <- 1
        
        for(i in 1:vars){
          if(i==j){
            tomin1 <- cbind(tomin1,matrix(FALSE,nrow=nrow(tomin),ncol=combined_df[i]))
            tomin2 <- cbind(tomin2,matrix(TRUE,nrow=nrow(tomin),ncol=combined_df[i]))
          }
          if(i != j){
            tomin1 <- cbind(tomin1,matrix(rep(tomin[,the.col],combined_df[i]),nrow=nrow(tomin),ncol=combined_df[i],byrow=FALSE))
            tomin2 <- cbind(tomin2,matrix(rep(tomin[,the.col],combined_df[i]),nrow=nrow(tomin),ncol=combined_df[i],byrow=FALSE))
            the.col <- the.col + 1
          }
        }
        
        for(k in 1:nrow(tomin)){
          d1 <- d
          coltochange1 <- (startcol:(startcol+sum(combined_df)-1))[as.logical(tomin1[k,])]
          d2 <- d
          coltochange2 <- (startcol:(startcol+sum(combined_df)-1))[as.logical(tomin2[k,])]
          d1[,coltochange1]=matrix(rep(minvalues[coltochange1],nrow(d)),nrow=nrow(d),byrow=TRUE)
          d2[,coltochange2]=matrix(rep(minvalues[coltochange2],nrow(d)),nrow=nrow(d),byrow=TRUE)
          eta1 <- as.matrix(cbind(1,d1[,2:ncol(d1)]),nrow=nrow(d2))%*%as.vector(m$coefficients)
          eta1 <- exp(eta1)
          eta2 <- as.matrix(cbind(1,d2[,2:ncol(d2)]),nrow=nrow(d2))%*%as.vector(m$coefficients)
          eta2 <- exp(eta2)
          pred.cases1 <- sum(w*eta1/(1+eta1))  ###  add the weights back in here
          pred.cases2 <- sum(w*eta2/(1+eta2))  ###  add the weights back in here
          SAF[k] <- (pred.cases1 - pred.cases2)/obs.cases
        }
        AFvec[j] <- weighted.mean(SAF,w=the.w)
      }
      the.ests[l,] <- c(AFvec,sum(AFvec))  
    }
    the.sd <- apply(the.ests,2,sd)
    the.df <- nsample_var - 1
    critval <- qt(1-(1-conf_level)/2,the.df)
    if(quantile_est) return(rbind(PARF,apply(the.ests,2,function(x){quantile(x,(1-conf_level)/2)}),apply(the.ests,2,function(x){quantile(x,1-(1-conf_level)/2)})))
    return(rbind(PARF,PARF-critval*the.sd,PARF + critval*the.sd))
    
    
  }  
  ### in the case of not using all the permutations - just sample some of them
  if(!allperm){
    
    if(!is.na(approx_error))
    {
      ## test
      est_sd <- affun_temp(the.form,the.data=the.data,nsample_perm = 100,  prev=prev, allperm=allperm)
      nsample_perm <- max(ceiling(est_sd^2*1.96^2/approx_error^2))+1 ## so there will be at least 2
      nsample_perm <- max(nsample_perm,nsample_var*2) ##  so a variance can be calculated
      print(paste(nsample_perm, " permutations will be performed",sep=""))
      chunksize = nsample_perm/nsample_var
    }
    
    indices <- matrix(0,nrow=nsample_perm,ncol=vars)
    for(i in 1:nsample_perm){
      indices[i,] <- sample(startcol:(vars+startcol-1),vars,replace=FALSE)
    }
    
    
      
    if(sum(combined_df) > 0){
      the.dim <- ncol(indices)+sum(combined_df-1)
      indices_w <- matrix(0,nrow(indices),the.dim)
      
      
      
      for(i in 1:nrow(indices)){
        o_vec <- indices[i,]
        new_vec <- numeric(the.dim)
        break_p <- (1:length(o_vec))[o_vec %in% combined_indices]
        nbreak_p <- setdiff((1:length(o_vec)),break_p)
        
         the_order <- order(o_vec[break_p])  #what is the position of the breakpoint for x1,the breakpoint for x2 ....
        combined_df_ordered <- combined_df[rank(o_vec[break_p])]  #which df comes first, which second and which last
        break_p <- break_p[order(o_vec[break_p])]
        count_col <- 1
        addfactor <- 0
        
        for(k in 1:length(break_p)){
          
            if(the_order[k]>1) start_i <- sum(combined_df_ordered[1:(the_order[k]-1)]-1)+break_p[k]
          if(the_order[k]==1) start_i <- break_p[k]
          end_i <- sum(combined_df_ordered[1:the_order[k]]-1)+break_p[k]
          new_vec[start_i:end_i] <- (o_vec[break_p[k]]+addfactor):(o_vec[break_p[k]]+sum(combined_df[1:k] -1))
          addfactor <- sum(combined_df[1:k]-1)
        }
        new_vec[new_vec == 0] <- o_vec[nbreak_p]
        
        if((100*i/nrow(indices))%%10 == 0){flush.console()
                                           print(paste('permutations ', 100*i/nrow(indices),'% done',sep=''))
        }
        indices_w[i,] <- new_vec
      }
      
      indices <- indices_w
    }
    
    
    
    beta_to_use <- MASS::mvrnorm(n = nsample_perm, mu=mean_beta , Sigma=cov_beta, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    ####
    d1=d
    ##  an error if indices matrix loses its dimensionality
    if(is.null(nrow(indices))) indices <- as.matrix(indices,nrow=length(indices),ncol=1)
    pred.cases.m=matrix(NA,nrow=nrow(indices),ncol=ncol(indices))
    prev.cases.m=matrix(NA,nrow=nrow(indices),ncol=ncol(indices))
    
    ### set minimum values for all the variables
    
    for(i in 1:nrow(indices)){
      ## Keep the same coefficients when you are within a chunk used to calculate a variance
      if((i %%  chunksize) == 1) m$coefficients <- beta_to_use[i,]
      minvalues <- numeric(ncol(d))
      minvalues[1] <- NA
      addfactor <- 0
      
      if(length(df_ordinal)>0) addfactor <- sum(df_ordinal-1)
      
      for(k in 1:ncol(indices)){
        d[,indices[i,k]]=rep(minvalues[indices[i,k]],nrow(d))
        
        eta <- as.matrix(cbind(1,d[,2:ncol(d)]),nrow=nrow(d))%*%as.vector(m$coefficients)
        eta <- exp(eta)
        pred.cases.m[i,k]<- sum(w*eta/(1+eta))  ###  add the weights back in here
      }
      d=d1
      if((100*i/nrow(indices))%%10 == 0){flush.console()
                                         print(paste('calculations ', 100*i/nrow(indices),'% done',sep=''))
      }
      
    }
    ###  now collate multiple basis functions corresonding to a single categorical variable
    pred.cases.m=cbind(rep(obs.cases,nrow(indices)),pred.cases.m)
    for(i in 1:nrow(indices)){
      for(k in 1:ncol(indices)){ prev.cases.m[i,indices[i,k]-(startcol-1)]=pred.cases.m[i,k]-pred.cases.m[i,k+1]
      }
    }
    prev.cases.m.new <- as.data.frame(prev.cases.m)
    
    if(length(combined_indices)>0){
      if(length(nbreak_p)>0) prev.cases.m.new <- as.data.frame(cbind(prev.cases.m[,1:length(nbreak_p)],matrix(0,nrow(indices),ncol=length(break_p))))
      if(length(nbreak_p)==0) prev.cases.m.new <- as.data.frame(matrix(0,nrow(indices),ncol=length(break_p)))
      
      ord_colnames <- c()
      if(lo>=1) ord_colnames <- paste('o',1:lo,sep="")
      cont_colnames <- c()
      if(lc>=1) cont_colnames <- paste('c',1:lc,sep="")
      colnames(prev.cases.m.new)[1:length(break_p)] <- c(ord_colnames,cont_colnames)
      addfactor <-1
      for(k in 1:length(break_p)){
        if(combined_df[k] == 1) prev.cases.m.new[,length(nbreak_p)+k] <- prev.cases.m[,(length(nbreak_p) + addfactor):(length(nbreak_p) + addfactor)]
        if(combined_df[k] > 1) prev.cases.m.new[,length(nbreak_p)+k] <- apply(prev.cases.m[,(length(nbreak_p) + addfactor):(length(nbreak_p) + sum(combined_df[1:k]))],1,sum)
        addfactor <- sum(combined_df[1:k])+1
      }
    }
    prev.cases.m.new <- cbind(prev.cases.m.new,total=apply(prev.cases.m.new,1,sum))
    colnames(prev.cases.m.new)[ncol(prev.cases.m.new)] <- "total"
    the.ests <- matrix(0,nrow=nsample_var,ncol=ncol(prev.cases.m.new))
    var.ests <- matrix(0,nrow=nsample_var,ncol=ncol(prev.cases.m.new))
    
    for(i in 1:nsample_var){
      themat <- prev.cases.m.new[(1 +  chunksize*(i-1)):(chunksize*i),]
      if(is.null(nrow(themat))) themat <- as.matrix(themat,ncol=1)
      the.ests[i,] <- apply(themat,2,function(x){mean(x,na.rm=TRUE)})/obs.cases
      var.ests[i,] <- apply(themat,2,function(x){var(x,na.rm=TRUE)})/(obs.cases^2)
    }
    if(!is.na(approx_error)) # print(paste("estimated approximation error is within", 1.96*sqrt(apply(var.ests,2,mean))/(obs.cases*sqrt(nsample_perm))))
      if(!correction_factor | allperm) the.var <- apply(the.ests,2,var)
    if(correction_factor & !(allperm)) the.var <- apply(the.ests,2,var)-apply(var.ests,2,mean)/chunksize + apply(var.ests,2,mean)/nsample_perm    ### dangerous to correct for the bias if chunksize is small
    the.sd <- sqrt(the.var)
    the.mean.est=apply(prev.cases.m.new,2,mean)
    #
    if(sep_est) PARF=affun(the.form,the.data=the.data,nsample_perm = nsample_perm, prev=prev, allperm=allperm, approx_error=NA)
    if(!sep_est) PARF <- the.mean.est/obs.cases
    names(PARF)=colnames(prev.cases.m.new)
    the.df <- nsample_var - 1
    critval <- qt(1-(1-conf_level)/2,the.df)
    if(quantile_est) return(rbind(PARF,apply(the.ests,2,function(x){quantile(x,(1-conf_level)/2)}),apply(the.ests,2,function(x){quantile(x,1-(1-conf_level)/2)})))
    return(rbind(PARF,PARF-critval*the.sd,PARF + critval*the.sd))
    
  }
}

###  sampling deciding which risk factors to put in th positions to the  left and their weights.

create_frame <- function(n){
  if(n==1) return(matrix(c(1,1),nrow=1))
  n.row <- 2^n
  out <- matrix(0,nrow=n.row,ncol=0)
  for(i in 1:n){
    gap <- 2^(n-i)  
    out <- cbind(out,rep(c(rep(1,gap),rep(0,gap)),2^(i-1)))  
  }
  the.sum <- apply(out,1,sum)
  out <- out[order(the.sum),]
  the.sum <- apply(out,1,sum)
  weight.vec <- factorial(0:n)*factorial(n:0)
  the.weights <- weight.vec[the.sum+1] 
  return(cbind(the.weights,out))
}




AF_exact <- function(the.form=NULL,the.data=NULL,prev=0.015){
  
  d=model.frame(the.form,the.data)
  y <- d[,1]
  w <- rep(1,length(y))
  if(!is.na(prev)){ 
    
    thenas <- apply(d[,2:ncol(d)],1,function(x){sum(is.na(x))>0})
    w[y==0] <- (sum(y==1 & !thenas)/sum(y==0 & !thenas))*((1-prev)/prev)
    w[y==1] <-  1  
  }
  m=the.form
  
  m$coefficients <- m$coefficients[!is.na(m$coefficients)]
  d <- d[,apply(d,2,function(x){sum(abs(x))})!=0]
  if(colnames(d)[ncol(d)]=="(weights)") d <- d[,1:(ncol(d)-1)]
  
  
  startcol <- grep("(c|o)[0-9]_",colnames(d),perl=TRUE)[1]
  
  vars=length(unique(substring(colnames(d)[startcol:ncol(d)],1,3)))
  n=nrow(d)
  obs.cases=sum(d[,1])
  continuous_indices <-  (startcol-1) + grep("c[0-9]+_*",unique(substring(colnames(d)[startcol:ncol(d)],1,3)),perl=TRUE)
  lc <- length(continuous_indices)
  ordinal_indices <-  (startcol-1) + grep("o[0-9]+_*",unique(substring(colnames(d)[startcol:ncol(d)],1,3)),perl=TRUE)
  lo <- length(ordinal_indices)
  combined_indices <- c(ordinal_indices, continuous_indices)
  df_ordinal <- numeric(length(ordinal_indices))
  if(length(df_ordinal) > 0){
    for(j in 1:length(df_ordinal)) df_ordinal[j] <- sum(grep(paste("o",j,"_",sep=""),substring(colnames(d),1,nchar(paste("o",j,"_",sep=""))),perl=TRUE)>0)
  }
  combined_df <- c(df_ordinal)
  combined_df <- combined_df[combined_df > 0] ###  remove 0's
  AFvec <- length(vars)
  minvalues <- numeric(ncol(d))
  minvalues[1] <- NA
  addfactor <- 0
  
  if(length(df_ordinal)>0) addfactor <- sum(df_ordinal-1)
  if(length(continuous_indices) > 0){
    for(k in 1:length(continuous_indices)){
      i1 <- continuous_indices[k]+addfactor
      i2 <- continuous_indices[k]+sum(df_ordinal-1)
      the.string <- substring(colnames(d),1,3)[i1]
      mat <- d[,substring(colnames(d),1,3) == the.string ]
      coeff <- m$coefficients[substring( names(m$coefficients),1,3) == the.string ]
      eta <- as.matrix(mat) %*% as.vector(coeff)
      target_row <- sort(eta)[ceiling(length(eta)/100)]
      best_row <- (1:nrow(eta))[eta==target_row][1]
      minvalues[i1:i2] <- as.numeric(d[best_row,i1:i2])
      addfactor <- sum(df_ordinal-1)
    }
  }
  AFvec <- numeric(vars)
  for(j in 1:vars){
    theframe <- create_frame(vars-1)
    SAF <- as.numeric(nrow(theframe)) ## this stores the sequential estimates
    the.w <- theframe[,1]
    tomin <- matrix(as.logical(theframe[,2:ncol(theframe)]),nrow=nrow(theframe),byrow=FALSE)
      tempdf <- combined_df[-j]
    tomin1 <- matrix(0,nrow(tomin),ncol=0)
    tomin2 <- matrix(0,nrow(tomin),ncol=0)
    count <- 1
       the.col <- 1
    
    for(i in 1:vars){
      if(i==j){
        tomin1 <- cbind(tomin1,matrix(FALSE,nrow=nrow(tomin),ncol=combined_df[i]))
        tomin2 <- cbind(tomin2,matrix(TRUE,nrow=nrow(tomin),ncol=combined_df[i]))
      }
      if(i != j){
        tomin1 <- cbind(tomin1,matrix(rep(tomin[,the.col],combined_df[i]),nrow=nrow(tomin),ncol=combined_df[i],byrow=FALSE))
        tomin2 <- cbind(tomin2,matrix(rep(tomin[,the.col],combined_df[i]),nrow=nrow(tomin),ncol=combined_df[i],byrow=FALSE))
        the.col <- the.col + 1
      }
    }
    
    for(k in 1:nrow(tomin)){
      d1 <- d
      coltochange1 <- (startcol:(startcol+sum(combined_df)-1))[as.logical(tomin1[k,])]
      d2 <- d
      coltochange2 <- (startcol:(startcol+sum(combined_df)-1))[as.logical(tomin2[k,])]
      d1[,coltochange1]=matrix(rep(minvalues[coltochange1],nrow(d)),nrow=nrow(d),byrow=TRUE)
      d2[,coltochange2]=matrix(rep(minvalues[coltochange2],nrow(d)),nrow=nrow(d),byrow=TRUE)
      eta1 <- as.matrix(cbind(1,d1[,2:ncol(d1)]),nrow=nrow(d2))%*%as.vector(m$coefficients)
      eta1 <- exp(eta1)
      eta2 <- as.matrix(cbind(1,d2[,2:ncol(d2)]),nrow=nrow(d2))%*%as.vector(m$coefficients)
      eta2 <- exp(eta2)
      pred.cases1 <- sum(w*eta1/(1+eta1))  ###  add the weights back in here
      pred.cases2 <- sum(w*eta2/(1+eta2))  ###  add the weights back in here
      SAF[k] <- (pred.cases1 - pred.cases2)/obs.cases
    }
    AFvec[j] <- weighted.mean(SAF,w=the.w)
  }
  return(c(AFvec,sum(AFvec)))
}

## this function is not used anymore - better to use AF_exact
perm=function(from,to,vec){
  if (to == 1) return(matrix(vec,from,1))
  else if (from==1) matrix(vec,1,to)
  else{
    X=NULL
    for(i in 1:from){
      X=rbind(X,cbind(vec[i],Recall(from-1,to-1,vec[-i])))
    }
    return(X)
  }
}


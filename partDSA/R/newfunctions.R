add.bas.fx <- function(BFs, new.add, m) {
  ## add a bas fx to the current basis function configuration
  first.BF <- BFs[[m]]
  var.num <- new.add[[5]]
  k.num <- new.add[[6]]
  second.BF <- lapply(first.BF, function(x) x[k.num])
  first.BF[[var.num]][[k.num]][2] <- new.add[[1]]
  second.BF[[var.num]][[1]][1] <- new.add[[1]]
  BFs[[m]] <- first.BF
  BFs[[length(BFs) + 1]] <- second.BF
  return(BFs)
}

assess.risk.on.add <- function(poss.add, BF.list, dat, y, wt, opts,x.temp,is.num,missing) {
  gc()
  
  ## from 'add.bas.fx' find best addition for each bas.fx
  ## now see which of M additions minimizes the RSS
  add.risk <- NULL
  new.bf <- vector("list", length=length(poss.add))

  for (M in 1:length(poss.add)) {
    ## check if we can do the split
    if (!poss.add[[M]][[4]]) {
      new.bf[[M]] <- add.bas.fx(BFs=BF.list, new.add=poss.add[[M]], m=M)
      ###update x.temp based on this proposed addition move
      dat2 <- update.missing(x.temp=x.temp,x=dat,y=y,bas=cbind(poss.add[[M]]$new.lt,poss.add[[M]]$new.gt),is.num=is.num,       todo="a",missing=missing)
      bas.fx <- assign.obs.to.bf(dat2=dat2, n=nrow(dat), p=ncol(dat),
                                 BFs=new.bf[[M]])
      add.risk[M] <- risk.fx(bas.fx=bas.fx, y=y, wt=wt, opts=opts)
      
    }
    else {
      add.risk[M] <- NA
    }
  }

  ind.a <- which.min(add.risk)[1]

  ### added
  if(sum(!is.na(add.risk))>0){
  ###update x.temp to reflect the chosen best addition move
      dat2 <- update.missing(x.temp=x.temp,x=dat,y=y,bas=cbind(poss.add[[ind.a]]
               $new.lt,poss.add[[ind.a]]$new.gt),is.num=is.num, todo="a",missing=missing)
  } else{
        dat2=x.temp
  }

  list(new.bf=new.bf[[ind.a]],
       add.risk=if(is.na(ind.a)) add.risk[1] else add.risk[ind.a],x.temp=dat2)
}

newLIST <- function(p) {
  return(list(rep(list(list(c(-Inf, Inf))), p)))
}

assign.obs.to.bf <- function(dat2, n, p, BFs) {
  ## p = num of variables
  ## n = num of obs
  ## dat2 = matrix with variables and corresponding values (nxp)
 
  fun <- function(BF) {
    or <- integer(n)
    for(K in 1:length(BF[[1]])) {
      add <- integer(n)
      for(P in 1:p) {
        x <- as.numeric(dat2[,P])
        add <- add + (BF[[P]][[K]][1] < x & x <= BF[[P]][[K]][2])
      }
      or <- or + (add == p) 
    }
    as.numeric(or)  # XXX for compatibility
  }
  names(BFs) <- paste('X', seq(along=BFs), sep='')
  do.call('data.frame', lapply(BFs, fun))
}


risk.fx <- function(bas.fx, y, wt, opts){
  if(opts$outcome.class=="factor"){
    get.risk<-categorical.risk(bas.fx=bas.fx,y=y,wt=wt,opts=opts)
  } else if (opts$outcome.class=="survival"){
    get.risk <- survival.overall.risk(bas.fx=bas.fx, y=y, wt=wt, opts=opts)
  } else get.risk<-rss.risk(bas.fx,y,wt)
  return(get.risk)
}


rss.risk <- function(bas.fx, y, wt) {
  n <- length(y)
  get.loss <- NULL
  lm.xy <- lm(y~., data=bas.fx, weights=wt)
  coef <- ifelse(is.na(lm.xy$coefficients), 0, lm.xy$coefficients)
  new <- as.matrix(cbind(1, bas.fx))
  pred.val <- new %*% coef
  get.loss  <- sum(wt * (y - pred.val) ^ 2)/sum(wt)
  return(get.loss)
}

num.x <- function(fv, lv, var.x) {
  return(sum(var.x > fv & var.x < lv))
}

categorical.risk <- function(bas.fx, y, wt, opts) {
  sum.impurities = 0
  for(i in 1:dim(bas.fx)[2]) {
    filter = bas.fx[,i]
    outcome.bucket = y[filter == 1]
    sum.impurities = sum.impurities + (length(outcome.bucket) / length(y)) *
        gini.impurity(outcome.bucket, opts)
  }

  return(sum.impurities)
}

gini.impurity <- function(y, opts) {
  p = rle(sort(as.integer(y)))$lengths / length(y)
  if (opts$loss.fx == "gini")         # gini
    1.0 - sum(p * p)
  else                                # entropy
    sum(-p * log(p))
}

get.bfs <-function(BFs, dat) {
  cbind(INDEX=rep(1, nrow(dat)),
        assign.obs.to.bf(dat2=dat, n=nrow(dat), p=ncol(dat), BFs=BFs))
}

get.coefs <- function(BFs, dat, y, wt) {
  outs <- assign.obs.to.bf(dat2=dat, n=nrow(dat), p=ncol(dat), BFs=BFs)
  lm.xy <- lm(y~., data=outs, weights=wt)
  coef <- ifelse(is.na(lm.xy$coefficients), 0, lm.xy$coefficients)
  r2 <- summary(lm.xy)$r.squared
  return(list(coef=coef, r2=r2))
}

get.votes <- function(BFs, dat, y, wt) {
  outs <- assign.obs.to.bf(dat2=dat, n=nrow(dat), p=ncol(dat), BFs=BFs)
  part.pred <- numeric(ncol(outs))
  for (j in 1:ncol(outs)) {
    ty <- table(y[outs[,j] == 1])
    d <- dimnames(ty)[[1]][which.max(ty)]
    part.pred[j] <- match(d, levels(y))
  }

  list(vote=part.pred)
}


##This function takes the best.ind which corresponds to the maximum goodness
## value. If there is only one such value, this same input value will be
## returned. If there are multiple maximum goodness values, we choose the one 
##that is most central. Therefore this function is called every time we find 
##best.ind.

central.split <- function(best.ind,n){
 
  ## This line finds the absolute value of the difference between each best.ind   ##and the central value (n/2) and takes that best.ind value which has the
  ## minimum distance to the center.
  
  return(best.ind[which.min(abs(best.ind-n/2))[1]])
}


## x.temp is an imputed version of x with NA values filled in based on the average of that variable within that observation's node. X represents the original X matrix so that we can determine the original placement of NA values. Then we update the x.temp values corresponding to the NA values based on the new basis functions, bas. is.num is a binary vector of length equal to the number variables with a 1 representing a continuous variable and 0 representing a categorical variable. todo is equal to a when we are creating two new basis functions and it equals d when we are creating (combining) one new basis function
update.missing <- function (x.temp,x, y, bas,is.num, todo, missing){
  
  if(missing=="no"){
    
    return(x)
  }else if(missing=="impute.at.split"){
  
  ### this is the case when we have two new basis function so dim(bas)[2]=2
  if (todo=="a"){
  
    # when we have a categorical outcome, the imputed value should be based on 
    # those observations in the node with the same outcome value
    if (is.factor(y)){
    
      # work on each variable individually
      for (i in 1:dim(x.temp)[2]){
 
        #get position of all na values for variable i
        na.value <- is.na(x[,i])

        #for each outcome level of y
        for (j in 1:length(levels(y))){
           
          #make sure there is at least one observation in
          # this basis function and level of y
          if(length((x[bas[,1] & y==levels(y)[j],i]))>0){
 
               
                #all observations in this node are missing on this 
                #variable so get average from this variale overall
          	if (length(which(!is.na(x[bas[,1]==1,i])))==0){
                	replace1 <- ifelse(is.num[i]==1,mean(x[,i],na.rm=TRUE),
                	as.numeric(names(which.max(table(x[,i])))))
               
               		 warning(paste("All values of variable",names(x.temp)[i],"are NA in node. Using overall average of this variable."))

          	#all observations in this node for this level are missing
                # on this variables so get average from node
                }else if (length(which(!is.na(x[bas[,1]==1&
                y==levels(y)[j],i])))==0){
			
                        replace1 <- (ifelse(is.num[i]==1,mean(x[bas[,1]==1,i],na.rm=TRUE),as.numeric(names(which.max(table(x[bas[,1]==1,i]))))))
                	warning(paste("All values of variable",names(x.temp)[i],"for level",levels(y)[j], "are NA in this node. Using overall average in the node"))

                #otherwise:no warnings messages, just get the imputed value from the mean or majority vote of those observations in proper level and proper node
          	}else{
                	replace1 <-  ifelse(is.num[i]==1,
                	mean(x[bas[,1]==1& y==levels(y)    [j],i],na.rm=TRUE),
                	as.numeric(names(which.max(table(x[bas[,1]==1&y==levels(y)[j],i])))))
         	}
         
                
                #update x.temp matrix based on above

          	x.temp[na.value & bas[,1]==1 & y==levels(y)[j],i] <- replace1
          }


          #do the exact same thing as above but this time with basis function 2
          if(length(x[bas[,2]==1 & y==levels(y)[j],i])>0){
          
		if (length(which(!is.na(x[bas[,2]==1,i])))==0){
                	replace2 <- ifelse(is.num[i]==1,mean(x[,i],na.rm=TRUE),
                	as.numeric(names(which.max(table(x[,i])))))
                	warning(paste("All values of variable",names(x.temp)[i],"are NA in node. Using overall average of this variable. "))
          
		}else if (length(which(!is.na(x[bas[,2]==1&y==levels(y)[j],i])))==0){
                
			replace2 <- ifelse(is.num[i]==1,mean(x[bas[,2]==1,i],na.rm=TRUE),
                	as.numeric(names(which.max(table(x[bas[,2]==1,i])))))
               		warning(paste("All values of variable",names(x.temp)[i],"for level",levels(y)[j], "are NA in this node. Using overall average in the node."))
          	}else{
                	replace2 <-  ifelse(is.num[i]==1,
                	mean(x[bas[,2]==1& y==levels(y)[j],i],na.rm=TRUE),
                	as.numeric(names(which.max(table(x[bas[,2]==1&y==levels(y)[j],i])))))
         	}
         
   
        
          	#do the same for the second basis function
          	x.temp[na.value & bas[,2]==1 & y==levels(y)[j],i] <- replace2
       
        }
      }
      }
  
    # this is the case for continuous outcomes, just take an mean or mode of that variable in the node
    }else{
      

      #for each variable
      for (i in 1:dim(x.temp)[2]){
        
      	#find location of  na values
        na.value <- is.na(x[,i])

        #if all observations in the node are NA for this variale,
        # take overall average
        if (length(which(!is.na(x[bas[,1]==1,i])))==0){
                replace1 <- ifelse(is.num[i]==1,mean(x[,i],na.rm=TRUE),
                as.numeric(names(which.max(table(x[,i])))))
                warning(paste("All values of variable",names(x.temp)[i],
                "are NA in node. Using overall average of this variable."))
        #otherwise take average within this node
        }else{
                replace1 <- ifelse (is.num[i]==1,
                mean(x[bas[,1]==1,i],na.rm=TRUE),
                as.numeric(names(which.max(table(x[bas[,1]==1,i])))))
        }
     
        #do the same as above, except don't worry about outcome class of y
        x.temp[na.value & bas[,1]==1,i] <- replace1

        #if all observaitons in the node are NA for this variable,
        # take overall average
	if (length(which(!is.na(x[bas[,2]==1,i])))==0){
                replace2 <- ifelse(is.num[i]==1,mean(x[,i],na.rm=TRUE),
                as.numeric(names(which.max(table(x[,i])))))
                warning(paste("All values of variable",names(x.temp)[i],
                "are NA in node. Using overall average of this variable."))
        #otherwise take average within this node
        }else{
                replace2 <- ifelse (is.num[i]==1,
                mean(x[bas[,2]==1,i],na.rm=TRUE),
                as.numeric(names(which.max(table(x[bas[,2]==1,i])))))
        }
        
      
        x.temp[na.value & bas[,2]==1,i] <- replace2
    }

    }
  }  


  #Now we only have one basis function, so bas is a vector.
  if (todo=="d"){
    
    # in the case of a categorical outcome, update based on the observations 
    # in that node with the same outcome class
    if (is.factor(y)){
    
      #for each variable
      for (i in 1:dim(x.temp)[2]){
        
        #get which observations are na on this variable
        na.value <- is.na(x[,i])
   
      
        
        for (j in 1:length(levels(y))){
 
          if(length((x[bas==1&y==levels(y)[j],i]))>0){

          	if (length(which(!is.na(x[bas==1,i])))==0){
                	replace <- ifelse(is.num[i]==1,mean(x[,i],na.rm=TRUE),
                	as.numeric(names(which.max(table(x[,i])))))
                        warning(paste("All values of variable",names(x.temp)[i],"are NA in node. Using overall average of this variable."))
          	}else if (length(which(!is.na(x[bas==1&y==levels(y)[j],i])))==0){
			replace <- ifelse(is.num[i]==1,mean(x[bas==1,i],na.rm=TRUE),as.numeric(names(which.max(table(x[bas==1,i])))))
                	warning(paste("All values of variable",names(x.temp)[i],"for level",levels(y)[j], "are NA in each node.
                         Using overall average in the node."))
          	}else{
               		replace <-  ifelse(is.num[i]==1,
                	mean(x[bas==1& y==levels(y)[j],i],na.rm=TRUE),
                	as.numeric(names(which.max(table(x[bas==1&y==levels(y)[j],i])))))       

         	}
          
          
          	#depending on whether variable is continuous or categorical,
          	# update with the mean or mode
         
        
          	x.temp[na.value & bas==1 & y==levels(y)[j],i] <- replace
       
     
           }
	}
      
      }


    }else{
     #in the case of a continuous outcome
      for (i in 1:dim(x.temp)[2]){
        na.value <- is.na(x[,i])
        if (length(which(!is.na(x[bas==1,i])))==0){
                replace <- ifelse(is.num[i]==1,mean(x[,i],na.rm=TRUE),
                as.numeric(names(which.max(table(x[,i])))))
                warning(paste("All values of variable",names(x.temp)[i],"are NA in node. Using overall average of this variable."))
        }else{
                replace <- ifelse (is.num[i]==1,mean(x[bas==1,i],na.rm=TRUE),
                as.numeric(names(which.max(table(x[bas==1,i])))))
        }
        x.temp[na.value& bas==1,i] <- replace
      }
    }

  }
  
  


  return(x.temp)
  }
}



impute.test <- function(x,y,x.test,missing){

 if (missing=="no"){
    return(x.test)
 }
 else if (missing=="impute.at.split"){
 ## impute the test matrix based on means and mode from training set
  for (k in 1:dim(x.test)[2]){
 
    #for categorical outcome, update based on the outcome class
    # as well as the observation's node
  #  if(is.factor(y)){
  #   for(i in levels(y)){
  #      x.test[(is.na(x.test[,k])&y.test==i),k] <- ifelse(is.factor(x[,k]),
  #      names(which.max(table(x[y==i,k]))),mean(x[y==i,k],na.rm=TRUE))
  #      }
    
   # }else {   
      x.test[(is.na(x.test[,k])),k] <- ifelse(is.factor(x[,k]),
      names(which.max(table(x[,k]))),mean(x[,k],na.rm=TRUE))
   # }
  }

 return(x.test)
 }
}

#This function returns a list with the attack and defend coefficients plus any added on predictors
#length of attack and defend coef must be the same
coef.fbRanks = function(object, ...){
  x=object
  the.fits=x$fit
  
  coef.vector=coef.list=list()
  
  for(clus in 1:length(the.fits)){
    fit=the.fits[[clus]]
    team.names = fit$xlevels$attack
    nteams.in.fit = length(team.names)
    
    if(!any(c("glm","speedglm","glmnet") %in% class(fit))){
      stop("fit needs to be of class glm, speedglm or glmnet")      
    }
           
    if("glmnet" %in% class(fit)){
      coef.vals = coef(fit)[-1,]
    }
    if("speedglm" %in% class(fit) | "glm" %in% class(fit)){
      #speedglm coef is sorted alpha by name
      #starts with attack+team name for all teams
      #next are defend+team name except first team is missing (set to 0)
      #the first defense name is always set to 0
      coef.vals = c(fit$coef[1:nteams.in.fit],0,fit$coef[(nteams.in.fit+1):(nteams.in.fit*2-1)])
      coef.names = c(paste("attack",team.names,sep=""),paste("defend",team.names,sep=""))
      add.col = nteams.in.fit*2 #START of next set of coef
      for(i in attr(fit$terms,"term.labels")[c(-1,-2)]){
        if(attr(fit$terms,"dataClasses")[[i]]=="factor"){ #factor
          if(length(fit$xlevels[[i]])>1){
            #first factor always set to 0
            coef.vals = c(coef.vals,0,fit$coef[add.col:(add.col+length(fit$xlevels[[i]])-2)])         
            add.col = add.col+length(fit$xlevels[[i]])-1 #START of next set of coef
          }else{
            coef.vals = c(coef.vals,0)           
          }
          coef.names = c(coef.names,paste(fit$xlevels[[i]],i,sep=""))
        }else{ #numeric
          coef.vals = c(coef.vals,fit$coef[add.col:add.col])
          add.col = add.col+1
          coef.names = c(coef.names,i)
        }      
      }
      names(coef.vals)=coef.names
    }
coefs=list()
    coefs$attack = coef.vals[1:nteams.in.fit]
      names(coefs$attack)=team.names
      #the first defense name is always set to 0
      coefs$defend  = coef.vals[(nteams.in.fit+1):(nteams.in.fit*2)]
      names(coefs$defend)=team.names
      add.col = nteams.in.fit*2+1 #START of next set of coef
      for(i in attr(fit$terms,"term.labels")[c(-1,-2)]){
        if(attr(fit$terms,"dataClasses")[[i]]=="factor"){ #factor
           #first factor always set to 0
           coefs[[i]] = coef.vals[add.col:(add.col+length(fit$xlevels[[i]])-1)]         
           add.col = add.col+length(fit$xlevels[[i]]) #START of next set of coef
          names(coefs[[i]])=fit$xlevels[[i]]
        }else{ #numeric
          coefs[[i]] = coef.vals[add.col:add.col] 
          add.col = add.col+1
        }      
    }
    coef.vector[[paste("cluster.",clus,sep="")]]=coef.vals
    coef.list[[paste("cluster.",clus,sep="")]]=coefs
  }
 return(list(coef.vector=coef.vector, coef.list=coef.list))
}
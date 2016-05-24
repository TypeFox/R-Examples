tree_Rasch <-
function(y,
                       DM_kov,
                       npersons,
                       nitems,
                       nvar,
                       ordered_values,
                       n_levels,
                       n_s,
                       alpha,
                       nperm,
                       trace,
                       penalize
  
){
  # design of rasch model 
  pp_design <- diag(npersons)                   # persons, person P reference 
  pp_design <- pp_design[rep(1:nrow(pp_design),each=nitems),]
  pp_design <- pp_design[,-npersons]
  
  ip_design <- -1*diag(nitems)                  # item parameter  
  ip_design <- ip_design[rep(1:nrow(ip_design),times=npersons),]
  
  dm_rasch <- cbind(pp_design,ip_design)               
  
  names_rasch    <- c(paste("theta",1:(npersons-1),sep=""),paste("beta",1:nitems,sep=""))
  colnames(dm_rasch) <- names_rasch      
  
  # functions to build design 
  thresholds <- lapply(1:nvar, function(j) ordered_values[[j]][-length(ordered_values[[j]])])
  
  v <- lapply(1:nvar,function(j) 1:(n_levels[j]-1)) 
  w <- lapply(1:nvar, function(j) rep(paste0("s",j),n_s[j]))
  
  design_one  <- function(x,threshold,upper){
    if(upper){
      ret <- ifelse(x > threshold,1,0)
    } else{
      ret <- ifelse(x > threshold,0,1)
    }
    return(ret)
  }
  
  design <- function(x,thresholds,upper){
    ret <- sapply(thresholds, function(j) design_one(x,j,upper))
    return(ret)
  }
  
  whole_design <- function(X,var,item,thresholds,upper=TRUE){
    design_tree <- matrix(0,nrow=nitems*npersons,ncol=length(thresholds[[var]]))
    rows        <- seq(item,(nitems*npersons),by=nitems)
    design_tree[rows,] <- design(X[,var],thresholds[[var]],upper)
    z <- rep(paste0("_u",item),length(thresholds[[var]]))
    colnames(design_tree) <- paste0(w[[var]],v[[var]],z)
    return(design_tree)
  }
  
  designlists <- function(X,thresholds,upper=TRUE){
    ret <- lapply(1:nitems, function(j){
      lapply(1:nvar, function(var){
        whole_design(X,var,j,thresholds,upper)
      })
    })
    
    return(ret)
  }
  
  #########################################################################################
  
  # initializations 
  mod_potential <- list()
  devs          <- c()
  crits         <- c()
  splits        <- c()
  pvalues       <- c() 
  ip            <- list()
  vars_evtl     <- list()
  splits_evtl   <- list()
  which_obs     <- list()
  numbers       <- list()
  count <- 1
  
  pp     <- paste("theta",1:(npersons-1),sep="")
  help_p <- paste0(pp,collapse="+")
  
  numbers[[1]]     <- lapply(1:nitems,function(j) 1)
  which_obs[[1]]   <- lapply(1:nitems,function(j) matrix(1:npersons,nrow=1))
  splits_evtl[[1]] <- lapply(1:nitems,function(j) lapply(1:nvar, function(var) matrix(1:n_s[var],nrow=1)))
  vars_evtl[[1]]   <- lapply(1:nitems,function(j) nvar)
  ip[[1]]          <- lapply(1:nitems,function(j) paste0("beta",j))
  help0            <- formula(paste("y~",help_p,"+",paste0(unlist(ip[[1]]),collapse="+"),"-1"))
  dat0             <- data.frame(y,dm_rasch)
  mod0             <- mod_potential[[1]] <- glm(help0,family=binomial(link="logit"),data=dat0)
  start            <- predict(mod0)
  
  design_upper <- designlists(DM_kov,thresholds)
  design_lower <- designlists(DM_kov,thresholds,upper=FALSE)
  sig   <- TRUE
  anysplit <- TRUE
  
  # function to compute all models in one knot
  allmodels <- function(i,var,kn,design_lower,design_upper){
    
    deviances <- rep(0,n_s[var])
    help_kn <- ip[[count]][[i]][kn]
    help1   <- paste0(unlist(ip[[count]])[-which(unlist(ip[[count]])==help_kn)],collapse="+")
    splits_aktuell <- splits_evtl[[count]][[i]][[var]][kn,] 
    splits_aktuell <- splits_aktuell[!is.na(splits_aktuell)]
    
    param_knot <- coef(mod0)[help_kn]
    obs_aktuell <- which_obs[[count]][[i]][kn,]
    obs_aktuell <- obs_aktuell[!is.na(obs_aktuell)]
    rows_obs   <- seq(i,(nitems*npersons),by=nitems)[obs_aktuell]
    start[rows_obs] <- start[rows_obs]+param_knot
    
    if(length(splits_aktuell)>0){
      
      for(j in splits_aktuell){
        dat   <- data.frame(dat0,design_lower[[i]][[var]][,j,drop=FALSE],design_upper[[i]][[var]][,j,drop=FALSE])
        help2 <- paste(ip[[count]][[i]][kn],c(colnames(design_lower[[i]][[var]])[j],colnames(design_upper[[i]][[var]])[j]),sep=":")
        help3 <- paste(help2,collapse="+")
        help4 <- formula(paste("y~",help3,"-1"))
        suppressWarnings(
          mod <- glm(help4,family=binomial(link="logit"),data=dat,offset=start)
        )
        deviances[j] <-  deviance(mod0)-deviance(mod)
      }
    }
    return(deviances)
  }  
  
  # estimate tree 
  while(sig & anysplit){
    
    # compute all models
    dv <- lapply(1:nvar,function(var) {
      lapply(1:nitems,function(i) {
        n_knots   <- length(ip[[count]][[i]])
        deviances <- matrix(rep(0,n_s[var]*n_knots),ncol=n_knots)
        for(kn in 1:n_knots){
          deviances[,kn] <- allmodels(i,var,kn,design_lower,design_upper)
        }
        return(deviances)
      })
    })
    
    # select optimum
    variable <- which.max(lapply(1:nvar,function(j) max(unlist(dv[[j]]))))
    item     <- which.max(lapply(1:nitems, function(j) max(dv[[variable]][[j]])))
    split    <- as.numeric(which(dv[[variable]][[item]]==max(dv[[variable]][[item]]),arr.ind=TRUE)[,1])
    knoten   <- as.numeric(which(dv[[variable]][[item]]==max(dv[[variable]][[item]]),arr.ind=TRUE)[,2])
    if(length(split)>1){
      split  <- split[1]
      knoten <- knoten[1]
      warning(paste("Maximum in iteration ",count," not uniquely defined"))
    }
    ip_old   <- ip[[count]][[item]][knoten]
    level    <- length(strsplit(ip_old,":")[[1]])
    number   <- numbers[[count]][[item]][knoten]
    left     <- max(numbers[[count]][[item]])+1
    right    <- max(numbers[[count]][[item]])+2
    
    # compute permutation test 
    dev <- rep(NA,nperm)
    
    for(perm in 1:nperm){
      dv_perm <- rep(0,n_s[variable])
      obs_aktuell <- which_obs[[count]][[item]][knoten,]
      obs_aktuell <- obs_aktuell[!is.na(obs_aktuell)]
      DM_kov_perm <- DM_kov
      DM_kov_perm[obs_aktuell,variable] <- sample(DM_kov_perm[obs_aktuell,variable],length(obs_aktuell))
      design_upper_perm      <- design_upper
      design_upper_perm[[item]][[variable]] <- whole_design(DM_kov_perm,variable,item,thresholds)
      design_lower_perm      <- design_lower
      design_lower_perm[[item]][[variable]] <- whole_design(DM_kov_perm,variable,item,thresholds,upper=FALSE)
      dv_perm <- allmodels(item,variable,knoten,design_lower_perm,design_upper_perm)
      dev[perm] <- max(dv_perm)
      if(trace){
        cat(".")
      }
    }
    
    # test decision 
    crit_val <- quantile(dev,1-(alpha/vars_evtl[[count]][[item]][knoten]))
    proof    <- max(dv[[variable]][[item]]) > crit_val
    devs[count]   <- max(dv[[variable]][[item]])
    crits[count]   <- crit_val
    pvalues[count] <- length(which(dev>max(dv[[variable]][[item]])))/nperm
    
    if(proof){
      
      # get new formula 
      help_kn2 <- ip[[count]][[item]][knoten]
      help5 <- paste0(unlist(ip[[count]])[-which(unlist(ip[[count]])==help_kn2)],collapse="+")
      help6 <- paste(ip[[count]][[item]][knoten],c(colnames(design_lower[[item]][[variable]])[split],colnames(design_upper[[item]][[variable]])[split]),sep=":")
      help7 <- paste(help6,collapse="+")
      help8 <- formula(paste("y~",help_p,"+",help5,"+",help7,"-1"))
      
      ######################
      if(level>1){
        help_kn4 <- lu(c(),1,level-1,c())
        help_kn5 <- unlist(strsplit(help_kn2,""))
        help_kn6 <- paste0(help_kn5[which(help_kn5=="_")+1],collapse="")
        knoten2  <- which(help_kn4==help_kn6)
      } else{
        knoten2 <- knoten
      }
      ######################
      
      splits <- rbind(splits,c(variable,item,split,level,knoten2,number,left,right))   
      
      
      # fit new model 
      dat <- dat0 <- data.frame(dat0,design_lower[[item]][[variable]][,split,drop=FALSE],design_upper[[item]][[variable]][,split,drop=FALSE])
      suppressWarnings(
        mod0  <- mod_potential[[count+1]] <- glm(help8,family=binomial(link="logit"),data=dat,etastart=start)
      )
      start <- predict(mod0)
      
      # generiere neue itemparameter 
      ip[[count+1]]                             <- ip[[count]]
      ip[[count+1]][[item]]                     <- rep("",length(ip[[count]][[item]])+1)
      ip[[count+1]][[item]][c(knoten,knoten+1)] <- help6
      ip[[count+1]][[item]][-c(knoten,knoten+1)]<- ip[[count]][[item]][-knoten]
      
      # passe splits_evtl an
      n_knots                                                       <- length(ip[[count+1]][[item]])
      splits_evtl[[count+1]]                                        <- splits_evtl[[count]]
      for(var in 1:nvar){
        splits_evtl[[count+1]][[item]][[var]]                       <- matrix(0,nrow=n_knots,ncol=n_s[var])
        splits_evtl[[count+1]][[item]][[var]][c(knoten,knoten+1),]  <- matrix(rep(splits_evtl[[count]][[item]][[var]][knoten,],2),nrow=2,byrow=T)
        splits_evtl[[count+1]][[item]][[var]][-c(knoten,knoten+1),] <- splits_evtl[[count]][[item]][[var]][-knoten,]
      }
      splits_evtl[[count+1]][[item]][[variable]][knoten,splits_evtl[[count+1]][[item]][[variable]][knoten,]>=split] <- NA 
      splits_evtl[[count+1]][[item]][[variable]][(knoten+1),splits_evtl[[count+1]][[item]][[variable]][(knoten+1),]<=split] <- NA
      
      # any split? 
      anysplit <- !all(is.na(unlist(splits_evtl[[count+1]])))
      
      # passe vars_evtl an 
      vars_evtl[[count+1]]                             <- vars_evtl[[count]]
      vars_evtl[[count+1]][[item]]                     <- rep(0,n_knots)
      vars_evtl[[count+1]][[item]][c(knoten,knoten+1)] <- rep(vars_evtl[[count]][[item]][knoten],2)
      vars_evtl[[count+1]][[item]][-c(knoten,knoten+1)]<- vars_evtl[[count]][[item]][-knoten]
      
      if(length(which(!is.na(splits_evtl[[count+1]][[item]][[variable]][knoten,])))==0){ 
        vars_evtl[[count+1]][[item]][knoten] <- vars_evtl[[count+1]][[item]][knoten]-1 
      }
      if(length(which(!is.na(splits_evtl[[count+1]][[item]][[variable]][knoten+1,])))==0){ 
        vars_evtl[[count+1]][[item]][knoten+1] <- vars_evtl[[count+1]][[item]][knoten+1]-1 
      }
      
      # passe which_obs an 
      which_obs[[count+1]]                               <- which_obs[[count]]
      which_obs[[count+1]][[item]]                       <- matrix(0,nrow=n_knots,ncol=npersons)
      which_obs[[count+1]][[item]][c(knoten,knoten+1),]  <- matrix(rep(which_obs[[count]][[item]][knoten,],2),nrow=2,byrow=T)
      which_obs[[count+1]][[item]][-c(knoten,knoten+1),] <- which_obs[[count]][[item]][-knoten,]
      thresh <- ordered_values[[variable]][1:n_s[variable]][split]
      which_obs[[count+1]][[item]][knoten,DM_kov[,variable]>thresh] <- NA
      which_obs[[count+1]][[item]][(knoten+1),DM_kov[,variable]<=thresh] <- NA
      
      # passe numbers an 
      numbers[[count+1]]                              <- numbers[[count]]
      numbers[[count+1]][[item]]                      <- numeric(length=n_knots)
      numbers[[count+1]][[item]][c(knoten,knoten+1)]  <- c(left,right)
      numbers[[count+1]][[item]][-c(knoten,knoten+1)] <- numbers[[count]][[item]][-knoten] 
      
      # trace
      if(trace){
        cat(paste0("\n Split"," ",count,";"," ","Item"," ",item,"\n"))
      }
      
      # erhoehe counter
      count <- count+1 
      
    } else{
      sig <- FALSE
      if(penalize){
        if(count>1){
          help9 <- formula(paste("~",0,"+",help_p))
          help10 <- formula(paste("~",help5,"+",help7))
          mod_potential[[count]] <- penalized(y,penalized=help10,unpenalized=help9,lambda2=1e-3,data=dat0,trace=FALSE)
        } else{
          help9  <- formula(paste("~",0,"+",help_p))
          help10 <- formula(paste("~",paste0(unlist(ip[[1]]),collapse="+")))
          mod_potential[[count]] <- penalized(y,penalized=help10,unpenalized=help9,lambda2=1e-3,data=dat0,trace=FALSE)
        }
      }
    }
  }
  
  ################################################################################### 
  
  # prettify results 
  mod_opt     <- mod_potential[[count]]
  ip_opt      <- ip[[count]]
  theta_hat   <- c(coef(mod_opt)[1:(npersons-1)],0)
  beta_hat    <- coef(mod_opt)[npersons:length(coef(mod_opt))]
  
  if(count>1){
    
    dif_items   <- unique(splits[,2])
    nodif_items <- c(1:nitems)[-dif_items]
    
    beta_hat_nodif <- sapply(nodif_items,function(j) beta_hat[ip_opt[[j]]])
    beta_hat_dif   <- lapply(dif_items, function(j)  beta_hat[ip_opt[[j]]])
    names(beta_hat_dif) <- dif_items
    
    help9 <- cumsum(c(0,(n_levels-1)))
    colnames(splits) <- c("var","item","split","level","node","number","left","right")
    splits <- data.frame(cbind(splits[,1:5,drop=FALSE],"variable"=rep(NA,nrow(splits)),"threshold"=rep(NA,nrow(splits)),splits[,6:8,drop=FALSE]))
    for(i in 1:nrow(splits)){
      splits[i,6] <- colnames(DM_kov)[splits[i,1]]
      v2 <- lapply(1:nvar,function(j) ordered_values[[j]][-length(ordered_values[[j]])])
      splits[i,7] <- v2[[splits[i,1]]][splits[i,3]]
    }
    splits <- splits[,-1]
    
    for(i in dif_items){
      info <- splits[splits[,"item"]==i,]
      endnodes <- get_endnodes(info)
      names(beta_hat_dif[[paste(i)]]) <- endnodes 
    }
  
  } else{
    
    beta_hat_nodif <- beta_hat
    beta_hat_dif   <- c()
    
  }
  
  to_return <- list("splits"=splits,
                    "thetas"=theta_hat,
                    "betas_nodif"=beta_hat_nodif,
                    "betas_dif"=beta_hat_dif,
                    "pvalues"=pvalues,
                    "devs"=devs,
                    "crits"=crits)
  
  return(to_return)
  
}

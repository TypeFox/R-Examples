tree_dif <-
function(Y,
                     DM_kov,
                     npersons,
                     nitems,
                     nvar,
                     sumscore,
                     ordered_values,
                     n_levels,
                     n_s,
                     alpha,
                     nperm,
                     trace)
{
  # initializations 
  mod_potential    <- list() 
  splits           <- list(c(),c())
  devs             <- c()
  crits            <- c()
  pvalues          <- c() 
  
  splits_evtl      <- list() 
  splits_evtl[[1]] <- lapply(1:nitems,function(j) replicate(2, lapply(1:nvar, function(var) matrix(1:n_s[var],nrow=1)),simplify=FALSE))
  
  vars_evtl        <- list() 
  vars_evtl[[1]]   <- lapply(1:nitems,function(j) replicate(2, nvar, simplify=FALSE))
  
  params           <- list() 
  params[[1]]      <- lapply(1:nitems, function(j) list("int","sumscore"))
  
  which_obs        <- list() 
  which_obs[[1]]   <- lapply(1:nitems,function(j) replicate(2, matrix(1:npersons,nrow=1), simplify=FALSE))
  
  numbers          <- list()
  numbers[[1]]     <- lapply(1:nitems,function(j) replicate(2, 1, simplify=FALSE))
  
  count            <- 1 
  
  design_lower <- designlists_logistic(DM_kov,nvar,ordered_values,n_levels,n_s)[[1]]
  design_upper <- designlists_logistic(DM_kov,nvar,ordered_values,n_levels,n_s)[[2]]
  sig <- TRUE
  anysplit <- TRUE
  
  # compute models without DIF 
  mod_potential[[count]] <- mod0 <- lapply(1:nitems, function(i){
    dat0 <- as.data.frame(cbind("y"=Y[,i],"int"=rep(1,npersons),sumscore))
    suppressWarnings(
      glm(y~int+sumscore-1,family=binomial(link="logit"),data=dat0)
    )
  })
  
  # function to compute all possible models 
  allmodels <- function(i,var,co,kn,design_lower,design_upper){
    
    deviances <- rep(0,n_s[var])
    dat   <- as.data.frame(cbind("y"=Y[,i],"int"=rep(1,npersons),sumscore,do.call(cbind,design_lower),do.call(cbind,design_upper)))
    splits_aktuell <- splits_evtl[[count]][[i]][[co]][[var]][kn,]
    splits_aktuell <- splits_aktuell[!is.na(splits_aktuell)]
    
    if(length(splits_aktuell)>0){
      
      for(j in splits_aktuell){
        
        help1 <- params[[count]][[i]][[co]]
        help2 <- params[[count]][[i]][[-co]]
        help3 <- help1[-kn]
        help4 <- paste(help1[kn],c(colnames(design_lower[[var]])[j],colnames(design_upper[[var]])[j]),sep=":")
        help5 <- paste(c(help2,help3,help4), collapse="+")
        help6 <- formula(paste("y~",help5,"-1"))
        suppressWarnings(
          mod            <- glm(help6,family=binomial(link="logit"),data=dat)
        )
        deviances[j] <- deviance(mod0[[i]]) - deviance(mod)
      }
    }
    return(deviances)
  }
  
  # estimate tree 
  while(sig & anysplit){
    
    # compute all models 
    dv <- lapply(1:nvar,function(var) {
      lapply(1:nitems,function(i) {
        lapply(1:2, function(co){
          n_knots   <- length(params[[count]][[i]][[co]])
          deviances <- matrix(rep(0,n_s[var]*n_knots),ncol=n_knots)
          for(kn in 1:n_knots){
            deviances[,kn] <- allmodels(i,var,co,kn,design_lower,design_upper)
          }
          return(deviances)
        })
      })
    })
    
    # select optimum
    variable    <- which.max(lapply(1:nvar,function(j) max(unlist(dv[[j]]))))
    item        <- which.max(lapply(1:nitems, function(j) max(unlist(dv[[variable]][[j]]))))
    comp        <- which.max(lapply(1:2, function(j) max(dv[[variable]][[item]][[j]])))
    
    split       <- as.numeric(which(dv[[variable]][[item]][[comp]]==max(dv[[variable]][[item]][[comp]]),arr.ind=TRUE)[,1])
    knoten      <- as.numeric(which(dv[[variable]][[item]][[comp]]==max(dv[[variable]][[item]][[comp]]),arr.ind=TRUE)[,2])
    
    level       <- length(strsplit(params[[count]][[item]][[comp]][[knoten]],":")[[1]])
    
    n_nodes     <- length(params[[count]][[item]][[comp]])
    number      <- numbers[[count]][[item]][[comp]][knoten]
    left        <- max(numbers[[count]][[item]][[comp]])+1
    right       <- max(numbers[[count]][[item]][[comp]])+2
    
    # compute permutation test 
    dev <- rep(NA,nperm)
    
    for(perm in 1:nperm){
      dv_perm <- rep(0,n_s[variable])
      obs_aktuell <- which_obs[[count]][[item]][[comp]][knoten,]
      obs_aktuell <- obs_aktuell[!is.na(obs_aktuell)]
      DM_kov_perm <- DM_kov
      DM_kov_perm[obs_aktuell,variable] <- sample(DM_kov_perm[obs_aktuell,variable],length(obs_aktuell))
      design_upper_perm      <- designlists_logistic(DM_kov_perm,nvar,ordered_values,n_levels,n_s)[[1]]
      design_lower_perm      <- designlists_logistic(DM_kov_perm,nvar,ordered_values,n_levels,n_s)[[2]]
      dv_perm <- allmodels(item,variable,comp,knoten,design_lower_perm,design_upper_perm)
      dev[perm] <- max(dv_perm)
      if(trace){
        cat(".")
      }
    }
    
    n_comp <- sum(sapply(1:2, function(j) !all(is.na(unlist(splits_evtl[[count]][[item]][[j]])))))
    
    # test decision 
    crit_val       <- quantile(dev,1-(alpha/n_comp/vars_evtl[[count]][[item]][[comp]][knoten]))
    proof          <- max(dv[[variable]][[item]][[comp]]) > crit_val
    devs[count]    <- max(dv[[variable]][[item]][[comp]])
    crits[count]   <- crit_val
    pvalues[count] <- length(which(dev>max(dv[[variable]][[item]][[comp]])))/nperm
    
    
    if(proof){
      
      # get new formula 
      help1 <- params[[count]][[item]][[comp]]
      help2 <- params[[count]][[item]][[-comp]]
      help3 <- help1[-knoten]
      help4 <- paste(help1[knoten],c(colnames(design_lower[[variable]])[split],colnames(design_upper[[variable]])[split]),sep=":")
      help5 <- paste(c(help2,help3,help4), collapse="+")
      help6 <- formula(paste("y~",help5,"-1"))
      
      
      ######################
      if(level>1){
        help_kn1 <- lu(c(),1,level-1,c())
        help_kn2 <- unlist(strsplit(params[[count]][[item]][[comp]][[knoten]],""))
        help_kn3 <- paste0(help_kn2[which(help_kn2=="_")+1],collapse="")
        knoten2  <- which(help_kn1==help_kn3)
      } else{
        knoten2 <- knoten
      }
      ######################     
      
      splits[[comp]] <- rbind(splits[[comp]],c(variable,item,split,level,knoten2,number,left,right))   
      
      # fit new model 
      dat   <- as.data.frame(cbind("y"=Y[,item],"int"=rep(1,npersons),sumscore,do.call(cbind,design_lower),do.call(cbind,design_upper)))
      mod_potential[[count+1]] <- mod_potential[[count]]
      suppressWarnings(
        mod_potential[[count+1]][[item]] <- glm(help6,family=binomial(link="logit"),data=dat)
      )
      mod0  <- mod_potential[[count+1]]
      
      # update params
      params[[count+1]]                 <- params[[count]]
      params[[count+1]][[item]][[comp]] <- rep("",n_nodes+1)
      params[[count+1]][[item]][[comp]][c(knoten,knoten+1)] <- help4
      params[[count+1]][[item]][[comp]][-c(knoten,knoten+1)]<- params[[count]][[item]][[comp]][-knoten]
      
      # passe splits_evtl an
      splits_evtl[[count+1]]                                        <- splits_evtl[[count]]
      for(var in 1:nvar){
        splits_evtl[[count+1]][[item]][[comp]][[var]]                       <- matrix(0,nrow=n_nodes+1,ncol=n_s[var])
        splits_evtl[[count+1]][[item]][[comp]][[var]][c(knoten,knoten+1),]  <- matrix(rep(splits_evtl[[count]][[item]][[comp]][[var]][knoten,],2),nrow=2,byrow=T)
        splits_evtl[[count+1]][[item]][[comp]][[var]][-c(knoten,knoten+1),] <- splits_evtl[[count]][[item]][[comp]][[var]][-knoten,]
      }
      splits_evtl[[count+1]][[item]][[comp]][[variable]][knoten,splits_evtl[[count+1]][[item]][[comp]][[variable]][knoten,]>=split] <- NA 
      splits_evtl[[count+1]][[item]][[comp]][[variable]][(knoten+1),splits_evtl[[count+1]][[item]][[comp]][[variable]][(knoten+1),]<=split] <- NA
      
      # any split? 
      anysplit <- !all(is.na(unlist(splits_evtl[[count+1]])))
      
      # passe vars_evtl an 
      vars_evtl[[count+1]]                                     <- vars_evtl[[count]]
      vars_evtl[[count+1]][[item]][[comp]]                     <- rep(0,n_nodes+1)
      vars_evtl[[count+1]][[item]][[comp]][c(knoten,knoten+1)] <- rep(vars_evtl[[count]][[item]][[comp]][knoten],2)
      vars_evtl[[count+1]][[item]][[comp]][-c(knoten,knoten+1)]<- vars_evtl[[count]][[item]][[comp]][-knoten]
      
      if(length(which(!is.na(splits_evtl[[count+1]][[item]][[comp]][[variable]][knoten,])))==0){ 
        vars_evtl[[count+1]][[item]][[comp]][knoten] <- vars_evtl[[count+1]][[item]][[comp]][knoten]-1 
      }
      if(length(which(!is.na(splits_evtl[[count+1]][[item]][[comp]][[variable]][knoten+1,])))==0){ 
        vars_evtl[[count+1]][[item]][[comp]][knoten+1] <- vars_evtl[[count+1]][[item]][[comp]][knoten+1]-1 
      }
      
      # passe which_obs an 
      which_obs[[count+1]]                                       <- which_obs[[count]]
      which_obs[[count+1]][[item]][[comp]]                       <- matrix(0,nrow=n_nodes+1,ncol=npersons)
      which_obs[[count+1]][[item]][[comp]][c(knoten,knoten+1),]  <- matrix(rep(which_obs[[count]][[item]][[comp]][knoten,],2),nrow=2,byrow=T)
      which_obs[[count+1]][[item]][[comp]][-c(knoten,knoten+1),] <- which_obs[[count]][[item]][[comp]][-knoten,]
      thresh <- ordered_values[[variable]][1:n_s[variable]][split]
      which_obs[[count+1]][[item]][[comp]][knoten,DM_kov[,variable]>thresh] <- NA
      which_obs[[count+1]][[item]][[comp]][(knoten+1),DM_kov[,variable]<=thresh] <- NA
      
      # passe numbers an 
      numbers[[count+1]]                              <- numbers[[count]]
      numbers[[count+1]][[item]][[comp]]                      <- numeric(length=n_nodes+1)
      numbers[[count+1]][[item]][[comp]][c(knoten,knoten+1)]  <- c(left,right)
      numbers[[count+1]][[item]][[comp]][-c(knoten,knoten+1)] <- numbers[[count]][[item]][[comp]][-knoten] 
      
      # trace
      if(trace){
        cat(paste0("\n Split"," ",count,";"," ","Item"," ",item,"\n"))
      }
      
      # counter
      count <- count+1 
    } else{
      sig <- FALSE
    }
  }
  
  
  ################################################################################### 
  
  mod_opt       <- mod_potential[[count]]
  params_opt    <- params[[count]]
  indices       <- lapply(1:nitems, function(i) lapply(1:2, function(j) numeric(length(params_opt[[i]][[j]]))))
  
  for(i in 1:nitems){
    for(j in 1:2){
      cn <- strsplit(params_opt[[i]][[j]],":")
      wn <- strsplit(names(coef(mod_opt[[i]])),":")
      for(k in 1:length(cn)){
        indices[[i]][[j]][k] <- which(sapply(1:length(wn), function(z) all(cn[[k]] %in% wn[[z]])))
      }
    }
  }
  
  if(count>1){
    
    dif_items   <- unique(c(splits[[1]][,2],splits[[2]][,2]))
    nodif_items <- c(1:nitems)[-dif_items]
    
    # parameters for non DIF items 
    if(length(nodif_items)!=0){
      gammas_hat_nodif <- sapply(nodif_items, function(j) coef(mod_opt[[j]])["int"])
      names(gammas_hat_nodif) <- paste0("gamma",nodif_items)
      alphas_hat_nodif <- sapply(nodif_items, function(j) coef(mod_opt[[j]])["sumscore"])
      names(alphas_hat_nodif) <- paste0("alpha",nodif_items)
    } else{
      gammas_hat_nodif <- c()
      alphas_hat_nodif <- c()
    }
    
    # parameters for DIF items 
    gammas_hat_dif   <- lapply(dif_items, function(j) {
      g_act <- coef(mod_opt[[j]])[indices[[j]][[1]]]
      names(g_act) <- params_opt[[j]][[1]]
      g_act
    })
    names(gammas_hat_dif) <- dif_items
    alphas_hat_dif   <- lapply(dif_items, function(j) {
      a_act <- coef(mod_opt[[j]])[indices[[j]][[2]]]
      names(a_act) <- params_opt[[j]][[2]]
      a_act
    })
    names(alphas_hat_dif) <- dif_items
    
    help9 <- cumsum(c(0,(n_levels-1)))
    
    modif_splits <- function(splits){
      
      colnames(splits) <- c("var","item","split","level","node","number","left","right")
      splits <- data.frame(cbind(splits[,1:5,drop=FALSE],"variable"=rep(NA,nrow(splits)),"threshold"=rep(NA,nrow(splits)),splits[,6:8,drop=FALSE]))
      for(i in 1:nrow(splits)){
        splits[i,6] <- colnames(DM_kov)[splits[i,1]]
        v2 <- lapply(1:nvar,function(j) ordered_values[[j]][-length(ordered_values[[j]])])
        splits[i,7] <- v2[[splits[i,1]]][splits[i,3]]
      }
      splits <- splits[,-1]
      return(splits)
      
    }
    
    for(i in 1:2){
      if(!is.null(splits[[i]])) {
        splits[[i]] <- modif_splits(splits[[i]])
      }
    }
    names(splits) <- c("gammas","alphas")
    
    for(i in dif_items){
      info1 <- splits[[1]][splits[[1]][,"item"]==i,]
      if(is.null(info1) || nrow(info1)==0){
        endnodes1 <- 1 
      } else{
        endnodes1 <- get_endnodes(info1)
      }
      names(gammas_hat_dif[[paste(i)]]) <- endnodes1   
      
      info2 <- splits[[2]][splits[[2]][,"item"]==i,]
        if(is.null(info2) || nrow(info2)==0){
          endnodes2 <- 1 
        } else{
          endnodes2 <- get_endnodes(info2) 
        }
      names(alphas_hat_dif[[paste(i)]]) <- endnodes2 
    }
    
    
  } else{
    
    gammas_hat_nodif  <- sapply(1:nitems, function(j) coef(mod_opt[[j]])["int"])
    names(gammas_hat_nodif) <- paste0("gamma",1:nitems)
    gammas_hat_dif    <- c()
    alphas_hat_nodif  <- sapply(1:nitems, function(j) coef(mod_opt[[j]])["sumscore"])
    names(alphas_hat_nodif) <- paste0("alpha",1:nitems)
    alphas_hat_dif    <- c() 
    
  }
  
  to_return <- list("splits"=splits,
                    "gammas_nodif"=gammas_hat_nodif,
                    "gammas_dif"=gammas_hat_dif,
                    "alphas_nodif"=alphas_hat_nodif,
                    "alphas_dif"=alphas_hat_dif,
                    "pvalues"=pvalues,
                    "devs"=devs,
                    "crits"=crits)
  
  return(to_return)
  
}

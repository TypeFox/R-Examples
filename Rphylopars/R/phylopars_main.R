convert_to_means <- function(df,index_col = 1,sort_vec,FUN = function(X) mean(X,na.rm=TRUE))
{
  ret <- matrix(NA,ncol = ncol(df)-1,nrow = length(unique(df[,index_col])))
  index <- df[,index_col]
  col_counter <- 1
  cols <- colnames(df)[-index_col]
  for(i in (1:ncol(df))[-index_col])
  {
    ret[,col_counter] <- sapply(split(df[,i],index),FUN = FUN)
    col_counter <- col_counter + 1
  }
  rownames(ret) <- names(split(df[,i],index))
  colnames(ret) <- cols
  if(missing(sort_vec)) return(ret) else return(ret[sort_vec,,drop=F])
}

simtraits <- function(ntaxa=15,ntraits=4,nreps=1,nmissing=0,tree,v,anc,intraspecific,model="BM",parameters,nsim=1,return.type="data.frame")
{
  if(nmissing>(ntaxa*ntraits*nreps)) nmissing <- round(runif(1,ntaxa*ntraits*nreps-1))
  if(missing(tree))
  {
    tree <- pbtree(n=ntaxa)
  } else ntaxa <- length(tree$tip.label)
  tree <- reorder(tree,"postorder")
  perm_tree <- tree
  if(model!="BM") tree <- transf.branch.lengths(phy = tree,model = model,parameters = parameters)$tree
  if(missing(v))
  {
    v <- matrix(0,ntraits,ntraits)
    npars <- length(v[upper.tri(v,TRUE)])
    pars <- rnorm(npars)
    v[upper.tri(v,TRUE)] <- pars
    v <- t(v)%*%v
  } else ntraits = length(diag(v))
  if(missing(anc))
  {
    anc <- rep(0,ntraits)
  } else if(length(anc)==1) anc <- rep(anc,ntraits)
  anc <- as.double(anc)
  
  if(missing(intraspecific)) intraspecific <- 0.1
  if(length(intraspecific)==1)
  {
    opt <- 1
    intraspecific <- matrix(rep(intraspecific,ntraits*ntaxa),nrow = ntaxa,ncol = ntraits)
  } else if(length(intraspecific)==ntraits)
  {
    opt <- 2
    intraspecific <- t(matrix(rep(intraspecific,ntaxa),nrow=ntraits))
  } else opt <- 3
  anc_mat <- matrix(1,ntaxa) %*% anc
  Xall <- sim.char(phy = tree,par = v,nsim = nsim)
  colnames(Xall) <- paste("V",1:ntraits,sep="")
  if(nreps==1 & nmissing==0 & nsim==1)
  {
    if(return.type=="matrix") return(list(trait_data=Xall[,,1],tree=perm_tree,sim_tree=tree)) else
      return(list(trait_data=data.frame(species=rownames(Xall[,,1]),Xall[,,1]),tree=perm_tree,sim_tree=tree))
  } else if(nreps==1 & nmissing==0) 
  {
    if(return.type=="matrix")
    {
      return(list(trait_data=lapply(apply(Xall,3,function(X) list(X)),function(X) X[[1]]),tree=perm_tree,sim_tree=tree))
    } else
      return(list(trait_data=lapply(apply(Xall,3,function(X) list(X)),function(X) data.frame(species=rownames(X[[1]]),X[[1]])),tree=perm_tree,sim_tree=tree))
  }
  
  X <- original_X <- rep(list(matrix(0,ntaxa*nreps,ntraits)),nsim)
  for(j in 1:nsim)
  {
    Xall[,,j] <- Xall[,,j] + anc_mat
    if(nreps==1)
    {
      X[[j]][1:(ntraits*ntaxa)] <- original_X[[j]][1:(ntraits*ntaxa)] <- Xall[,,j]
    } else
    {
      for(jj in 1:nreps)
      {
        original_X[[j]] <- Xall[,,j]
        X[[j]][1:ntaxa + (jj-1)*(ntaxa),] <- rnorm(n = ntraits*ntaxa,mean = Xall[,,j],sd = intraspecific)
      }
    }
    X[[j]][sample(1:length(X[[j]]),nmissing)] <- NA
    colnames(X[[j]]) <- paste("V",1:ncol(X[[j]]),sep = "")
    species <- rep(rownames(Xall[,,j]),nreps)
    rownames(X[[j]]) <- 1:nrow(X[[j]])
    X[[j]] <- data.frame(species=species,X[[j]])
    if(nreps==1) rownames(X[[j]]) <- species
  }
  if(nsim==1) list(trait_data=X[[1]],tree=perm_tree,sim_tree=tree,original_X=original_X[[1]]) else
    list(trait_data=X,tree=perm_tree,sim_tree=tree,original_X=original_X)
}

write.phylopars <- function(trait_data,tree,data_file,tree_file,species_identifier="species")
{ 
  rnms <- rownames(trait_data)
  rownames(trait_data) <- NULL
  trait_data <- as.data.frame(trait_data)
  if(is.null(colnames(trait_data))) colnames(trait_data) <- paste("V",1:ncol(trait_data),sep="")
  featurenames <- colnames(trait_data)
  if(length(featurenames[which(featurenames==species_identifier)])==0)
  {
    if(is.null(rnms))
    {
      warning("No tip labels and no row names. Assuming trait_data are in order of tips with no within-species replicates.")
      rnms <- tree$tip.label
    }
    trait_data[[species_identifier]] <- rnms
    featurenames <- colnames(trait_data)
  }
  featurenames[which(featurenames==species_identifier)] <- "species"
  colnames(trait_data) <- featurenames
  featurenames <- featurenames[featurenames!="species"]
  namecheck <- name.check(tree,data.names=trait_data$species)
  if((namecheck!="OK")[[1]])
  {
    if(length(namecheck$data_not_tree)>0)
    {
      warning(paste("trait_data for",paste(unique(namecheck$data_not_tree),collapse=", "),"not in tree. Pruning from trait_dataset."))
      drop <- which(!is.na(match(trait_data$species,namecheck$data_not_tree)))
      trait_data <- trait_data[-drop,]
    }
    if(length(namecheck$tree_not_data)>0)
    {
      missing_species <- namecheck$tree_not_data
      for(i in 1:length(missing_species))
      {
        trait_data <- rbind(trait_data,NA)
        trait_data$species[nrow(trait_data)] <- missing_species[i]
      }
    }
  }
  if(ncol(trait_data)==2)
  {
    trait_data[,-match("species",colnames(trait_data))] <- as.double(trait_data[,-match("species",colnames(trait_data))])
  } else
  {
    trait_data[,-match("species",colnames(trait_data))] <- apply(trait_data[,-match("species",colnames(trait_data)),drop=F],2,as.double)
  }
  species <- trait_data$species
  spcs <- unique(species)
  trait_data <- trait_data[,-which(colnames(trait_data)=="species"),drop=FALSE]
  out <- matrix(nrow=length(unique(species)),ncol=ncol(trait_data))
  for(i in 1:ncol(trait_data))
  {
    for(j in 1:length(unique(species)))
    {
      out[j,i] <- paste(trait_data[,i][species==spcs[j]],collapse = ";")
      out[j,i] <- gsub("NA","",out[j,i])
      semicolons <- which(is.finite(match(strsplit(out[j,i],"")[[1]],";")))
      if(length(semicolons)>0)
      {
        if(semicolons[length(semicolons)]==length(strsplit(out[j,i],"")[[1]]))
        {
          strng <- strsplit(out[j,i],"")[[1]]
          out[j,i] <- paste(strng[1:(length(strng)-1)],collapse = "")
        }
        if(semicolons[1]==1) out[j,i] <- sub(";","",out[j,i])
      }
      while(grepl(";;",out[j,i]))
      {
        out[j,i] <- sub(";;",";",out[j,i])
      }
      semicolons <- which(is.finite(match(strsplit(out[j,i],"")[[1]],";")))
      if(length(semicolons)>0)
      {
        if(semicolons[length(semicolons)]==length(strsplit(out[j,i],"")[[1]]))
        {
          strng <- strsplit(out[j,i],"")[[1]]
          out[j,i] <- paste(strng[1:(length(strng)-1)],collapse = "")
        }
        if(semicolons[1]==1) out[j,i] <- sub(";","",out[j,i])
      }
    }
  }
  rownames(out) <- unique(species)
  write.table(out,data_file,sep="\t",quote=FALSE,na = "",col.names=NA)
  write.tree(tree,tree_file)
}
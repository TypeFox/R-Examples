ultraFastAnc <- function(phy,x,vars=FALSE,CI=FALSE)
{
  original_tree <- reorder(phy,"postorder")
  if(!is.binary.tree(original_tree)) phy <- reorder(multi2di(phy),"postorder") else
    phy <- original_tree
  prep <- prep_multipic(as.matrix(x),phy,rescaled.tree=TRUE,pic_recon=TRUE)
  pics <- do.call(multipic,prep)
  nn <- phy$Nnode
  len_vec <- phy$edge.length
  nspecies <- length(phy$tip.label)
  nedge <- length(len_vec)
  pic_len_vec <- pics$pic_len
  E <- phy$edge
  m <- match(phy$edge[,2],phy$edge[,1],nomatch = 0)
  L <- cbind(m,m)
  L[m>0,2] <- L[m>0,2] + 1
  
  pic_ace <- pics$pic_recon
  mse <- mean(pics$contrasts^2)
  ace_hat <- var_hat <- pic_ace
  var_hat[nspecies+1] <- prod(pic_len_vec[(nedge-1):nedge])/sum(pic_len_vec[(nedge-1):nedge])
  ret <- C_ultrafast(nedge,L,E,pic_ace,ace_hat,var_hat,len_vec,pic_len_vec)
  ace_hat <- ret[,1]
  var_hat <- ret[,2]
  ret_seq <- (nspecies+1):length(ace_hat)
  ace_hat <- ace_hat[ret_seq]
  var_hat <- var_hat[ret_seq]
  names(ace_hat) <- names(var_hat) <- as.character(1:length(var_hat)+nspecies)
  var_hat <- mse*var_hat
  if(!is.binary.tree(original_tree))
  {
    ancNames <- matchNodes(original_tree, phy)
    ace_hat <- ace_hat[as.character(ancNames[,2])]
    names(ace_hat) <- ancNames[, 1]
    var_hat <- var_hat[as.character(ancNames[,2])]
    names(var_hat) <- ancNames[,1]
  }
  CI95 <- sqrt(var_hat)*1.96
  ace_hat_CI95 <- cbind(ace_hat-CI95,ace_hat+CI95)
  
  if(!vars & !CI) return(ace_hat)
  if(CI & !vars) return(list(ace=ace_hat,CI95=ace_hat_CI95))
  if(vars & !CI) return(list(ace=ace_hat,var=var_hat))
  return(list(ace=ace_hat,var=var_hat,CI95=ace_hat_CI95))
}

prep_ultraFastAnc <- function(phy,x,vars=FALSE,CI=FALSE)
{
  original_tree <- reorder(phy,"postorder")
  if(!is.binary.tree(original_tree)) phy <- reorder(multi2di(phy),"postorder") else
    phy <- original_tree
  prep <- prep_multipic(as.matrix(x),phy,rescaled.tree=TRUE,pic_recon=TRUE)
  nn <- phy$Nnode
  len_vec <- phy$edge.length
  nspecies <- length(phy$tip.label)
  nedge <- length(len_vec)
  pic_len_vec <- len_vec #pics$pic_len
  E <- phy$edge
  m <- match(phy$edge[,2],phy$edge[,1],nomatch = 0)
  L <- cbind(m,m)
  L[m>0,2] <- L[m>0,2] + 1
  
  pic_ace <- matrix(0,nrow=nspecies+phy$Nnode,ncol=1)
  mse <- 0
  ace_hat <- var_hat <- pic_ace
  var_hat[nspecies+1] <- 0
  multipic_args <- prep
  ultraFastAnc_args <- list(nedge=nedge,L=L,E=E,pic_ace=pic_ace,ace_hat=ace_hat,var_hat=var_hat,len_vec=len_vec,pic_len_vec=pic_len_vec)
  return(list(multipic_args=multipic_args,ultraFastAnc_args=ultraFastAnc_args,vars=vars,CI=CI,original_tree=original_tree,phy=phy))
}

update_ultraFastAnc <- function(new_x,args)
{
  original_tree <- args$original_tree
  phy <- args$phy
  nedge <- length(args$multipic_args$edge_len)
  nspecies <- args$multipic_args$ntip
  vars <- args$vars
  CI <- args$CI
  args[[1]]$phe[1:args[[1]]$ntip] <- new_x
  pics <- do.call(multipic,args[[1]])
  args[[2]]$pic_ace <- pics$pic_recon
  mse <- mean(pics$contrasts^2)
  ace_hat <- var_hat <- pic_ace <- args[[2]]$pic_ace
  pic_len_vec <- pics$pic_len
  args[[2]]$var_hat[nspecies+1] <- prod(pic_len_vec[(nedge-1):nedge])/sum(pic_len_vec[(nedge-1):nedge])
  args[[2]]$ace_hat <- pics$pic_recon
  args[[2]]$pic_len_vec <- pic_len_vec
  args[[2]]$pic_ace <- pics$pic_recon
  ret <- do.call(C_ultrafast,args[[2]])
  ace_hat <- ret[,1]
  var_hat <- ret[,2]
  ret_seq <- (nspecies+1):length(ace_hat)
  ace_hat <- ace_hat[ret_seq]
  var_hat <- var_hat[ret_seq]
  names(ace_hat) <- names(var_hat) <- as.character(1:length(var_hat)+nspecies)
  var_hat <- mse*var_hat
  if(!is.binary.tree(original_tree))
  {
    ancNames <- matchNodes(original_tree, phy)
    ace_hat <- ace_hat[as.character(ancNames[,2])]
    names(ace_hat) <- ancNames[, 1]
    var_hat <- var_hat[as.character(ancNames[,2])]
    names(var_hat) <- ancNames[,1]
  }
  CI95 <- sqrt(var_hat)*1.96
  ace_hat_CI95 <- cbind(ace_hat-CI95,ace_hat+CI95)
  
  if(!vars & !CI) return(ace_hat)
  if(CI & !vars) return(list(ace=ace_hat,CI95=ace_hat_CI95))
  if(vars & !CI) return(list(ace=ace_hat,var=var_hat))
  return(list(ace=ace_hat,var=var_hat,CI95=ace_hat_CI95))
}

prep_multipic <- function(x, phy, scaled = TRUE, var.contrasts = FALSE, rescaled.tree = FALSE, pic_recon = FALSE)
{
  if (!inherits(phy, "phylo")) 
    stop("object 'phy' is not of class \"phylo\"")
  if (is.null(phy$edge.length)) 
    stop("your tree has no branch lengths")
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  if (nb.node != nb.tip - 1) 
    stop("'phy' is not rooted and fully dichotomous")
  if (nrow(x) != nb.tip) 
    stop("length of phenotypic and of phylogenetic data do not match")
  if (any(is.na(x))) 
    stop("missing data in 'x': you may consider removing the species with missing data from your tree with the function 'drop.tip'.")
  phy <- reorder(phy, "postorder")
  nvar <- ncol(x)
  phenotype <- matrix(0,nb.tip + nb.node,nvar)
  if (is.null(rownames(x))) {
    phenotype[1:nb.tip,] <- x
  } else {
    if (all(rownames(x) %in% phy$tip.label)) 
      phenotype[1:nb.tip,] <- x[phy$tip.label,,drop=FALSE]
    else {
      phenotype[1:nb.tip,] <- x
      warning("the names of argument 'x' and the tip labels of the tree did not match: the former were ignored in the analysis.")
    }
  }
  contr <- matrix(0,phy$Nnode,nvar)
  list(ntip=as.integer(nb.tip), nnode=as.integer(nb.node), 
       edge1=as.integer(phy$edge[, 1]), edge2=as.integer(phy$edge[, 2]), 
       edge_len=as.double(phy$edge.length), phe=phenotype, contr=contr, 
       var_contr=double(phy$Nnode), scaled=as.integer(scaled), pic_len=as.integer(rescaled.tree), pic_recon=as.integer(pic_recon))
}

prep_multipic2 <- function(x, phy, edge_len_mat,scaled = TRUE, var.contrasts = FALSE, rescaled.tree = FALSE, pic_recon=FALSE)
{
  if (!inherits(phy, "phylo")) 
    stop("object 'phy' is not of class \"phylo\"")
  if (is.null(phy$edge.length)) 
    stop("your tree has no branch lengths")
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  if (nb.node != nb.tip - 1) 
    stop("'phy' is not rooted and fully dichotomous")
  if (nrow(x) != nb.tip) 
    stop("length of phenotypic and of phylogenetic data do not match")
  if (any(is.na(x))) 
    stop("missing data in 'x': you may consider removing the species with missing data from your tree with the function 'drop.tip'.")
  phy <- reorder(phy, "postorder")
  nvar <- ncol(x)
  phenotype <- matrix(0,nb.tip + nb.node,nvar)
  if (is.null(rownames(x))) {
    phenotype[1:nb.tip,] <- x
  } else {
    if (all(rownames(x) %in% phy$tip.label)) 
      phenotype[1:nb.tip,] <- x[phy$tip.label,,drop=FALSE]
    else {
      phenotype[1:nb.tip,] <- x
      warning("the names of argument 'x' and the tip labels of the tree did not match: the former were ignored in the analysis.")
    }
  }
  nedge <- nrow(phy$edge)
  if((dim(edge_len_mat)[1]==nedge | dim(edge_len_mat)[1]==(nedge+1)) & dim(edge_len_mat)[2]==nvar) edge_len_mat <- t(edge_len_mat) else
    if((dim(edge_len_mat)[2]!=nedge | dim(edge_len_mat)[2]!=(nedge+1)) | dim(edge_len_mat)[1]!=nvar) stop("edge_len_mat is not of dimension nvar x nedge")
  contr <- matrix(0,phy$Nnode,nvar)
  list(ntip=as.integer(nb.tip), nnode=as.integer(nb.node), 
       edge1=as.integer(phy$edge[, 1]), edge2=as.integer(phy$edge[, 2]), 
       edge_len=edge_len_mat, phe=phenotype, contr=contr, 
       var_contr=matrix(0,phy$Nnode,nvar), scaled=as.integer(scaled), pic_len=as.integer(rescaled.tree), pic_recon=as.integer(pic_recon))
}

paint.edges <- function(tree,species.groups,average.nodes = TRUE,root.edge = TRUE)
{
  tree <- reorder(tree,"postorder")
  r <- ace(species.groups,tree,type = "discrete")
  nspecies <- length(tree$tip.label)
  nedge <- nrow(tree$edge)
  marg <- matrix(0,nspecies+tree$Nnode,nlevels(species.groups))
  marg[1:nspecies,] <- sapply(levels(species.groups),function(X) as.integer(species.groups==X),USE.NAMES = TRUE)
  marg[(nspecies+1):(nspecies+tree$Nnode),] <- r$lik.anc
  rownames(marg) <- 1:nrow(marg)
  if(average.nodes)
  {
    if(root.edge) m <- rbind((marg[tree$edge[,2],] + marg[tree$edge[,1],])/2,marg[nspecies+1,]) else m <- (marg[tree$edge[,2],] + marg[tree$edge[,1],])/2
    rownames(m) <- 1:nrow(m)
  } else
  {
    if(root.edge) m <- rbind(marg[tree$edge[,2],],marg[nspecies+1,]) else   m <- marg[tree$edge[,2],]
    rownames(m) <- 1:nrow(m)
  }
  colnames(m) <- levels(species.groups)
  m
}

Rinv5 <- function(nspecies,nedge,edge.length,anc,des,L,R,painted.edges,phylocov_list)
{
  nvar <- ncol(phylocov_list[[1]])
  L1 <- LL <- R1 <- RR <- LR <- NA
  available_mats <- 0 # 0 indicates L and R are missing
  
  # make sure R is a matrix
  if(!missing(R))
  {
    if(!is.matrix(R))
    {
      if(is.data.frame(R))
      {
        R <- as.matrix(R)
      } else if(is.list(R))
      {
        stop("R should be in matrix or vector form.")
      } else if(length(R)==(nspecies*nvar))
      {
        R <- matrix(R,ncol=1)
      } else stop("R should be in matrix or vector form.")
    } else R <- matrix(R,ncol=1)
    mR <- ncol(R)
    RR <- rep(list(matrix(0,mR,mR)),nedge+1)
    R1 <- rep(list(matrix(0,nvar,mR)),nedge+1)
    available_mats <- available_mats + 2 # 2 indicates only R is present; 3 indicates both
  }
  
  if(!missing(L))
  {
    if(!is.matrix(L))
    {
      if(is.data.frame(L))
      {
        L <- as.matrix(L)
      } else if(is.list(L))
      {
        stop("L should be in matrix or vector form.")
      } else if(length(L)==nspecies)
      {
        L <- matrix(L,ncol=1)
      } else stop("L should be in matrix or vector form.")
    }
    mL <- ncol(L)
    LL <- rep(list(matrix(0,mL*nvar,mL*nvar)),nedge+1)
    L1 <- rep(list(matrix(0,nvar*mL,nvar)),nedge+1)
    available_mats <- available_mats + 1 # 1 indicates only L is present
  }
  
  if(available_mats==3) LR <- rep(list(matrix(0,mL*nvar,mR)),nedge+1)
  
  logd <- vector("double",nedge+1)
  p <- pA <- rep(list(matrix(0,nvar,nvar)),nedge+1)
  pA_des_value <- matrix(0,nvar,nvar)
  i <- des_node <- anc_node <- 0
  ngroups <- ncol(painted.edges)
  len <- matrix(0,nvar,nvar)
  for(i in 1:nedge)
  {
    des_node <- des[i]
    anc_node <- anc[i]
    len <- len*0
    for(z in 1:ngroups) len <- len + painted.edges[i,z]*edge.length[i]*phylocov_list[[z]] 
    if(des_node<=nspecies)
    {
      logd[des_node] <- determinant(len,logarithm = TRUE)$modulus[[1]]
      logd[anc_node] <- logd[anc_node] + logd[des_node]
      p[[des_node]] <- solve(len)
      pA[[anc_node]] <- pA[[anc_node]] + p[[des_node]]
      index <- des_node + nspecies*(1:nvar-1)
      if(available_mats>0 & available_mats!=2)
      {
        L1[[des_node]] <- p[[des_node]] %x% t(L[des_node,,drop=F])
        LL[[des_node]] <- p[[des_node]] %x% t(L[des_node,,drop=F]) %x% L[des_node,,drop=F]
        
        LL[[anc_node]] <- LL[[anc_node]] + LL[[des_node]]
        L1[[anc_node]] <- L1[[anc_node]] + L1[[des_node]]
      }
      if(available_mats>1)
      {
        R1[[des_node]] <- p[[des_node]] %*% R[index,,drop=FALSE]
        RR[[des_node]] <- t(R[index,,drop=FALSE]) %*% p[[des_node]] %*% R[index,,drop=FALSE] # Q
        
        RR[[anc_node]] <- RR[[anc_node]] + RR[[des_node]]
        R1[[anc_node]] <- R1[[anc_node]] + R1[[des_node]]
      }
      if(available_mats>2)
      {
        LR[[des_node]] <- p[[des_node]] %x% t(L[des_node,,drop=FALSE]) %*% R[index,,drop=FALSE]
        LR[[anc_node]] <- LR[[anc_node]] + LR[[des_node]]
      }
    } else
    {
      pA_des_value <- pA[[des_node]]
      itpa <- diag(nvar)+len%*%pA_des_value
      itpainv <- solve(itpa)
      
      logd[des_node] <- logd[des_node] + determinant(itpa,logarithm = TRUE)$modulus[[1]]
      logd[anc_node] <- logd[anc_node] + logd[des_node]
      
      p[[des_node]] <- pA_des_value %*% itpainv
      pA[[anc_node]] <- pA[[anc_node]] + p[[des_node]]
      
      if(available_mats>2)
      {
        LR[[des_node]] <- LR[[des_node]] - L1[[des_node]] %*% itpainv %*% len %*% (R1[[des_node]])
        LR[[anc_node]] <- LR[[anc_node]] + LR[[des_node]]
      }
      if(available_mats>0 & available_mats!=2)
      {
        LL[[des_node]] <- LL[[des_node]] - L1[[des_node]] %*% itpainv %*% len %*% t(L1[[des_node]])
        LL[[anc_node]] <- LL[[anc_node]] + LL[[des_node]]
        
        L1[[des_node]] <- L1[[des_node]] %*% itpainv
        L1[[anc_node]] <- L1[[anc_node]] + L1[[des_node]]
      }
      if(available_mats>1)
      {
        RR[[des_node]] <- RR[[des_node]] - t(R1[[des_node]]) %*% itpainv %*% len %*% (R1[[des_node]])
        RR[[anc_node]] <- RR[[anc_node]] + RR[[des_node]]
        
        R1[[des_node]] <- t(t(R1[[des_node]]) %*% itpainv)
        R1[[anc_node]] <- R1[[anc_node]] + R1[[des_node]]
      }
      
    }
  }
  i <- nedge+1
  len <- len*0
  for(z in 1:ngroups) len <- len + painted.edges[i,z]*edge.length[i]*phylocov_list[[z]] 
  
  pA_des_value <- pA[[anc_node]]
  itpa <- diag(nvar)+len%*%pA_des_value
  itpainv <- solve(itpa)
  logd[nedge+1] <- logd[nspecies+1] + determinant(itpa,logarithm = TRUE)[[1]]
  p[[nedge+1]] <- pA_des_value %*% itpainv
  
  ret <- list(p=p[[nedge+1]],logd=logd[nedge+1])
  
  if(available_mats>2)
  {
    LR[[nedge+1]] <- LR[[anc_node]] - L1[[anc_node]] %*% itpainv %*% len %*% (R1[[anc_node]])
    ret$LR <- LR[[nedge+1]]
  }  
  if(available_mats>0 & available_mats!=2)
  {
    LL[[nedge+1]] <- LL[[anc_node]] - L1[[anc_node]] %*% itpainv %*% len %*% t(L1[[anc_node]])    
    L1[[nedge+1]] <- L1[[anc_node]] %*% itpainv
    
    ret$LL <- LL[[nedge+1]]
    ret$L1 <- L1[[nedge+1]]
  }
  if(available_mats>1)
  {
    RR[[nedge+1]] <- RR[[anc_node]] - t(R1[[anc_node]]) %*% itpainv %*% len %*% (R1[[anc_node]])
    R1[[nedge+1]] <- t(t(R1[[anc_node]]) %*% itpainv)
    
    ret$R1 <- R1[[nedge+1]]
    ret$RR <- RR[[nedge+1]]
    if(available_mats==2)
    {
      ret$theta <- solve(ret$p,ret$R1)
      ret$LogLik <- -(nspecies*nvar*log(2*pi) + ret$logd + ret$RR - 2*t(ret$theta) %*% ret$R1 + t(ret$theta) %*% ret$p %*% ret$theta)/2
      ret$ReLogLik <- -((nspecies*nvar-nvar)*log(2*pi) + determinant(ret$p)$modulus[[1]] + ret$logd + ret$RR - 2*t(ret$theta) %*% ret$R1 + t(ret$theta) %*% ret$p %*% ret$theta)/2
    } else
    {
      ret$theta <- solve(ret$LL,ret$LR)
      ret$LogLik <- -(nspecies*nvar*log(2*pi) + ret$logd + ret$RR - 2*t(ret$theta) %*% ret$LR + t(ret$theta) %*% ret$LL %*% ret$theta)/2
      ret$ReLogLik <- -((nspecies*nvar-nvar)*log(2*pi) + determinant(ret$LL)$modulus[[1]] + ret$logd + ret$RR - 2*t(ret$theta) %*% ret$LR + t(ret$theta) %*% ret$LL %*% ret$theta)/2
    }
  }
  ret
}

print.evo.model <- function(x,...)
{
  cat("\nEvolutionary rate(s) (sigma2mult):\n")
  print(simplify2array(x$sigma2.mult))
  ndims <- ncol(x$evo.model.args$Y)
  if(x$evo.model.args$multirate) cat("Evolutionary rates constrained to be equal.\n")
  if(class(x$evo.model.args$species.groups)=="factor") cat("Species groups:",levels(x$evo.model.args$species.groups),"\n")
  if(x$evo.model.args$diag.phylocov) cat("Traits assumed to be independent (zero evolutionary covariance).")
  if(length(x$evo.model.args$force.zero.phylocov)>0) cat("Covariance of",length(x$evo.model.args$force.zero.phylocov),"traits with remaining traits constrained to zero.\n")
  if(length(x$evo.model.args$fixed.par)>0) cat(names(x$model.par)," constrained to equal ",x$evo.model.args$fixed.par,".\n",sep="")
  if(class(x$evo.model.args$fixed.effects)=="matrix") cat("Model fit with ", ncol(x$evo.model.args$fixed.effects)," fixed effects.\n",sep="")
  if(suppressWarnings(round(exp(lfactorial(ndims) - (log(2) + lfactorial(ndims-2)))))>x$evo.model.args$max.combn) cat("Pairwise log-likelihood approximated via Monte Carlo simulation.\n")
  cat("\nLog-likelihood: ",x$logL,"\n")
  cat("Method: ",x$evo.model.args$method,"\n")
  
  cat("\nEvolutionary model: ")
  if(x$evo.model.args$model[1]!="BM")
  {
    cat("\n")
    print(data.frame(Parameter=names(x$model.par),Value=x$model.par,row.names = x$evo.model.args$model))
  } else cat("BM")
  cat("\n")
}


evo.model <- function(tree, Y, fixed.effects = NA, species.groups, trait.groups,
                      model="BM", diag.phylocov=FALSE, method = "Pairwise REML", 
                      force.zero.phylocov=character(), species.id = "species", max.combn = 10000, painted.edges,
                      ret.level = 2, plot.LL.surface = TRUE, par.init.iters = 50, fixed.par = numeric(),
                      multirate = FALSE,subset=TRUE)
  # ret.level: 1=logL only, 2=multi rates, 3=full rates
{
  tree <- reorder(tree,"postorder")
  if(length(tree$edge.length)>nrow(tree$edge)) tree$edge.length <- tree$edge.length[1:nrow(tree$edge)]
  if(model=="OU") model <- "OUfixedRoot"
  if(missing(species.groups))
  {
    nrates.species <- 1
    species.group.levels <- NA
  } else
  {
    species.groups <- species.groups[match(tree$tip.label,names(species.groups))]
    if(missing(painted.edges))
    {
      painted.edges <- paint.edges(tree,species.groups)
    }
    nrates.species <- nlevels(species.groups)
    species.group.levels <- levels(species.groups)
  }
  if(missing(trait.groups))
  {
    nrates.traits <- 1
    trait.group.levels <- NA
  } else
  {
    nrates.traits <- nlevels(trait.groups)
    trait.group.levels <- levels(trait.groups)
  }
  if(multirate==FALSE & nrates.traits>1)
  {
    multirate <- TRUE
    warning("Setting multirate to TRUE because trait.groups was specified.",immediate. = TRUE)
  }
  if(method=="Pairwise REML" | method=="Full REML") REML <- 1 else REML <- 0
  evo.model.args <- as.list(environment())
  evo.model.args$REML <- evo.model.args$trait.group.levels <- evo.model.args$nrates.traits <- evo.model.args$species.group.levels <-
    evo.model.args$nrates.species <- NULL
  nspecies <- length(tree$tip.label)
  nedge <- nrow(tree$edge)
  
  intraspecific <- FALSE # this option is fixed for now; future updates will allow for within-species observations
  
  if(class(fixed.effects)=="logical")
  {
    X <- numeric(0)
    X1 <- matrix(1,nspecies)
    rownames(X1) <- tree$tip.label
  } else
  {
    X <- fixed.effects
    if(class(X)=="data.frame")
    {
      species_col <- match(species.id,colnames(X))
      if(is.na(species_col)) stop("Name of species.id much match a column in X (if supplying a data frame).")
      trait_names <- colnames(X)
      if(trait_names[1]!=species.id)
      {
        X <- data.frame(X[,species_col],X[,-species_col,drop=FALSE])
      }
      colnames(X) <- c("species",trait_names[-species_col])
      if(!intraspecific)
      {
        new_X <- as.matrix(X[,2:ncol(X),drop=FALSE])
        rownames(new_X) <- as.character(X[,1])
        X <- new_X
      }
    } else if(class(X)=="matrix")
    {
      if(is.null(rownames(X)))
      {
        rownames(X) <- tree$tip.label
        warning("Row names of X do not match tip labels of tree. Assuming data are in order of tips.",immediate. = TRUE)
      }
      nms <- name.check(phy = tree,data.names = rownames(X))
      if(class(nms)=="list")
      {
        rownames(X) <- tree$tip.label
        warning("Row names of X do not match tip labels of tree. Assuming data are in order of tips.",immediate. = TRUE)
      }
      if(is.null(colnames(X))) colnames(X) <- paste("X",1:ncol(X),sep="")
    } else if(class(X)=="numeric")
    {
      if(is.null(names(X))) 
      {
        warning("No species names provided for X. Assuming data are in order of tips.",immediate. = TRUE)
        names(X) <- tree$tip.label
      } else
      {
        nms <- name.check(phy = tree,data.names = names(X))
        if(class(nms)=="list")
        {
          warning("Names of X do not match tip labels of tree. Assuming data are in order of tips.",immediate. = TRUE)
          names(X) <- tree$tip.label
        }
        X <- matrix(data = X,nrow = nspecies,ncol = 1,dimnames=list(names(X),"X1"))
      }
    }
    X <- X[tree$tip.label,,drop=FALSE]
    evo.model.args$fixed.effects <- X
    X1 <- cbind(1,X)
  }
  
  # tree should be in postorder
  # Y should be in matrix or data frame format; if traits are not named, they will be named V1, V2, ..., Vm
  # species.groups should be a vector of designating species group membership, named by species
  # trait.groups should be a vector of factors designating trait group membership, named by traits
  # model can be "BM" "OUfixedRoot" "OUrandomRoot" "lambda" "kappa" "delta"
  # diag.phylocov = FALSE (default); set to TRUE to assume traits are evolutionarily independent
  # method = "Pairwise REML" by default; also can use "Pairwise ML" or "Full ML"
  # force.zero.phylocov allows one to test for covariance between two groups of traits (similar to phylogenetic partial least squares)
  #### should be a vector of trait names for inclusion in Y1 (Y2 will be assumed to be the remainder of Y)
  # species.id is only used if Y is a data frame, and corresponds to the column name containing species names
  # max.combn is the maximum number of allowed pairwise combinations of Y. If ncombn>max.combn, Monte Carlo estimation is used
  # supply.regimes 
  # supply.full.cov
  # supply.pairwise.cov
  # supply.model
  
  # make sure data (Y) is properly formatted
  # if intraspecific==FALSE, matrices will be used
  if(class(Y)=="data.frame")
  {
    species_col <- match(species.id,colnames(Y))
    if(is.na(species_col)) stop("Name of species.id much match a column in Y (if supplying a data frame).")
    trait_names <- colnames(Y)
    if(trait_names[1]!=species.id)
    {
      Y <- data.frame(Y[,species_col],Y[,-species_col,drop=FALSE])
    }
    colnames(Y) <- c("species",trait_names[-species_col])
    if(!intraspecific)
    {
      new_Y <- as.matrix(Y[,2:ncol(Y),drop=FALSE])
      rownames(new_Y) <- as.character(Y[,1])
      Y <- new_Y
    }
  } else if(class(Y)=="matrix")
  {
    if(is.null(rownames(Y)))
    {
      rownames(Y) <- tree$tip.label
      warning("Row names of Y do not match tip labels of tree. Assuming data are in order of tips.",immediate. = TRUE)
    }
    nms <- name.check(phy = tree,data.names = rownames(Y))
    if(class(nms)=="list")
    {
      rownames(Y) <- tree$tip.label
      warning("Row names of Y do not match tip labels of tree. Assuming data are in order of tips.",immediate. = TRUE)
    }
    if(is.null(colnames(Y))) colnames(Y) <- paste("V",1:ncol(Y),sep="")
  } else if(class(Y)=="numeric")
  {
    if(is.null(names(Y))) 
    {
      warning("No species names provided for Y. Assuming data are in order of tips.",immediate. = TRUE)
      names(Y) <- tree$tip.label
    } else
    {
      nms <- name.check(phy = tree,data.names = names(Y))
      if(class(nms)=="list")
      {
        warning("Names of Y do not match tip labels of tree. Assuming data are in order of tips.",immediate. = TRUE)
        names(Y) <- tree$tip.label
      }
      Y <- matrix(data = Y,nrow = nspecies,ncol = 1,dimnames=list(names(Y),"V1"))
    }
  }
  if(all(complete.cases(Y))) missing_data <- FALSE
  if(ncol(Y)==1) if(method=="Pairwise ML") method <- "Full ML" else if(method=="Pairwise REML") method <- "Full REML"
  if(!intraspecific & !missing_data)
  {
    Y <- Y[tree$tip.label,,drop=FALSE]
    evo.model.args$Y <- Y
    ndims <- ncol(Y)
    diagndims <- diag(ndims)
    diag2 <- diag(2)
    if(method=="Full ML" | method=="Full REML")
    {
      MC <- FALSE
    } else
    {
      ncombn <- round(exp(lfactorial(ndims) - (log(2) + lfactorial(ndims-2)))) # number of combinations for pairwise composite likelihood
      if(ncombn<=max.combn)
      {
        MC <- FALSE # use Monte Carlo sampling for likelihood if ncombn is too high
        subsets <- matrix(0,2,ncombn)
        counter <- 1
        for(r in 1:(ndims-1))
        {
          for(s in (r+1):ndims)
          {
            subsets[,counter] <- c(r,s)
            counter <- counter + 1
          }
        }
      } else
      {
        MC <- TRUE
        subsets <- apply(matrix(0,max.combn,2),1,function(X) sample(x = ndims,size = 2,replace = FALSE))
      }
    }

    if(nrates.species==1) # no need for distance transformations (costly)
    {
      simple_BM <- function(temp_tree,ditree,ret_pars=FALSE)
      {
        pY$edge_len <- ditree$edge.length
        a_results <- do.call(multipic,pY)
        suminvV_and_logdV <- c(a_results[[2]],a_results[[3]])
        a <- a_results[[1]]
        if(length(X)>0)
        {
          pX$edge_len <- ditree$edge.length
          x_results <- do.call(multipic,pX)
          x <- x_results[[1]]
          XANC <- x_results$root[1,]
          XX <- crossprod(cbind(0,x)) + tcrossprod(c(1,XANC))*suminvV_and_logdV[1]
          betas <- solve(crossprod(x),crossprod(x,a))
          C <- crossprod(a-x%*%betas) / (nspecies-REML)
        } else C <- crossprod(a) / (nspecies-REML)
        colnames(C) <- rownames(C) <- colnames(Y)
        if(ret_pars)
          if(length(X>0))
          {
            XY <- crossprod(cbind(0,x),a) + tcrossprod(c(1,XANC),a_results$root[1,])*suminvV_and_logdV[1]
            BETAS <- solve(XX,XY)
            predicted <- X1%*%BETAS
            multipic.results <- list(X.multipic=x_results,Y.multipic=a_results)
          } else
          {
            predicted <- matrix(1,nspecies)%*%a_results$root[1,]
            multipic.results <- list(Y.multipic=a_results)            
          }
        if(diag.phylocov) C[upper.tri(C)] <- C[lower.tri(C)] <- 0
        diagC <- diag(C)
        if(nrates.traits>1)
        {
          temp_C <- C
          for(i in 1:nrates.traits)
          {
            group_i <- which(trait.groups==trait.group.levels[i])
            diag(temp_C)[group_i] <- mean(diagC[group_i]) * if(subset) 1 else length(group_i)
            diagC[group_i] <- diag(temp_C)[group_i]
          }
          temp_C <- matrix(nearPD(x = temp_C,corr = FALSE,keepDiag = TRUE)$mat,ncol(C),ncol(C))
          C <- temp_C
        } else if(multirate)
        {
          diag(C) <- mean(diagC) * if(subset) 1 else ndims
          C <- matrix(nearPD(x = C,corr = FALSE,keepDiag = TRUE)$mat,ncol(C),ncol(C))
          diagC <- diag(C)
        }
        if(length(force.zero.phylocov)>0)
        {
          other <- (1:ncol(C))[-match(force.zero.phylocov,colnames(C))]
          C[force.zero.phylocov,other] <- C[other,force.zero.phylocov] <- 0
        }
        
        if(method=="Full REML" | method=="Full ML")
        {
          const <- (nspecies*ndims-ndims*REML)*(log(2*pi)) # constant term for full log likelihood
          ndims_logd <- ndims*suminvV_and_logdV[2] # ndims * log determimant
          ndims_invV <- ndims*log(suminvV_and_logdV[1]) # ndims * sum of inverse of vcv(tree)
          log_Csub <- determinant(C,logarithm = TRUE)$modulus[[1]]
          if(length(X)>0)
          {
            YY <- sum(tcrossprod(a-x%*%betas,backsolve(chol(C),diagndims,transpose=TRUE))^2)
            if(class(YY)=="try-error") YY <- sum(diag(solve(C)%*%crossprod(a-x%*%betas)))
            logdetXWX <- determinant(XX)$modulus[[1]]*ncol(C) - determinant(C)$modulus[[1]]*(ncol(X)+1)
          } else 
          {
            YY <- try(sum(tcrossprod(a,backsolve(chol(C),diagndims,transpose=TRUE))^2),silent=TRUE)
            if(class(YY)=="try-error") YY <- sum(diag(solve(C)%*%crossprod(a)))
            logdetXWX <- (-log_Csub + ndims_invV)
          }
          logL <- -(const + (nspecies*log_Csub + ndims_logd) + REML*logdetXWX + YY)/2
        } else if(method=="Pairwise REML" | method=="Pairwise ML")
        {
          const <- (nspecies*2-2*REML)*(log(2*pi)) # constant term for pairwise log likelihood
          ndims_logd <- 2*suminvV_and_logdV[2] # 2 * log determimant
          ndims_invV <- 2*log(suminvV_and_logdV[1]) # 2 * sum of inverse of vcv(tree)
          logL <- 0
          if(MC) subsets <- apply(matrix(0,max.combn,2),1,function(X) sample(x = ndims,size = 2,replace = FALSE))
          for(i in 1:ncol(subsets))
          {
            r <- subsets[1,i]
            s <- subsets[2,i]
            Csub <- C[c(r,s),c(r,s)]
            log_Csub <- log(Csub[1,1]*Csub[2,2]-Csub[1,2]^2)
            if(length(X)>0)
            {
              try(YY <- sum(tcrossprod(a[,c(r,s)]-x%*%betas[,c(r,s)],backsolve(chol(Csub),diag2,transpose=TRUE))^2),silent=TRUE)
              if(class(YY)=="try-error") YY <- sum(diag(solve(Csub)%*%crossprod(a[,c(r,s)]-x%*%betas[,c(r,s)])))
              logdetXWX <- determinant(XX)$modulus[[1]]*2 - log_Csub*(ncol(X)+1)
            } else 
            {
              try(YY <- sum(tcrossprod(a[,c(r,s)],backsolve(chol(Csub),diag2,transpose=TRUE))^2),silent=TRUE)
              if(class(YY)=="try-error") YY <- sum(diag(solve(Csub)%*%crossprod(a[,c(r,s)])))
              logdetXWX <- (-log_Csub + ndims_invV)
            }
            logL <- logL - (const + (nspecies*log_Csub + ndims_logd) + REML*logdetXWX + YY)
          }
          logL <- logL/i*ncombn/2
        }
        if(ret_pars) return(list(phylocov=C,logL=logL,predicted=predicted,transf_tree=temp_tree,multipic.results=multipic.results)) else return(logL)
      }
    } else # multiple groups -- requires transformation of data
    {
      ind <- 0L:1L
      tempY <- matrix(0,nspecies*2,1)
      tempY_span <- 1:(nspecies*2)
      simple_BM <- function(temp_tree,ditree,ret_pars=FALSE)
      {
        pY$edge_len <- ditree$edge.length
        a_results <- do.call(multipic,pY)
        suminvV_and_logdV <- c(a_results[[2]],a_results[[3]])
        a <- a_results[[1]]
        Y_anc <- a_results$root[1,]
        if(length(X)>0)
        {
          pX$edge_len <- ditree$edge.length
          x_results <- do.call(multipic,pX)
          x <- x_results[[1]]
          XANC <- x_results$root[1,]
          YANC <- a_results$root[1,]
          XX <- crossprod(cbind(0,x)) + tcrossprod(c(1,XANC))*suminvV_and_logdV[1]
          betas <- solve(crossprod(x),crossprod(x,a))
          BETAS <- solve(XX,XY)
          XY <- crossprod(cbind(0,x),a) + tcrossprod(c(1,XANC),YANC)*suminvV_and_logdV[1]
          C <- crossprod(a-x%*%betas) / (nspecies-REML)
          anc <- X1%*%BETAS
        } else
        {
          anc <- matrix(1,nspecies)%*%a_results$root[1,]
          C <- crossprod(a) / (nspecies-REML)
        }
        colnames(C) <- rownames(C) <- colnames(Y)
        if(ret_pars)
          if(length(X>0))
          {
            predicted <- anc
            multipic.results <- list(X.multipic=x_results,Y.multipic=a_results)
          } else
          {
            predicted <- anc
            multipic.results <- list(Y.multipic=a_results)            
          }
        
        resid <- fast_transform(vcv(temp_tree)) %*% (Y-anc)
        rownames(resid) <- rownames(Y)
        species.subsets <- Clist <- diagClist <- vector("list",nrates.traits)
        
        for(j in 1:nrates.species)
        {
          species.subsets[[j]] <- which(!is.na(match(species.groups,species.group.levels[j])))
          Clist[[j]] <- t(resid[species.subsets[[j]],]) %*% resid[species.subsets[[j]],] / (length(species.subsets[[j]]) - REML)
          if(diag.phylocov) Clist[[j]][upper.tri(Clist[[j]])] <- Clist[[j]][lower.tri(Clist[[j]])] <- 0
          diagClist[[j]] <- diag(Clist[[j]])
          if(nrates.traits>1)
          {
            temp_C <- Clist[[j]]
            for(i in 1:nrates.traits)
            {
              group_i <- which(trait.groups==trait.group.levels[i])
              diag(temp_C)[group_i] <- mean(diagClist[[j]][group_i]) * if(subset) 1 else length(group_i)
            }
            temp_C <- matrix(nearPD(x = temp_C,corr = FALSE,keepDiag = TRUE)$mat,ncol(Clist[[j]]),ncol(Clist[[j]]))
            Clist[[j]] <- temp_C
            diagClist[[j]] <- diag(Clist[[j]])
          } else if(multirate)
          {
            diag(Clist[[j]]) <- mean(diagClist[[j]]) * if(subset) 1 else length(diagClist[[j]])
            Clist[[j]] <- matrix(nearPD(x = Clist[[j]],corr = FALSE,keepDiag = TRUE)$mat,ncol(Clist[[j]]),ncol(Clist[[j]]))
            diagClist[[j]] <- diag(Clist[[j]])
          }
          if(length(force.zero.phylocov)>0)
          {
            other <- which(is.na(match(colnames(Clist[[j]]),force.zero.phylocov)))
            Clist[[j]][force.zero.phylocov,other] <- Clist[[j]][other,force.zero.phylocov] <- 0
          }
          if(is.null(temp_tree$root.edge)) temp_tree$edge.length <- c(temp_tree$edge.length[1:nedge],0) else temp_tree$edge.length <- c(temp_tree$edge.length[1:nedge],temp_tree$root.edge)
        }
        
        if(method=="Full REML" | method=="Full ML")
        {
          if(length(X)>0)
          {
            ret <- Rinv5(nspecies = nspecies,nedge = nedge,edge.length = temp_tree$edge.length,anc = 
                           tree$edge[,1],des = tree$edge[,2],L = X1,R = Y,painted.edges = painted.edges,phylocov_list = Clist)
          } else
          {
            ret <- Rinv3(nspecies = nspecies,nedge = nedge,edge.length = temp_tree$edge.length,anc = 
                           tree$edge[,1],des = tree$edge[,2],R = Y,painted.edges = painted.edges,phylocov_list = Clist)
          }
          if(method=="Full REML") logL <- ret$ReLogLik else logL <- ret$LogLik
        } else if(method=="Pairwise REML" | method=="Pairwise ML")
        {
          logL <- 0
          if(MC) subsets <- apply(matrix(0,max.combn,2),1,function(X) sample(x = ndims,size = 2,replace = FALSE))
          anc_nodes <- temp_tree$edge[,1]-1
          des_nodes <- temp_tree$edge[,2]-1
          len_vec <- temp_tree$edge.length
          for(i in 1:ncol(subsets))
          {
            r <- subsets[1,i]
            s <- subsets[2,i]
            tempC <- lapply(Clist,function(X) X[c(r,s),c(r,s)])
            tempY[tempY_span] <- Y[,c(r,s)]
            tempCmat <- matrix(0,2*nrates.species,2)
            for(zz in 1:length(tempC))
            {
              tempCmat[1:2 + 2*(zz-1),] <- tempC[[zz]]
            }
            if(length(X)>0)
            {
              ret <- Rinv6(nspecies = nspecies,nedge = nedge,ngroups=nrates.species,ind=ind,len_vec = len_vec,anc = 
                                             anc_nodes,des = des_nodes,L = X1,R = tempY,painted_edges = painted.edges,phylocovs = tempCmat,REML=REML)
            } else
            {
              ret <- Rinv4(nspecies = nspecies,nedge = nedge,ngroups=nrates.species,ind=ind,len_vec = len_vec,anc = 
                                             anc_nodes,des = des_nodes,R = tempY,painted_edges = painted.edges,phylocovs = tempCmat,REML=REML)
            }
            
            logL <- logL+ret
          }
          logL <- logL/i*ncombn
        }
        if(ret_pars)
        {
          theta <- matrix(0,ncol(X1),ndims)
          for(i in 1:ndims)
            theta[,i] <- theta_Rinv6(nspecies = nspecies,nedge = nedge,ngroups=nrates.species,ind=0L,len_vec = temp_tree$edge.length,anc = 
                                                       tree$edge[,1]-1,des = tree$edge[,2]-1,L = X1,R = Y[,i,drop=FALSE],painted_edges = painted.edges,phylocovs = matrix(sapply(Clist,function(X) X[i,i])),REML=REML)
          predicted <- X1 %*% theta
        }
        if(ret_pars) return(list(phylocov=Clist,logL=logL,predicted=predicted,transf_tree=temp_tree)) else return(logL)
      }
    }
  } else
  {
    stop("Missing data and intraspecific observations are not yet supported.")
  }
  
  if(!is.binary.tree(tree))
  {
    binary <- TRUE
    ditree <- multi2di(tree,random = FALSE)
  } else
  {
    binary <- FALSE
    ditree <- tree
  }
  pY <- prep_multipic(Y,tree)
  if(length(X)>0) pX <- prep_multipic(X,tree)
  if(model!="BM")
  {
    ## Calculate bounds for evolutionary model parameters (adapted from phylolm package)
    if(model!="BM") if((model=="OUfixedRoot" | model=="OUrandomRoot") & !is.ultrametric(tree) & nrates.species>1) stop("Non-ultrametric trees not currently supported for OU model with multiple species groups.") 
    dis <- pruningwise.distFromRoot(reorder(tree,"pruningwise"))[1:nspecies]
    Tmax <- mean(dis)
    bounds.default <- matrix(c(1e-7/Tmax,50/Tmax,1e-7/Tmax,50/Tmax,1e-7,1,1e-6,1,1e-5,3,-3/Tmax,0), ncol=2, byrow=TRUE)
    rownames(bounds.default) <- c("alpha","alpha","lambda","kappa","delta","rate")
    colnames(bounds.default) <- c("min","max")
    starting.values.default <- c(0.5/Tmax,0.5/Tmax,0.5,0.5,0.5,-1/Tmax)
    names(starting.values.default) <- c("OUrandomRoot","OUfixedRoot","lambda","kappa","delta","EB")
    model_i <- match(model,names(starting.values.default))
    bounds.default <- bounds.default[model_i,,drop=FALSE]
    starting.values.default <- starting.values.default[model_i]
    names(starting.values.default) <- rownames(bounds.default)
    if(length(fixed.par)>0)
    {
      plot.LL.surface <- FALSE
      bounds.default[1:2] <- rep(fixed.par,2)
      starting.values.default[1] <- fixed.par
    }
    if(model!="EB" & model!="BM")
    {
      starting.values.default <- log(starting.values.default)
      bounds.default <- log(bounds.default)
    }
    if(model[1]=="BM") starting.values.default <- bounds.default <- NULL
    
    geiger_args <- list(x=tree,model=if(model=="OUfixedRoot" | model=="OUrandomRoot") "OU" else model,par=5)
    names(geiger_args)[3] <- names(starting.values.default)
    if(model=="EB") names(geiger_args)[3] <- "a"
    if(model=="OUrandomRoot" | model=="OUfixedRoot" | model=="OU")
    {
      temp_args <- list(phy = tree,model = model,parameters = list(pars=5))
      names(temp_args$parameters) <- names(starting.values.default)
    }
    model_LL <- function(par,ret_pars = FALSE)
    {
      if(model!="EB") par <- exp(par)
      if(model=="OUrandomRoot" | model=="OUfixedRoot" | model=="OU")
      {
        temp_args$parameters[[1]] <- par
        temp_tree <- do.call(transf.branch.lengths,temp_args)$tree
      } else
      {
        geiger_args[[3]] <- par
        temp_tree <- do.call(rescale,geiger_args)
      }
      if(!binary) temp_ditree <- multi2di(temp_tree,random=FALSE) else temp_ditree <- temp_tree
      ret <- try(simple_BM(temp_tree = temp_tree,ditree = temp_ditree,ret_pars = ret_pars),silent=TRUE)
      if(class(ret)=="try-error") return(-.Machine$double.xmax) else return(ret)
    }
    
    if(MC) try_length <- max(50,par.init.iters) else try_length <- par.init.iters
    if(length(fixed.par)>0) try_length <- 1
    if(model!="EB") try_theta <- log(seq(exp(min(bounds.default)),exp(max(bounds.default)),length=try_length)) else
      try_theta <- (seq((min(bounds.default)),(max(bounds.default)),length=try_length))
    try_LL <- double(length(try_theta))
    for(i in 1:length(try_theta))
    {
      try_LL[i] <- model_LL(try_theta[i])
    }
    
    if(plot.LL.surface & !MC) try(
      {
        plot(if(model!="EB") exp(try_theta) else try_theta,try_LL)
      },silent=TRUE)
    starting.values.default <- try_theta[which.max(try_LL)]
    
    if(MC & length(fixed.par)==0)
    {
      if(model!="EB") try_theta <- exp(try_theta)
      if(model!="EB") bounds.default <- exp(bounds.default)
      best <- double(round(length(try_theta)/2))
      for(i in 1:length(best))
      {
        best[i] <- BIC(lm(try_LL~poly(try_theta,degree = i)))
      }
      degree <- which.min(best)
      LL_surface <- lm(try_LL~poly(try_theta,degree = degree))
      theta_start <- try_theta[which(predict(LL_surface)==max(predict(LL_surface)))][1]
      o <- optim(par = theta_start,fn = function(X) predict.lm(LL_surface,newdata = data.frame(try_theta=X)),control=list(fnscale=-1),
                 method = "L-BFGS-B",lower = bounds.default[1],upper = bounds.default[2])
      if(plot.LL.surface) try(
        {
          plot(try_theta,try_LL)
          points(try_theta,predict.lm(LL_surface),col="blue")
          points(o$par,o$value,col="red",pch=16)
        },silent=TRUE)
      ret <- model_LL(par = if(model!="EB") log(o$par) else o$par,ret_pars = ret.level>1)
      ret$model.par[[rownames(bounds.default)]] <- o$par
      if(ret.level>1 & (model=="OUfixedRoot" | model=="OUrandomRoot" | model=="OU"))
      {
        if(nrates.species>1)
        {
          for(i in 1:nrates.species) ret$phylocov[[i]] <- 2 * ret$phylocov[[i]] * ret$alpha
        } else ret$phylocov <- 2 * ret$phylocov * ret$alpha
      }
    } else
    {
      if(length(fixed.par)>0)
      {
        o <- list(par=starting.values.default,value=model_LL(starting.values.default))
      } else  o <- optim(starting.values.default,model_LL,control=list(fnscale=-1),
                         method = "L-BFGS-B",lower = bounds.default[1],upper = bounds.default[2])
      ret <- model_LL(par = o$par,ret_pars = ret.level>1)
      if(model!="EB") o$par <- exp(o$par)
      if(ret.level>1)
      {
        ret$model.par[[rownames(bounds.default)]] <- o$par
        if(ret.level>1 & (model=="OUfixedRoot" | model=="OUrandomRoot" | model=="OU"))
        {
          if(nrates.species>1)
          {
            for(i in 1:nrates.species) ret$phylocov[[i]] <- 2 * ret$phylocov[[i]] * ret$model.par[[1]]
          } else ret$phylocov <- 2 * ret$phylocov * ret$model.par[[1]]
        }
      }
    }
    
  } else ret <- simple_BM(temp_tree = tree,ditree = ditree,ret_pars = ret.level>1)
  
  if(ret.level<=1) return(ret)
  if(ret.level>1)
  {
    ret$method <- method
    sigma2.mult <- vector("list",nrates.species)
    if(nrates.species>1)
    {
      for(i in 1:nrates.species)
      {
        if(nrates.traits<2) sigma2.mult[[i]] <- mean(diag(ret$phylocov[[i]])) else
        {
          for(j in 1:nrates.traits)
          {
            group_j <- which(trait.groups==trait.group.levels[j])
            sigma2.mult[[i]] <- c(sigma2.mult[[i]],mean(diag(ret$phylocov[[i]])[group_j]) * if(subset) 1 else length(group_j))
          }
          names(sigma2.mult[[i]]) <- levels(trait.groups)
        }
      }
      names(sigma2.mult) <- levels(species.groups)
    } else
    {
      if(nrates.traits<2) sigma2.mult <- mean(diag(ret$phylocov)) else
      {
        sigma2.mult <- double()
        for(j in 1:nrates.traits)
        {
          group_j <- which(trait.groups==trait.group.levels[j])
          sigma2.mult <- c(sigma2.mult,mean(diag(ret$phylocov)[group_j]) * if(subset) 1 else length(group_j))
        }
        names(sigma2.mult) <- levels(trait.groups)
      }
    }
    if(ret.level<3) ret$phylocov <- NULL
    ret <- c(ret,list(sigma2.mult=sigma2.mult,evo.model.args=evo.model.args))
  }
  class(ret) <- "evo.model"
  return(ret)
}

Rinv3 <- function(nspecies,nedge,edge.length,anc,des,L,R,painted.edges,phylocov_list)
{
  nvar <- ncol(phylocov_list[[1]])
  L1 <- LL <- R1 <- RR <- LR <- NA
  available_mats <- 0 # 0 indicates L and R are missing
  
  if(!missing(L))
  {
    if(!is.matrix(L))
    {
      if(is.data.frame(L))
      {
        L <- as.matrix(L)
      } else if(is.list(L))
      {
        stop("L should be in matrix or vector form.")
      } else if(length(L)==nspecies)
      {
        L <- matrix(L,ncol=1)
      } else stop("L should be in matrix or vector form.")
    }
    mL <- ncol(L)
    LL <- rep(list(matrix(0,mL,mL)),nedge+1)
    L1 <- rep(list(matrix(0,mL,nvar)),nedge+1)
    available_mats <- 1 # 1 indicates only L is present
  }
  
  
  # make sure R is a matrix
  if(!missing(R))
  {
    if(!is.matrix(R))
    {
      if(is.data.frame(R))
      {
        R <- as.matrix(R)
      } else if(is.list(R))
      {
        stop("R should be in matrix or vector form.")
      } else if(length(R)==(nspecies*nvar))
      {
        R <- matrix(R,ncol=1)
      } else stop("R should be in matrix or vector form.")
    } else R <- matrix(R,ncol=1)
    mR <- ncol(R)
    RR <- rep(list(matrix(0,mR,mR)),nedge+1)
    R1 <- rep(list(matrix(0,nvar,mR)),nedge+1)
    available_mats <- available_mats + 2 # 2 indicates only R is present; 3 indicates both
    if(available_mats==3) LR <- rep(list(matrix(0,mL,mR)),nedge+1)
  }
  
  
  logd <- vector("double",nedge+1)
  p <- pA <- rep(list(matrix(0,nvar,nvar)),nedge+1)
  pA_des_value <- matrix(0,nvar,nvar)
  i <- des_node <- anc_node <- 0
  ngroups <- ncol(painted.edges)
  len <- matrix(0,nvar,nvar)
  for(i in 1:nedge)
  {
    des_node <- des[i]
    anc_node <- anc[i]
    len <- len*0
    for(z in 1:ngroups) len <- len + painted.edges[i,z]*edge.length[i]*phylocov_list[[z]] 
    
    if(des_node<=nspecies)
    {      
      logd[des_node] <- determinant(len,logarithm = TRUE)$modulus[[1]]
      logd[anc_node] <- logd[anc_node] + logd[des_node]
      p[[des_node]] <- solve(len)
      pA[[anc_node]] <- pA[[anc_node]] + p[[des_node]]
      index <- des_node + nspecies*(1:nvar-1)
      if(available_mats>0 & available_mats!=2)
      {
        L1[[des_node]] <- t(L[des_node,,drop=FALSE]) %*% p[[des_node]]
        LL[[des_node]] <- t(L[des_node,,drop=FALSE]) %*%  p[[des_node]] %*% L[des_node,,drop=FALSE]
        
        LL[[anc_node]] <- LL[[anc_node]] + LL[[des_node]]
        L1[[anc_node]] <- L1[[anc_node]] + L1[[des_node]]
      }
      if(available_mats>1)
      {
        R1[[des_node]] <- p[[des_node]] %*% R[index,,drop=FALSE]
        RR[[des_node]] <- t(R[index,,drop=FALSE]) %*% p[[des_node]] %*% R[index,,drop=FALSE] # Q
        
        RR[[anc_node]] <- RR[[anc_node]] + RR[[des_node]]
        R1[[anc_node]] <- R1[[anc_node]] + R1[[des_node]]
      }
      if(available_mats>2)
      {
        LR[[des_node]] <- t(L[index,,drop=FALSE]) %*% p[[des_node]] %*% R[index,,drop=FALSE]
        LR[[anc_node]] <- LR[[anc_node]] + LR[[des_node]]
      }
    } else
    {
      pA_des_value <- pA[[des_node]]
      itpa <- diag(nvar)+len%*%pA_des_value
      itpainv <- solve(itpa)
      
      logd[des_node] <- logd[des_node] + determinant(itpa,logarithm = TRUE)$modulus[[1]]
      logd[anc_node] <- logd[anc_node] + logd[des_node]
      
      p[[des_node]] <- pA_des_value %*% itpainv
      pA[[anc_node]] <- pA[[anc_node]] + p[[des_node]]
      
      if(available_mats>2)
      {
        LR[[des_node]] <- LR[[des_node]] - L1[[des_node]] %*% itpainv %*% len %*% (R1[[des_node]])
        LR[[anc_node]] <- LR[[anc_node]] + LR[[des_node]]
      }
      if(available_mats>0 & available_mats!=2)
      {
        LL[[des_node]] <- LL[[des_node]] - L1[[des_node]] %*% itpainv %*% len %*% t(L1[[des_node]])
        LL[[anc_node]] <- LL[[anc_node]] + LL[[des_node]]
        
        L1[[des_node]] <- L1[[des_node]] %*% itpainv
        L1[[anc_node]] <- L1[[anc_node]] + L1[[des_node]]
      }
      if(available_mats>1)
      {
        RR[[des_node]] <- RR[[des_node]] - t(R1[[des_node]]) %*% itpainv %*% len %*% (R1[[des_node]])
        RR[[anc_node]] <- RR[[anc_node]] + RR[[des_node]]
        
        R1[[des_node]] <- t(t(R1[[des_node]]) %*% itpainv)
        R1[[anc_node]] <- R1[[anc_node]] + R1[[des_node]] #t(t(R1[[des_node]]) %*% pAinv)
      }
      
    }
  }
  i <- nedge+1
  len <- len*0
  for(z in 1:ngroups) len <- len + painted.edges[i,z]*edge.length[i]*phylocov_list[[z]] 
  
  pA_des_value <- pA[[anc_node]]
  itpa <- diag(nvar)+len%*%pA_des_value
  itpainv <- solve(itpa)
  logd[nedge+1] <- logd[nspecies+1] + determinant(itpa,logarithm = TRUE)[[1]]
  p[[nedge+1]] <- pA_des_value %*% itpainv
  
  ret <- list(p=p[[nedge+1]],logd=logd[nedge+1])
  
  if(available_mats>2)
  {
    LR[[nedge+1]] <- LR[[anc_node]] - L1[[anc_node]] %*% itpainv %*% len %*% (R1[[anc_node]])
    ret$LR <- LR[[nedge+1]]
  }  
  if(available_mats>0 & available_mats!=2)
  {
    LL[[nedge+1]] <- LL[[anc_node]] - L1[[anc_node]] %*% itpainv %*% len %*% t(L1[[anc_node]])    
    L1[[nedge+1]] <- L1[[anc_node]] %*% itpainv
    
    ret$LL <- LL[[nedge+1]]
    ret$L1 <- L1[[nedge+1]]
  }
  if(available_mats>1)
  {
    RR[[nedge+1]] <- RR[[anc_node]] - t(R1[[anc_node]]) %*% itpainv %*% len %*% (R1[[anc_node]])
    R1[[nedge+1]] <- t(t(R1[[anc_node]]) %*% itpainv)
    
    ret$R1 <- R1[[nedge+1]]
    ret$RR <- RR[[nedge+1]]
    if(available_mats==2)
    {
      ret$theta <- solve(ret$p,ret$R1)
      ret$LogLik <- -(nspecies*nvar*log(2*pi) + ret$logd + ret$RR - 2*t(ret$theta) %*% ret$R1 + t(ret$theta) %*% ret$p %*% ret$theta)/2
      ret$ReLogLik <- -((nspecies*nvar-nvar)*log(2*pi) + determinant(ret$p)$modulus[[1]] + ret$logd + ret$RR - 2*t(ret$theta) %*% ret$R1 + t(ret$theta) %*% ret$p %*% ret$theta)/2
    } else
    {
      ret$theta <- solve(ret$LL,ret$LR)
      ret$LogLik <- -(nspecies*nvar*log(2*pi) + ret$logd + ret$RR - 2*t(ret$theta) %*% ret$LR + t(ret$theta) %*% ret$LL %*% ret$theta)/2
      ret$ReLogLik <- -((nspecies*nvar-nvar)*log(2*pi) + determinant(ret$LL)$modulus[[1]] + ret$logd + ret$RR - 2*t(ret$theta) %*% ret$LR + t(ret$theta) %*% ret$LL %*% ret$theta)/2
    }
  }
  ret
}

pre_parallel <- function()
{
  if (!is.na(Sys.getenv()["R_GUI_APP_VERSION"])[[1]]) {
    return(FALSE)
  }
  tmp <- rownames(installed.packages())
  if (("doParallel" %in% tmp | "doSNOW" %in% tmp) & "foreach" %in% tmp)
  {
    if(Sys.info()["sysname"] == "Windows" & "doSNOW" %in% tmp)
    {
      requireNamespace("doSNOW")
      return(TRUE)
    } else if("doParallel" %in% tmp)
    {
      requireNamespace("doParallel")
      return(TRUE)
    }
  } else return(FALSE)
}

multi.distFromRoot <- function (phy,edge.lengths)
{
  nt = length(phy$tip.label)
  nvar <- nrow(edge.lengths)
  xx <- matrix(0,nvar,phy$Nnode + nt)
  for (i in nrow(phy$edge):1) xx[,phy$edge[i, 2]] <- xx[,phy$edge[i,1]] + edge.lengths[,i]
  colnames(xx) <- if (is.null(phy$node.label)) 
    1:(nt + phy$Nnode) else c(phy$tip.label, phy$node.label)
  if(any(edge.lengths[,nrow(phy$edge)+1]!=0)) xx <- xx + edge.lengths[,nrow(phy$edge)+1] %*% matrix(0,ncol=dim(xx)[2])
  return(xx[,1:nt])
}

K.mult <- function(model,nsim=1000,plot=TRUE)
{
  estimate_power <- TRUE
  if(!is.binary.tree(model$evo.model.args$tree)) FLAG <- TRUE else FLAG <- FALSE
  nvar <- ncol(model$evo.model.args$Y)
  null_model <- model
  null_model$evo.model.args$tree <- rescale(null_model$evo.model.args$tree,model = "lambda",lambda=0)
  sim_null <- sim.model(model = null_model,nsim = nsim)
  sim_alt <- sim.model(model = model,nsim = nsim)
  
  args <- model$evo.model.args
  args$ret.level <- 3
  if(args$model!="BM") args$fixed.par <- model1$model.par
  model1 <- do.call(what = evo.model,args = args)
  
  BMtree <- model$evo.model.args$tree
  tree <- model1$transf_tree
  nspecies <- length(tree$tip.label)
  nedge <- nrow(tree$edge)
  if(class(model$evo.model.args$species.groups)=="factor") gps <- TRUE else gps <- FALSE
  if(length(tree$edge.length)==nrow(tree$edge))
    if(!is.null(tree$root.edge)) tree$edge.length <- c(tree$edge.length[1:nedge],tree$root.edge) else tree$edge.length <- c(tree$edge.length[1:nedge],0)
  model$evo.model.args$ret.level <- 3
  new_edge <- matrix(0,nvar,nrow(tree$edge)+1)
  if(nrow(new_edge)>1) new_trees <- rep(tree,nrow(new_edge)) else new_trees <- list(tree)
  original_heights <- pruningwise.distFromRoot(reorder(tree,"pruningwise"))[1:nspecies]
  new_heights <- matrix(0,nvar,nspecies)
  sum_original_heights <- sum(original_heights)
  des_nodes <- tree$edge[,2]-1
  anc_nodes <- tree$edge[,1]-1
  if(model$evo.model.args$model=="BM") model$evo.model.args$max.combn <- 1
  pY <- prep_multipic2(model$evo.model.args$Y,phy = model$evo.model.args$tree,edge_len_mat = t(new_edge))
  if(class(model$evo.model.args$fixed.effects)=="matrix")
  {
    pX <- prep_multipic(model$evo.model.args$fixed.effects,phy = model$evo.model.args$tree)
  }
  
  calc_denom <- function(Y,fixed_effects)
  {
    model$evo.model.args$Y <- Y
    if(class(model$evo.model.args$fixed.effects)=="matrix") model$evo.model.args$fixed.effects <- fixed_effects
    model <- do.call(evo.model,model$evo.model.args)
    
    if(gps)
    {
      for(i in 1:length(model$phylocov))
      {
        new_edge <- new_edge + diag(model$phylocov[[i]]) %*% (t(tree$edge.length)*model$evo.model.args$painted.edges[,i])
      }
    } else new_edge <- diag(model$phylocov) %*% t(tree$edge.length)
    new_heights <- multi.distFromRoot(tree,new_edge)
    new_edge <- new_edge * t(matrix(1,dim(new_edge)[2]) %*% (mean(original_heights) / rowMeans(new_heights)))
    new_heights <- multi.distFromRoot(tree,new_edge)
    pY$edge_len <- t(new_edge)
    pY$phe[1:nspecies,] <- Y
    pY_results <- do.call(multipic2,pY)
    
    if(class(model$evo.model.args$fixed.effects)=="matrix")
    {
      pX$phe[1:nspecies,] <- model$evo.model.args$fixed.effects
      pX_results <- apply(new_edge,1,function(X) 
      {
        pX$edge_len <- X
        do.call(multipic,pX)
      })
    }
    for(i in 1:nvar)
    {
      new_trees[[i]]$edge.length <- new_edge[i,1:nrow(tree$edge)]
    }
    
    #scale_res <- apply(new_heights,1,function(X) X/mean(X))
    scale_res <- apply(new_heights,1,function(X) X/original_heights)
    num <- sum(diag(crossprod((Y-model$predicted)/scale_res,Y-model$predicted)))
    sum_new_heights <- rowSums(new_heights)
    invsums <- pY_results$sum_invV[1,]
    if(class(model$evo.model.args$fixed.effects)!="matrix")
    {
      ratio <- num / sum(diag(crossprod(pY_results$contrasts)))
    } else
    {
      ret_list <- vector("list")
      for(i in 1:nvar) ret_list[[i]] <- diag(crossprod(pY_results$contrasts[,i] - pX_results[[i]]$contrasts %*% solve(crossprod(pX_results[[i]]$contrasts),crossprod(pX_results[[i]]$contrasts,pY_results$contrasts[,i]))))
      denom <- sum(simplify2array(ret_list))
      ratio <- num/denom
    }
    return(ratio)
  }
  
  denom_ratio <- numeric(nsim)
  
  for(i in 1:nsim)
  {
    denom_ratio[i] <- if(class(model$evo.model.args$fixed.effects)!="matrix")
      calc_denom(Y = sim_alt$trait_data[[i]]) else
        calc_denom(Y = sim_alt$trait_data[[i]],fixed_effects = sim_alt$fixed_effectes[[i]])
  }
  
  e_ratio <- mean(denom_ratio)
  
  get_K <- function(Y,fixed_effects)
  {
    model$evo.model.args$Y <- Y
    if(class(model$evo.model.args$fixed.effects)=="matrix") model$evo.model.args$fixed.effects <- fixed_effects
    model <- do.call(evo.model,model$evo.model.args)
    
    if(gps)
    {
      for(i in 1:length(model$phylocov))
      {
        new_edge <- new_edge + diag(model$phylocov[[i]]) %*% (t(tree$edge.length)*model$evo.model.args$painted.edges[,i])
      }
    } else new_edge <- diag(model$phylocov) %*% t(tree$edge.length)
    new_heights <- multi.distFromRoot(tree,new_edge)
    new_edge <- new_edge * t(matrix(1,dim(new_edge)[2]) %*% (mean(original_heights) / rowMeans(new_heights)))
    new_heights <- multi.distFromRoot(tree,new_edge)
    pY$edge_len <- t(new_edge)
    pY$phe[1:nspecies,] <- Y
    pY_results <- do.call(multipic2,pY)
    
    if(class(model$evo.model.args$fixed.effects)=="matrix")
    {
      pX$phe[1:nspecies,] <- model$evo.model.args$fixed.effects
      pX_results <- apply(new_edge,1,function(X) 
      {
        pX$edge_len <- X
        do.call(multipic,pX)
      })
    }
    for(i in 1:nvar)
    {
      new_trees[[i]]$edge.length <- new_edge[i,1:nrow(tree$edge)]
    }
    
    #scale_res <- apply(new_heights,1,function(X) X/mean(X))
    scale_res <- apply(new_heights,1,function(X) X/original_heights)
    num <- sum(diag(crossprod((Y-model$predicted)/scale_res,Y-model$predicted)))
    sum_new_heights <- rowSums(new_heights)
    invsums <- pY_results$sum_invV[1,]
    if(class(model$evo.model.args$fixed.effects)!="matrix")
    {
      K <- (num / sum(diag(crossprod(pY_results$contrasts)))) / e_ratio
    } else
    {
      ret_list <- vector("list")
      for(i in 1:nvar) ret_list[[i]] <- diag(crossprod(pY_results$contrasts[,i] - pX_results[[i]]$contrasts %*% solve(crossprod(pX_results[[i]]$contrasts),crossprod(pX_results[[i]]$contrasts,pY_results$contrasts[,i]))))
      denom <- sum(simplify2array(ret_list))
      K <- (num/denom) / e_ratio
    }
    K
  }
  nullK <- double(nsim)
  altK <- double(nsim)
  
  
    suppressWarnings
    {
    cat("\nBootstrapping under null model.\n")
      for(i in 1:nsim) nullK[i] <- if(class(model$evo.model.args$fixed.effects)=="matrix") get_K(sim_null$trait_data[[i]],sim_null$fixed_effects[[i]]) else get_K(sim_null$trait_data[[i]])
      altK <- denom_ratio / e_ratio
    }
  
  K <- if(class(model$evo.model.args$fixed.effects)=="matrix") get_K(model$evo.model.args$Y,model$evo.model.args$fixed.effects) else
    get_K(model$evo.model.args$Y)
  critical_K <- sort(nullK)[min(round(nsim*.95)+1,nsim)]
  if(plot)
  {
    if(estimate_power)
    {
      xlim <- range(c(0,1.2,min(max(altK),quantile(ecdf(altK),.99)[[1]]),min(max(altK),quantile(ecdf(nullK),.99)[[1]]),K,critical_K))
      d1 <- density(altK)
      d2 <- density(nullK)
      ylim <- c(range(0,range(d1$y),range(d2$y)))
      #plot(d1,xlim=xlim,ylim=ylim,xlab="K",main="Null vs Alt. Distribution")
      #polygon(d1,col=rgb(0,0,1,.5))
      #lines(d2)
      #polygon(d2,col=rgb(.1,.1,.1,.5))
    } else
    {
      xlim <- range(c(0,1.2,c(range(nullK)),K,critical_K))
      d2 <- density(nullK)
      #plot(d2,xlab="K",main="Null Distribution",xlim=xlim)
      #polygon(d2,col=rgb(.1,.1,.1,.5))
    }
    #abline(v=K,lwd=5,lty=3)
    #abline(v=critical_K,col="red",lwd=5,lty=3)
  }
  ret <- list(K = K,Pval=length(which(K<=nullK))/nsim)
  if(estimate_power) ret <- c(ret,power=length(which(altK>=critical_K)) / nsim,list(K.expectation = mean(altK)))
  plot.tools <- list(test_statistic=K,critical_test_statistic=critical_K,test_statistic_name="K",null_sim_test_statistic=nullK)
  if(estimate_power) plot.tools <- c(plot.tools,list(alt_sim_test_statistic=altK))
  ret$plot.tools <- plot.tools
  class(ret) <- "compare.model"
  if(plot) plot(ret)
  ret
}

sim.model <- function(model,nsim=1000,return.type="matrix")
{
  FLAG <- FALSE
  if(nsim==1)
  {
    FLAG <- TRUE
    nsim <- 2
  }
  model1 <- model
  model1SAVE <- model1
  perm_fixed_par_1 <- model1$evo.model.args
  nvar <- ncol(model1$evo.model.args$Y)
  if(is.null(model1$phylocov))
  {
    args <- model1$evo.model.args
    args$ret.level <- 3
    if(args$model!="BM") args$fixed.par <- model1$model.par
    if(class(perm_fixed_par_1$fixed.effects)=="matrix")
    {
      args$Y <- cbind(args$Y,args$fixed.effects)
      args$fixed.effects <- NA
      colnames(args$Y)[(nvar+1):ncol(args$Y)] <- paste("temp_fixed_effect_",(nvar+1):ncol(args$Y),sep="")
      model1 <- do.call(what = evo.model,args = args)
    } else model1 <- do.call(what = evo.model,args = args)
  }
  tree <- model1$evo.model.args$tree
  fixed_effect_1 <- rep(list(NA),nsim)
  if(is.list(model1$phylocov))
  {
    sim1 <- sim.groups(tree = tree,groups = model1$evo.model.args$species.groups,
                                       painted.edges = model1$evo.model.args$painted.edges,model = model1$evo.model.args$model,
                                       parameters=if(model1$evo.model.args$model=="BM") list("BM") else 
                                         as.list(model1$evo.model.args$fixed.par),
                                       phylocov = model1$phylocov,nsim = nsim,return.type=return.type)
  } else
  {
    sim1 <- sim.traits(nreps=1,nmissing=0,tree=tree,v=model1$phylocov,
                                      model=model1$evo.model.args$model,
                                      parameters=if(model1$evo.model.args$model=="BM") list("BM") else 
                                        as.list(model1$evo.model.args$fixed.par),nsim=nsim,return.type=return.type)
  }
  
  if(class(perm_fixed_par_1$fixed.effects)=="matrix")
  {
    for(i in 1:nsim)
    {
      fixed_effect_1[[i]] <- as.matrix(sim1$trait_data[[i]][,(nvar+1):ncol(args$Y)+is.data.frame(sim1$trait_data[[i]])])
      sim1$trait_data[[i]] <- sim1$trait_data[[i]][,1:(nvar+is.data.frame(sim1$trait_data[[i]]))]
      rownames(fixed_effect_1[[i]]) <- tree$tip.label
      colnames(fixed_effect_1[[i]]) <- colnames(model1SAVE$evo.model.args$fixed.effects)
    }
    if(FLAG) fixed_effect_1 <- fixed_effect_1[1]
  }
  if(return.type=="matrix") 
  {
    for(i in 1:nsim) 
      colnames(sim1$trait_data[[i]]) <- colnames(model1SAVE$evo.model.args$Y)
  } else
  {
    for(i in 1:nsim) 
      colnames(sim1$trait_data[[i]]) <- c(model1SAVE$evo.model.args$species.id,colnames(model1SAVE$evo.model.args$Y))
  }
  
  if(FLAG) sim1$trait_data <- sim1$trait_data[1]
  if(class(perm_fixed_par_1$fixed.effects)=="matrix")
  {
    return(list(trait_data=sim1$trait_data,fixed_effects=fixed_effect_1))
  } else
  {
    return(list(trait_data=sim1$trait_data))
  }
}

compare.models <- function(model1,model2,nsim=1000,plot=TRUE,estimate_power=TRUE,parallel=TRUE,conf.int=TRUE)
{
  model1SAVE <- model1
  model2SAVE <- model2
  perm_fixed_par_1 <- model1$evo.model.args
  perm_fixed_par_2 <- model2$evo.model.args
  conf_int <- (conf.int & (!is.null(model2$model.par) & length(model2$evo.model.args$fixed.par)==0))
  nvar <- ncol(model1$evo.model.args$Y)
  if(conf_int) estimate_power <- TRUE
  if(estimate_power & !is.null(model2$model.par)) conf_int <- TRUE
  if(is.null(model1$phylocov))
  {
    args <- model1$evo.model.args
    args$ret.level <- 3
    if(args$model!="BM") args$fixed.par <- model1$model.par
    if(class(perm_fixed_par_1$fixed.effects)=="matrix")
    {
      args$Y <- cbind(args$Y,args$fixed.effects)
      args$fixed.effects <- NA
      colnames(args$Y)[(nvar):ncol(args$Y)] <- paste("temp_fixed_effect_",(nvar+1):ncol(args$Y),sep="")
      model1 <- do.call(what = evo.model,args = args)
    } else model1 <- do.call(what = evo.model,args = args)
  }
  if(is.null(model2$phylocov))
  {
    args <- model2$evo.model.args
    args$ret.level <- 3
    if(args$model!="BM") args$fixed.par <- model2$model.par
    if(class(perm_fixed_par_2$fixed.effects)=="matrix")
    {
      args$Y <- cbind(args$Y,args$fixed.effects)
      args$fixed.effects <- NA
      colnames(args$Y)[(nvar):ncol(args$Y)] <- paste("temp_fixed_effect_",(nvar+1):ncol(args$Y),sep="")
      model2 <- do.call(what = evo.model,args = args)
    } else model2 <- do.call(what = evo.model,args = args)
  }
  tree <- model1$evo.model.args$tree
  fixed_effect_1 <- fixed_effect_2 <- rep(list(NA),nsim)
  if(is.list(model1$phylocov))
  {
    sim1 <- sim.groups(tree = tree,groups = model1$evo.model.args$species.groups,
                                       painted.edges = model1$evo.model.args$painted.edges,model = model1$evo.model.args$model,
                                       parameters=if(model1$evo.model.args$model=="BM") list("BM") else 
                                         as.list(model1$evo.model.args$fixed.par),
                                       phylocov = model1$phylocov,nsim = nsim)
  } else
  {
    sim1 <- sim.traits(nreps=1,nmissing=0,tree=tree,v=model1$phylocov,
                                      model=model1$evo.model.args$model,
                                      parameters=if(model1$evo.model.args$model=="BM") list("BM") else 
                                        as.list(model1$evo.model.args$fixed.par),nsim=nsim)
  }
  
  if(class(perm_fixed_par_1$fixed.effects)=="matrix")
  {
    for(i in 1:nsim)
    {
      fixed_effect_1[[i]] <- as.matrix(sim1$trait_data[[i]][,(nvar+1):ncol(args$Y)+is.data.frame(sim1$trait_data[[i]])])
      sim1$trait_data[[i]] <- sim1$trait_data[[i]][,1:(nvar+is.data.frame(sim1$trait_data[[i]]))]
      rownames(fixed_effect_1[[i]]) <- tree$tip.label
    }
  }
  
  if(estimate_power)
  {
    if(is.list(model2$phylocov))
    {
      sim2 <- sim.groups(tree = tree,groups = model2$evo.model.args$species.groups,
                                         painted.edges = model2$evo.model.args$painted.edges,model = model2$evo.model.args$model,
                                         parameters=if(model2$evo.model.args$model=="BM") list("BM") else as.list(model2$evo.model.args$fixed.par),
                                         phylocov = model2$phylocov,nsim = nsim)
    } else
    {
      sim2 <- sim.traits(nreps=1,nmissing=0,tree=tree,v=model2$phylocov,
                                        model=model2$evo.model.args$model,
                                        parameters=if(model2$evo.model.args$model=="BM") list("BM") else as.list(model2$evo.model.args$fixed.par),nsim=nsim)
      
    }
    if(class(perm_fixed_par_2$fixed.effects)=="matrix")
    {
      for(i in 1:nsim)
      {
        fixed_effect_2[[i]] <- as.matrix(sim2$trait_data[[i]][,(nvar+1):ncol(args$Y)+is.data.frame(sim2$trait_data[[i]])])
        sim2$trait_data[[i]] <- sim2$trait_data[[i]][,1:(nvar+is.data.frame(sim2$trait_data[[i]]))]
        rownames(fixed_effect_2[[i]]) <- tree$tip.label
      }
    }
  } else sim2 <- sim1
  
  for(i in 1:nsim) colnames(sim1$trait_data[[i]]) <- colnames(sim2$trait_data[[i]]) <- colnames(model1SAVE$evo.model.args$Y)
  args1 <- perm_fixed_par_1
  args2 <- perm_fixed_par_2
  args1$ret.level <- args2$ret.level <- 1
  args1$fixed.par <- perm_fixed_par_1$fixed.par
  args2$fixed.par <- perm_fixed_par_2$fixed.par
  args1$plot.LL.surface <- args2$plot.LL.surface <- FALSE
  par_try <- character()
  class(par_try) <- "try-error"
  if(parallel)
  {
    par_try <- try(
      {
        if(pre_parallel())
        {
          ncores <- detectCores()
          cl <- if(Sys.info()["sysname"] == "Windows") makeCluster(ncores) else makeCluster(ncores)
          cat("registering...")
          if(Sys.info()["sysname"] == "Windows") registerDoSNOW(cl) else registerDoParallel(cl)
          cat("\nBootstrapping under null model.\n")
          null_sim_LR <- simplify2array(foreach(i=1:nsim) %dopar%
          {
            args1$Y <- sim1$trait_data[[i]]
            args2$Y <- sim1$trait_data[[i]]
            args1$fixed.effects <- fixed_effect_1[[i]]
            args2$fixed.effects <- if(class(fixed_effect_1[[i]])=="matrix" & class(fixed_effect_2[[i]])=="matrix") fixed_effect_1[[i]] else perm_fixed_par_2$fixed.effects
            -2*(do.call(phylocurve::evo.model,args1) - do.call(phylocurve::evo.model,args2))
          })
          if(estimate_power)
          {
            cat("Bootstrapping under alternative model.\n\n")
            if(conf_int)
            {
              alt_sim <- foreach(i=1:nsim) %dopar%
              {
                args1$Y <- sim2$trait_data[[i]]
                args2$Y <- sim2$trait_data[[i]]
                args1$fixed.effects <- if(class(fixed_effect_1[[i]])=="matrix" & class(fixed_effect_2[[i]])=="matrix") fixed_effect_2[[i]] else perm_fixed_par_1$fixed.effects
                args2$fixed.effects <- fixed_effect_2[[i]]
                args2$ret.level <- 2
                null_LL <- do.call(phylocurve::evo.model,args1)
                alt_LL <- do.call(phylocurve::evo.model,args2)
                list(D=-2*(null_LL - alt_LL$logL),model.par.sim=alt_LL$model.par[[1]])
              }
              alt_sim_LR <- sapply(alt_sim,function(X) X[[1]])
              model.par.sim <- sapply(alt_sim,function(X) X[[2]])
            } else
            {
              alt_sim_LR <- simplify2array(foreach(i=1:nsim) %dopar%
              {
                args1$Y <- sim2$trait_data[[i]]
                args2$Y <- sim2$trait_data[[i]]
                args1$fixed.effects <- if(class(fixed_effect_1[[i]])=="matrix" & class(fixed_effect_2[[i]])=="matrix") fixed_effect_2[[i]] else perm_fixed_par_1$fixed.effects
                args2$fixed.effects <- fixed_effect_2[[i]]
                -2*(do.call(phylocurve::evo.model,args1) - do.call(phylocurve::evo.model,args2))
              })
            }
          }
        }
      },silent=TRUE)
    try(if(Sys.info()["sysname"] == "Windows") stopCluster(cl) else stopCluster(cl),silent=TRUE)
  }
  if(!parallel | class(par_try)=="try-error")
  {
    cat("\nBootstrapping under null model.\n")
    null_sim_LR <- double(nsim)
    for(i in 1:nsim)
    {
      args1$Y <- sim1$trait_data[[i]]
      args2$Y <- sim1$trait_data[[i]]
      args1$fixed.effects <- fixed_effect_1[[i]]
      args2$fixed.effects <- if(class(fixed_effect_1[[i]])=="matrix" & class(fixed_effect_2[[i]])=="matrix") fixed_effect_1[[i]] else perm_fixed_par_2$fixed.effects
      null_sim_LR[i] <- -2*(do.call(evo.model,args1) - do.call(evo.model,args2))
    }
    if(estimate_power)
    {
      cat("Bootstrapping under alternative model.\n\n")
      if(conf_int)
      {
        alt_sim <- vector("list",nsim)
        suppressWarnings(
          {
            for(i in 1:nsim)
            {
              args1$Y <- sim2$trait_data[[i]]
              args2$Y <- sim2$trait_data[[i]]
              args2$ret.level <- 2
              args1$fixed.effects <- if(class(fixed_effect_1[[i]])=="matrix" & class(fixed_effect_2[[i]])=="matrix") fixed_effect_2[[i]] else perm_fixed_par_1$fixed.effects
              args2$fixed.effects <- fixed_effect_2[[i]]
              null_LL <- do.call(evo.model,args1)
              alt_LL <- do.call(evo.model,args2)
              alt_sim[[i]] <- list(D=-2*(null_LL - alt_LL$logL),model.par.sim=alt_LL$model.par[[1]])
            }
          })
        alt_sim_LR <- sapply(alt_sim,function(X) X[[1]])
        model.par.sim <- sapply(alt_sim,function(X) X[[2]])
      } else
      {
        alt_sim_LR <- double(nsim)
        suppressWarnings(
          {
            for(i in 1:nsim)
            {
              args1$Y <- sim2$trait_data[[i]]
              args2$Y <- sim2$trait_data[[i]]
              args1$fixed.effects <- if(class(fixed_effect_1[[i]])=="matrix" & class(fixed_effect_2[[i]])=="matrix") fixed_effect_2[[i]] else perm_fixed_par_1$fixed.effects
              args2$fixed.effects <- fixed_effect_2[[i]]
              alt_sim_LR[i] <- -2*(do.call(evo.model,args1) - do.call(evo.model,args2))
            }
          })
      }
    }
  }
  delta <- as.double(-2*(model1SAVE$logL-model2SAVE$logL))
  critical_delta <- sort(null_sim_LR)[min(round(nsim*.95)+1,nsim)]
  
  if(plot)
  {
    if(estimate_power)
    {
      xlim <- range(c(range(alt_sim_LR),c(range(null_sim_LR)),delta,critical_delta))
      h1 <- hist(alt_sim_LR,plot=FALSE)
      h2 <- hist(null_sim_LR,plot=FALSE)
      ylim <- range(c(0,range(h1$counts),range(h2$counts)))
      plot(h1,xlim=xlim,ylim=ylim,col=rgb(0,0,1,0.5),xlab = "Likelihood ratio",main = "Null vs Alt. Distribution")
      plot(h2,col=rgb(.1,.1,.1,.5),add=TRUE)
    } else hist(null_sim_LR,col=rgb(.1,.1,.1,.5),xlim=range(c(range(null_sim_LR),delta,critical_delta)),xlab = "Likelihood ratio",main = "Null Distribution")
    abline(v=delta,lwd=5,lty=3)
    abline(v=critical_delta,col="red",lwd=5,lty=3)
  }
  ret <- list(delta=delta,critical_delta=critical_delta,Pval=length(which(null_sim_LR>delta)) / nsim)
  if(estimate_power)
  {
    ret$power <- length(which(alt_sim_LR>=critical_delta)) / nsim
  }
  if(conf_int)
  {
    ret$model.par.sim <- model.par.sim
    ret$model.par.CI <- sort(model.par.sim)[c(round(nsim*(.025)),round(nsim*(1-.025)))]
  }
  plot.tools <- list(test_statistic=delta,critical_test_statistic=critical_delta,test_statistic_name="Likelihood Ratio",
                     null_sim_test_statistic=null_sim_LR)
  if(estimate_power) plot.tools <- c(plot.tools,list(alt_sim_test_statistic=alt_sim_LR))
  ret$plot.tools <- plot.tools
  
  ret$model1 <- model1SAVE
  ret$model2 <- model2SAVE
  class(ret) <- "compare.model"
  return(ret)
}

print.compare.model <- function(x,...)
{
  if(x$plot.tools$test_statistic_name!="K")
  {
    cat("\n**********Null model details**********")
    print(x$model1)
    cat("\n**********Alternative model details**********")
    print(x$model2)
  }  
  mat <- matrix(0,2+as.integer(!is.null(x$power)),1)
  mat[1] <- x$plot.tools$test_statistic
  mat[2] <- x$plot.tools$critical_test_statistic
  i <- 2
  if(!is.null(x$power))
  {
    i <- i+1
    mat[i] <- x$power
  }
  rownames(mat) <- c(paste("Test statistic (",x$plot.tools$test_statistic_name,")",sep=""),"Critical test statistic","Estimated Power")[1:(2+as.integer(!is.null(x$power)))]
  colnames(mat) <- ""
  cat("\n**********Simulation results**********")
  print(mat)
  cat("\n")
  if(!is.null(x$model.par.CI)) cat(names(x$model2$model.par)," 95% confidence interval: (",x$model.par.CI[1],",",x$model.par.CI[2],")\n",sep="")
  cat("P-value:",x$Pval)
  cat("\n")
}

plot.compare.model <- function(x,hist=FALSE,plot.alt=TRUE,
                               plot.test.statistic=TRUE,plot.critical.test.statistic=TRUE,xlim,ylim,
                               null.area.col=rgb(.1,.1,.1,.5),alt.area.col=rgb(0,0,1,0.5),
                               test.statistic.col="black",critical.test.statistic.col="red",
                               test.statistic.lwd=5,critical.test.statistic.lwd=5,
                               test.statistic.lty=3,critical.test.statistic.lty=3,...)
{
  estimate_power <- !is.null(x$plot.tools$alt_sim_test_statistic)
  if(estimate_power) estimate_power <- plot.alt
  test_statistic <- x$plot.tools$test_statistic
  critical_test_statistic <- x$plot.tools$critical_test_statistic
  null_sim_test_statistic <- x$plot.tools$null_sim_test_statistic
  test_statistic_name <- x$plot.tools$test_statistic_name
  if(estimate_power) alt_sim_test_statistic <- x$plot.tools$alt_sim_test_statistic
  if(estimate_power)
  {
    if(missing(xlim))
      xlim <- range(c(max(min(alt_sim_test_statistic),quantile(ecdf(alt_sim_test_statistic),.01)[[1]]),
                      max(min(null_sim_test_statistic),quantile(ecdf(null_sim_test_statistic),.01)[[1]]),
                      min(max(alt_sim_test_statistic),quantile(ecdf(alt_sim_test_statistic),.99)[[1]]),
                      min(max(null_sim_test_statistic),quantile(ecdf(null_sim_test_statistic),.99)[[1]]),
                      test_statistic,critical_test_statistic))
    if(hist)
    {
      h1 <- hist(alt_sim_test_statistic,plot=FALSE)
      h2 <- hist(null_sim_test_statistic,plot=FALSE)
      if(missing(ylim)) ylim <- range(c(0,range(h1$counts),range(h2$counts)))
      plot(h1,xlim=xlim,ylim=ylim,col=alt.area.col,xlab = test_statistic_name,main = "Null vs Alt. Distribution")
      plot(h2,col=null.area.col,add=TRUE)
    } else
    {
      d1 <- density(alt_sim_test_statistic)
      d2 <- density(null_sim_test_statistic)
      if(missing(ylim)) ylim <- c(range(0,range(d1$y),range(d2$y)))
      plot(d1,xlim=xlim,ylim=ylim,xlab=test_statistic_name,main="Null vs Alt. Distribution")
      polygon(d1,col=alt.area.col)
      lines(d2)
      polygon(d2,col=null.area.col)
    }
  } else
  {
    if(missing(xlim))
      xlim <- range(c(min(max(null_sim_test_statistic),quantile(ecdf(null_sim_test_statistic),.99)[[1]]),
                      max(min(null_sim_test_statistic),quantile(ecdf(null_sim_test_statistic),.01)[[1]]),
                      test_statistic,critical_test_statistic))
    if(hist)
    {
      if(missing(ylim))
        hist(null_sim_test_statistic,col=null.area.col,xlim=xlim,xlab = test_statistic_name,main = "Null Distribution") else
          hist(null_sim_test_statistic,col=null.area.col,xlim=xlim,xlab = test_statistic_name,main = "Null Distribution",ylim=ylim)
    } else
    {
      d2 <- density(null_sim_test_statistic)
      if(missing(ylim))
        plot(d2,xlab=test_statistic_name,main="Null Distribution",xlim=xlim) else
          plot(d2,xlab=test_statistic_name,main="Null Distribution",xlim=xlim,ylim=ylim)
      polygon(d2,col=null.area.col)
    }
  }
  if(plot.test.statistic) abline(v=test_statistic,lwd=test.statistic.lwd,lty=test.statistic.lty,col=test.statistic.col)
  if(plot.critical.test.statistic) abline(v=critical_test_statistic,col=critical.test.statistic.col,
                                          lwd=critical.test.statistic.lwd,lty=critical.test.statistic.lty)
}

sim.groups <- function(tree,groups,painted.edges,model="BM",parameters=list(),phylocov,nsim=100,fixed.effects,return.type="matrix")
{
  tree <- reorder(tree,"postorder")
  
  ngroups <- nlevels(groups)
  group_levels <- levels(groups)
  if(!is.list(phylocov)) stop("phylocov must be a list equal in length to the number of groups.")
  if(length(phylocov)==1) phylocov <- rep(phylocov,ngroups)
  if(length(model)==1) model <- rep(model,ngroups)
  if(length(parameters)==1) parameters <- rep(parameters,ngroups)
  ntraits <- ncol(phylocov[[1]])
  Y <- rep(list(matrix(0,length(tree$tip.label),ntraits,dimnames=list(tree$tip.label,paste("V",1:ntraits,sep="")))),nsim)
  
  nspecies <- length(tree$tip.label)
  nedge <- nrow(tree$edge)
  m <- painted.edges
  
  tree_list <- rep(tree,ngroups)
  if(ngroups==1) tree_list <- list(tree)
  sims <- vector("list",ngroups)
  for(i in 1:ngroups)
  {
    tree_list[[i]]$edge.length <- tree$edge.length*m[1:nedge,i]
    tree_list[[i]] <- reorder(tree_list[[i]],"pruningwise")
    temp_tr <- transf.branch.lengths(phy = tree_list[[i]],model = model[i],parameters = parameters[i])$tree
    temp_tr$edge.length[tree_list[[i]]$edge.length==0] <- 0
    tree_list[[i]] <- temp_tr
    sims[[i]] <- sim.char(phy = tree_list[[i]],par = phylocov[[i]],nsim = nsim,root = 0)
    for(j in 1:nsim)
    {
      Y[[j]] <- Y[[j]] + sims[[i]][,,j]
    }
  }
  if(return.type!="matrix") for(j in 1:nsim) Y[[j]] <- data.frame(species=rownames(Y[[j]]),Y[[j]])
  tree <- reorder(tree,"postorder")
  if(nsim==1) Y <- Y[[j]]
  list(trait_data=Y,tree=tree)
  
}

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

phylocurve <- function(formula,tree,data,ymin=.01,ymax=.99,ylength=30,tip.coefficients,species.identifier="species",verbose=FALSE)
{
  if(missing(tip.coefficients))
  {
    if(!(species.identifier %in% colnames(data))) stop("Add species names column to data.")
    tip.coefficients <- get.tip.coefficients(formula = formula,tree = tree,data = data,ymin = ymin,ymax = ymax,ylength = ylength,species.identifier = species.identifier,verbose = verbose)
  } else if(verbose) cat("Phase 1: Tip estimates already provided\n")
  pgls_curve(tree = tree,tip.coefficients = tip.coefficients,ymin = ymin,ymax = ymax,ylength = ylength,verbose = verbose)
}

phylocurve.generalized <- function(tree,X,Y)
{
  nspecies <- length(tree$tip.label)
  p <- fast_anc_hand(x = X,Y = Y,tree = tree)
  tree1 <- multi2di(tree)
  ancx <- apply(t(matrix(p$aligned_coordinates$x,ncol = nspecies,dimnames = list(NULL,tree1$tip.label))),2,function(X) ace(X,tree1,method="pic")$ace[1])
  ancy <- apply(t(matrix(p$aligned_coordinates$y,ncol = nspecies,dimnames = list(NULL,tree1$tip.label))),2,function(X) ace(X,tree1,method="pic")$ace[1])
  p$nr <- nrow(p$aligned_coordinates)/nspecies
  
  ret <- c(p,list(anc_X=ancx,anc_Y=ancy,tree=tree))
  ret
}

polynomial.fit <- function(data, x_variable, y_variable, method = "BIC",
                           nterms = 2,min_x = -Inf,max_x = Inf,min_y = 0, max_y = Inf,eval_length=30)
{
  if(is.null(colnames(data))) stop("data variable must have column names")
  speciesInd <- match("species", colnames(data))
  xInd <- match(x_variable, colnames(data))
  yInd <- match(y_variable, colnames(data))
  if(is.na(speciesInd)) stop("data variable must have column named 'species'")
  if(is.na(xInd)) stop("data variable must have a column name matching x_variable")
  if(is.na(yInd)) stop("data variable must have a column name matching y_variable")
  dat <- data[,c(speciesInd,xInd,yInd)]
  colnames(dat) <- c("species","x","y")
  dat <- dat[dat$x>=min_x & dat$x<=max_x,]
  new_x <- unique(dat[,2])
  new_x <- seq(min(new_x),max(new_x),length=max(length(new_x),eval_length))
  species <- unique(as.character(dat[,1]))
  nspecies <- length(species)
  
  mod_list <- vector("list",nspecies)
  names(mod_list) <- unique(as.character(dat$species))
  
  ret <- matrix(NA,nrow = nspecies,ncol = length(new_x))
  colnames(ret) <- new_x
  rownames(ret) <- species
  mod <- "y ~ x"
  for(i in 1:nterms)
  {
    mod <- paste(mod," + I(x^",i,") ",sep="",collapse="")
  }
  mod <- formula(mod)
  
  for(i in 1:nspecies)
  {
    temp_species <- species[i]
    x <- dat[dat[,1]==temp_species,2]
    y <- dat[dat[,1]==temp_species,3]
    complete_cases_y <- which(complete.cases(y))
    x <- x[complete_cases_y]
    y <- y[complete_cases_y]
    n <- length(complete_cases_y)
    temp_lm <- lm(mod)
    k <- if(method=="AIC") 2 else log(n)
    mod_list[[i]] <- step(temp_lm,direction = "both",k = k,trace = 0)
    if(length(coef(mod_list[[i]]))==1) mod_list[[i]] <- lm(y~x)
    pred_y <- predict.lm(mod_list[[i]],newdata = data.frame(x=new_x))
    pred_y[pred_y<min_y] <- min_y
    pred_y[pred_y>max_y] <- max_y
    ret[i,] <- pred_y
  }
  ret <- ret[complete.cases(ret),]
  list(X=new_x,Y=ret)
}

nonlinear.fit <- function(data, x_variable, y_variable, fct = LL2.3(),
                          min_x = -Inf,max_x = Inf,min_y = 0, max_y = Inf,eval_length = 30,...)
{
  if(is.null(colnames(data))) stop("data variable must have column names")
  speciesInd <- match("species", colnames(data))
  xInd <- match(x_variable, colnames(data))
  yInd <- match(y_variable, colnames(data))
  if(is.na(speciesInd)) stop("data variable must have column named 'species'")
  if(is.na(xInd)) stop("data variable must have a column name matching x_variable")
  if(is.na(yInd)) stop("data variable must have a column name matching y_variable")
  dat <- data[,c(speciesInd,xInd,yInd)]
  colnames(dat) <- c("species","x","y")
  dat <- dat[dat$x>=min_x & dat$x<=max_x,]
  new_x <- unique(dat[,2])
  new_x <- seq(min(new_x),max(new_x),length=max(length(new_x),eval_length))
  species <- unique(as.character(dat[,1]))
  nspecies <- length(species)
  
  mod_list <- vector("list",nspecies)
  names(mod_list) <- unique(as.character(dat$species))
  
  ret <- matrix(NA,nrow = nspecies,ncol = length(new_x))
  colnames(ret) <- new_x
  rownames(ret) <- species
  
  for(i in 1:nspecies)
  {
    temp_species <- species[i]
    x <- dat[dat[,1]==temp_species,2]
    y <- dat[dat[,1]==temp_species,3]
    complete_cases_y <- which(complete.cases(y))
    x <- x[complete_cases_y]
    y <- y[complete_cases_y]
    temp_lm <- try(drm(y~x,fct = fct,...),silent=TRUE)
    if(class(temp_lm)=="try-error")
    {
      warning("Convergence failed for ",temp_species,". Using simple linear regression.",immediate. = TRUE)
      pred_y <- predict(lm(y~x),newdata = data.frame(x=new_x))
    } else
    {
      mod_list[[i]] <- temp_lm
      pred_y <- predict(mod_list[[i]],newdata = data.frame(x=new_x))
    }
    pred_y[pred_y<min_y] <- min_y
    pred_y[pred_y>max_y] <- max_y
    ret[i,] <- pred_y
  }
  
  ret <- ret[complete.cases(ret),]
  list(X=new_x,Y=ret)
}

GP.fit <- function(data, x_variable, y_variable,
                   min_x = -Inf,max_x = Inf,min_y = 0, max_y = Inf,eval_length = 30,...)
{
  if(is.null(colnames(data))) stop("data variable must have column names")
  speciesInd <- match("species", colnames(data))
  xInd <- match(x_variable, colnames(data))
  yInd <- match(y_variable, colnames(data))
  if(is.na(speciesInd)) stop("data variable must have column named 'species'")
  if(is.na(xInd)) stop("data variable must have a column name matching x_variable")
  if(is.na(yInd)) stop("data variable must have a column name matching y_variable")
  dat <- data[,c(speciesInd,xInd,yInd)]
  colnames(dat) <- c("species","x","y")
  dat <- dat[dat$x>=min_x & dat$x<=max_x,]
  new_x <- unique(dat[,2])
  new_x <- seq(min(new_x),max(new_x),length=max(length(new_x),eval_length))
  species <- unique(as.character(dat[,1]))
  nspecies <- length(species)
  
  mod_list <- vector("list",nspecies)
  names(mod_list) <- unique(as.character(dat$species))
  
  ret <- matrix(NA,nrow = nspecies,ncol = length(new_x))
  colnames(ret) <- new_x
  rownames(ret) <- species
  
  for(i in 1:nspecies)
  {
    temp_species <- species[i]
    x <- dat[dat[,1]==temp_species,2]
    y <- dat[dat[,1]==temp_species,3]
    complete_cases_y <- which(complete.cases(y))
    x <- x[complete_cases_y]
    y <- y[complete_cases_y]
    temp_lm <- try(GP_fit(X = normalize_to_01(c(range(new_x),x))[-(1:2)],Y = y,...),silent=TRUE)
    if(class(temp_lm)=="try-error")
    {
      warning("Convergence failed for ",temp_species,". Using simple linear regression.",immediate. = TRUE)
      pred_y <- predict(lm(y~x),newdata = data.frame(x=new_x))
    } else
    {
      mod_list[[i]] <- temp_lm
      pred_y <- predict.GP(mod_list[[i]],xnew = normalize_to_01(new_x))[[1]]
    }
    pred_y[pred_y<min_y] <- min_y
    pred_y[pred_y>max_y] <- max_y
    ret[i,] <- pred_y
  }
  ret <- ret[complete.cases(ret),]
  list(X=new_x,Y=ret)
}

phylocurve.trim <- function(phylocurve.generalized,min_Y = -Inf,max_Y = Inf,min_X = -Inf,max_X = Inf)
{
  p <- phylocurve.generalized
  tree <- p$tree
  ancx <- p$anc_X
  ancy <- p$anc_Y
  range_x <- which(ancx>min_X & ancx<max_X)
  range_y <- which(ancy>min_Y & ancy<max_Y)
  range <- intersect(range_x,range_y)
  nspecies <- length(tree$tip.label)
  nr <- nrow(p$aligned_coordinates)/nspecies
  p$aligned_coordinates <- p$aligned_coordinates[as.double(sapply((1:nspecies-1)*max(nr),function(X) X + range)),]
  p$aligned_data <- p$aligned_data[,c(1,range+1,range+nr+1)]
  p$aligned_X <- p$aligned_X[,range]
  p$aligned_Y <- p$aligned_Y[,range]
  p$anc_X <- p$anc_X[range]
  p$anc_Y <- p$anc_Y[range]
  p$nr <- nr
  p
}

get.tip.coefficients <- function(formula,tree,data,ymin=.01,ymax=.99,ylength=30,species.identifier="species",verbose=FALSE)
{
  if(!(species.identifier %in% colnames(data))) stop("Add species names column to data.")
  dat <- model.frame(formula,data=data)
  dat$species <- data$species
  data <- dat
  nspecies <- length(tree$tip.label)
  y_name <- colnames(data)[1]
  x_name <- colnames(data)[2]
  taxa <- tree$tip.label
  tip.coefficients <- matrix(NA,nrow = nspecies,ncol = 2,dimnames = list(taxa,c("(Intercept)",x_name)))
  if(verbose) cat("Phase 1: Fitting species curves: ")
  counter <- 0
  for(i in 1:nspecies)
  {
    if(verbose)
    {
      if((ceiling(i/nspecies*100)>counter))
      {
        cat(paste(counter,"%   ",sep=""))
        counter <- counter + 10
      }
      if(i==nspecies) cat(paste(100,"%   ",sep=""))
      
    }
    tip.coefficients[taxa[i],] <- coef(glm(formula,family = quasibinomial("logit"),data = data[data$species==taxa[i],]))
  }
  if(verbose) cat("\n")
  tip.coefficients
}

get.aligned.function.data <- function(tip.coefficients,ylength=30,ymin=.01,ymax=.99)
{
  Y <- t(apply(tip.coefficients,1,logit_inv,y=seq(ymin,ymax,length=ylength)))
  data.frame(species=rownames(Y),Y)
}

pgls_curve <- function(tree,tip.coefficients,varAY,vals_only=FALSE,ymin,ymax,ylength,verbose)
{
  if(is.null(rownames(tip.coefficients)))
  {
    warning("No row names for tip.coefficients. Assuming rows match up with tree tips.")
    rownames(tip.coefficients) <- tree$tip.label
  }
  if(missing(varAY))
  {
    D <- dist.nodes(tree)
    nspecies <- length(tree$tip.label)
    Nnode <- tree$Nnode
    MRCA <- mrca(tree, full = TRUE)
    M <- D[as.character(nspecies + 1), MRCA]
    dim(M) <- rep(sqrt(length(M)), 2)
    varAY <- M[(nspecies+1):(nspecies+Nnode), 1:nspecies]
    varA <- M[(nspecies+1):(nspecies+Nnode), (nspecies+1):(nspecies+Nnode)]
    colnames(varAY) <- tree$tip.label
  }
  
  nspecies <- length(tree$tip.label)
  Nnode <- tree$Nnode
  X <- matrix(1,nspecies,1)
  y <- seq(ymin,ymax,length=ylength)
  rownames(X) <- tree$tip.label
  tip.coefficients <- tip.coefficients[tree$tip.label,]
  Yinv <- t(apply(tip.coefficients,1,function(X) logit_inv(X,y)))
  if(verbose) cat("Phase 2: Reconstructing root curve","\n")
  root_three_point <- three.point.compute(phy = tree,P = X,Q = Yinv)
  root <- solve(root_three_point$PP) %*% root_three_point$QP
  root <- (matrix(rep(root,nspecies),length(root),nspecies))
  if(verbose) cat("Phase 3: Reconstructing internal node curves: ")
  P <- t(varAY)
  Q <- (Yinv-t(root))
  nvar <- ncol(P)
  var_BM <- apply(Q,2,function(X) t(three.point.compute(tree,X)$PP)/(nspecies-1))
  anc_three_point <- matrix(NA,ncol(Q),ncol(P))
  DD <- D1 <- matrix(NA,nvar,1)
  for(i in 1:ceiling(nvar/50))
  {
    if(verbose) cat(paste(round(i/ceiling(nvar/50)*100),"%   ",sep=""))
    temp <- three.point.compute(phy = tree,P = P[,(1+((i-1)*50)):min(50+(i-1)*50,nvar)],Q = Q)
    DD[(1+((i-1)*50)):min(50+(i-1)*50,nvar)] <- diag(temp$PP)
    D1[(1+((i-1)*50)):min(50+(i-1)*50,nvar)] <- temp$P1
    anc_three_point[,(1+((i-1)*50)):min(50+(i-1)*50,nvar)] <- temp$QP
  }
  vec11 <- temp$vec11
  var_nodes <- t((diag(varA) - DD + diag((1 - D1) %*% solve(vec11) %*% t(1 - D1))))
  var_nodes <- matrix(rep(var_nodes,ylength),ncol=ylength)
  var_BM <- matrix(rep(var_BM,tree$Nnode),ncol=tree$Nnode)
  CI <- t(var_nodes)*(var_BM)*1.96
  anc_vals <- (anc_three_point)+root[,1:Nnode]
  
  if(vals_only) return(t(anc_vals))
  if(verbose) cat("\nPhase 4: Fitting ancesral curves")
  ret <- t(apply(anc_vals,2,function(X) coef(glm(y~X,family=quasibinomial("logit")))))
  rownames(ret) <- (nspecies+1):(nspecies+Nnode)
  rownames(CI) <- rownames(anc_vals) <- seq(ymin,ymax,length=ylength)
  colnames(CI) <- colnames(anc_vals) <- (nspecies+1):(nspecies+tree$Nnode)
  lower_CI <- anc_vals - CI
  upper_CI <- anc_vals + CI
  return(list(node_coefficients=ret,fitted_x=anc_vals,lower_CI_x=lower_CI,upper_CI_x=upper_CI,y_vals=seq(ymin,ymax,length=ylength),tip.coefficients=tip.coefficients,tip_X=Yinv))
}

# logit function
# first parameter is the glm intercept
# second parameter the glm slope
logit_fx <- function(theta,x)
{
  b <- theta[1]
  m <- theta[2]
  exp(m*x+b)/(1+exp(m*x+b))
}

# inverse logit function
logit_inv <- function(theta,y)
{
  b <- theta[1]
  m <- theta[2]
  ((log(-y/(y-1))-b)/m)
}

# transform logit glm parameters to uncorrelated slope parameter
logit_slope_func <- function(theta)
{
  b <- theta[1]
  m <- theta[2]
  m/4
}

# transform logit glm parameters to uncorrelated EC50 parameter
logit_ec50_func <- function(theta)
{
  b <- theta[1]
  m <- theta[2]
  -b/m
}

# back-calculates glm intercept parameter
logit_b_func <- function(slope,ec50)
{
  m <- slope*4
  b <- -ec50 * m
  return(b)
}

# back-calculates glm slope parameter
logit_m_func <- function(slope,ec50)
{
  m <- slope*4
  return(m)
}

sim.curves <- function(nspecies = 30,x_length=20,startree=FALSE,lambda=1,seed)
{
  if(!missing(seed)) set.seed(seed)
  x <- seq(0,10,length=x_length) # environmental gradient
  tree <- pbtree(n=nspecies)
  if(!missing(lambda)) simtree <- rescale(tree,"lambda",lambda=lambda) else simtree <- tree
  if(startree) simtree <- starTree(tree$tip.label,rep(1,length(tree$tip.label)))
  
  logit_ec50_true <- fastBM(simtree,5,bounds=c(2,20),internal = TRUE)
  logit_slope_true <- fastBM(simtree,.5,bounds=c(.2,1),internal = TRUE)
  
  logit_ec50 <- logit_ec50_true[1:nspecies]
  logit_slope <- logit_slope_true[1:nspecies]
  
  logit_b <- logit_b_func(slope=logit_slope,ec50=logit_ec50)
  logit_m <- logit_m_func(slope=logit_slope,ec50=logit_ec50)
  logit_coefs <- cbind(logit_b,logit_m)
  logit_sim_coefs <- logit_coefs
  
  logit_y <- data.frame(species =   as.character(t(matrix(rep(simtree$tip.label,length(x)),nrow=nspecies))),
                        x =  rep(x,nspecies),y=rep(0,nspecies*length(x)))
  for(i in 1:nspecies)
  {
    logit_y[1:length(x)+length(x)*(i-1),3] <- logit_fx(c(logit_b[i],logit_m[i]),x)
  }
  ret <- list(data = logit_y,tree = simtree,
              true_coefs = data.frame(Intercept=logit_b_func(logit_slope_true,logit_ec50_true),slope=logit_m_func(logit_slope_true,logit_ec50_true),row.names = names(logit_ec50_true))
  )
  ret
}

# calculates Felsenstein's ancestral state reconstruction values
ace_hand <- function(x,Y,tree,gpr_fit=TRUE,gp_list)
{
  gp_missing <- missing(gp_list)
  nspecies <- length(tree$tip.label)
  if(gp_missing)
  {
    gp_list <- vector("list",nspecies)
    names(gp_list) <- tree$tip.label
  }
  if(gpr_fit) x <- normalize_to_01((x-min(x))/(max(x)-min(x)))
  
  Y <- Y[tree$tip.label,]
  Y <- rbind(Y,matrix(0,nrow=tree$Nnode,ncol=ncol(Y)))
  v1_lookup <- v2_lookup <- n1 <- n2 <- d1 <- d2 <- double(1)
  corrected_d <- tree$edge.length
  for(i in (tree$Nnode+nspecies):(nspecies+1))
  {
    v1_lookup <- which(tree$edge[,1]==i)[1]
    v2_lookup <- which(tree$edge[,1]==i)[2]
    n1 <- tree$edge[v1_lookup,2] # node 1 number
    n2 <- tree$edge[v2_lookup,2] # node 2 number
    d1 <- (1-(corrected_d[v1_lookup]/(corrected_d[v1_lookup]+corrected_d[v2_lookup])))
    d2 <- (1-(corrected_d[v2_lookup]/(corrected_d[v1_lookup]+corrected_d[v2_lookup])))
    if(gpr_fit)
    {
      if(!(gp_missing)) gp1 <- gp_list[[n1]] else gp1 <- GP_fit(x,Y[n1,])
      if(!(gp_missing)) gp2 <- gp_list[[n2]] else gp2 <- GP_fit(x,Y[n2,])
      if(n1<=nspecies) gp_list[[n1]] <- gp1
      if(n2<=nspecies) gp_list[[n2]] <- gp2
      time_warp <- dtw(predict(gp1,x)$Y_hat,predict(gp2,x)$Y_hat)
    } else
    {
      gp1 <- coef(glm(Y[n1,]~x,family=quasibinomial("logit")))
      gp2 <- coef(glm(Y[n2,]~x,family=quasibinomial("logit")))
      time_warp <- dtw(logit_fx(gp1,x),logit_fx(gp2,x))
    }
    
    new_x1 <- (((max(x)-min(x))*(time_warp$index1-min(time_warp$index1))) / (max(time_warp$index1) - min(time_warp$index1))) + min(x)
    new_x2 <- (((max(x)-min(x))*(time_warp$index2-min(time_warp$index2))) / (max(time_warp$index2) - min(time_warp$index2))) + min(x)
    if(gpr_fit)
    {
      new_y1 <- predict(gp1,new_x1)$Y_hat
      new_y2 <- predict(gp2,new_x2)$Y_hat
    } else
    {
      new_y1 <- logit_fx(gp1,new_x1)
      new_y2 <- logit_fx(gp2,new_x2)
    }
    anc_x <- (d1*new_x1) + (d2*new_x2)
    anc_y <- (d1*new_y1) + (d2*new_y2)
    if(gpr_fit)
    {
      anc_gp <- GP_fit(anc_x,anc_y)
      Y[i,] <- predict(anc_gp,x)$Y_hat
    } else
    {
      anc_gp <- coef(glm(anc_y~anc_x,family=quasibinomial("logit")))
      Y[i,] <- logit_fx(anc_gp,x)
    }
    corrected_d[which(tree$edge[,2]==i)] <- corrected_d[which(tree$edge[,2]==i)]+((corrected_d[v1_lookup]*corrected_d[v2_lookup])/(corrected_d[v1_lookup]+corrected_d[v2_lookup]))
  }
  anc_y <- Y[(length(tree$tip.label)+1):(tree$Nnode+length(tree$tip.label)),]
  return(list(anc_y=anc_y,gp_list=gp_list,anc_gp=anc_gp))
}

# calculates ancestral state reconstruction with pairwise dynamic time warping and 
# the fast PIC method (adapted from fastAnc in phytools)
fast_anc_hand <- function(x,Y,tree,root_only=TRUE,gpr_fit=TRUE)
{
  if(is.null(rownames(Y)))
  {
    rownames(Y) <- tree$tip.label
    warning("No taxa labels attached to coefficients. Assuming order matches tips on tree.")
  } else
  {
    temp <- match(tree$tip.label,rownames(Y))
    Y <- Y[temp,]
  }
  
  tip_gp <- vector("list",length = nrow(Y))
  
  M <- tree$Nnode
  N <- length(tree$tip.label)
  a <- tree
  node_vals <- matrix(0,M,ncol(Y))
  # calculate root state
  for(i in 1:M + N)
  {
    a <- multi2di(root(tree,node=(i)))
    if(i==(1+N))
    {
      res <- ace_hand(x,Y,tree=a,gpr_fit = gpr_fit)
      gp_list <- res$gp_list
      anc_Y <- res[[1]][1,]
      anc_gp <- res$anc_gp
    } else
    {
      res <- ace_hand(x,Y,tree=a,gpr_fit = gpr_fit,gp_list = gp_list)
    }
    node_vals[i-N,] <- res[[1]][1,]
    if(root_only)
    {
      #return(node_vals[1,])
      break
    }
  }
  rownames(node_vals) <- 1:M + N
  
  alignments <- vector("list",nrow(Y))
  names(alignments) <- names(res$gp_list)
  max_dtw <- rep(1,length(x))
  for(i in 1:nrow(Y))
  {
    alignments[[i]] <- dtw(Y[i,],anc_Y)
    temp_max <- tapply(alignments[[i]]$index2,alignments[[i]]$index2,FUN = function(X) length(X))
    max_dtw[temp_max > max_dtw] <- temp_max[temp_max > max_dtw]
  }
  anc_align <- numeric()
  for(i in 1:length(x))
  {
    anc_align <- c(anc_align,rep(i,max_dtw[i]))
  }
  
  x_align <- matrix(0,nrow = nrow(Y),ncol = length(anc_align))
  
  for(j in 1:nrow(Y))
  {
    low <- 1
    for(i in 1:length(x))
    {
      temp_align <- alignments[[j]]$index1
      temp_ref <- alignments[[j]]$index2
      temp_align <- temp_align[temp_ref==i]
      temp_range <- (round(seq(min(temp_align),max(temp_align),length=max_dtw[i])))
      high <- low + max_dtw[i] - 1
      x_align[j,low:high] <-  temp_range
      low <- high + 1
    }
  }
  convert_x <- (x-min(x))/(max(x)-min(x))
  new_x <- apply(x_align,1,function(X) x[X])
  new_convert_x <- (new_x-min(x))/(max(x)-min(x))
  new_Y <- matrix(0,nrow = nrow(Y),ncol = nrow(new_x))
  for(i in 1:nrow(Y))
  {
    new_Y[i,] <- predict.GP(gp_list[[i]],xnew = new_convert_x[,i])$Y_hat
  }
  new_x <- t(new_x)
  rownames(new_Y) <- rownames(new_x) <- tree$tip.label
  aligned_data <- data.frame(species=tree$tip.label,new_x,new_Y,row.names=tree$tip.label)
  colnames(aligned_data) <- c("species",paste("x",1:ncol(new_x),sep=""),paste("y",1:ncol(new_Y),sep=""))
  aligned_coordinates <- data.frame(species = as.character(t(matrix(rep(tree$tip.label,ncol(new_x)),ncol=ncol(new_x)))),
                                    x = as.double(t(new_x)),
                                    y = as.double(t(new_Y)))
  return(list(aligned_data=aligned_data,aligned_coordinates=aligned_coordinates,aligned_X=new_x,aligned_Y=new_Y))
}

normalize_to_01 <- function(x)
{
  max_x <- max(x)
  min_x <- min(x)
  ret <- ((x-min_x) / (max_x - min_x))
  ret[ret<0] <- 0
  ret[ret>1] <- 1
  ret
}

sim.traits <- function(ntaxa=15,ntraits=4,nreps=1,nmissing=0,tree,v,anc,intraspecific,model="BM",parameters,nsim=1,return.type="matrix")
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

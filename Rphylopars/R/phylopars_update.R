# pheno_correlated <- TRUE
# pheno_error <- TRUE
# REML <- TRUE
# EM_Fels_limit <- 1e3
# model <- "lambda"
# repeat_optim_limit <- 1
# EM_missing_limit <- 50
# repeat_optim_tol <- 1e-2
# model_par_evals <- 10
phylopars <- function(trait_data,tree,model="BM",pheno_error,phylo_correlated=TRUE,pheno_correlated=TRUE,REML=TRUE,full_alpha=TRUE,phylocov_start,phenocov_start,model_par_start,phylocov_fixed,phenocov_fixed,model_par_fixed,skip_optim=FALSE,skip_EM=FALSE,EM_Fels_limit=1e3,repeat_optim_limit=1,EM_missing_limit=50,repeat_optim_tol = 1e-2,model_par_evals=10,max_delta=1e4,EM_verbose=FALSE,optim_verbose=FALSE,npd=FALSE,nested_optim=FALSE,usezscores=TRUE)
{
  tree <- reorder(tree,"postorder")
  trait_data[,1] <- as.character(trait_data[,1])
  f_args <- as.list(environment())
  if(model=="white" | model=="star") tree <- starTree(tree$tip.label,rep(1,length(tree$tip.label)))
  nvar <- ncol(trait_data)-1
  if(nvar==1 & model=="mvOU") model <- "OU"
  
  height <- mean(pruningwise.distFromRoot((tree)))
  if(model=="lambda")
    par_bounds <- seq(0,1,length=model_par_evals) else
      if(model=="kappa")
        par_bounds <- seq(0,3,length=model_par_evals) else
          if(model=="delta")
            par_bounds <- seq(10^(-7)/height,3/height,length=model_par_evals) else
              if(model=="OU")
                par_bounds <- seq(10^(-7)/height,50/height,length=model_par_evals) else
                  if(model=="EB")
                    par_bounds <- seq(-3/height,0,length=model_par_evals)
  
  evec <- function(bounds)
  {
    if((model=="OU" & bounds==0) | (model=="EB" & bounds==0) | (model=="lambda" & bounds==1) | (model=="delta" & bounds==1) | (model=="kappa" & bounds==1)) return(edge_vec)
    if(model=="OU") edge_vec <- c(reorder(rescale(tree,model="OU",alpha=bounds),"postorder")$edge.length,0)
    if(model=="EB") edge_vec <- c(reorder(rescale(tree,model="EB",a=bounds),"postorder")$edge.length,0)
    if(model=="lambda") edge_vec <- c(reorder(rescale(tree,model="lambda",lambda=bounds),"postorder")$edge.length,0)
    if(model=="kappa") edge_vec <- c(reorder(rescale(tree,model="kappa",kappa=bounds),"postorder")$edge.length,0)
    if(model=="delta") edge_vec <- c(reorder(rescale(tree,model="delta",delta=bounds),"postorder")$edge.length,0)
    edge_vec
  }
  
  if(!is.ultrametric(tree) & (model=="OU" | model=="mvOU")) stop("OU model not currently supported for non-ultrametric trees.")
  if(model=="mvOU" & FALSE)
  {
    pheno_error <- FALSE
    pheno_correlated <- FALSE
  }
  
  if(missing(phylocov_start)) phylocov_start <- matrix(NA) else if(is.na(phylocov_start)[[1]]) phylocov_start <- matrix(NA)
  if(missing(phenocov_start)) phenocov_start <- matrix(NA) else if(is.na(phenocov_start)[[1]]) phenocov_start <- matrix(NA)
  if(missing(model_par_start)) model_par_start <- matrix(NA) else if(is.na(model_par_start)[[1]]) model_par_start <- matrix(NA)
  
  
  zstand <- list(phylocov=diag(nvar),mu=matrix(rep(0,nvar)))
  if(usezscores)
  {
    for(i in 1:nvar)
    {
      dat <- setNames(trait_data[,i+1],as.character(trait_data$species))
      dat <- tapply(dat,names(dat),function(X) mean(X,na.rm=TRUE))
      #dat <- dat[!is.na(dat)]
      
      zstand$mu[i] <- mean(dat,na.rm=TRUE)
      zstand$phylocov[i,i] <- var(dat,na.rm=TRUE)*(length(tree$tip.label)-1)/(length(tree$tip.label)-REML)
      
      #drp <- name.check(phy = tree,data.names = names(dat))
      #if(length(drp)>1)
      #{
      #  if(length(drp[[1]])>0) temp_tree <- drop.tip(tree,drp$tree_not_data) else temp_tree <- tree
      #  if(length(drp[[2]])>0) dat <- dat[-match(drp$data_not_tree,names(dat))]
      #} else
      #{
      #  temp_tree <- tree
      #}
      #zstand$phylocov[i,i] <- sum(pic(x = dat,phy = multi2di(temp_tree,random = FALSE))^2) / (length(temp_tree$tip.label)-REML)
      #zstand$mu[i] <- ace(x = dat,phy = multi2di(temp_tree,random=FALSE),method="pic")$ace[[1]]
    }
  }
  
  if(missing(phylocov_fixed)) phylocov_fixed <- zphylocov_fixed <- matrix(NA) else if(!is.na(phylocov_fixed)[[1]]) zphylocov_fixed <- as.matrix(phylocov_fixed / (sqrt(diag(zstand$phylocov)) %*% t(sqrt(diag(zstand$phylocov))))) else phylocov_fixed <- zphylocov_fixed <- matrix(NA)
  if(missing(phenocov_fixed)) phenocov_fixed <- zphenocov_fixed <- matrix(NA) else if(!is.na(phenocov_fixed)[[1]]) zphenocov_fixed <- as.matrix(phenocov_fixed / (sqrt(diag(zstand$phylocov)) %*% t(sqrt(diag(zstand$phylocov))))) else phenocov_fixed <- zphenocov_fixed <- matrix(NA)
  if(missing(model_par_fixed)) model_par_fixed <- matrix(NA) else if(model=="mvOU") alpha <- model_par_fixed
  
  if(is.na(model_par_fixed)[[1]]) model_par_fixed <- matrix(NA)
  
  if(!is.na(phylocov_start)[[1]]) phylocov <- as.matrix(phylocov_start)
  if(!is.na(phylocov_fixed)[[1]]) phylocov <- as.matrix(phylocov_fixed) else phylocov <- diag(nvar)
  if(!is.na(phenocov_start)[[1]]) phenocov <- as.matrix(phenocov_start)
  if(!is.na(phenocov_fixed)[[1]]) phenocov <- as.matrix(phenocov_fixed) else phenocov <- diag(nvar)
  if(!is.na(model_par_start)[[1]]) model_par <- as.matrix(model_par_start)
  if(!is.na(model_par_fixed)[[1]]) model_par <- as.matrix(model_par_fixed) else if(model=="mvOU") model_par <- diag(nvar)
  
  if(!is.na(model_par_fixed[[1]])) nested_optim <- TRUE
  
  mu <- matrix(0)
  nspecies <- length(tree$tip.label)
  nedge <- length(tree$edge.length)
  if(!missing(pheno_error)) if(pheno_error==FALSE) pheno_correlated <- FALSE
  if(pheno_correlated) pheno_error <- TRUE
  if(max(table(as.character(trait_data$species)))>1 & if(!missing(pheno_error)) !pheno_error else FALSE)
  {
    temp_dat <- convert_to_means(trait_data,sort_vec=tree$tip.label)
    temp_dat <- data.frame(species=rownames(temp_dat),temp_dat)
    trait_data <- temp_dat
    rm(temp_dat)
  } else if(max(table(as.character(trait_data$species)))==1) pheno_error <- pheno_correlated <- FALSE else pheno_error <- TRUE
  
  if(pheno_correlated) pheno_error <- 2 else if(pheno_error) pheno_error <- 1 else pheno_error <- 0
  
  if(pheno_error==0)
  {
    phenocov_fixed <- phenocov_start <- matrix(NA)
  }
  
  # index for model parameters
  if(!is.na(phylocov_fixed)[[1]]) phylocov_pars <- numeric() else if(phylo_correlated) phylocov_pars <- 1:((nvar^2-nvar)/2+nvar) else phylocov_pars <- 1:nvar
  if(!is.na(phenocov_fixed)[[1]]) phenocov_pars <- numeric() else if(pheno_error>0) phenocov_pars <- if(pheno_correlated) length(phylocov_pars) + 1:((nvar^2-nvar)/2+nvar) else (length(phylocov_pars)+1):(length(phylocov_pars)+nvar)
  if(model=="mvOU")
  {
    if(!is.na(model_par_fixed)[[1]]) alpha_pars <- numeric() else alpha_pars <- if(full_alpha) length(phylocov_pars) + 1:((nvar^2-nvar)/2+nvar) else (length(phylocov_pars)+1):(length(phylocov_pars)+nvar)
    if(is.na(model_par_fixed)[[1]]) if(pheno_error>0) alpha_pars <- alpha_pars + length(phenocov_pars)
  } else alpha_pars <- numeric()
  
  if(!is.na(model_par_fixed[[1]])) perm_model <- model else perm_model <- NA
  
  # which variables are represented (by row) in trait_data
  species_subset <- lapply(apply(trait_data[,1:nvar+1,drop=FALSE],1,function(X) list(which(!is.na(X))-1)),function(X) X[[1]])
  
  # which variables are NOT represented (by row) in trait_data
  un_species_subset <- as.list(apply(trait_data[,1:nvar+1,drop=FALSE],1,function(X) which(is.na(X))-1))
  
  # what unique combinations of variables subsets exist
  subset_list <- unique(species_subset)
  
  # what combination of unique variable combinations do each row correspond to
  tip_combn <- match(species_subset,subset_list)-1
  
  # edge indices
  anc <- tree$edge[,1]-1
  des <- tree$edge[,2]-1
  
  # number of individuals (rows) observed
  nind <- nrow(trait_data)
  
  # is the data complete
  complete_data <- if(all(complete.cases(trait_data))) TRUE else FALSE
  if(complete_data) EM_missing_limit <- 1
  # which row of X and vectorized Y corresponds to a given observation (row)
  inds <- t((!is.na(trait_data[,1:nvar+1,drop=FALSE])) * 1:nind)
  inds[inds==0] <- NA
  inds[(1:length(inds))[!is.na(inds)]] <- 1:length(inds[(1:length(inds))[!is.na(inds)]])
  ind_list <- lapply(apply(inds,2,function(X) list(X[!is.na(X)]-1)),function(X) X[[1]])
  
  # design matrix for a vectorized Y (zeros and ones)
  X <- matrix(0,length(na.exclude(as.integer(inds))),nvar)
  
  # design matrix if X had complete data (for EM algorithm)
  X_complete <- matrix(0,nrow(trait_data)*nvar,nvar)
  for(i in 1:nind)
  {
    X_complete[1:nvar+(i-1)*nvar,] <- diag(nvar)
    if(length(species_subset[[i]])>1)
    {
      diag(X[inds[which(!is.na(inds[,i])),i],species_subset[[i]]+1]) <- 1
    } else X[inds[which(!is.na(inds[,i])),i],species_subset[[i]]+1] <- 1
  }
  
  # vectorized Y
  R <- matrix(na.exclude(as.double(t(trait_data[,1:nvar+1,drop=FALSE]))),ncol=1)
  
  # matrix form of Y
  Rmat <- trait_data[,1:nvar+1,drop=FALSE]
  
  # total number of observations (e.g., for complete data, nind*nvar)
  nob <- nrow(X)
  
  # edge lengths
  edge_vec <- tree$edge.length
  species_ind <- match(trait_data$species,tree$tip.label)
  
  # what edge number corrsponds to the descendant for a given individual (row)
  edge_ind <- match(species_ind,tree$edge[,2])-1
  
  # what individual (row number) corresponds to a given terminal edge
  ind_edge <- match(tree$edge[,2],species_ind)-1
  
  # map to link descendant nodes to parent nodes
  parent_edges <- match(tree$edge[,1],tree$edge[,2],nomatch = nrow(tree$edge)+1)
  edge_vec <- c(edge_vec,0)
  parent_edges <- c(parent_edges,parent_edges[length(parent_edges)])-1
  
  # does a given edge number correspond to a terminal edge
  is_edge_ind <- rep(FALSE,nedge+1)
  is_edge_ind[which(tree$edge[,2]<=nspecies)] <- TRUE
  
  # same as above, but if data were complete
  species_subset_complete <- rep(list(1:nvar-1),nind)
  un_species_subset_complete <- rep(list(),nind)
  subset_list_complete <- list(1:nvar-1)
  ind_list_complete <- lapply(apply(matrix(1:(nvar*nind)-1,nrow=nvar),2,function(X) list(X)),function(X) X[[1]])
  tip_combn_complete <- rep(0,nind)
  
  em_tol <- 1e-3
  
  # function for estimating phylogenetic/phenotypic covariance using EM and numerical optimization
  estim_pars <- function(edge_vec,do_optim=FALSE,phylocov=diag(nvar),phenocov=diag(nvar),mu=matrix(0))
  {
    tree$edge.length <- edge_vec[1:nedge]
    
    # function for numerical optimization
    # if the difference between the EM log-likelihood and the current log-likelihood is extremely high, reject (likely numerical precision issues)
    f <- function(pars,em_ll,R,Rmat,phylocov_fixed,phenocov_fixed)
    {
      if(!nested_optim & (model=="lambda" | model=="OU" | model=="EB" | model=="kappa" | model=="delta"))
      {
        edge_vec <- evec(((max(par_bounds) - min(par_bounds)) / (1+exp(-(pars[length(pars)])))) + min(par_bounds))
        tree$edge.length[1:nedge] <- edge_vec[1:nedge]
      }
      if(is.na(phylocov_fixed)[[1]]) phylocov <- pars_to_mat(pars[phylocov_pars],nvar,as.integer(!phylo_correlated)) else phylocov <- phylocov_fixed
      if(pheno_error>0) if(is.na(phenocov_fixed)[[1]]) phenocov <- pars_to_mat(pars[phenocov_pars],nvar,pheno_error) else phenocov <- phenocov_fixed
      pars <- numeric()
      if(is.na(phylocov_fixed)[[1]])
      {
        if(npd) phylocov <- try(as.matrix(nearPD(phylocov)$mat),silent=TRUE)
        if(class(phylocov)=="try-error") return(em_ll - max_delta)
        pars <- mat_to_pars(phylocov,nvar,as.integer(!phylo_correlated))
        if(class(pars)=="try-error") return(em_ll - max_delta)
      }
      if(is.na(phenocov_fixed)[[1]])
      {
        if(npd) phenocov <- try(as.matrix(nearPD(phenocov)$mat),silent=TRUE)
        if(class(phenocov)=="try-error") return(em_ll - max_delta)
        if(pheno_error>0) pars <- c(pars,mat_to_pars(phenocov,nvar,pheno_error))
        if(class(pars)=="try-error") return(em_ll - max_delta)
      }
      
      ll <- try(tp(L=X, R=R, Rmat = as.matrix(Rmat),mL=ncol(X), mR=1, pheno_error=pheno_error, edge_vec=edge_vec, 
                   edge_ind=edge_ind,ind_edge=ind_edge, parent_edges = parent_edges,pars=pars, nvar=nvar, 
                   phylocov_diag=as.integer(!phylo_correlated), nind=nind, nob=nob, nspecies=nspecies, nedge=nedge, anc=anc, des=des, REML=as.integer(REML), 
                   species_subset=species_subset, un_species_subset = un_species_subset,subset_list=subset_list,
                   ind_list=ind_list, tip_combn=tip_combn,is_edge_ind=is_edge_ind,fixed_mu=matrix(0),ret_level=1,
                   is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=phylocov_fixed,
                   is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=phenocov_fixed,OU_len=list())$logl[[1]],silent=TRUE)
      if(class(ll)=="try-error") ll <- em_ll-max_delta
      if(is.na(ll)) ll <- em_ll - max_delta
      if(abs(abs(em_ll)-abs(ll))>max_delta) ll <- em_ll - max_delta
      if(ll==(em_ll-max_delta)) return(ll)
      
      temp_pars <- numeric()
      if(is.na(phylocov_fixed)[[1]])
      {
        temp_pars <- mat_to_pars(phylocov+diag(nvar)*1e-6,nvar,as.integer(!phylo_correlated))
        if(class(pars)=="try-error") return(em_ll - max_delta)
      }
      if(is.na(phenocov_fixed)[[1]])
      {
        if(pheno_error>0) temp_pars <- c(temp_pars,mat_to_pars(phenocov+diag(nvar)*1e-6,nvar,pheno_error))
        if(class(pars)=="try-error") return(em_ll - max_delta)
      }
      
      temp_ll2 <- try(tp(L=X, R=R, Rmat = as.matrix(Rmat),mL=ncol(X), mR=1, pheno_error=pheno_error, edge_vec=edge_vec, 
                         edge_ind=edge_ind,ind_edge=ind_edge, parent_edges = parent_edges,pars=temp_pars, nvar=nvar, 
                         phylocov_diag=as.integer(!phylo_correlated), nind=nind, nob=nob, nspecies=nspecies, nedge=nedge, anc=anc, des=des, REML=as.integer(REML), 
                         species_subset=species_subset, un_species_subset = un_species_subset,subset_list=subset_list,
                         ind_list=ind_list, tip_combn=tip_combn,is_edge_ind=is_edge_ind,fixed_mu=matrix(0),ret_level=1,
                         is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=phylocov_fixed,
                         is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=phenocov_fixed,OU_len=list())$logl[[1]],silent=TRUE)[[1]]
      
      
      if(class(temp_ll2)!="try-error") if(!is.na(temp_ll2)) if(abs(temp_ll2-ll[[1]])>100)
      {
        ll <- em_ll-max_delta
      }
      
      if(optim_verbose) cat(c(ll,"\n"))
      ll
    }
    if(!skip_EM | complete_data)
    {
      # if estimating phenotypic covariance, use Felsenstein 2008 EM algorithm
      # if missing data, nest Felsenstein 2008 EM algorithm within missing data EM algorithm
      if(pheno_error!=0)
      {
        pars <- numeric()
        if(is.na(phylocov_fixed)[[1]]) pars <- mat_to_pars(phylocov,nvar,as.integer(!phylo_correlated))
        if(is.na(phenocov_fixed)[[1]]) pars <- c(pars,mat_to_pars(phenocov,nvar,pheno_error))
        imputed <- tp(L=X, R=R, Rmat = as.matrix(Rmat),mL=ncol(X), mR=1, pheno_error=pheno_error, edge_vec=edge_vec, 
                      edge_ind=edge_ind,ind_edge=ind_edge, parent_edges = parent_edges,pars=pars, nvar=nvar, 
                      phylocov_diag=as.integer(!phylo_correlated), nind=nind, nob=nob, nspecies=nspecies, nedge=nedge, anc=anc, des=des, REML=as.integer(REML), 
                      species_subset=species_subset, un_species_subset = un_species_subset,subset_list=subset_list,
                      ind_list=ind_list, tip_combn=tip_combn,is_edge_ind=is_edge_ind,fixed_mu=mu,ret_level=3,
                      is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=phylocov_fixed,
                      is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=phenocov_fixed,OU_len=list())
        imputed_dat <- data.frame(species=trait_data$species,imputed$recon_ind)
        
        args <- prep_em(imputed_dat,tree)
        AP <- EM_Fels2008(pics=args[[1]],vars=args[[2]],phylocov=phylocov,phenocov=phenocov,nvar=nvar,tol=em_tol,diag_pheno=as.integer(!pheno_correlated),EM_Fels_limit=EM_Fels_limit,REML=as.integer(REML),diag_phylo=as.integer(!phylo_correlated),
                          is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=phylocov_fixed,is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=phenocov_fixed)
        best_phylocov <- phylocov <- AP[1:nvar,1:nvar,drop=FALSE]
        best_phenocov <- phenocov <- AP[1:nvar+nvar,1:nvar,drop=FALSE]
        
        pars <- numeric()
        if(is.na(phylocov_fixed)[[1]]) pars <- mat_to_pars(phylocov,nvar,as.integer(!phylo_correlated))
        if(is.na(phenocov_fixed)[[1]]) pars <- c(pars,mat_to_pars(phenocov,nvar,pheno_error))
        
        ll2 <- tp(L=X, R=R, Rmat = as.matrix(Rmat),mL=ncol(X), mR=1, pheno_error=pheno_error, edge_vec=edge_vec, 
                  edge_ind=edge_ind,ind_edge=ind_edge, parent_edges = parent_edges,pars=pars, nvar=nvar, 
                  phylocov_diag=as.integer(!phylo_correlated), nind=nind, nob=nob, nspecies=nspecies, nedge=nedge, anc=anc, des=des, REML=as.integer(REML), 
                  species_subset=species_subset, un_species_subset = un_species_subset,subset_list=subset_list,
                  ind_list=ind_list, tip_combn=tip_combn,is_edge_ind=is_edge_ind,fixed_mu=matrix(0),ret_level=2,
                  is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=phylocov_fixed,
                  is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=phenocov_fixed,OU_len=list())
        
        mu <- ll2$theta
        ll <- ll2[[1]]
        
        
        for(i in 1:EM_missing_limit)
        {
          imputed <- tp(L=X, R=R, Rmat = as.matrix(Rmat),mL=ncol(X), mR=1, pheno_error=pheno_error, edge_vec=edge_vec, 
                        edge_ind=edge_ind,ind_edge=ind_edge, parent_edges = parent_edges,pars=pars, nvar=nvar, 
                        phylocov_diag=as.integer(!phylo_correlated), nind=nind, nob=nob, nspecies=nspecies, nedge=nedge, anc=anc, des=des, REML=as.integer(REML), 
                        species_subset=species_subset, un_species_subset = un_species_subset,subset_list=subset_list,
                        ind_list=ind_list, tip_combn=tip_combn,is_edge_ind=is_edge_ind,fixed_mu=mu,ret_level=3,
                        is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=phylocov_fixed,
                        is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=phenocov_fixed,OU_len=list())
          imputed_dat <- data.frame(species=trait_data$species,imputed$recon_ind)
          args <- prep_em(imputed_dat,tree)
          AP <- EM_Fels2008(pics=args[[1]],vars=args[[2]],phylocov=phylocov,phenocov=phenocov,nvar=nvar,tol=em_tol,diag_pheno=as.integer(!pheno_correlated),EM_Fels_limit=EM_Fels_limit,REML=as.integer(REML),diag_phylo=as.integer(!phylo_correlated),
                            is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=phylocov_fixed,is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=phenocov_fixed)
          AP[1:nvar+nvar,1:nvar] <- AP[1:nvar+nvar,1:nvar,drop=FALSE]
          phylocov <- AP[1:nvar,1:nvar,drop=FALSE]
          phenocov <- AP[1:nvar+nvar,1:nvar,drop=FALSE]
          pars <- numeric()
          if(is.na(phylocov_fixed)[[1]]) pars <- mat_to_pars(phylocov,nvar,as.integer(!phylo_correlated))
          if(is.na(phenocov_fixed)[[1]]) pars <- c(pars,mat_to_pars(phenocov,nvar,pheno_error))
          AP <- EM_Fels2008(pics=args[[1]],vars=args[[2]],phylocov=phylocov,phenocov=phenocov,nvar=nvar,tol=em_tol,diag_pheno=as.integer(!pheno_correlated),EM_Fels_limit=EM_Fels_limit,REML=as.integer(REML),diag_phylo=as.integer(!phylo_correlated),
                            is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=phylocov_fixed,is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=phenocov_fixed)
          AP[1:nvar+nvar,1:nvar] <- AP[1:nvar+nvar,1:nvar,drop=FALSE] + imputed$tip_uncertainty/(nind-nvar)
          phylocov <- AP[1:nvar,1:nvar,drop=FALSE]
          if(is.na(phenocov_fixed)[[1]]) phenocov <- AP[1:nvar+nvar,1:nvar,drop=FALSE]
          pars <- numeric()
          if(is.na(phylocov_fixed)[[1]]) pars <- mat_to_pars(phylocov,nvar,as.integer(!phylo_correlated))
          if(is.na(phenocov_fixed)[[1]]) pars <- c(pars,mat_to_pars(phenocov,nvar,pheno_error))
          mu <- tp(L=X_complete, R=matrix(as.double(t(imputed$recon_ind)),ncol=1), Rmat = imputed$recon_ind,mL=ncol(X), mR=1, pheno_error=pheno_error, edge_vec=edge_vec, 
                   edge_ind=edge_ind,ind_edge=ind_edge, parent_edges = parent_edges,pars=pars, nvar=nvar, 
                   phylocov_diag=as.integer(!phylo_correlated), nind=nind, nob=nob, nspecies=nspecies, nedge=nedge, anc=anc, des=des, REML=as.integer(REML), 
                   species_subset=species_subset_complete, un_species_subset = un_species_subset_complete,subset_list=subset_list_complete,
                   ind_list=ind_list_complete, tip_combn=tip_combn_complete,is_edge_ind=is_edge_ind,fixed_mu=matrix(0),ret_level=2,
                   is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=phylocov_fixed,
                   is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=phenocov_fixed,OU_len=list())$theta
          ll2 <- tp(L=X, R=R, Rmat = as.matrix(Rmat),mL=ncol(X), mR=1, pheno_error=pheno_error, edge_vec=edge_vec, 
                    edge_ind=edge_ind,ind_edge=ind_edge, parent_edges = parent_edges,pars=pars, nvar=nvar, 
                    phylocov_diag=as.integer(!phylo_correlated), nind=nind, nob=nob, nspecies=nspecies, nedge=nedge, anc=anc, des=des, REML=as.integer(REML), 
                    species_subset=species_subset, un_species_subset = un_species_subset,subset_list=subset_list,
                    ind_list=ind_list, tip_combn=tip_combn,is_edge_ind=is_edge_ind,fixed_mu=matrix(0),ret_level=1,
                    is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=phylocov_fixed,
                    is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=phenocov_fixed,OU_len=list())$logl
          if(ll2<ll)
          {
            phylocov <- best_phylocov
            phenocov <- best_phenocov
            ll2 <- ll
            break
          } else
          {
            best_phylocov <- phylocov
            best_phenocov <- phenocov
            ll <- ll2
            if(EM_verbose) cat(c(ll2,"\n"))
          }
        }
      } else if(is.na(phylocov_fixed)[[1]])
      {
        # EM algorithm for missing data (no phenotypic error)
        pars <- mat_to_pars(phylocov,nvar,as.integer(!phylo_correlated))
        mu <- matrix(colMeans(trait_data[,1:nvar+1,drop=FALSE],na.rm = TRUE))
        dat <- as.matrix(trait_data[,1:nvar+1,drop=FALSE])
        rownames(dat) <- rownames(trait_data$species)
        for(i in 1:ncol(Rmat)) dat[is.na(dat[,i]),i] <- mu[i]
        mu <- tp(L=X_complete, R=as.matrix(as.double(dat)), Rmat = as.matrix(dat),mL=ncol(X), mR=1, pheno_error=pheno_error, edge_vec=edge_vec, 
                 edge_ind=edge_ind,ind_edge=ind_edge, parent_edges = parent_edges,pars=pars, nvar=nvar, 
                 phylocov_diag=as.integer(!phylo_correlated), nind=nind, nob=nob, nspecies=nspecies, nedge=nedge, anc=anc, des=des, REML=as.integer(REML), 
                 species_subset=species_subset_complete, un_species_subset = un_species_subset_complete,subset_list=subset_list_complete,
                 ind_list=ind_list_complete, tip_combn=tip_combn_complete,is_edge_ind=is_edge_ind,fixed_mu=matrix(0),ret_level=2,
                 is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=phylocov_fixed,
                 is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=phenocov_fixed,OU_len=list())$theta
        best_phylocov <- phylocov <- crossprod(apply(dat,2,pic,phy=multi2di(tree,random=FALSE)))/(nspecies-REML)
        
        if(!phylo_correlated) best_phylocov <- phylocov <- diag(diag(phylocov))
        pars <- mat_to_pars(phylocov,nvar,as.integer(!phylo_correlated))
        ll <- ll2 <- tp(L=X, R=R, Rmat = as.matrix(Rmat),mL=ncol(X), mR=1, pheno_error=pheno_error, edge_vec=edge_vec, 
                        edge_ind=edge_ind,ind_edge=ind_edge, parent_edges = parent_edges,pars=pars, nvar=nvar, 
                        phylocov_diag=as.integer(!phylo_correlated), nind=nind, nob=nob, nspecies=nspecies, nedge=nedge, anc=anc, des=des, REML=as.integer(REML), 
                        species_subset=species_subset, un_species_subset = un_species_subset,subset_list=subset_list,
                        ind_list=ind_list, tip_combn=tip_combn,is_edge_ind=is_edge_ind,fixed_mu=matrix(0),ret_level=1,
                        is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=phylocov_fixed,
                        is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=phenocov_fixed,OU_len=list())$logl
        
        for(k in 1:EM_missing_limit)
        {
          imputed <- tp(L=X, R=R, Rmat = as.matrix(Rmat),mL=ncol(X), mR=1, pheno_error=pheno_error, edge_vec=edge_vec, 
                        edge_ind=edge_ind,ind_edge=ind_edge, parent_edges = parent_edges,pars=pars, nvar=nvar, 
                        phylocov_diag=as.integer(!phylo_correlated), nind=nind, nob=nob, nspecies=nspecies, nedge=nedge, anc=anc, des=des, REML=as.integer(REML), 
                        species_subset=species_subset, un_species_subset = un_species_subset,subset_list=subset_list,
                        ind_list=ind_list, tip_combn=tip_combn,is_edge_ind=is_edge_ind,fixed_mu=mu,ret_level=3,
                        is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=phylocov_fixed,
                        is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=phenocov_fixed,OU_len=list())
          dat <- imputed$recon_ind
          rownames(dat) <- rownames(trait_data$species)
          biased_sigma <- crossprod(apply(dat,2,pic,phy=multi2di(tree,random=FALSE)))/(nspecies-REML)
          biased_pars <- mat_to_pars(biased_sigma,nvar,as.integer(!phylo_correlated))
          tip_un <- tp(L=X, R=R, Rmat = as.matrix(Rmat),mL=ncol(X), mR=1, pheno_error=pheno_error, edge_vec=edge_vec, 
                       edge_ind=edge_ind,ind_edge=ind_edge, parent_edges = parent_edges,pars=biased_pars, nvar=nvar, 
                       phylocov_diag=as.integer(!phylo_correlated), nind=nind, nob=nob, nspecies=nspecies, nedge=nedge, anc=anc, des=des, REML=as.integer(REML), 
                       species_subset=species_subset, un_species_subset = un_species_subset,subset_list=subset_list,
                       ind_list=ind_list, tip_combn=tip_combn,is_edge_ind=is_edge_ind,fixed_mu=matrix(0),ret_level=3,
                       is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=phylocov_fixed,
                       is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=phenocov_fixed,OU_len=list())$tip_uncertainty
          mu <- as.matrix(apply(dat,2,function(X) ace(X,multi2di(tree,random=FALSE),method="pic")$ace[[1]]))
          phylocov <- biased_sigma + tip_un/(nspecies-REML)
          if(!phylo_correlated) phylocov <- diag(diag(phylocov))
          pars <- mat_to_pars(phylocov,nvar,as.integer(!phylo_correlated))
          ll2 <- tp(L=X, R=R, Rmat = as.matrix(Rmat),mL=ncol(X), mR=1, pheno_error=pheno_error, edge_vec=edge_vec, 
                    edge_ind=edge_ind,ind_edge=ind_edge, parent_edges = parent_edges,pars=pars, nvar=nvar, 
                    phylocov_diag=as.integer(!phylo_correlated), nind=nind, nob=nob, nspecies=nspecies, nedge=nedge, anc=anc, des=des, REML=as.integer(REML), 
                    species_subset=species_subset, un_species_subset = un_species_subset,subset_list=subset_list,
                    ind_list=ind_list, tip_combn=tip_combn,is_edge_ind=is_edge_ind,fixed_mu=matrix(0),ret_level=1,
                    is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=phylocov_fixed,
                    is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=phenocov_fixed,OU_len=list())$logl
          if(ll2<ll)
          {
            phylocov <- best_phylocov
            #phenocov <- best_phenocov
            ll2 <- ll
            break
          } else
          {
            best_phylocov <- phylocov
            #best_phenocov <- phenocov
            ll <- ll2
            if(EM_verbose) cat(c(ll2,"\n"))
          }
          if(EM_verbose) cat(c(ll2,"\n"))
        }
      }
    }
    if(model=="BM" | model=="white" | model=="star")
    {
      if(!is.na(phylocov_start)[[1]] & is.na(phylocov_fixed)[[1]]) phylocov <- phylocov_start
      if(!is.na(phenocov_start)[[1]] & is.na(phenocov_fixed)[[1]]) phenocov <- phenocov_start
    }
    
    pars <- numeric()
    if(is.na(phylocov_fixed)[[1]]) pars <- mat_to_pars(phylocov,nvar,as.integer(!phylo_correlated))
    if(is.na(phenocov_fixed)[[1]]) if(pheno_error!=0) pars <- c(pars,mat_to_pars(phenocov,nvar,pheno_error))
    
    ll2 <- tp(L=X, R=R, Rmat = as.matrix(Rmat),mL=ncol(X), mR=1, pheno_error=pheno_error, edge_vec=edge_vec, 
              edge_ind=edge_ind,ind_edge=ind_edge, parent_edges = parent_edges,pars=pars, nvar=nvar, 
              phylocov_diag=as.integer(!phylo_correlated), nind=nind, nob=nob, nspecies=nspecies, nedge=nedge, anc=anc, des=des, REML=as.integer(REML), 
              species_subset=species_subset, un_species_subset = un_species_subset,subset_list=subset_list,
              ind_list=ind_list, tip_combn=tip_combn,is_edge_ind=is_edge_ind,fixed_mu=matrix(0),ret_level=3,
              is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=phylocov_fixed,
              is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=phenocov_fixed,OU_len=list())
    mu <- ll2[[2]]
    if(!(!do_optim | skip_optim))
    {
      tedge_vec <- edge_vec
      # numerical optimization
      ll <- ll2[[1]]
      if(!nested_optim & (model=="lambda" | model=="OU" | model=="EB" | model=="kappa" | model=="delta")) MOD_PAR <- -log(-(-max(par_bounds)+mean(par_bounds))/(-min(par_bounds)+mean(par_bounds)))
      if((complete_data & (model!="BM" & model!="white" & model!="star" & is.na(model_par_fixed[[1]]))) | pheno_error | !complete_data)
        for(i in 1:repeat_optim_limit)
        {
          zRmat <- (trait_data[,1:nvar+1,drop=FALSE] - matrix(1,nind) %*% t(zstand$mu)) / sqrt(matrix(1,nind) %*% diag(zstand$phylocov))
          zR <- matrix(na.exclude(as.double(t(zRmat[,1:nvar]))),ncol=1)
          pars2 <- numeric()
          if(is.na(phylocov_fixed)[[1]]) pars2 <- mat_to_pars(phylocov / (sqrt(diag(zstand$phylocov)) %*% t(sqrt(diag(zstand$phylocov)))),nvar,as.integer(!phylo_correlated))
          if(is.na(phenocov_fixed)[[1]]) if(pheno_error!=0) pars2 <- c(pars2,mat_to_pars(phenocov / (sqrt(diag(zstand$phylocov)) %*% t(sqrt(diag(zstand$phylocov)))),nvar,pheno_error))
          
          z_em_ll <- tp(L=X, R=zR, Rmat = as.matrix(zRmat),mL=ncol(X), mR=1, pheno_error=pheno_error, edge_vec=edge_vec, 
                        edge_ind=edge_ind,ind_edge=ind_edge, parent_edges = parent_edges,pars=pars2, nvar=nvar, 
                        phylocov_diag=as.integer(!phylo_correlated), nind=nind, nob=nob, nspecies=nspecies, nedge=nedge, anc=anc, des=des, REML=as.integer(REML), 
                        species_subset=species_subset, un_species_subset = un_species_subset,subset_list=subset_list,
                        ind_list=ind_list, tip_combn=tip_combn,is_edge_ind=is_edge_ind,fixed_mu=matrix(0),ret_level=1,
                        is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=zphylocov_fixed,
                        is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=zphenocov_fixed,OU_len=list())$logl
          if(!nested_optim & (model=="lambda" | model=="OU" | model=="EB" | model=="kappa" | model=="delta")) pars2 <- c(pars2,MOD_PAR)
          pars2 <- optim(par = pars2,fn = f,control=list(fnscale=-1),method=if(i%%2 == 1) "BFGS" else "Nelder-Mead",em_ll=z_em_ll,R=zR,Rmat=zRmat,phylocov_fixed=zphylocov_fixed,phenocov_fixed=zphenocov_fixed)
          if(is.na(phylocov_fixed)[[1]]) phylocov2 <- pars_to_mat(pars2$par[phylocov_pars],nvar,as.integer(!phylo_correlated)) else phylocov2 <- zphylocov_fixed
          if(is.na(phenocov_fixed)[[1]])
          {
            if(pheno_error!=0) phenocov2 <- pars_to_mat(pars2$par[phenocov_pars],nvar,pheno_error) else phenocov2 <- zphenocov_fixed
          }
          if(!nested_optim & (model=="lambda" | model=="OU" | model=="EB" | model=="kappa" | model=="delta"))
          {
            MOD_PAR <- pars2$par[length(pars2$par)]
            tedge_vec <- evec(((max(par_bounds) - min(par_bounds)) / (1+exp(-(MOD_PAR)))) + min(par_bounds))
          } else tedge_vec <- edge_vec
          pars2 <- numeric()
          if(is.na(phylocov_fixed)[[1]]) pars2 <- mat_to_pars(phylocov2 * (sqrt(diag(zstand$phylocov)) %*% t(sqrt(diag(zstand$phylocov)))),nvar,as.integer(!phylo_correlated))
          if(is.na(phenocov_fixed)[[1]]) if(pheno_error!=0) pars2 <- c(pars2,mat_to_pars(phenocov2 * (sqrt(diag(zstand$phylocov)) %*% t(sqrt(diag(zstand$phylocov)))),nvar,pheno_error))
          
          ll2 <- tp(L=X, R=R, Rmat = as.matrix(Rmat),mL=ncol(X), mR=1, pheno_error=pheno_error, edge_vec=tedge_vec, 
                    edge_ind=edge_ind,ind_edge=ind_edge, parent_edges = parent_edges,pars=pars2, nvar=nvar, 
                    phylocov_diag=as.integer(!phylo_correlated), nind=nind, nob=nob, nspecies=nspecies, nedge=nedge, anc=anc, des=des, REML=as.integer(REML), 
                    species_subset=species_subset, un_species_subset = un_species_subset,subset_list=subset_list,
                    ind_list=ind_list, tip_combn=tip_combn,is_edge_ind=is_edge_ind,fixed_mu=matrix(0),ret_level=1,
                    is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=phylocov_fixed,
                    is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=phenocov_fixed,OU_len=list())
          if(is.na(phenocov_fixed)[[1]]) if(pheno_error!=0) phenocov <- phenocov2 * (sqrt(diag(zstand$phylocov)) %*% t(sqrt(diag(zstand$phylocov))))
          if(is.na(phylocov_fixed)[[1]]) phylocov <- phylocov2 * (sqrt(diag(zstand$phylocov)) %*% t(sqrt(diag(zstand$phylocov))))
          pars <- numeric()
          if(is.na(phylocov_fixed)[[1]]) pars <- mat_to_pars(phylocov,nvar,as.integer(!phylo_correlated))
          if(is.na(phenocov_fixed)[[1]]) if(pheno_error!=0) pars <- c(pars,mat_to_pars(phenocov,nvar,pheno_error))
          if(i>1) if(abs(ll-ll2$logl)<repeat_optim_tol) break else ll <- ll2$logl
        }
      
      ll2 <- tp(L=X, R=R, Rmat = as.matrix(Rmat),mL=ncol(X), mR=1, pheno_error=pheno_error, edge_vec=tedge_vec, 
                edge_ind=edge_ind,ind_edge=ind_edge, parent_edges = parent_edges,pars=pars, nvar=nvar, 
                phylocov_diag=as.integer(!phylo_correlated), nind=nind, nob=nob, nspecies=nspecies, nedge=nedge, anc=anc, des=des, REML=as.integer(REML), 
                species_subset=species_subset, un_species_subset = un_species_subset,subset_list=subset_list,
                ind_list=ind_list, tip_combn=tip_combn,is_edge_ind=is_edge_ind,fixed_mu=matrix(0),ret_level=3,
                is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=phylocov_fixed,
                is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=phenocov_fixed,OU_len=list())
    }
    ret <- list(ll=ll2$logl,phylocov=phylocov,phenocov=phenocov,mu=ll2[[2]],pars=pars,ll2=ll2)
    if(!nested_optim & (model=="lambda" | model=="OU" | model=="EB" | model=="kappa" | model=="delta"))
    {
      ret[model] <- ((max(par_bounds) - min(par_bounds)) / (1+exp(-(MOD_PAR)))) + min(par_bounds)
    }
    ret
  }
  if(nested_optim & model!="BM" & model!="white" & model!="star" & model!="mvOU")
  {
    max_i <- max_logl <- -Inf
    max_par <- mean(par_bounds)
    max_phylocov <- phylocov
    max_phenocov <- phenocov
    max_mu <- mu
    
    
    if(!is.na(model_par_fixed)[[1]])
    {
      tree$edge.length <- evec(model_par_fixed)[1:nedge]
      ret <- phylopars(trait_data=trait_data,tree=tree,pheno_error=pheno_error>0,pheno_correlated=pheno_error==2,phylo_correlated=phylo_correlated,
                       REML=REML,EM_Fels_limit=EM_Fels_limit,repeat_optim_limit=repeat_optim_limit,EM_missing_limit=EM_missing_limit,repeat_optim_tol=repeat_optim_tol,phylocov_start=phylocov_start,
                       phenocov_start=phenocov_start,phylocov_fixed=phylocov_fixed,phenocov_fixed=phenocov_fixed,skip_optim=skip_optim,skip_EM=skip_EM,EM_verbose=EM_verbose,optim_verbose=optim_verbose)
      if(!is.na(perm_model))
      {
        if(perm_model=="lambda" | perm_model=="OU" | perm_model=="EB" | perm_model=="kappa" | perm_model=="delta")
        {
          ret$model <- list(model=perm_model,mod_par=model_par_fixed)
          names(ret$model)[2] <- if(ret$model[[1]]=="OU") "alpha" else if(ret$model[[1]]=="EB") "rate" else ret$model[[1]]
          ret$npars <- ret$npars + 1
          
          if(perm_model=="OU")
          {
            stationary_cov <- ret$pars[[1]]/(2*ret$model[[2]])
            dimnames(stationary_cov) <-   list(colnames(trait_data)[1:nvar+1],colnames(trait_data)[1:nvar+1])
            ret$model <- c(ret$model,list(stationary_cov=stationary_cov))
          }
          return(ret)
        }
      }
    }
    
    temp_edge_vec <- edge_vec
    for(j in 1:1)
    {
      for(i in 1:model_par_evals)
      {
        temp_edge_vec <- evec(par_bounds[i])
        pars <- try(estim_pars(temp_edge_vec,do_optim = FALSE),silent=TRUE)
        if(class(pars)=="try-error")
          pars <- try(estim_pars(temp_edge_vec,do_optim = TRUE,skip_EM=TRUE),silent=TRUE)
        if(class(pars)=="try-error")
        {
          next
        }
        if(is.na(phylocov_fixed)[[1]]) phylocov <- pars$phylocov
        if(is.na(phenocov_fixed)[[1]]) phenocov <- pars$phenocov
        mu <- pars$mu
        if(i==1 | is.infinite(max_logl))
        {
          max_logl <- pars[[1]]
          max_i <- i
          max_par <- par_bounds[i]
          max_phylocov <- phylocov
          max_phenocov <- phenocov
          max_mu <- mu
        }
        if(pars[[1]] > max_logl)
        {
          max_logl <- pars[[1]]
          max_i <- i
          max_par <- par_bounds[i]
          max_phylocov <- phylocov
          max_phenocov <- phenocov
          max_mu <- mu
        }
      }
      max_i <- c(max_i-2,max_i+2)
      max_i[max_i<1] <- 1
      max_i[max_i>model_par_evals] <- model_par_evals
      if(max_i[1]==1 & max_i[2]==model_par_evals) break
      par_bounds <- seq(min(par_bounds[max_i]),max(par_bounds[max_i]),length=model_par_evals)
      if(max_par < min(par_bounds)) par_bounds[which.min(par_bounds)] <- max_par
      if(max_par > max(par_bounds)) par_bounds[which.max(par_bounds)] <- max_par
      par_bounds <- seq(min(par_bounds),max(par_bounds),length=model_par_evals)
    }
    if(!is.na(model_par_start)[[1]]) max_par <- model_par_start
    a <- -log(-(-max(par_bounds)+max_par)/(-min(par_bounds)+max_par))
    if(is.infinite((a))) a <- -log(-(-max(par_bounds)+mean(par_bounds))/(-min(par_bounds)+mean(par_bounds)))
    
    o1 <- optim(par = a,fn = function(X) estim_pars(evec(((max(par_bounds) - min(par_bounds)) / (1+exp(-(X)))) + min(par_bounds)),do_optim = TRUE,phylocov = max_phylocov,phenocov = max_phenocov,mu = max_mu)[[1]],control=list(fnscale=-1),method="BFGS")
    #o1 <- optimize(f = function(X) estim_pars(evec(X),do_optim = TRUE,phylocov = phylocov,phenocov = phenocov,mu = mu)[[1]],interval = c(min(par_bounds[max_i]),max(par_bounds[max_i])),maximum = TRUE)
    o2 <- estim_pars(evec(min(par_bounds[max_i])),do_optim = TRUE,phylocov=max_phylocov,phenocov=max_phenocov,mu=max_mu)
    o3 <- estim_pars(evec(max(par_bounds[max_i])),do_optim = TRUE,phylocov=max_phylocov,phenocov=max_phenocov,mu=max_mu)
    
    w <- which.max(c(o1$value,o2$ll,o3$ll))
    if(w==1)
    {
      model_par <- ((max(par_bounds) - min(par_bounds)) / (1+exp(-(o1$par)))) + min(par_bounds)
      pars <- estim_pars(evec(model_par),do_optim = TRUE,phylocov = max_phylocov,phenocov=max_phenocov,mu = max_mu)
      ll2 <- tp(L=X, R=R, Rmat = as.matrix(Rmat),mL=ncol(X), mR=1, pheno_error=pheno_error, edge_vec=evec(model_par), 
                edge_ind=edge_ind,ind_edge=ind_edge, parent_edges = parent_edges,pars=pars$pars, nvar=nvar, 
                phylocov_diag=as.integer(!phylo_correlated), nind=nind, nob=nob, nspecies=nspecies, nedge=nedge, anc=anc, des=des, REML=as.integer(REML), 
                species_subset=species_subset, un_species_subset = un_species_subset,subset_list=subset_list,
                ind_list=ind_list, tip_combn=tip_combn,is_edge_ind=is_edge_ind,fixed_mu=matrix(0),ret_level=3,
                is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=phylocov_fixed,
                is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=phenocov_fixed,OU_len=list())
    } else if(w==2)
    {
      model_par <- min(par_bounds[max_i])
      ll2 <- o2$ll2
    } else if(w==3)
    {
      model_par <- max(par_bounds[max_i])
      ll2 <- o3$ll2
    }
    ll2[model] <- model_par
    
    if((model=="OU" & round(model_par-0,6)==0) | (model=="EB" & round(model_par-0,6)==0) | (model=="lambda" & round(model_par-1,6)==0) | (model=="delta" & round(model_par-1,6)==0) | (model=="kappa" & round(model_par-1,6)==0)) warning(paste("Estimated",model,"parameter is at bounds and is indistinguishable from a Brownian Motion model."))
    ret <- list(ll=ll2$logl,phylocov=phylocov,phenocov=phenocov,mu=ll2[[2]],pars=pars,ll2=ll2)
    if(model=="lambda" | model=="OU" | model=="EB" | model=="kappa" | model=="delta")
    {
      ret[model] <- ((max(par_bounds) - min(par_bounds)) / (1+exp(-(model_par)))) + min(par_bounds)
    }
    ll2 <- ret
    
  } else if((model=="BM" | model=="white" | model=="star") | (!nested_optim & model!="mvOU")) ll2 <- estim_pars(edge_vec = edge_vec,do_optim = TRUE,phylocov = phylocov,phenocov = phenocov,mu = mu) else 
    if(model=="mvOU")
    {
      hts <- getHeights(tree)
      edge_mat <- postorder_tools(tree)
      des_order <- edge_mat[[2]]
      edge_mat <- edge_mat[[1]]
      
      
      OU_fun <- function(pars,R,Rmat,phylocov_fixed,phenocov_fixed,ret_level=1,BM_ll=NA)
      {
        if(is.na(phylocov_fixed)[[1]]) phylocov <- pars_to_mat(pars[phylocov_pars],nvar,as.integer(!phylo_correlated)) else phylocov <- phylocov_fixed
        if(pheno_error>0) if(is.na(phenocov_fixed)[[1]]) phenocov <- pars_to_mat(pars[phenocov_pars],nvar,pheno_error) else phenocov <- phenocov_fixed
        if(is.na(model_par_fixed)[[1]]) alpha <- pars_to_mat(pars[alpha_pars],nvar,diag = abs(full_alpha-1))
        
        if(is.na(phylocov_fixed)[[1]])
        {
          if(npd) phylocov <- try(as.matrix(nearPD(phylocov)$mat),silent=TRUE)
          if(class(phylocov)=="try-error") return(phylocov)
          pars <- mat_to_pars(phylocov,nvar,as.integer(!phylo_correlated))
          if(class(pars)=="try-error") return(pars)
        }
        if(pheno_error>0) if(is.na(phenocov_fixed)[[1]])
        {
          if(npd) phenocov <- try(as.matrix(nearPD(phenocov)$mat),silent=TRUE)
          if(class(phenocov)=="try-error") return(phenocov)
          pars <- c(pars,mat_to_pars(phenocov,nvar,pheno_error))
          if(class(pars)=="try-error") return(pars)
        }
        if(npd) alpha <- try(as.matrix(nearPD(alpha)$mat),silent=TRUE)
        if(class(alpha)=="try-error") return(alpha)
        
        #HOU <- calc_OU_heights(heights = hts,edge_mat = edge_mat,des_order = des_order,nedge = nedge,alpha = alpha,sigma = phylocov)
        #len_vec <- apply(HOU,2,function(X) convert_to_tree(heights = X,nspecies = nspecies,nedge = nedge,anc = tree$edge[,1],des = tree$edge[,2]))
        #len_vec <- lapply(apply(len_vec,1,function(X) list(matrix(X,nvar,nvar))),function(X) X[[1]])
        pars <- numeric()
        if(is.na(phylocov_fixed)[[1]]) pars <- mat_to_pars(phylocov,nvar,as.integer(!phylo_correlated))
        if(!is.na(BM_ll)) if(class(pars)=="try-error") return(BM_ll - max_delta)
        if(pheno_error>0) if(is.na(phenocov_fixed)[[1]]) pars <- c(pars,mat_to_pars(phenocov,nvar,pheno_error))
        if(is.na(model_par_fixed)[[1]]) pars <- c(pars,mat_to_pars(alpha,nvar,abs(full_alpha-1)))
        if(!is.na(BM_ll)) if(class(pars)=="try-error") return(BM_ll - max_delta)
        ea <- try(eigen(alpha),silent=TRUE)
        if(class(ea)=="try-error") return(ea)
        P <- ea$vectors
        lambda <- ea$values
        if(any(lambda<=1e-8))
        {
          class(alpha) <- "try-error"
          return(alpha)
        }
        
        #temp_phylocov <- phylocov + diag(diag(phylocov)*1e-12)
        temp_phylocov <- phylocov + diag(nvar)*1e-6
        temp_len_vec <- calc_OU_len(heights=hts,edge_mat=edge_mat,des_order=des_order,nedge=nedge,P=P,lambda=lambda,sigma=temp_phylocov,anc=tree$edge[,1],des=tree$edge[,2],nvar=nvar,nspecies=nspecies)
        temp_ll2 <- tp(L=X, R=R, Rmat = as.matrix(Rmat),mL=ncol(X), mR=1, pheno_error=pheno_error, edge_vec=edge_vec, 
                       edge_ind=edge_ind,ind_edge=ind_edge, parent_edges = parent_edges,pars=pars, nvar=nvar, 
                       phylocov_diag=as.integer(!phylo_correlated), nind=nind, nob=nob, nspecies=nspecies, nedge=nedge, anc=anc, des=des, REML=as.integer(REML), 
                       species_subset=species_subset, un_species_subset = un_species_subset,subset_list=subset_list,
                       ind_list=ind_list, tip_combn=tip_combn,is_edge_ind=is_edge_ind,fixed_mu=matrix(0),ret_level=1,OU_len=temp_len_vec,OU_par=1,use_LL=0,
                       is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=phylocov_fixed,
                       is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=phenocov_fixed)[[1]]
        
        
        len_vec <- calc_OU_len(heights=hts,edge_mat=edge_mat,des_order=des_order,nedge=nedge,P=P,lambda=lambda,sigma=phylocov,anc=tree$edge[,1],des=tree$edge[,2],nvar=nvar,nspecies=nspecies)
        ll2 <- tp(L=X, R=R, Rmat = as.matrix(Rmat),mL=ncol(X), mR=1, pheno_error=pheno_error, edge_vec=edge_vec, 
                  edge_ind=edge_ind,ind_edge=ind_edge, parent_edges = parent_edges,pars=pars, nvar=nvar, 
                  phylocov_diag=as.integer(!phylo_correlated), nind=nind, nob=nob, nspecies=nspecies, nedge=nedge, anc=anc, des=des, REML=as.integer(REML), 
                  species_subset=species_subset, un_species_subset = un_species_subset,subset_list=subset_list,
                  ind_list=ind_list, tip_combn=tip_combn,is_edge_ind=is_edge_ind,fixed_mu=matrix(0),ret_level=ret_level,OU_len=len_vec,OU_par=1,use_LL=0,
                  is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=phylocov_fixed,
                  is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=phenocov_fixed)
        if(is.na(ll2)[[1]])
        {
          class(alpha) <- "try-error"
          return(alpha)
        }
        
        if(!is.na(temp_ll2)) if(abs(temp_ll2-ll2[[1]])>100)
        {
          class(alpha) <- "try-error"
          return(alpha)
        }
        
        temp_alpha <- alpha + diag(nvar)*1e-6
        temp_ea <- try(eigen(temp_alpha),silent=TRUE)
        if(class(temp_ea)=="try-error") return(temp_ea)
        temp_P <- temp_ea$vectors
        temp_lambda <- temp_ea$values
        if(any(temp_lambda<=1e-8))
        {
          class(temp_alpha) <- "try-error"
          return(temp_alpha)
        }
        
        temp_len_vec <- calc_OU_len(heights=hts,edge_mat=edge_mat,des_order=des_order,nedge=nedge,P=temp_P,lambda=temp_lambda,sigma=phylocov,anc=tree$edge[,1],des=tree$edge[,2],nvar=nvar,nspecies=nspecies)
        temp_ll2 <- tp(L=X, R=R, Rmat = as.matrix(Rmat),mL=ncol(X), mR=1, pheno_error=pheno_error, edge_vec=edge_vec, 
                       edge_ind=edge_ind,ind_edge=ind_edge, parent_edges = parent_edges,pars=pars, nvar=nvar, 
                       phylocov_diag=as.integer(!phylo_correlated), nind=nind, nob=nob, nspecies=nspecies, nedge=nedge, anc=anc, des=des, REML=as.integer(REML), 
                       species_subset=species_subset, un_species_subset = un_species_subset,subset_list=subset_list,
                       ind_list=ind_list, tip_combn=tip_combn,is_edge_ind=is_edge_ind,fixed_mu=matrix(0),ret_level=1,OU_len=temp_len_vec,OU_par=1,use_LL=0,
                       is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=phylocov_fixed,
                       is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=phenocov_fixed)[[1]]
        
        if(!is.na(temp_ll2)) if(abs(temp_ll2-ll2[[1]])>100)
        {
          class(alpha) <- "try-error"
          return(alpha)
        }
        
        
        if(optim_verbose) cat(c(ll2[[1]],"\n"))
        if(ret_level>1) return(list(ll=ll2$logl,phylocov=phylocov,phenocov=phenocov,alpha=alpha,mu=ll2[[2]],pars=pars,ll2=ll2))
        list(ll=ll2[[1]])
      }
      
      ret <- phylopars(trait_data=trait_data,tree=tree,pheno_error=pheno_error>0,pheno_correlated=pheno_error==2,phylo_correlated=phylo_correlated,
                       REML=REML,EM_Fels_limit=EM_Fels_limit,repeat_optim_limit=repeat_optim_limit,EM_missing_limit=EM_missing_limit,repeat_optim_tol=repeat_optim_tol,phylocov_start=phylocov_start,
                       phenocov_start=phenocov_start,phylocov_fixed=phylocov_fixed,phenocov_fixed=phenocov_fixed,skip_optim=skip_optim,skip_EM=skip_EM,EM_verbose=EM_verbose,optim_verbose=optim_verbose)
      phylocov <- ret$pars[[1]]
      if(pheno_error!=0) phenocov <- ret$pars[[2]]
      mu <- ret$mu
      
      if(!is.na(phylocov_start)[[1]]) phylocov <- phylocov_start
      if(!is.na(phenocov_start)[[1]]) phenocov <- phenocov_start
      
      pars <- numeric()
      if(is.na(phylocov_fixed)[[1]]) pars <- mat_to_pars(phylocov,nvar,as.integer(!phylo_correlated))
      if(is.na(phenocov_fixed)[[1]]) if(pheno_error!=0) pars <- c(pars,mat_to_pars(phenocov,nvar,pheno_error))
      #pars[alpha_pars] <- mat_to_pars(diag(runif(nvar)),nvar,diag = abs(full_alpha-1))
      
      if(is.na(model_par_start)[[1]] & is.na(model_par_fixed)[[1]])
      {
        pars[alpha_pars] <- mat_to_pars(diag(nvar),nvar,diag = abs(full_alpha-1))
      } else if(!is.na(model_par_start)[[1]])
      {
        pars[alpha_pars] <- mat_to_pars(model_par_start,nvar,diag = abs(full_alpha-1))
      }
      
      zRmat <- (trait_data[,1:nvar+1,drop=FALSE] - matrix(1,nind) %*% t(zstand$mu)) / sqrt(matrix(1,nind) %*% diag(zstand$phylocov))
      zR <- matrix(na.exclude(as.double(t(zRmat[,1:nvar]))),ncol=1)
      pars2 <- numeric()
      if(is.na(phylocov_fixed)[[1]]) pars2 <- mat_to_pars(phylocov / (sqrt(diag(zstand$phylocov)) %*% t(sqrt(diag(zstand$phylocov)))),nvar,as.integer(!phylo_correlated))
      #pars2 <- mat_to_pars(phylocov / (sqrt(diag(zstand$phylocov)) %*% t(sqrt(diag(zstand$phylocov)))),nvar,as.integer(!phylo_correlated))
      if(is.na(phenocov_fixed)[[1]]) if(pheno_error!=0) pars2 <- c(pars2,mat_to_pars(phenocov / (sqrt(diag(zstand$phylocov)) %*% t(sqrt(diag(zstand$phylocov)))),nvar,pheno_error))
      if(is.na(model_par_fixed)[[1]]) pars2[alpha_pars] <- pars[alpha_pars]
      zll2 <- tp(L=X, R=zR, Rmat = as.matrix(zRmat),mL=ncol(X), mR=1, pheno_error=pheno_error, edge_vec=edge_vec, 
                 edge_ind=edge_ind,ind_edge=ind_edge, parent_edges = parent_edges,pars=pars2, nvar=nvar, 
                 phylocov_diag=as.integer(!phylo_correlated), nind=nind, nob=nob, nspecies=nspecies, nedge=nedge, anc=anc, des=des, REML=as.integer(REML), 
                 species_subset=species_subset, un_species_subset = un_species_subset,subset_list=subset_list,
                 ind_list=ind_list, tip_combn=tip_combn,is_edge_ind=is_edge_ind,fixed_mu=matrix(0),ret_level=3,OU_par=0,use_LL=0,
                 is_phylocov_fixed=as.integer(!is.na(phylocov_fixed)[[1]]),phylocov_fixed=zphylocov_fixed,
                 is_phenocov_fixed=as.integer(!is.na(phenocov_fixed)[[1]]),phenocov_fixed=zphenocov_fixed,OU_len=list())
      zBM_ll <- zll2[[1]]
      
      o <- optim(pars2,fn = function(X,BM_ll,R,Rmat,phylocov_fixed,phenocov_fixed)
      {
        ll <- try(OU_fun(X,R=R,Rmat=Rmat,phylocov_fixed=phylocov_fixed,phenocov_fixed=phenocov_fixed,BM_ll=BM_ll),silent=TRUE)
        if(class(ll)=="try-error") ll <- BM_ll-max_delta
        ll <- ll[[1]]
        if(is.na(ll)) ll <- BM_ll - max_delta
        if(abs(abs(BM_ll)-abs(ll[[1]]))>max_delta) ll <- BM_ll - max_delta
        ll[[1]]
      },control=list(fnscale=-1),method="BFGS",BM_ll=zBM_ll,R=zR,Rmat=zRmat,phylocov_fixed=zphylocov_fixed,phenocov_fixed=zphenocov_fixed)
      if(is.na(phylocov_fixed)[[1]]) phylocov2 <- pars_to_mat(o$par[phylocov_pars],nvar,as.integer(!phylo_correlated)) else phylocov2 <- zphylocov_fixed
      if(is.na(phenocov_fixed)[[1]]) if(pheno_error>0) phenocov2 <- pars_to_mat(o$par[phenocov_pars],nvar,pheno_error) else phenocov2 <- zphenocov_fixed
      if(is.na(model_par_fixed)[[1]]) alpha <- pars_to_mat(o$par[alpha_pars],nvar,diag = abs(full_alpha-1))
      if(npd) phylocov2 <- as.matrix(nearPD(phylocov2)$mat)
      if(npd) if(pheno_error>0) phenocov2 <- as.matrix(nearPD(phenocov2)$mat)
      if(npd) alpha <- as.matrix(nearPD(alpha)$mat)
      
      pars2 <- numeric()
      if(is.na(phylocov_fixed)[[1]]) pars2 <- mat_to_pars(phylocov2 * (sqrt(diag(zstand$phylocov)) %*% t(sqrt(diag(zstand$phylocov)))),nvar,as.integer(!phylo_correlated))
      if(is.na(phenocov_fixed)[[1]]) if(pheno_error!=0) pars2 <- c(pars2,mat_to_pars(phenocov2 * (sqrt(diag(zstand$phylocov)) %*% t(sqrt(diag(zstand$phylocov)))),nvar,pheno_error))
      if(is.na(model_par_fixed)[[1]]) pars2[alpha_pars] <- o$par[alpha_pars]
      
      if(is.na(phylocov_fixed)[[1]]) phylocov <- phylocov2 * (sqrt(diag(zstand$phylocov)) %*% t(sqrt(diag(zstand$phylocov))))
      if(is.na(phenocov_fixed)[[1]]) if(pheno_error!=0) phenocov <- phenocov2 * (sqrt(diag(zstand$phylocov)) %*% t(sqrt(diag(zstand$phylocov))))
      
      pars <- pars2
      ll2 <- OU_fun(pars,R=R,Rmat=Rmat,phylocov_fixed=phylocov_fixed,phenocov_fixed=phenocov_fixed,ret_level=3)
      o <- optim(pars,fn = function(X,BM_ll,R,Rmat,phylocov_fixed,phenocov_fixed) 
      {
        ll <- try(OU_fun(X,R=R,Rmat=Rmat,phylocov_fixed=phylocov_fixed,phenocov_fixed=phenocov_fixed,BM_ll=BM_ll),silent=TRUE)
        if(class(ll)=="try-error") ll <- BM_ll-max_delta
        ll <- ll[[1]]
        if(is.na(ll)) ll <- BM_ll - max_delta
        if(abs(abs(BM_ll)-abs(ll[[1]]))>max_delta) ll <- BM_ll - max_delta
        ll[[1]]
      },control=list(fnscale=-1),method="BFGS",BM_ll=ll2[[1]],R=R,Rmat=Rmat,phylocov_fixed=phylocov_fixed,phenocov_fixed=phenocov_fixed)
      if(is.na(phylocov_fixed)[[1]]) phylocov <- pars_to_mat(o$par[phylocov_pars],nvar,as.integer(!phylo_correlated))
      if(is.na(phenocov_fixed)[[1]]) if(pheno_error>0) phenocov <- pars_to_mat(o$par[phenocov_pars],nvar,pheno_error)
      if(is.na(model_par_fixed)[[1]]) alpha <- pars_to_mat(o$par[alpha_pars],nvar,diag = abs(full_alpha-1))
      ll2 <- OU_fun(o$par,R=R,Rmat=Rmat,phylocov_fixed=phylocov_fixed,phenocov_fixed=phenocov_fixed,ret_level=3)
      pars <- o$par
      if(usezscores)
      {
        f_args$usezscores <- FALSE
        raw_ll2 <- do.call(phylopars,f_args)
        if(raw_ll2[[1]]>ll2[[1]]) return(raw_ll2)
      }
      #ll2 <- list(ll=ll2$logl,phylocov=phylocov,phenocov=phenocov,mu=ll2[[2]],pars=pars,ll2=ll2,alpha=alpha)
    }
  npars <- 0
  dimnames(ll2$phylocov) <- dimnames(ll2$phenocov) <- list(colnames(trait_data)[1:nvar+1],colnames(trait_data)[1:nvar+1])
  if(phylo_correlated) npars <- ((nvar^2-nvar)/2+nvar) else npars <- nvar
  if(pheno_error==1) npars <- npars + nvar else if (pheno_error==2) npars <- npars + ((nvar^2-nvar)/2+nvar)
  if(pheno_error==0) pars <- list(phylocov=ll2$phylocov) else pars <- list(phylocov=ll2$phylocov,phenocov=ll2$phenocov)
  
  if(model=="mvOU")
  {
    dimnames(alpha) <- list(colnames(trait_data)[1:nvar+1],colnames(trait_data)[1:nvar+1])
    model <- list(model=model,alpha=alpha)
    npars <- npars + if(full_alpha) ((nvar^2-nvar)/2+nvar) else nvar
  } else if(model!="BM" & model!="white" & model!="star")
  {
    model <- list(model=model,mod_par=ll2[[model]])
    names(model)[2] <- if(model[[1]]=="OU") "alpha" else if(model[[1]]=="EB") "rate" else model[[1]]
    npars <- npars + 1
  } else
  {
    model <- list(model=model)
  }
  logLik <- ll2[[1]]
  anc_recon <- ll2$ll2$anc_recon
  anc_var <- matrix(0,nrow(anc_recon),ncol(anc_recon))
  anc_cov <- vector("list",nrow(anc_recon))
  colnames(anc_recon) <- colnames(anc_var) <- colnames(trait_data)[1:nvar+1]
  
  for(i in 1:nrow(anc_recon))
  {
    anc_cov[[i]] <- ll2$ll2$recon_var[1:nvar+(i-1)*nvar,,drop=FALSE]
    dimnames(anc_cov[[i]]) <- list(colnames(trait_data)[1:nvar+1],colnames(trait_data)[1:nvar+1])
    anc_var[i,] <- as.matrix(diag(anc_cov[[i]]))
  }
  
  rownames(anc_var) <- rownames(anc_recon) <- names(anc_cov) <- 1:length(anc_cov)
  rownames(anc_var)[1:nspecies] <- rownames(anc_recon)[1:nspecies] <- names(anc_cov)[1:nspecies] <- tree$tip.label
  rownames(anc_var)[(nspecies+1):nrow(anc_var)] <- rownames(anc_recon)[(nspecies+1):nrow(anc_var)] <- names(anc_cov)[(nspecies+1):nrow(anc_var)] <- (nspecies+1):nrow(anc_var)
  
  
  ret_list <- list(logLik=logLik,pars=pars,model=model,mu=anc_recon[nspecies+1,],npars=npars,anc_recon=anc_recon,anc_var=anc_var,anc_cov=anc_cov,tree=tree,trait_data=trait_data,REML=REML)
  if(model[[1]]=="mvOU" | model[[1]]=="OU")
  {
    if(model[[1]]=="OU") stationary_cov <- pars[[1]]/(2*model[[2]])
    else
    {
      ea <- eigen(model[[2]])
      P <- ea$vectors
      lambda <- ea$values
      lambda_mat <- matrix(0,nvar,nvar)
      for(i in 1:nvar)
        for(j in 1:i)
        {
          lambda_mat[i,j] <- lambda[i] + lambda[j]
        }
      lambda_mat[upper.tri(lambda_mat)] <- t(lambda_mat)[upper.tri(lambda_mat)]
      stationary_cov <- P %*% ((1/lambda_mat)*(solve(P)%*%pars[[1]]%*%t(solve(P)))) %*% t(P)
    }
    dimnames(stationary_cov) <-   list(colnames(trait_data)[1:nvar+1],colnames(trait_data)[1:nvar+1])
    ret_list$model <- c(ret_list$model,list(stationary_cov=stationary_cov))
  }
  
  
  class(ret_list) <- "phylopars"
  #g <- grad(f,pars2[[1]],em_ll=z_em_ll[[1]],R=R,Rmat=Rmat)
  #pars_to_mat(g[phylocov_pars],nvar,as.integer(!phylo_correlated))
  ret_list
}

prep_em <- function(trait_data,tree)
{
  nvar <- ncol(trait_data)-1
  p1 <- pic.ortho(x = setNames(lapply(tree$tip.label,function(X) trait_data[as.character(trait_data$species)==X,2]),tree$tip.label),phy = multi2di(tree,random=FALSE),var.contrasts = TRUE,intra = TRUE)
  pp1 <- cbind(c(unlist(attributes(p1)$intra),p1[,1]))
  vars <- c(unlist(attributes(p1)$intra)*0,p1[,2])
  if(ncol(trait_data)>2)
    for(i in 2:nvar)
    {
      p1 <- pic.ortho(x = setNames(lapply(tree$tip.label,function(X) trait_data[as.character(trait_data$species)==X,i+1]),tree$tip.label),phy = multi2di(tree,random=FALSE),intra = TRUE)
      pp1 <- cbind(pp1,c(unlist(attributes(p1)$intra),p1))
    }
  list(pics=pp1,vars=vars)
}

postorder_tools <- function(tree)
{
  # this function is linear in time with the number of species
  # this function returns a matrix called "edge_mat" and a vector called "des_order"
  ###
  ### edge_mat has 3 columns, constructed as follows:
  ####### the node number of the most recent
  ####### common ancestor of the species number in
  ####### column 1 and the species number in column 2
  ####### is found in column 3. The order is
  ####### in reverse postorder tree transversal
  ####### (from the root to the tip).
  #######
  ###
  ### des_order returns the vector of indexes that would put
  ### tree$edge[,2] in order, with the exception that a root edge
  ### is inserted at position (nspecies+1)
  
  tree <- reorder(tree,"postorder")
  nedge <- nrow(tree$edge)
  nspecies <- length(tree$tip.label)
  
  des <- tree$edge[,2] # descendant nodes for edge i
  anc <- tree$edge[,1] # ancestral nodes for edge i
  des_order <- order(des)
  
  # the next steps reorder the edge matrix so that
  # the descendant nodes are in order
  # A root edge is temporarily inserted at row (nspecies+1)
  new_edge <- tree$edge[des_order,1:2]
  new_edge <- rbind(new_edge[1:nspecies,],c(nspecies+1,nspecies+1),new_edge[(nspecies+1):nedge,])
  
  # 'edges' is the matrix to be returned
  ### first column is the low species number
  ### second column is the high species number
  ### third column is the descendant number
  edges <- cbind(new_edge*0,new_edge[,2])
  edges[,1] <- nspecies^2 # make sure the low species number initially contains a number >(nedge+1)
  edges[1:nspecies,1] <- edges[1:nspecies,2] <- 1:nspecies # set the mrca for species as themselves
  for(i in 1:nedge)
  {
    des_node <- des[i]
    anc_node <- anc[i]
    if(des_node<=nspecies)
    {
      edges[anc_node,1] <- min(edges[anc_node,1],des_node) # the lowest "low species" so far to represent anc[i]
      edges[anc_node,2] <- max(edges[anc_node,2],des_node) # the highest "high species" so far to represent anc[i]
    } else
    {
      # the lowest "low species" so far to represent anc[i]
      edges[anc_node,1] <- min(edges[des_node,1],edges[anc_node,1])
      # the highest "high species" so far to represent anc[i], but also make sure it's not >nspecies
      edges[anc_node,2] <- max(edges[des_node,2],ifelse(edges[anc_node,1]<=nspecies,edges[anc_node,1],0))
    }
  }
  
  new_des_order <- numeric(nedge+1)
  new_des_order[1:nspecies] <- des_order[1:nspecies]
  new_des_order[nspecies+1] <- nedge+1 # add in a root edge for des_order
  new_des_order[(nspecies+2):(nedge+1)] <- des_order[(nspecies+1):(nedge)]
  
  edges <- edges[tree$edge[,2],c(1,2,3)] # sort the edge_mat in postorder, and remove root edge
  list(edge_mat=edges,des_order=new_des_order)
}


getHeights <- function(tree) # identical to phylolm::pruningwise.distFromRoot
{
  tree <- reorder(tree,"postorder")
  nt = length(tree$tip.label)
  xx <- numeric(tree$Nnode + nt)
  for (i in length(tree$edge.length):1) xx[tree$edge[i, 2]] <- xx[tree$edge[i,1]] + tree$edge.length[i]
  names(xx) <- if (is.null(tree$node.label)) 
    1:(nt + tree$Nnode)
  else c(tree$tip.label, tree$node.label)
  return(xx)
}

convert_to_tree <- function(heights,nspecies,nedge,anc,des)
{
  new_len_vec <- numeric(nedge+1)
  edge_i <- nedge
  for(i in seq_along(des))
  {
    des_i <- des[edge_i]
    anc_i <- anc[edge_i]
    new_len_vec[edge_i] <- heights[des_i] - if(anc_i>(nspecies+1)) heights[anc_i] else 0
    edge_i <- edge_i - 1
  }
  new_len_vec[nedge+1] <- heights[nspecies+1] # root edge
  new_len_vec
}

calc_OU_heights <- function(heights,edge_mat,des_order,nedge,alpha,sigma)
{
  ntraits <- dim(alpha)[1]
  P <- eigen(alpha)$vectors
  tP <- t(P)
  invP <- solve(P)
  tinvP <- t(invP)
  len_mat <- matrix(0,nedge+1,ntraits^2)
  lambda <- eigen(alpha)$values
  lambda_mat <- matrix(0,ntraits,ntraits)
  for(i in 1:ntraits)
    for(j in 1:i)
    {
      lambda_mat[i,j] <- lambda[i] + lambda[j]
    }
  lambda_mat[upper.tri(lambda_mat)] <- t(lambda_mat)[upper.tri(lambda_mat)]
  for(i in 1:nedge)
  {
    temp_edge <- edge_mat[i,]
    t1 <- heights[temp_edge[1]]
    t2 <- heights[temp_edge[2]]
    t12 <- heights[temp_edge[3]]
    if(ntraits>1)
    {
      eAC1 <- P%*%diag(exp(-lambda*(t1-t12)))%*%invP
      teAC2 <- P%*%diag(exp(-lambda*(t2-t12)))%*%invP
    } else
    {
      eAC1 <- P*(exp(-lambda*(t1-t12)))*invP
      teAC2 <- P*(exp(-lambda*(t2-t12)))*invP
    }
    M <- P %*% (((1/lambda_mat)*(1-exp(-lambda_mat*t12)))*(invP%*%sigma%*%tinvP)) %*% tP
    temp <- eAC1 %*% M %*% teAC2
    len_mat[i,] <- temp[1:length(temp)]
  }
  
  #if(!is.ultrametric(tree))
  #{
  #  flag <- 1
  #  stop("Non-ultrametric trees not supported for OU model.")
  #externalEdge <- tree$edge[,2] <= nspecies
  #u <- c(mean(hts[1:nspecies]) - hts[1:nspecies],rep(0,nedge-nspecies+1))
  #tree$edge.length[externalEdge] <- tree$edge.length[externalEdge] + u[des[externalEdge]]
  #hts <- getHeights(tree)
  #HOU <- calc_OU_heights(heights = hts,edge_mat = edge_mat,des_order = des_order,nedge = nedge,alpha = alpha,sigma = sigma)
  #transf_heights <- 
  #  calc_OU_heights(heights = u,edge_mat = edge_mat,des_order = des_order,nedge = nedge,alpha = alpha,sigma = sigma)
  #HOU <- HOU - transf_heights
  #} else
  #{
  #  flag <- 0
  #  hts <- getHeights(tree)
  #  edge_mat <- postorder_tools(tree)
  #  des_order <- edge_mat[[2]]
  #  edge_mat <- edge_mat[[1]]
  #  HOU <- calc_OU_heights(heights = hts,edge_mat = edge_mat,des_order = des_order,nedge = nedge,alpha = alpha,sigma = sigma)
  #  len_vec <- apply(HOU,2,function(X) convert_to_tree(heights = X,nspecies = nspecies,nedge = nedge,anc = tree$edge[,1],des = tree$edge[,2]))
  #  
  #}
  #len_vec <- apply(HOU,2,function(X) convert_to_tree(heights = X,nspecies = nspecies,nedge = nedge,anc = tree$edge[,1],des = tree$edge[,2]))
  
  #if(flag==1)
  #{
  #ealpha <- eigen(alpha)
  #ea_vectors <- ealpha$vectors
  #inv_ea_vectors <- solve(ea_vectors)
  #ea_values <- ealpha$values
  
  #enalpha <- eigen(-alpha)
  #ena_vectors <- enalpha$vectors
  #inv_ena_vectors <- solve(ena_vectors)
  #ena_values <- enalpha$values
  #invL <- lapply(u[1:nspecies],function(X) ea_vectors %*% diag(exp(-ea_values*X)) %*% inv_ea_vectors)
  #L <- t(matrix(as.double(simplify2array(invL)),nrow=nvar))
  #R <- L %*% as.double(Y)
  #}
  
  len_mat[des_order,]
}

mat_to_pars <- function(M,nvar,diag,log_chol=TRUE,mod_chol=TRUE)
{
  len <- (nvar*nvar-nvar)/2+nvar
  if(diag==1)
  {
    ret <- if(log_chol) log(sqrt(diag(M))) else sqrt(diag(M))
    return(ret)
  } else
  {
    M <- try(chol(M),silent=TRUE)
    if(class(M)=="try-error") return(M)
    diag(M) <- if(log_chol) log(diag(M)) else diag(M)
  }
  if(diag!=1 & nvar>1)
  {
    if(log_chol & mod_chol) M[1,2] <- M[1,2]/exp(M[1,1])
  }
  
  vec_count = 1
  ret <- numeric(len)
  
  for(j in 1:nvar)
  {
    for(i in 1:j)
    {
      ret[vec_count] <- M[i,j]
      vec_count <- vec_count + 1
    }
  }
  ret
}

print.phylopars <- function(x, ...)
{
  PPE <- x
  
  cat("Phylogenetic trait variance-covariance\n")
  print(PPE$pars[[1]])
  cat("\n")
  if(length(x$pars)>1)
  {
    cat("Phenotypic trait variance-covariance\n")
    print(PPE$pars[[2]])
    
    nvar <- ncol(PPE$pars[[1]])
    xs <- as.list(PPE$trait_data[,1:nvar+1,drop=FALSE])
    finite <- lapply(xs,is.finite)
    xs <- lapply(xs,function(X) X[is.finite(X)])
    
    featurecount <- unlist(lapply(xs,length))
    featuresum <- unlist(lapply(xs,function(X) mean(sum(X),na.rm=TRUE)))
    featuremean <- featuresum/featurecount
    feature2sum <- unlist(lapply(xs,function(X) mean(sum(X^2),na.rm=TRUE)))
    featurevar <- feature2sum/featurecount - featuremean^2
    cat("\n% variance explained by phlogeny\n")
    percent_var <- (1-diag(PPE$pars[[2]])/featurevar)*100
    names(percent_var) <- colnames(PPE$pars[[1]])
    print(percent_var)
    
  }
  
  if(length(PPE$model)>1)
  {
    cat(paste("",PPE$model[[1]],names(PPE$model)[2]," = "))
    if(PPE$model[[1]]!="mvOU") cat((PPE$model[[2]])) else
    {
      cat("\n")
      print(PPE$model[[2]])
      cat("\n")
    }
    if(PPE$model[[1]]=="OU" | PPE$model[[1]]=="mvOU")
    {
      cat("\n")
      cat("Stationary covariance = ")
      cat("\n")
      print(PPE$model[[3]])
    }
  } else if(PPE$model=="BM") cat("Brownian motion model") else if(PPE$model=="white" | PPE$model=="star") cat("Star phylogeny model")
  cat("\n")
}

summary.phylopars <- function(object, ...)
{
  PPE <- object
  if(length(PPE$model)>1)
  {
    cat(paste("",PPE$model[[1]],names(PPE$model)[2]," = "))
    if(PPE$model[[1]]!="mvOU") cat((PPE$model[[2]])) else
    {
      cat("\n")
      print(PPE$model[[2]])
    }
  } else if(PPE$model=="BM") cat("Brownian motion model") else if(PPE$model=="white" | PPE$model=="star") cat("Star phylogeny model")
  cat("\n")
  
  ret <- matrix(0,nrow(PPE$pars[[1]]),4,dimnames = list(colnames(PPE$pars[[1]]),c("phylogenetic mean","phylogenetic sd","phenotypic sd",
                                                                                  "% variance explained by phylogeny")))
  ret[,1] <- PPE$mu
  ret[,2] <- round(sqrt(diag(PPE$pars[[1]])),4)
  ret[,3] <- round(sqrt(diag(PPE$pars[[2]])),4)
  nvar <- ncol(PPE$pars[[1]])
  xs <- as.list(PPE$trait_data[,1:nvar+1,drop=FALSE])
  finite <- lapply(xs,is.finite)
  xs <- lapply(xs,function(X) X[is.finite(X)])  
  featurecount <- unlist(lapply(xs,length))
  featuresum <- unlist(lapply(xs,function(X) mean(sum(X),na.rm=TRUE)))
  featuremean <- featuresum/featurecount
  feature2sum <- unlist(lapply(xs,function(X) mean(sum(X^2),na.rm=TRUE)))
  featurevar <- feature2sum/featurecount - featuremean^2
  percent_var <- (1-diag(PPE$pars[[2]])/featurevar)*100
  ret[,4] <- round(percent_var,2)
  ret
}

pval <- function(r,tree)
{
  N <- length(tree$tip.label)
  r <- abs(r)
  t <- r*sqrt((N-2)/(1-r^2))
  p <- 2*pt(t,N-2-1,lower.tail = FALSE)
  p[p>1] <- 1
  p
}

phylopars.lm <- function()
{
  args <- as.list(match.call())
  args <- args[3:length(args)]
  #trait_data$species <- factor(trait_data$species, levels=tree$tip.label)
  trait_data <- trait_data[,c(which(colnames(trait_data)=="species"),which(colnames(trait_data)!="species"))]
  original_data <- trait_data
  original_option <- getOption("na.action")
  options(na.action="na.pass")
  mod.mat <- model.matrix(object = formula,data = trait_data)
  intercept <- attr(terms(formula,data = trait_data),"intercept")==1
  if(!intercept) stop("Intercept-free PGLS not currently supported.")
  y_var <- model.frame(formula,data=trait_data)
  var_name <- colnames(y_var)[1]
  y_var <- y_var[,1,drop=FALSE]
  mod.mat <- cbind(mod.mat,y_var)
  colnames(mod.mat)[ncol(mod.mat)] <- var_name
  if(intercept)
  {
    trait_data <- data.frame(species=trait_data$species,mod.mat[,2:ncol(mod.mat)])
    colnames(trait_data) <- c("species",colnames(mod.mat)[2:ncol(mod.mat)])
  } else
  {
    trait_data <- data.frame(species=trait_data$species,mod.mat[,1:ncol(mod.mat)])
    colnames(trait_data) <- c("species",colnames(mod.mat)[1:ncol(mod.mat)])
  }
  options(na.action = original_option)
  args$trait_data <- trait_data
  PPE <- do.call(phylopars,args)
  n <- nspecies <- length(PPE$tree$tip.label)
  means <- PPE$anc_recon[nspecies+1,]
  trait_data <- PPE$trait_data
  df.int <- as.integer(intercept)
  k <- ncol(PPE$pars[[1]])
  rdf <- n - k
  if(PPE$model[[1]]!="mvOU" & PPE$model[[1]]!="OU") covX <- PPE$pars[[1]] else covX <- PPE$model[[3]]
  npred <- ncol(covX)-1
  y_pos <- ncol(covX)
  
  if(ncol(covX)==1 & intercept)
  {
    R2 <- 0
    ts <- ps <- SEs <- NA
    
  } else
  {
    coefs <- solve(covX[1:npred,1:npred,drop=FALSE])%*%covX[1:npred,y_pos,drop=FALSE]
    R2 <- as.double(sum(covX[1:npred,y_pos,drop=FALSE] * coefs) / covX[y_pos,y_pos,drop=FALSE])
  }
  R2adj <- 1-(1-R2)*(n-df.int)/(rdf)
  SST <- as.double(covX[y_pos,y_pos]) * (n-1)
  
  SSreg <- SST * R2
  SSres <- SST - SSreg
  MSres <- SSres / ((rdf))
  sigma <- sqrt(MSres)
  if(!(ncol(covX)==1 & intercept))
  {
    SEs <- sqrt(diag(solve((covX)[1:npred,1:npred,drop=FALSE]) * MSres / (n-1) ))
    ts <- coefs / SEs
    ps <- 2*(1-pt(abs(ts),rdf))
  }
  if(intercept==1)
  {
    if(ncol(covX)==1 & intercept) coefs <- setNames(means[y_pos],"(Intercept)") else
    {
      coefs <- as.double(c(Intercept=means[y_pos] - means[1:npred] %*% solve(covX[1:npred,1:npred,drop=FALSE])%*%covX[1:npred,y_pos,drop=FALSE],coefs))
      SEs <- c(NA,as.double(SEs))
      ts <- c(NA,as.double(ts))
      ps <- c(NA,as.double(ps))
      names(coefs) <- c("(Intercept)",colnames(covX)[1:npred])
    }
  } else names(coefs) <- colnames(covX)[1:npred]
  Fstat <- rdf / (k-df.int)*R2 / (1-R2)
  pval <- as.double(pf(Fstat,k-df.int,rdf,lower.tail = FALSE))
  logdet <- three.point.compute(tree,cbind(setNames(rep(1,n),tree$tip.label)))$logd
  ll <- -n/2 * log(2*pi) - n/2 * log((n-k) * MSres/n) - logdet/2 - n/2
  if(any(is.na(trait_data))) ll <- NA
  ret <- list(coefficients=coefs,SEs=SEs,ts=ts,ps=ps,R2=R2,R2adj=R2adj,sigma=sigma,Fstat=Fstat,pval=pval,df1=k,df2=rdf,dims=list(N=n,p=npred,REML=PPE$REML,df.int=df.int),model=formula,SST=SST,SSres=SSres,SSreg=SSreg,logLik=ll,PPE=PPE,original_data=original_data,covX=covX)
  class(ret) <- "phylopars.lm"
  ret
}
formals(phylopars.lm) <- c(alist(formula = ),formals(phylopars))

anova.phylopars.lm <- function(object,...)
{
  trait_data <- object$original_data
  covX <- object$covX
  SST <- object$SST
  tlabels <- attr(terms(object$model),"term.labels")
  k <- length(tlabels)
  n <- object$dims$N
  y_pos <- ncol(covX)
  NR <- length(tlabels) + 1
  rss <- resdf <- rep(NA, NR)
  rss[1] <- object$SST
  resdf[1] <- n - 1
  coefs <- coef(object)
  vars <- double()
  accum <- 0
  for (i in 1:(k)) {
    fmla <- as.formula(paste(colnames(object$PPE$pars[[1]])[ncol(object$PPE$pars[[1]])],"~",paste(tlabels[1:i],collapse="+")))
    rdf <-  length(colnames(model.matrix(fmla,data=trait_data))[-1]) - length(vars)
    accum <- accum + rdf
    vars <- colnames(model.matrix(fmla,data=trait_data))[-1]
    coefs <- solve(covX[1:accum,1:accum,drop=FALSE])%*%covX[1:accum,y_pos,drop=FALSE]
    R2 <- as.double(sum(covX[1:accum,y_pos,drop=FALSE] * coefs) / covX[y_pos,y_pos,drop=FALSE])
    SSreg <- SST * R2
    SSres <- SST - SSreg
    rss[i + 1] <- SSres
    resdf[i + 1] <- n - rdf
    resdf[i+1] <- resdf[i]-rdf
  }
  ss <- c(abs(diff(rss)), object$SSres)
  df <- c(abs(diff(resdf)), n - object$df1)
  ms <- ss/df
  fval <- ms/ms[NR]
  P <- pf(fval, df, df[NR], lower.tail = FALSE)
  table <- data.frame(df, ss, ms, f = fval, P)
  table[length(P), 4:5] <- NA
  dimnames(table) <- list(c(tlabels, "Residuals"), c("Df", 
                                                     "Sum Sq", "Mean Sq", "F value", "Pr(>F)"))
  structure(table, heading = c("Analysis of Variance Table\nSequential SS",
                               paste("Response:", deparse(object$model[[2L]]))), 
            class = c("anova", "data.frame"))
}

coef.phylopars.lm <- function(object,...)
{
  object$coefficients
}

print.phylopars.lm <- function(x,...)
{
  PPE <- x
  dd <- x$dims
  mCall <- x$call
  cat("Generalized least squares fit by ")
  cat(ifelse(dd$REML == 1, "REML\n", "maximum likelihood\n"))
  cat("  Model:", deparse(x$model), "\n\n")
  if(length(PPE$model)>1)
  {
    cat(paste("",PPE$model[[1]],names(PPE$model)[2]," = "))
    if(PPE$model[[1]]!="mvOU") cat((PPE$model[[2]])) else
    {
      cat("\n")
      print(PPE$model[[2]])
      cat("\n")
    }
    if(PPE$model[[1]]=="OU" | PPE$model[[1]]=="mvOU")
    {
      cat("\n")
      cat("Stationary covariance = ")
      cat("\n")
      print(PPE$model[[3]])
    }
  } else if(PPE$model=="BM") cat("Brownian motion model") else if(PPE$model=="white" | PPE$model=="star") cat("Star phylogeny model")
  cat("\n")
  cat("  Log-", ifelse(dd$REML == 1, "restricted-", 
                       ""), "likelihood: ", (x$logLik), "\n", sep = "")
  cat("\nCoefficients:\n")
  print(coef(x))
  cat("\n")
  
  cat("Degrees of freedom:", dd[["N"]] - dd[["p"]]-dd[["df.int"]], "\n")
  cat("Residual standard error:", (x$sigma), "\n")
  invisible(x)
}

logLik.phylopars <- function(object,...)
{
  val <- object$logLik
  attr(val,"nobs") <- attr(val,"nall") <- length(na.exclude(as.double(unlist(object$trait_data[,2:ncol(object$trait_data)]))))
  attr(val,"df") <- object$npars
  class(val) <- "logLik"
  val
}

logLik.phylopars.lm <- function(object,...)
{
  val <- object$logLik
  attr(val,"nall") <- object$dims$N
  attr(val,"nobs") <- object$dims$N
  attr(val,"df") <- object$df1 + 1
  class(val) <- "logLik"
  val
}

summary.phylopars.lm <- function(object,...)
{
  PPE <- object
  tTable <- data.frame(coef(object), object$SEs, object$ts, object$ps)
  dimnames(tTable) <- list(names(coef(object)), c("Value", "Std.Error", 
                                                  "t-value", "p-value"))
  dd <- object$dims
  mCall <- object$call
  cat("Generalized least squares fit by ")
  cat(ifelse(dd$REML == 1, "REML\n", "maximum likelihood\n"))
  cat("  Model:", deparse(object$model), "\n")
  
  if(length(PPE$model)>1)
  {
    cat(paste("",PPE$model[[1]],names(PPE$model)[2]," = "))
    if(PPE$model[[1]]!="mvOU") cat((PPE$model[[2]])) else
    {
      cat("\n")
      print(PPE$model[[2]])
      cat("\n")
    }
    if(PPE$model[[1]]=="OU" | PPE$model[[1]]=="mvOU")
    {
      cat("\n")
      cat("Stationary covariance = ")
      cat("\n")
      print(PPE$model[[3]])
    }
  } else if(PPE$model=="BM") cat("Brownian motion model") else if(PPE$model=="white" | PPE$model=="star") cat("Star phylogeny model")
  cat("\n")
  
  aux <- logLik(object)
  object$BIC <- BIC(aux)
  object$AIC <- AIC(aux)
  
  print(data.frame(AIC = object$AIC, BIC = object$BIC, logLik = as.vector(object$logLik), 
                   row.names = " "))
  cat("\nCoefficients:\n")
  colnames(tTable) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
  printCoefmat(tTable,P.values = TRUE,has.Pvalue = TRUE,na.print = "")
  
  cat("\n")
  cat("Residual standard error:", (object$sigma), "\n")
  cat("Degrees of freedom:", dd[["N"]] - 
        dd[["p"]] - dd[["df.int"]], "\n")
  cat("Multiple R-squared:", (object$R2), "Adjusted R-squared:", (object$R2adj), "\n")
  cat("F-statistic:", (object$Fstat), "on", (object$df1-1), "and", (object$df2), "DF, p-value:", (object$pval), "\n")
  invisible(object)
}


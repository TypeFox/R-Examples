#########################################################################################
log.likelihood <- function(cross,DG,markers,phenotypes,genotypes,addcov=NULL,intcov=NULL)
{
  aux1 <- 0
  aux2 <- 0
  n.phe <- length(phenotypes)
  for(i in 1:n.phe){
    covs <- get.cov.loglik(resp=phenotypes[i],DG=DG)
    aux <- log.likelihood.s(cross=cross, 
                            node=phenotypes[i], markers=markers[[phenotypes[i]]], 
                            cov.node=covs, genotypes=genotypes, 
                            addcov=addcov, intcov=intcov)
    aux1 <- aux1 + aux[[1]]
    aux2 <- aux2 + aux[[2]]
  }
  list(aux1,aux2)
}
###################################
get.cov.loglik <- function(resp,DG)
{
  nDG <- length(DG[,1])
  cov1 <- c()
  for(i in 1:nDG){	
    if((DG[i,1] == resp) & (DG[i,2] == "<----")){
      cov1 <- c(cov1,DG[i,3])
    }
    if((DG[i,3] == resp) & (DG[i,2] == "---->")){
      cov1 <- c(cov1,DG[i,1])
    }
  }
  return(cov1)
}
#################################################################################################
log.likelihood.s <- function(cross, node, markers, cov.node, genotypes, addcov=NULL, intcov=NULL)
{
  node.col <- find.pheno(cross, pheno=node)
  y <- data.frame(cross$pheno[,node.col])
  names(y) <- "y"
  nQ <- length(markers)
  if(nQ > 0){
    mygeno <- data.frame(genotypes[,markers])
    nQ <- length(mygeno[1,])
    for(i in 1:nQ){
      mygeno[,i] <- as.factor(mygeno[,i])
    }
    names(mygeno) <- paste("Q",1:nQ,sep="")
    if(!is.null(cov.node)){
      covar.node <- data.frame(cross$pheno[,c(intcov,addcov)],cross$pheno[,find.pheno(cross, pheno=cov.node)])
      names(covar.node) <- c(intcov,addcov,cov.node)
      aux <- s.formula(addcov=names(covar.node), intcov=intcov, nQ=nQ, dat=covar.node)
      mylm <- lm(aux[[1]],data.frame(y,aux[[2]],mygeno))
      mylist <- list(logLik(mylm)[[1]],AIC(mylm,k=log(length(y[,1]))))
    }
    else{
      covar.node <- data.frame(cross$pheno[,c(intcov,addcov)])
      names(covar.node) <- c(intcov,addcov)
      aux <- s.formula(addcov=names(covar.node), intcov=intcov, nQ=nQ, dat=covar.node)
      if(!is.null(aux[[2]])) mylm <- lm(aux[[1]],data.frame(y,aux[[2]],mygeno))
      if(is.null(aux[[2]])) mylm <- lm(aux[[1]],data.frame(y,mygeno))
      mylist <- list(logLik(mylm)[[1]],AIC(mylm,k=log(length(y[,1]))))
    }
  }
  else{
    if(!is.null(cov.node)){
      covar.node <- data.frame(cross$pheno[,c(intcov,addcov)],cross$pheno[,find.pheno(cross, pheno=cov.node)])
      names(covar.node) <- c(intcov,addcov,cov.node)
      aux <- s.formula(addcov=names(covar.node), intcov=intcov, nQ=0, dat=covar.node)
      mylm <- lm(aux[[1]],data.frame(y,aux[[2]]))
      mylist <- list(logLik(mylm)[[1]],AIC(mylm,k=log(length(y[,1]))))
    }
    else{
      covar.node <- data.frame(cross$pheno[,c(intcov,addcov)])
      names(covar.node) <- c(intcov,addcov)
      aux <- s.formula(addcov=names(covar.node), intcov=intcov, nQ=0, dat=covar.node)
      if(!is.null(aux[[2]])) mylm <- lm(aux[[1]],data.frame(y,aux[[2]]))
      if(is.null(aux[[2]])) mylm <- lm(aux[[1]],y)
      mylist <- list(logLik(mylm)[[1]],AIC(mylm,k=log(length(y[,1]))))
    }
  }
  mylist
}
#################################################################################################
s.formula <- function(addcov=NULL, intcov=NULL, nQ, dat)
{
  if(!is.null(addcov) & is.null(intcov)){ 
    mycovs <- dat[,addcov]
    form <- as.formula(paste(" ~ ", paste(addcov,collapse="+")))
    mycovs <- data.frame(model.matrix(form,dat)[,-1])
    addcov <- names(mycovs)
    if(nQ > 0){
      Qnames <- paste("Q", 1:nQ, sep = "")
      myform <- as.formula(paste("y ~ ", paste(c(addcov,Qnames), collapse = "+")))
    }
    else{
      myform <- as.formula(paste("y ~ ", paste(addcov, collapse = "+")))
    }
  }
  if(!is.null(intcov)){
    le <- length(intcov)	
    intaddcov <- unique(c(intcov,addcov))
    mycovs <- dat[,c(intaddcov)]
    form <- as.formula(paste(" ~ ", paste(intaddcov,collapse="+")))
    mycovs <- data.frame(model.matrix(form,dat)[,-1])
    intaddcov <- names(mycovs)
    if(nQ > 0){
      Qnames <- paste("Q", 1:nQ, sep = "")
      intQnames <- c()
      for(i in 1:le){
        intQnames <- c(intQnames,paste(intaddcov[i], Qnames, sep=":"))
      }
      myform <- as.formula(paste("y ~ ", paste(c(intaddcov,Qnames,intQnames), collapse = "+")))
    }
    else{
      myform <- as.formula(paste("y ~ ", paste(intaddcov, collapse = "+")))
    }
  }
  if(is.null(addcov) & is.null(intcov)){
    if(nQ > 0){
      Qnames <- paste("Q", 1:nQ, sep = "")
      myform <- as.formula(paste("y ~ ", paste(Qnames, collapse = "+")))
      mycovs <- NULL
    }
    else{
      myform <- as.formula("y ~ 1")
      mycovs <- NULL
    }
  }	
  list(myform,mycovs)
}

###############################################################################################################
get.all.solutions <- function(DG, rc, n.shuffles, cross, QTLs, markers, phenotypes, genotypes,
                              addcov = NULL, intcov = NULL)
{
  mylist <- list()
  rc <- order.as(rc)
  mylist[[1]] <- rc
  myloglik <- c()
  myBIC <- c()
  n.arrows <- length(DG[,1])
  aux <- log.likelihood(cross = cross, DG = rc, markers = markers, phenotypes = phenotypes,
                        genotypes = genotypes, addcov = addcov, intcov = intcov) 
  myloglik[1] <- aux[[1]]
  myBIC[1] <- aux[[2]]
  for(i in 2:n.shuffles){
    DG <- shuffle.DG(DG)
    rc <- recheck.directions(cross=cross,QTLs=QTLs,oldDG=DG,addcov=addcov,intcov=intcov)
    if(attr(rc, "message") == "algorithm converged"){
      rc <- order.as(rc)
      mylist[[i]] <- rc
      aux <- log.likelihood(cross = cross, DG = rc, markers = markers, phenotypes = phenotypes,
                            genotypes = genotypes, addcov = addcov, intcov = intcov)	
      myloglik[i] <- aux[[1]]
      myBIC[i] <- aux[[2]]
    }
  }
  ed <- unique(myloglik)
  le <- length(ed)
  newlist <- list()
  loglikelihood <- bic <- rep(0,le)
  for(i in 1:le){
    aux.pos <- which(myloglik==ed[i])
    pos <- aux.pos[1]
    newlist[[i]] <- mylist[[pos]]
    loglikelihood[i] <- myloglik[pos]
    bic[i] <- myBIC[pos]
  }
  outlist <- list(newlist,loglikelihood,bic)
  names(outlist) <- c("solutions","loglikelihood","BIC")
  outlist
}
########################
order.as <- function(as)
{
  aux1 <- row.names(as)
  n <- length(aux1)
  ordered <- data.frame(matrix(0,n,4))
  for(i in 1:n){
    aux2 <- which(aux1==i)
    ordered[i,] <- as[aux2,]
  }
  names(ordered) <- c("node1","direction","node2","lod")
  ordered
}

##########################################################################
recheck.directions <- function(cross,QTLs,oldDG,addcov=NULL,intcov=NULL)
{
  DG1 <- DG2 <- oldDG
  aux1 <- "not equal"
  cont <- 0
  le <- length(DG1[,1])
  while((aux1 == "not equal") & (cont < 30)){
    for(i in 1:le){
      DG2 <- orient.graph.edges.cov(cross=cross, QTLs=QTLs, oldDG=DG2, i=i, addcov=addcov, intcov=intcov)
    }
    aux1 <- check.DG(newDG=DG2,oldDG=DG1)
    DG1 <- DG2 
    cont <- cont + 1
  }
  if(aux1 == "equal"){
    newDG <- DG1
    message <- "algorithm converged"
  }
  if(aux1 != "equal"){
    newDG <- oldDG
    message <- "algorithm didn't converge"
  }
  attr(newDG, "edgemode") <- "directed"
  attr(newDG, "message") <- message
  attr(newDG, "cont") <- cont
  newDG
}
#################################
check.DG <- function(newDG,oldDG)
{
  aux <- all.equal(newDG[,2],oldDG[,2])
  return( ifelse(aux == TRUE, "equal", "not equal") )
}
###################################
get.covariates <- function(pair,DG)
{
  nDG <- length(DG[,1])
  node1 <- DG[pair,1]
  node2 <- DG[pair,3]
  cov1 <- c()
  cov2 <- c()
  for(i in 1:nDG) {
    if((i != pair) & (DG[i,1] == node1) & (DG[i,2] == "<----"))
      cov1 <- c(cov1,DG[i,3])
    if((i != pair) & (DG[i,3] == node1) & (DG[i,2] == "---->"))
      cov1 <- c(cov1,DG[i,1])
    if((i != pair) & (DG[i,1] == node2) & (DG[i,2] == "<----"))
      cov2 <- c(cov2,DG[i,3])
    if((i != pair) & (DG[i,3] == node2) & (DG[i,2] == "---->"))
      cov2 <- c(cov2,DG[i,1])
  }
  return(list(cov1,cov2))
}
##########################
shuffle.DG <- function(DG)
{
  le <- length(DG[,1])
  aux <- sample(c(1:le),le,replace=FALSE)
  return(DG[aux,]) 
}
###################################
pull.geno.argmax <- function(cross)
{
  le <- length(cross$geno)
  all <- cross$geno[[1]]$argmax
  for(i in 2:le){
    aux <- cross$geno[[i]]$argmax
    all <- cbind(all,aux)
  }
  all
}
##############################################################################
orient.graph.edges <- function(cross,UDG,QTLs,addcov=NULL,intcov=NULL)
{
  UDG <- subset(UDG,UDG[,3]==1)
  le <- length(UDG[,1])
  DG <- data.frame(matrix(0,le,4))
  for(i in 1:le){
    node1 <- DG[i,1] <- UDG[i,1]
    node2 <- DG[i,3] <- UDG[i,2]
    s <- lod.score(cross=cross, node1=node1, node2=node2, 
                   qtl.node1=QTLs[[node1]], qtl.node2=QTLs[[node2]], 
                   cov.node1=addcov, cov.node2=addcov, intcov=intcov, artfact.qtl=QTLs[[1]])
    DG[i,2] <- ifelse(s >= 0, "---->", "<----")
    DG[i,4] <- s
  }
  names(DG) <- c("node1", "direction", "node2", "lod score")
  class(DG) <- c("qdg", "data.frame")
  attr(DG, "edgemode") <- "directed"
  DG
}
##############################################################################
orient.graph.edges.cov <- function(cross,QTLs,oldDG,i,addcov=NULL,intcov=NULL)
{
  newDG <- oldDG
  cov <- get.covariates(pair=i,DG=oldDG)
  node1 <- oldDG[i,1]
  node2 <- oldDG[i,3]
  ls <- lod.score(cross=cross, node1=node1, node2=node2, 
                  qtl.node1=QTLs[[node1]], qtl.node2=QTLs[[node2]], 
                  cov.node1=c(cov[[1]],addcov), cov.node2=c(cov[[2]],addcov), 
                  intcov=intcov, artfact.qtl=QTLs[[1]])
  if(ls >= 0){newDG[i,2] <- "---->"}
  else{newDG[i,2] <- "<----"}
  newDG[i,4] <- ls
  names(newDG) <- c("node1", "direction", "node2", "lod score")
  return(newDG)
}

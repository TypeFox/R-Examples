selectRho <- function(simMat,rhovec=NULL) {

  require(WGCNA)
  
  if (any(is.na(simMat))) stop("Similarity Matrix cannot contain missing values. \n")
  if (any(simMat < 0)) stop("Similarity Matrix cannot contain negative values. \n")
  if (any(diag(simMat)>0)) {
    cat("GOGANPA ignore self-links, diagonal entries of simMat set to 0. \n")
    diag(simMat) <- 0
  }
  if (is.null(rhovec)) {
    maxSim <- max(simMat)
    rhovec <- seq(maxSim*0.1,maxSim*0.9,length=10)
  }
  
  crit <- c()
  for (i in 1:length(rhovec)) {
    cat('checking rho =',rhovec[i])
    adjMat <- simMat>=rhovec[i]
    conn <- rowSums(adjMat)
    crit <- rbind(crit,scaleFreeFitIndex(conn))
    cat(' Rsq =',crit[i,1], 'Slope =',crit[i,2], '\n')
  }

  crit <- cbind(rhovec,crit)
  colnames(crit) <- c('rho','Rsq','slope','truncatedExpRsq')
  goodCrit <- crit[crit$slope<0,]
  bestrho <- goodCrit$rho[which.max(goodCrit$Rsq)]
  list(criterion=crit,bestrho=bestrho)

}

getGNET <- function(simMat,rho) {

  if (any(is.na(simMat))) stop("Similarity Matrix cannot contain missing values. \n")
  if (any(simMat < 0)) stop("Similarity Matrix cannot contain negative values. \n")
  if (any(diag(simMat)>0)) {
    cat("GOGANPA ignore self-links, diagonal entries of simMat set to 0. \n")
    diag(simMat) <- 0
  }
    
  adjMat <- simMat>=rho
  conn <- rowSums(adjMat)
  if (any(conn==0)) adjMat <- adjMat[-which(conn==0),-which(conn==0)]
  gNET <- vector('list',nrow(adjMat))
  geneNames <- colnames(adjMat)
  for (i in 1:nrow(adjMat)) gNET[[i]] <- geneNames[adjMat[i,]==1]
  names(gNET) <- geneNames

  gNET
  
}

GOGANPA <- function(gExprs.obj, gsets, gNET=NULL, simMat=NULL, rho=NULL,
                    msp.groups, check.exprs=TRUE, msp.correction=TRUE,
                    size.min=15, size.max=500, permN=2000, randN=30,
                    permFDR.cutoff=0.15,output.label="GOGANPAResult") {

  require(GANPA)
  
  if (!is.null(gNET) & !is.null(simMat)) stop('exactly one of gNET and simMat must be NULL. \n')
  if (is.null(gNET)) {
    if (is.null(rho)) {
      cat('finding best rho based on scale-free-criterion. \n')
      rho <- selectRho(simMat)$bestrho
      cat('rho selected as',rho,'. \n')
    }
    gNET <- getGNET(simMat,rho)
  }
  
  GSE.Test.Main(gExprs.obj=gExprs.obj, gsets=gsets, gNET=gNET,
                  msp.groups=msp.groups, check.exprs=check.exprs, size.min=size.min, size.max=size.max,
                  permN=permN, randN=randN, permFDR.cutoff=permFDR.cutoff,
                  output.label=output.label, msp.correction = msp.correction)

}

##Function: expression matrix check
.checkExpMat <- function( expMat1, sampleVec1, expMat2,  sampleVec2 ){
  
  #check input data
  if( is.null( rownames(expMat1)) | is.null(rownames(expMat2)) )
    stop("Error: no rownames for expMat1 or expMat2")
  if( length( which( (rownames(expMat1) == rownames(expMat2)) == FALSE ) ) > 0 )
    stop("Error: confict rownames between expMat1 and expMat2")
  if( length(sampleVec1) != ncol(expMat1) )
    stop("Error: conflict between expMat1 and sampleVec")
  if( length(sampleVec2) != ncol(expMat2) )
    stop("Error: conflict between expMat2 and sampleVec")
  
  if( length( which( (rownames(expMat1) == rownames(expMat2)) == FALSE )) > 0 ) 
    stop("Error:conflicted gene order in expMat1 and expMat2")

  sample.un1 <- unique(sampleVec1)
  sample.un2 <- unique(sampleVec2)
  if( length(sample.un1) != length(sample.un2) )
    stop("Error: different time_points/tissues/unique samples between expMat1 and expMat2.\n")
  
  if( length( which( (rownames(sample.un1) == rownames(sample.un2)) == FALSE ) ) > 0 )
    stop("Error: conflict time_points/tissues/unique samples between expMat1 and expMat2.\n")

}





.expMean <- function( expMat, sampleVec, samplelabel = NULL, logTransformed = TRUE, base = 2 ) {
  
   sample.un <- unique(sampleVec)
   if( length(sample.un) == length(sampleVec) ) {
      return(expMat)
   }

   #get expression data before log-transformed
   expMat_org <- expMat
   if( logTransformed == TRUE ) {
      expMat_org <- .expFun( expMat, base )
    }

    #generate matrix to record final result
    res <- matrix(0, nrow = nrow(expMat), ncol = length(sample.un) )
    rownames(res) <- rownames(expMat)
    if( is.null(samplelabel) ) {
      colnames(res) <- paste( "sample", sample.un, sep = "" )
    }else {
      colnames(res) <- samplelabel
    }   
    

    #average duplicates
    for( ii in 1:length(sample.un) ) {
       idx <- which( sampleVec == sample.un[ii] )
       if( length(idx) == 1 ) {
          res[,ii] <- expMat_org[,idx]
       }else {
          res[,ii] <- apply(expMat_org[,idx], 1, mean )
       }#end if-else
    }##end for ii

    if( logTransformed == TRUE ) 
      res <- log( res, base )

    return(res)

}






#Function: gene expression comparision for two different biologically condtions
##base: numeric or EXP
.difExp <- function( expMat1, sampleVec1, expMat2, sampleVec2, logTransformed = TRUE, base = 2 ){

  .checkExpMat( expMat1, sampleVec1, expMat2, sampleVec2)
  
  sample.un1 <- unique(sampleVec1)
  
  result <- matrix( 0, nrow = nrow(expMat1), ncol= length(sample.un1))
  rownames(result) <- rownames(expMat1)
  colnames(result) <- sample.un1

  ##exp mean
  expMat1.mean <- .expMean( expMat = expMat1, sampleVec = sampleVec1, samplelabel = NULL, logTransformed = logTransformed, base = base )
  expMat2.mean <- .expMean( expMat = expMat2, sampleVec = sampleVec2, samplelabel = NULL, logTransformed = logTransformed, base = base )

  
  expMat1.mean_org <- expMat1.mean
  expMat2.mean_org <- expMat2.mean
  if( logTransformed == TRUE ) {
     expMat1.mean_org <- .expFun( expMat1.mean, base )
     expMat2.mean_org <- .expFun( expMat2.mean, base )
  }

  ##exp z-score
  zScore1 <- .matZScore( expMat1.mean_org )
  zScore2 <- .matZScore( expMat2.mean_org )

  ##fold change  
  fcMat <- NULL
  if( logTransformed == TRUE ) {
     fcMat <- expMat2.mean - expMat1.mean
  }else {
     fcMat <- log( expMat2.mean_org/expMat1.mean_org, base)
  }
  colnames(fcMat) <- paste( "fc", sample.un1, sep = "" )
  
  
  return( list(expMat1.mean = expMat1.mean, expMat2.mean = expMat2.mean, zScore1 = zScore1, zScore2 = zScore2, fcMat = fcMat ) )
}



##Functions: generating express-based features 
expFeatureMatrix <- function( expMat1, sampleVec1, expMat2, sampleVec2, logTransformed = TRUE, base = 2, 
                              features = c("zscore", "foldchange", "cv", "expression")) {
    
  res <- .difExp( expMat1, sampleVec1, expMat2, sampleVec2, logTransformed, base )

  sampleUnique <- unique( sampleVec1 )
  
  featureMat <- matrix( 0, nrow = nrow(expMat1), ncol = 1)
  colnames(featureMat) <- "tmp"
  for( ii in 1:length(features) ){
    tmp <- NULL
    if( features[ii] == "zscore" ) {
      tmp <- c( colnames(featureMat), paste( "zScore1", sampleUnique, sep = "_"), paste( "zScore2", sampleUnique, sep = "_") )
      featureMat <- cbind(featureMat, res$zScore1, res$zScore2)
    }else if( features[ii] == "foldchange"  ) {
      tmp <- c( colnames(featureMat), paste( "foldchange", sampleUnique, sep = "_"))
      featureMat <- cbind( featureMat, res$fcMat )
    }else if( features[ii] == "cv" ){
      cv1 <- apply( res$expMat1.mean, 1, .cv )
      cv2 <- apply( res$expMat2.mean, 1, .cv )
      tmp <- c( colnames(featureMat), "cv1", "cv2" )
      featureMat <- cbind( featureMat, cv1, cv2 )
    }else if( features[ii] == "expression" ) {
      tmp <- c( colnames(featureMat), paste( "expression1", sampleUnique, sep = "_"), paste( "expression2", sampleUnique, sep = "_") )
      featureMat <- cbind( featureMat, res$expMat1.mean, res$expMat2.mean )  
    }else {
      stop("Error: undefined expression-based features")
    }
    colnames(featureMat) <- tmp
  }#end for ii
  featureMat <- featureMat[,-1]
  featureMat
}




##gene selection with different statistical methods
geneRanker <- function( expmat1, expmat2, genes, rankers = c("ttest", "SAM", "Limma"), verbose = FALSE ) {

 #if( !require(GeneSelector) ) {
 #   source("http://bioconductor.org/biocLite.R")
 #   biocLite("GeneSelector")
 #   require(GeneSelector)
 #}
  
  
  dataMatrix <- cbind( expmat1[genes,], expmat2[genes,] )
  classes <- c( rep(0, ncol(expmat1)), rep(1, ncol(expmat2)) )
  
  res <- list() 
  for( ranker in rankers ) {
    if( verbose )
      cat( ranker, '\n')
    if(ranker == "ttest") ranked <- RankingTstat(dataMatrix, classes, type="paired", gene.names = genes)
    if(ranker == "FC") ranked <- RankingFC(dataMatrix, classes, type="paired")
    if(ranker == "Wilcoxon") ranked <- RankingWilcoxon(dataMatrix, classes, type="paired")
    if(ranker == "BaldiLong") ranked <- RankingBaldiLong(dataMatrix, classes, type="paired")
    if(ranker == "ShrinkageT" ) ranked <- RankingShrinkageT(dataMatrix, classes, type="paired")
    if(ranker == "SAM") ranked <- RankingSam(dataMatrix, classes, type="paired")
    if(ranker == "SoftthresholdT") ranked <- RankingSoftthresholdT(dataMatrix, classes, type="paired")
    if(ranker == "Limma") ranked <- RankingLimma(dataMatrix, classes, type="paired")
    
    
    ranked.mat <- toplist( ranked, top = length(genes), show =FALSE )
    
    tt <- ranked.mat[,1]
    idx <- tt[!is.na(tt)]
    other.idx <- 1:length(genes)
    other.idx <- setdiff( other.idx, idx )
    rownames(ranked.mat) <- genes[c(idx, other.idx)]
    
    ranked.mat[other.idx, "statistic"] <- 0
    ranked.mat[other.idx, "pval"] <- 1
    
    ranked.mat <- ranked.mat[genes,]
    
    res[[ranker]] <- ranked.mat
  }
  
  
  #for pvalue matrix
  pvalMat <- matrix(0, nrow = length(genes), ncol = length(rankers) ) 
  rownames(pvalMat) <- genes
  colnames(pvalMat) <- rankers
  for( i in 1:length(rankers) ){
    pvalMat[,i] <- res[[rankers[i]]][genes,"pval"]
  }
  res[["pvalMat"]] <- pvalMat
  
  
  res
}





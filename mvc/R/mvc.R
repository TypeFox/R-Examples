source("R/mvc-utils.R")


#' Multi-View Clustering using Spherical k-Means for categorical data.
#' See: S. Bickel, T. Scheffer: Multi-View Clustering, ICDM 04.
#' Hierachical clustering used to determine the initial centers for k-Means
#'
#' @param view1 View number one, a data frame with the same row names as view2. All columns numeric. Entries are natural numbers, starting from 1.
#' @param view2 View number two, a data frame with the same row names as view1. All columns numeric. Entries are natural numbers, starting from 1.
#' @param k The maximum number of clusters to create
#' @param startView The view on which to perform the initial E step, one of "view1", "view2"
#' @param nthresh The number of iterations to run without improvement of the objective function
#' @param doOutput Whether output to the console should be done
#' @param doDebug Whether debug output to the console should be done (implies normal output)
#' @param plotFile File name where the hierarchical clustering plot should be stored
#' @return A list reporting the final clustering, with names finalIndices, agreementRate, indicesSv, indicesOv. They designate final cluster indices as a vector, as well as agreement rate of the two views, and the individual indices given by the two views over the course of iterations as data frames.
#' @examples {
#' # Demo program, showing how to run Multi-
#' # View Clustering using Spherical k-Means 
#' # AM, 2011
#'
#' # load toy data 'toyView1' and 'toyView2'
#' data(toyViews)
#' 
#' mvcsph(
#'    view1=toyView1,
#'    view2=toyView2,
#'    nthresh=20,
#'    k=4,
#'    startView="view1"
#' )
#' 
#' }


mvcsph <- function(
   view1=NULL,
   view2=NULL,
   k=Inf,
   startView="view1",
   nthresh=20,
   doOutput=F,
   doDebug=F,
   plotFile="Rplots.pdf"
  ) {

  if (doDebug) doOutput=T
  if (doOutput) cat(paste("\nMVC (Spherical k-Means)."))
  if (doOutput) cat(paste("\nk",k))
  if (doOutput) cat(paste("\nnthresh",nthresh))
  if (doOutput) cat("\n")

  if (startView == "view1") {
    startView = view1
    otherView = view2
  } 
  else {
    if (startView == "view2") {
      startView = view2
      otherView = view1
    }
    else {
      stop("startView argument wrong")
    }
  }
  checkViews(startView, otherView)
  uniqueVals <- viewsClasses(startView, otherView)
  if (doDebug) {
    cat("\nUnique values per view:\n")
    print(uniqueVals)
  }
  nClassesSV = length(uniqueVals$view1)
  nClassesOV = length(uniqueVals$view2)
  if (sum(uniqueVals$view1 != c(1:nClassesSV))) stop(paste("Start view classes must be in 1 ...",nClassesSV))
  if (sum(uniqueVals$view2 != c(1:nClassesOV))) stop(paste("Other view classes must be in 1 ...",nClassesOV))

    
  # # # Start MVCSph # # #

  # First thing: convert rows to unit length
  for (i in 1:NROW(startView)) {
    startView[i,]=UL(unlist(matrix(startView[i,]))) # ugly to do it like that (need vectors)
  }
  for (i in 1:NROW(otherView)) {
    otherView[i,]=UL(unlist(matrix(otherView[i,])))
  }

  # Get k centers from startView
  if (doOutput) cat("\n\nHierarchical Clustering (for centers to start with)\n")
  pdf(file=plotFile)
  heatmap(as.matrix(startView),revC=T,margins=c(10,10)) # with normal distance function (since we have already converted to UL)
  #Sys.sleep(2)
  dev.off()
  hclustering = hclust( dist( as.matrix(startView ) ) ) # with normal distance function (since we have already converted to UL)
  if (k>NROW(startView)) { 
    k = NROW(startView)
    cat(paste("\nNote: Truncated k to",k))
  }
  initialCenters<-centers.hclust( as.matrix(startView), hclustering, nclust=k, use.median=FALSE )


  # perform initial E-Step, i.e. estimate cluster for each data point on start view
  if (doOutput) cat("\nInitial E-Step on start view\n")
  CIdx = conceptIndicesSkm(startView, initialCenters, doOutput)
  if (doDebug) {
    print(initialCenters)
    str(CIdx)
  }


  # Main loop until objective function yields termination criterion
  # Uses lists to exchange views
  # Starts with processing otherView
  views = list(CV = otherView, HV = startView)
  conVs = list(CCV = NULL, HCV = initialCenters)
  conIs = list(CCI = NULL, HCI = CIdx)
  maxima = list(mStart = -Inf, mOther = -Inf)
  rounds = list(rStart = 0, rOther = 0)

  itCount = 0
  t = 0
  agreementRateDf = NULL
  indicesSVDf = NULL
  indicesOVDf = NULL
  repeat {
    t = t + 1
    itCount = itCount + 1

    # M-Step: cluster centers
    if (doOutput) cat(paste("\nM-Step:",itCount,"\n"))
    conVs$CCV = conceptVectorsSkm(views$CV,conIs$HCI, doOutput)
    if (doDebug) print(conVs$CCV)

    # E-Step: cluster assignment for data
    if (doOutput) cat(paste("\nE-Step:",itCount,"\n"))
    conIs$CCI = conceptIndicesSkm(views$CV,conVs$CCV, doOutput)
    if (doDebug) str(conIs$CCI)

    # Objective Function after each iteration of BOTH views
    if (itCount == 2) {
      itCount = 0
      rounds$rStart = rounds$rStart + 1
      rounds$rOther = rounds$rOther + 1
      ofStart = oFSkm(views$CV, conVs$CCV, conIs$CCI)
      ofOther = oFSkm(views$HV,  conVs$HCV, conIs$HCI)
      if (ofStart > maxima$mStart) {
        maxima$mStart = ofStart
        rounds$rStart = 0
      }
      if (ofOther > maxima$mOther) {
        maxima$mOther = ofOther
        rounds$rOther = 0
      }

      if (doOutput) {
        cat("\nObjective Function:\n")
        str(maxima)
        cat("\nRounds without improvement:\n")
        str(rounds)
      }

      agreementRate = ( ( sum( conIs$CCI == conIs$HCI ) )  / length( row.names(view1) ) )
      timeStamp = sub(" ", "-", Sys.time())
      newRow = data.frame(t=t,AgreementRate=sprintf("%.2f", agreementRate),TimeStamp=timeStamp)

      if (is.null(agreementRateDf)) agreementRateDf = newRow
      else agreementRateDf = rbind(agreementRateDf,newRow)

      if (is.null(indicesSVDf)) { indicesSVDf = data.frame(matrix(conIs$CCI,1,length(conIs$CCI))); names(indicesSVDf)=rownames(view1) }
      else indicesSVDf = rbind(indicesSVDf, conIs$CCI)

      if (is.null(indicesOVDf)) { indicesOVDf = data.frame(matrix(conIs$HCI,1,length(conIs$HCI))); names(indicesOVDf)=rownames(view1) }
      else indicesOVDf = rbind(indicesOVDf, conIs$HCI)
      
    }

    if ( (rounds$rStart >= nthresh) && (rounds$rOther >= nthresh) ) {
      break
    }

    views = list(CV = views$HV, HV = views$CV)
    conVs = list(CCV = conVs$HCV, HCV = conVs$CCV)
    conIs = list(CCI = conIs$HCI, HCI = conIs$CCI)

  }

  if (doOutput) {
    cat("\nTerminating EM step")
    cat("\nConsensus Means per Cluster and View\n")
  }
  consensusMeans = consensusMeansPerClVSkm( views$CV, views$HV, conIs$CCI, conIs$HCI )

  cat("\nFinal Cluster Assignment\n")
  finalIndices = assignFinIdxPerClSkm( views$CV, views$HV, consensusMeans )

  timeStamp = sub(" ", "-", Sys.time())
  #paste( t, "F", sprintf("%.2f", agreementRate), paste(finalIndices,sep="",collapse=":"), timeStamp , sep=",")

  #print(consensusMeans)
  #if (doOutput) str(finalIndices)

  list(finalIndices=finalIndices, agreementRate=agreementRateDf, indicesSv=indicesSVDf, indicesOv=indicesOVDf)

}


#' Multi-View Clustering using mixture of categoricals EM.
#' See: S. Bickel, T. Scheffer: Multi-View Clustering, ICDM 04.
#'
#' @param view1 View number one, a data frame with the same row names as view2. All columns numeric. Entries are natural numbers, starting from 1.
#' @param view2 View number two, a data frame with the same row names as view1. All columns numeric. Entries are natural numbers, starting from 1.
#' @param k The maximum number of clusters to create
#' @param startView String designating the view on which to perform the initial E step, one of "view1", "view2"
#' @param nthresh The number of iterations to run without improvement of the objective function
#' @param doOutput Whether output to the console should be done
#' @param doDebug Whether debug output to the console should be done (implies normal output)
#' @return A list reporting the final clustering, with names finalIndices, agreementRate, indicesSv, indicesOv. They designate final cluster indices as a vector, as well as agreement rate of the two views, and the individual indices given by the two views over the course of iterations as data frames.
#' @examples {
#' # Demo program, showing how to run Multi-
#' # View Clustering using Mixture of Binomials EM.
#' # AM, 2011
#'
#' # load toy data 'toyView1' and 'toyView2'
#' data(toyViews)
#' 
#' mvcmb(
#'   view1=toyView1,
#'   view2=toyView2,
#'   nthresh=20,
#'   k=4,
#'   startView="view1"
#' )
#' 
#' }


mvcmb <- function(
   view1=NULL,
   view2=NULL,
   k=Inf,
   startView="view1",
   nthresh=20,
   doOutput=F,
   doDebug=F
  ) {

  if (doDebug) doOutput=T
  if (doOutput) cat(paste("\nMVC (MB-EM)."))
  if (doOutput) cat(paste("\nk",k))
  if (doOutput) cat(paste("\nnthresh",nthresh))
  if (doOutput) cat("\n")

  if (startView == "view1") {
    startView = view1
    otherView = view2
  } 
  else {
    if (startView == "view2") {
      startView = view2
      otherView = view1
    }
    else {
      stop("startView argument wrong")
    }
  }
  checkViews(startView, otherView)
  uniqueVals <- viewsClasses(startView, otherView)
  if (doDebug) {
    cat("\nUnique values per view:\n")
    print(uniqueVals)
  }
  nClassesSV = length(uniqueVals$view1)
  nClassesOV = length(uniqueVals$view2)
  if (sum(uniqueVals$view1 != c(1:nClassesSV))) stop(paste("Start view classes must be in 1 ...",nClassesSV))
  if (sum(uniqueVals$view2 != c(1:nClassesOV))) stop(paste("Other view classes must be in 1 ...",nClassesOV))


  # # # Start MVCmb # # #

  # initialize variables m = #Documents
  m = NROW(startView)

  # Word Probabilities
  nWordsSV = NCOL(startView)
  nWordsOV = NCOL(otherView)
  PwSV = log( array( runif(k*nWordsSV*nClassesSV, min=0, max=1), c(k,nWordsSV,nClassesSV) ) ) # random init
  PwOV = log( array( runif(k*nWordsOV*nClassesOV, min=0, max=1), c(k,nWordsOV,nClassesOV) ) )
  if (doDebug) {
    cat("\nStart view Word Probabilities (by cluster):\n")
    print(exp(PwSV))
  }

  # Prior Probabilities
  alphaSV = log(rep( 1/k, k ))
 
  if (doOutput) cat("\nInitial E-Step on start view\n")
  # perform initial E-Step: estimate likelihood on start view

  # j-th row: document probs for cluster j
  Px = matrix ( rep(0, k*m), k, m, byrow=T ) 
  for (j in 1:k) {
    Px[j,] = estLogPxCatGthetaJ( startView, PwSV[j,,] )
  }
  if (doDebug) {
    cat("\nStart view Document Probabilities (by cluster):\n")
    print(exp(Px))
    cat("\nCol Margins:\n")
    print(exp(apply(Px,2,function(logx) logsum(logx))))
    #Sys.sleep(2)
  }
  
  # i-th row: cluster probs for document i
  DPS = Px + alphaSV # multiply each row with its alpha
  Pj = t(DPS) - apply( DPS, 2, function(logx) logsum(logx) ) # normalize column-wise (to obtain distribution on clusters)
  Pj[is.nan(Pj)] = 0
  if (doDebug) {
    cat("\nStart view Cluster Probabilities (by document):\n")
    print(exp(Pj))
    cat("\nRow Margins:\n")
    print(exp(apply(Pj,1,function(logx) logsum(logx))))
    #Sys.sleep(2)
  }

  # data
  views =  list ( CV = otherView, HV = startView )

  # distributions. To inherit correct matrices, init CVs to HVs, where applicable (Px, Pj).
  alpha = list ( CV = NULL, HV = alphaSV )   # cluster prior (vector)
  Pw = list ( CV = PwOV, HV = PwSV )         # word likelihood (matrix by cluster)
  Px = list ( CV = Px, HV = Px )             # document likelihood (matrix by cluster)
  Pj = list ( CV = Pj, HV = Pj )             # cluster posterior (matrix by document)

  # helper
  maxima = list ( mStart = -Inf, mOther = -Inf )
  rounds = list ( rStart = 0, rOther = 0 )
  DPS = list ( CV = NULL, HV = DPS)

  itCount = 0
  t = 0

  agreementRateDf = NULL
  indicesSVDf = NULL
  indicesOVDf = NULL
  repeat {

    t = t + 1
    itCount = itCount + 1

    # M-Step: Maximum Likelihood Parameters
    if (doOutput) cat(paste("\nM-Step:",itCount,"\n"))

    # word probabilities
    for (c in 1:dim(Pw$CV)[3]) {
      binInd = views$CV==c
      for (j in 1:k) {
        ViewClusterWgt = binInd * exp(Pj$HV[,j])  # view weighted by cluster (document-wise)
        Pw$CV[j,,c] = as.vector(colSums(ViewClusterWgt))    # words weighted by cluster (sum across documents)
      }
    }
    for (j in 1:k) {
      Pw$CV[j,,] = log(Pw$CV[j,,] / rowSums(Pw$CV[j,,]))
      Pw$CV[is.nan(Pw$CV)] = 0
    }

    if (doDebug) {
      cat("\nWord Probabilities (by cluster):\n")
      for (j in 1:k) {
        print(exp(Pw$CV[j,,]))
        cat("\nRow Margins:\n")
        print(exp(apply( Pw$CV[j,,], 1, function(logx) logsum(logx))))
      }
      #Sys.sleep(2)
    }

    # prior probabilities
    alpha$CV = apply(Pj$HV,2,function(logx) logsum(logx)) - log(m)            # mean of cluster probabilities across documents
    alpha$CV[is.nan(alpha$CV)] = 0
    if (doDebug) {
      cat("\nPrior Probabilities:\n")
      print(exp(alpha$CV))
      #Sys.sleep(2)
    }

    # E-Step: Likelihood
    if (doOutput) cat(paste("\nE-Step:",itCount,"\n"))
    for (j in 1:k) {
      Px$CV[j,] = estLogPxCatGthetaJ( views$CV, Pw$CV[j,,] )
    }
    
    if (doDebug) {
      cat("\nDocument Probabilities (by cluster):\n")
      print(exp(Px$CV))
      #Sys.sleep(2)
    }

    DPS$CV = Px$CV + alpha$CV

    Pj$CV = t(DPS$CV) - apply( DPS$CV, 2, function(logx) logsum(logx) ) # normalize column-wise (to obtain distribution on clusters)

    Pj$CV[is.nan(Pj$CV)] = 0
    if (doDebug) {
      cat("\nCluster Probabilities (by document):\n")
      print(exp(Pj$CV))
      cat("\nRow Margins:\n")
      print(exp(apply(Pj$CV,1,function(logx) logsum(logx))))
      #Sys.sleep(2)
    }
    

    # Objective Function after each iteration of BOTH views
    if (itCount == 2) {
      itCount = 0
      rounds$rStart = rounds$rStart + 1
      rounds$rOther = rounds$rOther + 1
      ofStart = oFMixBinEM(DPS$CV)
      ofOther = oFMixBinEM(DPS$HV)
      if (ofStart > maxima$mStart) {
        maxima$mStart = ofStart
        rounds$rStart = 0
      }
      if (ofOther > maxima$mOther) {
        maxima$mOther = ofOther
        rounds$rOther = 0
      }

      if (doOutput) {
       cat("\nObjective Function:\n")
       str(maxima)
       cat("\nRounds without improvement:\n")
       str(rounds)
      }

      agreementRate = agreementRateBinM( Pj$CV, Pj$HV )
      indicesSV = assignIdxPerClMBinEM( Pj$CV, Pj$CV )
      indicesOV = assignIdxPerClMBinEM( Pj$HV, Pj$HV )
      timeStamp = sub(" ", "-", Sys.time())
      newRow = data.frame(t=t,AgreementRate=sprintf("%.2f", agreementRate),TimeStamp=timeStamp)

      if (is.null(agreementRateDf)) agreementRateDf = newRow
      else agreementRateDf = rbind(agreementRateDf,newRow)

      if (is.null(indicesSVDf)) { indicesSVDf = data.frame(matrix(indicesSV,1,length(indicesSV))); names(indicesSVDf)=rownames(view1) }
      else indicesSVDf = rbind(indicesSVDf, indicesSV)

      if (is.null(indicesOVDf)) { indicesOVDf = data.frame(matrix(indicesOV,1,length(indicesOV))); names(indicesOVDf)=rownames(view1) }
      else indicesOVDf = rbind(indicesOVDf, indicesOV)

    }

    if ( (rounds$rStart >= nthresh) && (rounds$rOther >= nthresh) ) {
      break
    }

    views = list ( CV = views$HV, HV = views$CV )
    alpha = list ( CV = alpha$HV, HV = alpha$CV )
    Pw = list ( CV = Pw$HV, HV = Pw$CV )
    Px = list ( CV = Px$HV, HV = Px$CV )
    Pj = list ( CV = Pj$HV, HV = Pj$CV )
    DPS = list ( CV = DPS$HV, HV = DPS$CV )

  }
  
  if (doOutput) {
    cat("\nTerminating EM step")
    cat("\nFinal Cluster Assignment\n")
  }
  finalIndices = assignIdxPerClMBinEM( Pj$CV, Pj$HV )

  timeStamp = sub(" ", "-", Sys.time())
  #paste( t, "F", sprintf("%.2f", agreementRate), paste(finalIndices,sep="",collapse=":"), timeStamp , sep=",")
  #if (doOutput) str(finalIndices)

  list(finalIndices=finalIndices, agreementRate=agreementRateDf, indicesSv=indicesSVDf, indicesOv=indicesOVDf)

}


barplotMixCluster <-
function(combiRes, data = combiRes@mixmodOutput@data, nbCluster = combiRes@mixmodOutput@bestResult@nbCluster, x = combiRes@mixmodOutput@bestResult, variables=colnames(data), main=paste("Barplot of",variables), ...){
  # check the options
  if ( !is(combiRes,"MixmodCombi") )
    stop("'combiRes' must be a combiRes object!")
  if ( !is.matrix(data) & !is.data.frame(data) & !is.vector(data) )
    stop("'data' must be a vector, a data.frame or a matrix object!")
  if (!(is.numeric(nbCluster) & length(nbCluster) == 1 & nbCluster > 0 & nbCluster <= combiRes@mixmodOutput@bestResult@nbCluster) )
  	stop("'nbCluster' must be a positive numeric value between 1 and combiRes@mixmodOutput@bestResult@nbCluster")
  if ( !is(x,"MixmodResults") )
    stop("'x' must be a MixmodResults object!")
  if ( (length(variables)==0) & (ncol(data)>1) )
    stop("'variables' is empty!")
  if ( length(variables)>ncol(data) )
    stop("List of variables too long!")
  if ( sum(!(variables %in% colnames(data))) )
    stop("At least one variable is unknown!") 
  
  # get old par 
  op <- par(no.readonly = TRUE) # the whole list of settable par's.
  
  # get the indices of variables
  if ( ncol(data)==1 ){ indices<-1 }
  else { indices<-which(colnames(data) %in% variables) }
  nvar<-length(indices)
  
  if ( is(x@parameters, "GaussianParameter") ){
    stop("x must contain multinomial parameters. See hist() to plot Gaussian parameters.")
  }
  else if ( is(x@parameters, "MultinomialParameter") ){
    if ( !isQualitative(data) )
      stop("data must contain only qualitative variables!")      
    if ( is.factor(data) ) data<-as.integer(data)
    else if ( !is.vector(data) ){
      # loop over columns to check whether type is factor
      for ( j in 1:ncol(data) ){
        if ( is.factor(data[,j]) ) data[,j] <- as.integer(data[,j])
      }
    }
    # split the layout
    #par( mfrow = c( nvar, x@nbCluster+1 ) )
    #par(mar = par("mar")*.75)
    # split the layout
    if( nvar < 4 & nvar > 1) par( mfrow = c( 1, nvar ) )
    else if ( nvar >= 4 ){
      nrow<-round(sqrt(nvar))
      if (is.wholenumber(sqrt(nvar))) ncol<-sqrt(nvar)
      else ncol<-sqrt(nvar)+1
      par( mfrow = c( nrow, ncol ) ) 
    }
    # get number of observations
    nobs<-nrow(data)
    i<-1
    # loop over variables
    for (j in indices ){
      f<-x@parameters@factor[j]
      t<-colSums(data[ ,j] == matrix(rep(1:f, each = nobs), ncol = f)) # Unconditional frequencies * nobs
      names(t)<-1:f
      # freq<-barplot(as.vector(t), names.arg=names(t), xlab=xlab[i], main="", ...)
      proba<-matrix(0,nrow=nbCluster,ncol=f)
      for ( k in 1:nbCluster){
      	# Conditional frequencies of each factor value, conditional on the class k, for dimension j
        proba[k,] <- colSums(data[mixmodMap_V2M(combiRes@hierarchy[[nbCluster]]@partition)[,k] == 1, j] == matrix(rep(1:f, each = sum(mixmodMap_V2M(combiRes@hierarchy[[nbCluster]]@partition)[,k] == 1)), ncol = f))/sum(mixmodMap_V2M(combiRes@hierarchy[[nbCluster]]@partition)[,k] == 1)
        #prob<-barplot(proba_k[1:f], names.arg=names(t), main="", col=k+1, ylim=c(0,1), ...)
        #title(paste("Cluster",k),cex.main=1)
        #text(prob,proba_k[1:f]+(max(proba_k)/10),ifelse(proba_k[1:f]<0.01,format(proba_k[1:f],scientific=TRUE,digits=2),round(proba_k[1:f],2)),xpd=TRUE,font=2)
      }
      # Do the barplot and save the bar midpoints
      mp<-barplot(proba, beside = TRUE, axisnames = FALSE, names.arg=names(t), main="", ylab="Conditional frequency", ylim=c(0,1))
      # add unconditional frequencies
      for ( k in 1:f){
        lines(c(min(mp[,k])-1,max(mp[,k])+1),rep(t[k]/nobs,2),lty=2)
      }
      # add unconditional frequency legend only for the first variable
      if ( i == 1 ) text(min(mp[,1])-.5,t[1]/nobs,labels="Unconditional frequency",cex=1, pos=4, adj=c(0,0), font=2) 
      # add title
      title(main[i],cex.main=1)
      # Add the individual bar labels
      mtext(1, at = mp, text = paste("C",1:nbCluster), line = 0, cex = 0.5)
      # Add the group labels for each pair
      #mtext(1, at = rbind(mp[1,],colMeans(mp[-1,])), text = rep(c("freq", "proba"), f), line = 1, cex = 0.75)
      # Add the labels for each group
      mtext(1, at = colMeans(mp), text = names(t), line = 2)
      i<-i+1
    }
  }
  else{
    stop("Uknown type of parameters!")
  }

	currTitle <- 	if (nbCluster == combiRes@mixmodOutput@bestResult@nbCluster & !nbCluster == combiRes@ICLNbCluster) paste( "BIC solution (",as.character(nbCluster)," clusters)", sep = "") else if (nbCluster == combiRes@ICLNbCluster & !nbCluster == combiRes@mixmodOutput@bestResult@nbCluster) paste( "Combined solution with ",as.character(nbCluster)," clusters (Number of clusters selected with ICL)", sep = "") else if (nbCluster == combiRes@ICLNbCluster & nbCluster == combiRes@mixmodOutput@bestResult@nbCluster) paste( "BIC and ICL solution (",as.character(nbCluster)," clusters)", sep = "") else	paste( "Combined solution with ",as.character(nbCluster)," clusters", sep = "")
  	par(oma=c(0,0,1.5,0))#, mar = c(5,4,0,2)+0.1)
	title(currTitle, outer = TRUE)

  par(op)
}

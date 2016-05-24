scatterPlot <- function(model.obj,data.obj=NULL,xvar=NULL,yvar=NULL){

  if ( class(model.obj)[[1]] != "Rcpp_bei" ){
    stop("model.obj should be 'Rcpp_bei'")
  }
  
  logD_input = log(model.obj$Y.input)
  logD_edited = log(model.obj$Y.edited)
  n_record = dim(logD_input)[[1]]
  n_var = dim(logD_input)[[2]]
  whichFaulty = model.obj$FaultyRecordID
  whichPass = c(1:n_record)[-whichFaulty] 
  
  logRatio_U = array(NA,c(n_var,n_var))
  if (!is.null(data.obj)){
    Ratio = data.obj$ratio
    if (!is.null(Ratio)){
      n_ratio = dim(Ratio)[[1]]
      for (i_row in 1:n_ratio){
        if (sum(Ratio[i_row,]==1)>0){
          var_NUM = which(Ratio[i_row,]==1)
          var_DEN = which( (Ratio[i_row,]!=1)&(Ratio[i_row,]!=0) )
          logRatio_U[var_NUM,var_DEN] = log(-Ratio[i_row,var_DEN])
        }
      }
    } # if (!is.null(data.obj))
  } # if (!is.null(ratio)) 
    
  if ( is.null(xvar) & is.null(yvar) ){
    
    for (i_var in 1:(n_var-1)){
      minX = min(logD_edited[,i_var])-2 ; maxX = max(logD_edited[,i_var])+2
      
      for (j_var in (i_var+1):n_var){
      
        minY = min(logD_edited[,j_var])-2 ; maxY = max(logD_edited[,j_var])+2
        
        par( mfrow=c(1,2), ask=TRUE )
        
        plot( logD_input[whichPass,i_var], logD_input[whichPass,j_var], pch=20, xlim=c(minX,maxX), ylim=c(minY,maxY), xlab=paste("log (Var ",i_var,")",sep=""), ylab=paste("log (Var ",j_var,")",sep=""), col="skyblue", cex = 1.0, main="log (Y.input)" )
        points( logD_input[whichFaulty,i_var], logD_input[whichFaulty,j_var], pch=20, col="blue", cex = 0.4 )
        if (!is.na(logRatio_U[j_var,i_var])){
          points( c(-100,200),c(-100+logRatio_U[j_var,i_var],200+logRatio_U[j_var,i_var]),col="red",type="l",lty="dotted",lwd=1.5) # upper bound 
        }
        if (!is.na(logRatio_U[i_var,j_var])){
          points( c(-100,200),c(-100-logRatio_U[i_var,j_var],200-logRatio_U[i_var,j_var]),col="red",type="l",lty="dotted",lwd=1.5)
        }
        
        plot( logD_edited[whichPass,i_var], logD_edited[whichPass,j_var], pch=20, xlim=c(minX,maxX), ylim=c(minY,maxY), xlab=paste("log (Var ",i_var,")",sep=""), ylab=paste("log (Var ",j_var,")",sep=""), col="skyblue", cex = 1.0, main="log (Y.edited)" )
        points( logD_edited[whichFaulty,i_var], logD_edited[whichFaulty,j_var], pch=20, col="blue", cex = 0.4 )
        if (!is.na(logRatio_U[j_var,i_var])){
          points( c(-100,200),c(-100+logRatio_U[j_var,i_var],200+logRatio_U[j_var,i_var]),col="red",type="l",lty="dotted",lwd=1.5) # upper bound 
        }
        if (!is.na(logRatio_U[i_var,j_var])){
          points( c(-100,200),c(-100-logRatio_U[i_var,j_var],200-logRatio_U[i_var,j_var]),col="red",type="l",lty="dotted",lwd=1.5)
        }
        
      } # for (j_var) 
    } # for (i_var)
  
  } # if ( is.null(xvar) & is.null(yvar) )
  
  if ( (!is.null(xvar)) & is.null(yvar) ){
    
    if ( (xvar<1)|(xvar>n_var) ) stop("xvar should be an interger between 1 and p") ;
    
      i_var = xvar
    
      minX = min(logD_edited[,i_var])-2 ; maxX = max(logD_edited[,i_var])+2
      
      for (j_var in c(1:n_var)[-i_var]){
        
        minY = min(logD_edited[,j_var])-2 ; maxY = max(logD_edited[,j_var])+2
        
        par( mfrow=c(1,2), ask=TRUE )
        
        plot( logD_input[whichPass,i_var], logD_input[whichPass,j_var], pch=20, xlim=c(minX,maxX), ylim=c(minY,maxY), xlab=paste("log (Var ",i_var,")",sep=""), ylab=paste("log (Var ",j_var,")",sep=""), col="skyblue", cex = 1.0, main="log (Y.input)" )
        points( logD_input[whichFaulty,i_var], logD_input[whichFaulty,j_var], pch=20, col="blue", cex = 0.4 )
        if (!is.na(logRatio_U[j_var,i_var])){
          points( c(-100,200),c(-100+logRatio_U[j_var,i_var],200+logRatio_U[j_var,i_var]),col="red",type="l",lty="dotted",lwd=1.5) # upper bound 
        }
        if (!is.na(logRatio_U[i_var,j_var])){
          points( c(-100,200),c(-100-logRatio_U[i_var,j_var],200-logRatio_U[i_var,j_var]),col="red",type="l",lty="dotted",lwd=1.5)
        }
        
        plot( logD_edited[whichPass,i_var], logD_edited[whichPass,j_var], pch=20, xlim=c(minX,maxX), ylim=c(minY,maxY), xlab=paste("log (Var ",i_var,")",sep=""), ylab=paste("log (Var ",j_var,")",sep=""), col="skyblue", cex = 1.0, main="log (Y.edited)" )
        points( logD_edited[whichFaulty,i_var], logD_edited[whichFaulty,j_var], pch=20, col="blue", cex = 0.4 )
        if (!is.na(logRatio_U[j_var,i_var])){
          points( c(-100,200),c(-100+logRatio_U[j_var,i_var],200+logRatio_U[j_var,i_var]),col="red",type="l",lty="dotted",lwd=1.5) # upper bound 
        }
        if (!is.na(logRatio_U[i_var,j_var])){
          points( c(-100,200),c(-100-logRatio_U[i_var,j_var],200-logRatio_U[i_var,j_var]),col="red",type="l",lty="dotted",lwd=1.5)
        }
        
      } # for (j_var) 
    
  } # if ( (!is.null(xvar)) & is.null(yvar) )
  
  if ( is.null(xvar) & (!is.null(yvar)) ){
    
    if ( (yvar<1)|(yvar>n_var) ) stop("yvar should be an interger between 1 and p") ;
    
    j_var = yvar
    
    for (i_var in c(1:n_var)[-j_var]){

        minX = min(logD_edited[,i_var])-2 ; maxX = max(logD_edited[,i_var])+2
        
        minY = min(logD_edited[,j_var])-2 ; maxY = max(logD_edited[,j_var])+2
        
        par( mfrow=c(1,2), ask=TRUE )
        
        plot( logD_input[whichPass,i_var], logD_input[whichPass,j_var], pch=20, xlim=c(minX,maxX), ylim=c(minY,maxY), xlab=paste("log (Var ",i_var,")",sep=""), ylab=paste("log (Var ",j_var,")",sep=""), col="skyblue", cex = 1.0, main="log (Y.input)" )
        points( logD_input[whichFaulty,i_var], logD_input[whichFaulty,j_var], pch=20, col="blue", cex = 0.4 )
        if (!is.na(logRatio_U[j_var,i_var])){
          points( c(-100,200),c(-100+logRatio_U[j_var,i_var],200+logRatio_U[j_var,i_var]),col="red",type="l",lty="dotted",lwd=1.5) # upper bound 
        }
        if (!is.na(logRatio_U[i_var,j_var])){
          points( c(-100,200),c(-100-logRatio_U[i_var,j_var],200-logRatio_U[i_var,j_var]),col="red",type="l",lty="dotted",lwd=1.5)
        }
        
        plot( logD_edited[whichPass,i_var], logD_edited[whichPass,j_var], pch=20, xlim=c(minX,maxX), ylim=c(minY,maxY), xlab=paste("log (Var ",i_var,")",sep=""), ylab=paste("log (Var ",j_var,")",sep=""), col="skyblue", cex = 1.0, main="log (Y.edited)" )
        points( logD_edited[whichFaulty,i_var], logD_edited[whichFaulty,j_var], pch=20, col="blue", cex = 0.4 )
        if (!is.na(logRatio_U[j_var,i_var])){
          points( c(-100,200),c(-100+logRatio_U[j_var,i_var],200+logRatio_U[j_var,i_var]),col="red",type="l",lty="dotted",lwd=1.5) # upper bound 
        }
        if (!is.na(logRatio_U[i_var,j_var])){
          points( c(-100,200),c(-100-logRatio_U[i_var,j_var],200-logRatio_U[i_var,j_var]),col="red",type="l",lty="dotted",lwd=1.5)
        }
        
    } # for (i_var)
    
  } # if ( is.null(xvar) & (!is.null(yvar)) )
  
  if ( (!is.null(xvar)) & (!is.null(yvar)) ){
    
    if ( (xvar<1)|(xvar>n_var) ) stop("xvar should be an interger between 1 and p") ;
    if ( (yvar<1)|(yvar>n_var) ) stop("yvar should be an interger between 1 and p") ;
    
    i_var = xvar ; j_var = yvar ; 
      
      minX = min(logD_edited[,i_var])-2 ; maxX = max(logD_edited[,i_var])+2
      
      minY = min(logD_edited[,j_var])-2 ; maxY = max(logD_edited[,j_var])+2
      
      par( mfrow=c(1,2), ask=FALSE )
      
      plot( logD_input[whichPass,i_var], logD_input[whichPass,j_var], pch=20, xlim=c(minX,maxX), ylim=c(minY,maxY), xlab=paste("log (Var ",i_var,")",sep=""), ylab=paste("log (Var ",j_var,")",sep=""), col="skyblue", cex = 1.0, main="log (Y.input)" )
      points( logD_input[whichFaulty,i_var], logD_input[whichFaulty,j_var], pch=20, col="blue", cex = 0.4 )
      if (!is.na(logRatio_U[j_var,i_var])){
        points( c(-100,200),c(-100+logRatio_U[j_var,i_var],200+logRatio_U[j_var,i_var]),col="red",type="l",lty="dotted",lwd=1.5) # upper bound 
      }
      if (!is.na(logRatio_U[i_var,j_var])){
        points( c(-100,200),c(-100-logRatio_U[i_var,j_var],200-logRatio_U[i_var,j_var]),col="red",type="l",lty="dotted",lwd=1.5)
      }
      
      plot( logD_edited[whichPass,i_var], logD_edited[whichPass,j_var], pch=20, xlim=c(minX,maxX), ylim=c(minY,maxY), xlab=paste("log (Var ",i_var,")",sep=""), ylab=paste("log (Var ",j_var,")",sep=""), col="skyblue", cex = 1.0, main="log (Y.edited)" )
      points( logD_edited[whichFaulty,i_var], logD_edited[whichFaulty,j_var], pch=20, col="blue", cex = 0.4 )
      if (!is.na(logRatio_U[j_var,i_var])){
        points( c(-100,200),c(-100+logRatio_U[j_var,i_var],200+logRatio_U[j_var,i_var]),col="red",type="l",lty="dotted",lwd=1.5) # upper bound 
      }
      if (!is.na(logRatio_U[i_var,j_var])){
        points( c(-100,200),c(-100-logRatio_U[i_var,j_var],200-logRatio_U[i_var,j_var]),col="red",type="l",lty="dotted",lwd=1.5)
      }
    
  } # if ( (!is.null(xvar)) & (!is.null(yvar)) )
  
} # scatterPlot
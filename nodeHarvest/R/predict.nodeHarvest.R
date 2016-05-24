predict.nodeHarvest <-
function(object, newdata=NULL, explain=NULL, maxshow=5, weight=sapply(object[["nodes"]],attr,"weight"), ...){
  if(is.null(newdata)){
    predtest <-  object$predicted 
  }else{

      if(is.null(nrow(newdata)))  newdata <- as.data.frame(matrix(unlist(newdata),nrow=1))

      if(is.data.frame(newdata)){
      for (k in 1:ncol(newdata)){
        if(!class(newdata[,k])%in%c("numeric","factor")){
          newdata[,k] <- as.numeric(newdata[,k])
        }
      }
    }
        
    Z <- object[["nodes"]]
    
      weight <- sapply(Z,attr,"weight")
      ITEST <- getI(Z, newdata,Y=NULL, mode="mean")$I
      ISIGN <- getI(Z, newdata,Y=NULL, mode="member")$I
    predtest <- as.numeric( ITEST %*% weight ) / as.numeric( ISIGN %*% weight)

    if(!is.null(object[["bias"]])){
      predtest <- object[["bias"]][1]+object[["bias"]][2]*predtest
    }else{
    
      if(!is.null(explain) & maxshow>0){
        
        for (k in explain){
          whichnodes <- which(ITEST[k,]!=0)
          ord <- order(sapply(Z,attr,"weight")[whichnodes],decreasing=TRUE)
          lw <- length(whichnodes)
          cat("\n \t Observation ",k," has a predicted value ", signif(predtest[k],3) ,"\n \t  since this is the weighted average response across the ", lw , if(lw>1) " nodes" else " node"," it is a member of:",sep="")
          if( lw > maxshow){
            cat("\n \t\t (showing the most important maxshow=",maxshow,"nodes only)")
            lw <- maxshow
          }
          for (kc in 1:lw){
            
            nodenumber <- (whichnodes[ord])[kc]
            cat("\n\n \t \t ",kc,") Node ",nodenumber,", containing ",attr(Z[[nodenumber]],"n")," training observations, with node mean ", signif(attr(Z[[nodenumber]],"mean"),3)," and weight ", signif( attr(Z[[nodenumber]],"weight"),3)," :",sep="")
            if(attr(Z[[nodenumber]],"depth")>0){
              node <- drawtext(Z,nodenumber,varnames=object[["varnames"]],plot=FALSE)
              for (lc in 1:length(node)){
                cat("\n\t\t \t  ", node[[lc]])
              }
            }else{
              cat("\n\t\t\t ROOT NODE")
            }
          }
          cat("\n")
          
        }
      }
    }
    
  }
  return(predtest)

}


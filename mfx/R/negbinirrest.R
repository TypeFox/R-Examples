negbinirrest <-
function(formula, data, robust = FALSE, clustervar1 = NULL, 
                        clustervar2 = NULL, start = NULL, control = glm.control()){
  
  if(is.null(formula)){
    stop("formula is missing")
  }
  if(!is.data.frame(data)){
    stop("data arguement must contain data.frame object")
  }
  # cluster sort part
  if(is.null(clustervar1) & !is.null(clustervar2)){
    stop("use clustervar1 arguement before clustervar2 arguement")
  }    
  if(!is.null(clustervar1)){
    if(is.null(clustervar2)){
      if(!(clustervar1 %in% names(data))){
        stop("clustervar1 not in data.frame object")
      }    
      data = data.frame(model.frame(formula, data, na.action=NULL),data[,clustervar1])
      names(data)[dim(data)[2]] = clustervar1
      data=na.omit(data)
    }
    if(!is.null(clustervar2)){
      if(!(clustervar1 %in% names(data))){
        stop("clustervar1 not in data.frame object")
      }    
      if(!(clustervar2 %in% names(data))){
        stop("clustervar2 not in data.frame object")
      }    
      data = data.frame(model.frame(formula, data, na.action=NULL),
                        data[,c(clustervar1,clustervar2)])
      names(data)[c(dim(data)[2]-1):dim(data)[2]] = c(clustervar1,clustervar2)
      data=na.omit(data)
    }
  }
  fit = glm.nb(formula, data=data, x = T, start = NULL, control = control)    
  x1 = model.matrix(fit)
  if (any(alias <- is.na(coef(fit)))) {
    x1 <- x1[, !alias, drop = FALSE]
  }
  
  # get variances
  vcv = vcov(fit) 
  
  if(robust){
    if(is.null(clustervar1)){
      # white correction
      vcv = vcovHC(fit, type = "HC0")
    } else {
      if(is.null(clustervar2)){
        vcv = clusterVCV(data=data, fm=fit, cluster1=clustervar1,cluster2=NULL)
      } else {
        vcv = clusterVCV(data=data, fm=fit, cluster1=clustervar1,cluster2=clustervar2)
      }
    }
  }
  
  if(robust==FALSE & is.null(clustervar1)==FALSE){
    if(is.null(clustervar2)){
      vcv = clusterVCV(data=data, fm=fit, cluster1=clustervar1,cluster2=NULL)
    } else {
      vcv = clusterVCV(data=data, fm=fit, cluster1=clustervar1,cluster2=clustervar2)
    }
  }
  
  irr = data.frame(irr=exp(na.omit(coef(fit))), se=NA)
  # get standard errors
  gr = diag(exp(na.omit(coef(fit))))
  
  irr$se = sqrt(diag(gr %*% vcv %*% t(gr)))
  
  # pick out constant and remove from mfx table
  temp1 = apply(x1,2,function(x)length(table(x))==1)
  const = names(temp1[temp1==TRUE])
  irr = irr[row.names(irr)!=const,]
  
  output = list(fit=fit, irr=irr)
  return(output)
}

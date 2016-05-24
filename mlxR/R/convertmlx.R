
convertmlx <- function(data, dataIn,trt,iop.group,id.out=FALSE,id.ori=NULL,gr.ori=NULL){
  
  g <- dataIn$group
  
  #   if (is.null(g))
  #     iop.group <- 0
  #   else
  #     iop.group <- 1
  
  iop.gout <- 0
  N <- 0
  for(k in seq(1,length(g))){
    g[[k]]$output=NULL
    gk.size <- prod(g[[k]]$size)
    N <- N + gk.size
    if (gk.size > 1)
      iop.gout <- 1
  }
  
  cv <- dataIn$catvar
  var <- dataIn$variability
  if (length(unique(var$id))==1)
    var$id <- NULL
  
  if (length(g)>1){
    gr=numeric(0)
    for(k in seq(1,length(g))){
      pgk <- prod(g[[k]]$size)
      gr=c(gr,rep(k,pgk))
    }
  }else{
    iop.group <- 0
  }  
  
  
  gr=numeric(0)
  for(k in seq(1,length(g))){
    gr=c(gr,rep(k,prod(g[[k]]$size)))
  } 
  
  dd=list()
  if (sum(gr)==1){
    df <- data.frame()
  }else{
    df <- NULL
  }
  for(k in seq(1,length(data))){
    ak=data[[k]]
    if (length(unlist(ak$value))>0)
    {
      nk =length(ak$value)
      vk=numeric(0)
      idk=numeric(0)
      tk=numeric(0)
      gk=numeric(0)
      for(i in seq(1,nk)){
        vki = ak$value[[i]]
        vk=c(vk, vki)
        nki=length(vki)
        idk=c(idk, rep(i,nki))
        if(iop.group==1)
          gk=c(gk, rep(gr[i],nki))
        
        tki = ak$time[[i]]
        tk=c(tk, tki)
      }
      ick <- which(ak$name==cv$name)
      if (length(ick)>0){
        vk <- cv$categories[[ick]][vk]
      }else{
        if (isfield(ak,"categories")){
          vk <- ak$categories[vk]
        }
      }
      if(length(tk)>0){
        iop.tk=1
      }else{
        iop.tk=0
      }
      # if(length(unique(idk))>1){
      if(N>1){
        iop.id=1
      }else{
        iop.id=0
      }
      if(iop.id==0){
        if(iop.tk==1){
          dk=data.frame(time=tk, value=vk)
        }else{
          dk=data.frame(value=vk)
        }
      }else{
        if(iop.group==0){
          if(iop.tk==1){
            dk=data.frame(id=factor(idk), time=tk, value=vk)
          }else{
            dk=data.frame(id=factor(idk), value=vk)
          }
        }else{
          if(iop.tk==1){
            dk=data.frame(id=factor(idk), group=factor(gk), time=tk, value=vk)
          }else{
            dk=data.frame(id=factor(idk), group=factor(gk), value=vk)
          }
        }
      }
      names(dk)[names(dk)=="value"] <- ak$name
      
      if (id.out==TRUE){
        if (is.null(dk$id)){
          dk$id <- 1
          nk <- length(dk)
          dk <- dk[,c(nk,(1:(nk-1)))]
        }
        if (is.null(dk$group)){
          dk$group <- 1
          nk <- length(dk)
          dk <- dk[,c(1,nk,(2:(nk-1)))]
        }
      }
      if (iop.tk==0){
        if(iop.id==0){
          df <- c(df,dk)
        }else{
          if (is.null(df)){
            df <- dk
          }else{
            df <- cbind(df,dk)
            
            j1 <- which(names(df)=="id")
            if (length(j1>1))
              j1 <- j1[-1]
            
            j2 <- which(names(df)=="group")
            if (length(j2>1))
              j2 <- j2[-1]
            df <- df[-c(j1,j2)]
            
          }
        }
      }
      if(iop.id==0)
        df <- data.frame(df)
      
      attr(dk,"name")=ak$name
      attr(dk,"type")=ak$label
      dd[[ak$name]] = dk
    }
  }
  
  #   if (length(df)>0){
  #     if (length(df[[1]])>1)
  #       dd$parameter = df
  #   }
  # if (length(df)>1){
  if (!is.null(df)  && nrow(df)>0)
  {
    attr(df,"type") <- "parameter"
    dd$parameter = df
    dd[names(dd$parameter)] <- NULL
  }
  
  if (!is.null(var)){
    v <- data.frame(var$value)
    names(v) <- var$colNames
    if (N>1){
      id0 <- 0
      vv <- NULL
      for (j in (1:length(g))){
        vj <- v[which(v$id==j),]
        dj <- nrow(vj)
        gj.size <- prod(g[[j]]$size)
        vj <- do.call("rbind", replicate(gj.size, vj, simplify = FALSE))
        vj$id <- rep((1:gj.size),each=dj) +id0
        id0 <- id0 + gj.size
        if (length(g)>1 & iop.group==1)
          vj$group <- j
        vv <- rbind(vv,vj)
      } 
    }else{
      vv <- v
      vv$id <- NULL    
    }
    attr(vv,"type") <- "varlevel"      
    dd$varlevel <- vv
    
    vv$group <- NULL
    vv$id <- NULL
    for(k in seq(1,length(dd))){
      if (is.null(dd[[k]]$time)){
        vdk <- cbind(vv, dd[[k]])
        j=which(names(vdk)=="group")
        if (length(j)>0){
          u=(1:dim(vdk)[2])
          vdk <- vdk[,c(j,u[-j])]
        }
        j=which(names(vdk)=="id")
        if (length(j)>0){
          u=(1:dim(vdk)[2])
          vdk <- vdk[,c(j,u[-j])]
        }
        dd[[k]] <- vdk
      }
    }   
  }
  
  #   if (!is.null(id.ori)){
  #     for(k in seq(1,length(dd))){
  #       if (!is.null(dd[[k]]$id)){
  #         dd[[k]]$id <- id.ori[dd[[k]]$id]
  #       }
  #     }
  #   }
  
  if (iop.gout==1)
    dd$group=g
  
  if (!is.null(dataIn$regressor)){
    reg <- data.frame(dataIn$regressor$value)
    names(reg) <- dataIn$regressor$colNames
    nreg <- ncol(reg)-2
    for (k in (1:nreg)){
      xk <- reg[k+2]
      nk <- names(xk)
      idk <- which(!is.na(xk))
      regk <- reg[idk,c(1,2,k+2)]
      if (N>1){
        id0 <- 0
        reg.gk <- NULL
        for (j in (1:length(g))){
          regkj <- regk[which(regk$id==j),]
          dj <- nrow(regkj)
          gj.size <- prod(g[[j]]$size)
          regkj <- do.call("rbind", replicate(gj.size, regkj, simplify = FALSE))
          regkj$id <- rep((1:gj.size),each=dj) +id0
          if (length(g)>1 & iop.group==1)
            regkj$group <- j
          id0 <- id0 + gj.size
          reg.gk <- rbind(reg.gk,regkj)
        }        
      }else{
        regk$id <- NULL
        reg.gk <- regk
      }
      attr(reg.gk,"type") <- "regressor"      
      attr(reg.gk,"name") <- nk
      dd[[nk]] <- reg.gk
    }
  }
  
  if (!is.null(trt))
  {
    if (N>1)
    {
      ng <- length(trt)
      id0 <- 0
      treatment <- NULL
      for (k in (1:ng))
      {
        trtk <- as.data.frame(trt[[k]])
        dk <- nrow(trtk)
        nk <- prod(dataIn$group[[k]]$size)
        if (dk>0)
        {
          trtk <- trtk[,c("time","amount","rate","type")]
          trtk <- do.call("rbind", replicate(nk, trtk, simplify = FALSE))
          if (ng>1 & iop.group==1)
            trtk <- cbind(list(group=k),trtk)
          trtk <- cbind(list(id=rep(((1:nk)+id0),each=dk)),trtk)
          treatment <- rbind(treatment,trtk)
        }
        id0 <- id0 + nk
      }
    }
    else
    {
      treatment <- as.data.frame(trt)
      treatment <- treatment[,c("time","amount","rate","type")]
      treatment$id <- NULL
    }
    if (all(unique(treatment$type)==1))
      treatment$type <- NULL
    if (all(unique(treatment$rate)==Inf))
      treatment$rate <- NULL
    attr(treatment,"type") <- "treatment"
    dd$treatment <- treatment
  }
  
  if (!is.null(gr.ori)){
    for (k in (1:length(dd))){
      ddk <- dd[[k]]
      if (!is.null(ddk$id)){
        ddk$group <- gr.ori[as.numeric(as.character(ddk$id))]
        dd[[k]] <- ddk
      }
    }
  }
  return(dd)
  
}
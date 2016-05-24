hformat  <-  function(list.input)
{
  #hformat 
  treatment=list.input$treatment
  parameter=list.input$parameter
  output=list.input$output
  regressor=list.input$regressor
  varlevel=list.input$varlevel
  group=list.input$group
  id=list.input$id
  if (isfield(id,"newId"))
    id <- id$newId
  nv=c('treatment','parameter','output','regressor','varlevel')
  lv=list(treatment,parameter,output,regressor,varlevel)
  
  #------------------------------------
  if (!is.null(group))
  {
    if (!is.null(names(group))) 
      group <- list(group) 
    N=length(group)
    
    model <- list.input$model
    
    if (identical(file_ext(model),"R")) 
      Rfile <- TRUE 
    else 
      Rfile <- FALSE
    if ( !is.null(model) && exists(model, mode="function") )
      Rsource <- TRUE
    else 
      Rsource <- FALSE
    
    if (Rsource || Rfile)
    {
      for (k in (1:N))
      {
        if (is.null(group[[k]]$size))
          group[[k]]$size <- 1
        if (length(group$size)>1)
          stop("Define group$size as a scalar instead of a vector.")
        if (!is.null(group[[k]]$level))
          warning("'level' defined in 'group' is not used with a R model.")
      }
    }
    else
    {    
      model.info <- parse.model(model)
      for (k in (1:N))
      {
        if (is.null(group[[k]]$size))
        {
          group[[k]]$size <- rep(1, length(model.info$level))
          group[[k]]$level <- model.info$level
        }
        else
        {
          if (is.null(group[[k]]$level))
          {
            if (length(group[[k]]$size)>1)
              stop("levels of randomization associated to the sizes have not been defined in 'group'", call.=FALSE)      
            if ("population" %in% model.info$level)
              group[[k]]$level <- "population"
            else if ("covariate" %in% model.info$level)
              group[[k]]$level <- "covariate"
            else if ("individual" %in% model.info$level)
              group[[k]]$level <- "individual"
            else
              group[[k]]$level <- "longitudinal"
            
            if (group[[k]]$level != "longitudinal")
              warning(paste0("level of randomization has not been defined in group ", k,": '", 
                             group[[k]]$level, "' is used"), call.=FALSE)      
            
          }
          else
          {
            if (!is.null(group[[k]]$size) & (length(group[[k]]$size)!=length(group[[k]]$level)))
              stop("'level' and 'size' defined in 'group' have different lengths.", call.=FALSE)      
            gk <- group[[k]]$level
            sk <- rep(1, length(model.info$level))
            for (jk in (1:length(gk)))
            {
              rk <- grep(gk[jk],model.info$level)
              if (length(rk)>0)
                sk[rk] <- group[[k]]$size[jk]
            }
            group[[k]]$size <- sk
            group[[k]]$level <- model.info$level
          }
        }
      }
    }
  }
  else
    N <- nrow(list.input$id)
  # N=NULL
  
  Nid<-NULL
  for (k in seq(1,length(lv))) 
  {
    lvk <- lv[[k]]
    if (!is.null(names(lvk)))
    {  
      lvk <- list(lvk) 
      lv[[k]] <- lvk
    }
    if (length(lvk)>0)
    {
      for (j in seq(1,length(lvk))) 
      {
        lvkj <- lvk[[j]]
        # dkj <- dim(lvkj)
        #     if (!is.null(dkj)){
        #   N=c(N,dkj[1])
        if (isfield(lvkj,"id"))
        {
          # N <- c(N,length(unique(lvkj$id)))
          # lu <- unique(lvkj$id)
          lvkj$id <- match(lvkj$id, id)
          Nid <- c(Nid,lvkj$id)
        }
        else if(isfield(lvkj,"design"))
        {
          # lu <- unique(lvkj$design$id)
          lvkj$design$id <- match(lvkj$design$id, id)
          Nid <- c(Nid,lvkj$design$id)
          #             names(lvkj)[names(lvkj)=="design"]<-"time"
          lv[[k]][[j]]=lvkj
        }
        else if(isfield(lvkj,"time"))
        {
          if (isfield(lvkj$time,"id"))
          {
            # lu <- unique(lvkj$time$id)
            lvkj$time$id <- match(lvkj$time$id, id)
            Nid <- c(Nid,lvkj$time$id)
          }
        }
        else if(isfield(lvkj,'header'))
        {
          warning("deprecated syntax:  use 'colNames' instead of 'header'",immediate.=TRUE)
          names(lvkj)[names(lvkj)=="header"]<-"colNames"
        }
        
        if (isfield(lvkj,'colNames'))
        {
          #N=c(N,length(unique(lvkj$value[,1]))) #assuming that first column = id
          Nid=c(Nid,lvkj$value[,1]) #assuming that first column = id
        }
        lv[[k]][[j]]=lvkj
      }
    }
  }    
  Nid = length(unique(Nid))
  
  if (is.null(N))
    N <- Nid 
  
  if (is.null(N)||N==0)
    N=1
  else
    N = unique(N)
  
  if (is.null(group))
  {
    group <- vector("list",N)
    for (i in seq(1,N))
      group[[i]] <- list(size=1)
  }
  
  id.ori=list()
  #---  parameters
  iv <- which(nv=="parameter")
  parameter <- vector("list", N)
  if (!is.null(lv[[iv]])){
    r <- format.parameter(parameter,lv[[iv]],seq(1,N))
    if (length(r$id)>0) {id.ori=c(id.ori,r$id)}
    parameter <- r$parameter
  }
  for (i in seq(1,N))
  {
    gi <- group[[i]]
    if (isfield(gi,'parameter')) 
    {
      r <- format.parameter(parameter,gi$parameter,i)
      parameter <- r$parameter
      group[[i]]$parameter <- NULL
    }
  }
  parameter <- merge.parameter(parameter)
  list.output=list(parameter=parameter)  
  
  #---  outputs
  output <- list(individual=vector("list", N),group=vector("list", N), id=id)
  for (i in seq(1,N))
  {
    output$individual[[i]]=list(name=NULL,time=NULL)
    output$group[[i]]=list(name=NULL)
  }
  iv <- which(nv=="output")
  r <- format.output(output,lv[[iv]],seq(1,N))
  if (length(r$id)>0) {id.ori=c(id.ori,r$id)}
  output <- r$output
  for (i in seq(1,N))
  {
    gi <- group[[i]]
    if (isfield(gi,'output')) 
    {
      r <- format.output(output,gi$output,i)
      output <- r$output
    }
    group[[i]]$output <- output$group[[i]]$name
  }
  output <- output$individual
  list.output$output <- output
  
  #---  treatments
  iv <- which(nv=="treatment")
  if (!is.null(lv[[iv]]))
  {
    r <- format.treatment(lv[[iv]],seq(1,N))
    if (length(r$id)>0) 
      id.ori=c(id.ori,r$id)
    treatment <- r$treatment
  }
  for (i in seq(1,N))
  {
    gi <- group[[i]]
    if (isfield(gi,'treatment')) 
    {
      r <- format.treatment(gi$treatment,i)
      pgi <- r$treatment
      treatment <- c(treatment, pgi)
      group[[i]]$treatment <- NULL
    }
  }
  if (!is.null(treatment))
  {
    r <- merge.treatment(treatment,N)
    list.output <- c(list.output, r)
  }
  
  #---  regressor
  iv <- which(nv=="regressor")
  if (!is.null(lv[[iv]]))
  {
    r <- format.regressor(lv[[iv]],seq(1,N))
    if (length(r$id)>0) 
      id.ori=c(id.ori,r$id)
    regressor <- r$regressor
  }
  for (i in seq(1,N))
  {
    gi <- group[[i]]
    if (isfield(gi,'regressor')) 
    {
      r <- format.regressor(gi$regressor,i)
      pgi <- r$regressor
      regressor <- c(regressor, pgi)
      group[[i]]$regressor <- NULL
    }
  }
  if (!is.null(regressor))
    list.output$regressor=merge.regressor(regressor,N)
  
  #---  varlevel
  iv <- which(nv=="varlevel")
  if (!is.null(lv[[iv]]))
  {
    r <- format.varlevel(lv[[iv]],seq(1,N))
    if (length(r$id)>0) 
      id.ori=c(id.ori,r$id)
    varlevel <- r$varlevel
  }
  for (i in seq(1,N)){
    gi <- group[[i]]
    if (isfield(gi,'varlevel')) 
    {
      r <- format.varlevel(gi$varlevel,i)
      pgi <- r$varlevel
      varlevel <- c(varlevel, pgi)
      group[[i]]$varlevel <- NULL
    }
  }
  if (!is.null(varlevel))
    list.output$varlevel <- merge.varlevel(varlevel,N)
  
  #--------------------
  if (length(unique(id.ori))>0)
  {
    if (length(unique(id.ori))>1)
    {
      #stop("\n\nid's should be the same in the different input arguments (parameters, treatment, output,...)\n")
    }
    else
      list.output$id.ori <- unique(id.ori)#[[1]]
  }
  else
    list.output$id.ori <- NULL
  
  #--------------------
  list.output$group=group
  return(list.output)
}

#-----------------------------
format.treatment <- function(treatment,uN) 
{
  N <- length(uN)
  if (!is.null(names(treatment))) 
    treatment <- list(treatment) 
  id.ori <- NULL
  for (k in seq(1,length(treatment)))
  {
    trtk <- treatment[[k]]
    if (isfield(trtk,'colNames'))
    {
      pp <- as.data.frame(trtk$value)
      names(pp) <- trtk$colNames
      trtk=pp
    }
    else if (!is.data.frame(trtk))
    {
      trtk <- as.data.frame(trtk,stringsAsFactors =FALSE)
      n <- nrow(trtk)
      trtk <- trtk[rep(1:n,each=N),] 
      trtk$id <- rep(uN,n)
    }
    else
    {
      idk <- sort(unique(trtk$id))
      # if (any(idk != uN)) 
      #id.ori[[length(id.ori)+1]]<-idk
      # id.ori<-c(id.ori,idk)
      # trtk$id <- match(trtk$id, idk)
      
    }
    names(trtk)[names(trtk)=="amt"] <- "amount"
    names(trtk)[names(trtk)=="adm"] <- "type"
    
    if (!is.null(trtk$rate))
    {
      trtkrate=rep(Inf,length(trtk$rate))
      irate <- (trtk$rate!=".")
      trtkrate[irate]=as.numeric(as.character(trtk$rate[irate]))
      trtk$rate <- trtkrate
    }
    if (!is.null(trtk$tinf))
    {
      names(trtk)[names(trtk)=="tinf"] <- "rate"
      trtkrate=rep(Inf,length(trtk$rate))
      irate <- trtk$rate!="."
      trtkrate[irate]=trtk$amount[irate]/as.numeric(as.character(trtk$rate[irate]))
      trtk$rate <- trtkrate
    }
    if (is.null(trtk$rate))
      trtk$rate <- Inf
    treatment[[k]] <- trtk   
  }
  #r <- list(treatment=treatment, id=id.ori)
  r <- list(treatment=treatment, id=unique(id.ori))
  return(r)
}
#--------------------------
merge.treatment <- function(treatment,N)  
{
  p <- vector("list",N)
  ptr <- treatment[[1]]
  np <- length(treatment)
  if (np>1)
  {
    for (j in seq(2,np))
      ptr=rbind(ptr,treatment[[j]])
  }
  
  if (isfield(ptr,"target"))
  {
    ptr$target = as.factor(ptr$target)
    list.target <- levels(ptr$target)
    nt <- length(list.target)
    depot <- vector("list", nt)
    for (k in seq(1,nt))
      depot[[k]] <- list(type=k, target=list.target[k])
    ptr$type    <- as.numeric(ptr$target)
    ptr$target  <- NULL
  }
  else
  {
    depot<-NULL
  }
  
  for (i in seq(1,N))
  {
    pi <- NULL
    ij <- which(ptr$id==i)
    if (length(ij)>0) 
      pi <- ptr[ij,]
    pi$id=NULL
    
    if (!isfield(pi,"type"))
      pi$type <- 1
    
    if (!is.null(pi$time))
    {
      pi <- pi[with(pi,order(time)),]
      p[[i]] <- as.list(pi)
    }
  }
  r <- list(treatment=p, depot=depot)
  return(r)
}

#-----------------------------
format.parameter <- function(parameter,param,uN) 
{
  N <- length(uN)
  if (!is.null(names(param)))
    param=list(param) 
  param <- param[!unlist(lapply(param, is.null))]
  #id.ori <- list()
  id.ori <- NULL
  for (k in seq(1,length(param)))
  {
    paramk <- param[[k]]
    if (isfield(paramk,'colNames'))
    {
      pp=as.data.frame(paramk$value)
      names(pp)=paramk$colNames
      paramk=pp
    }
    id <- uN
    if (!is.data.frame(paramk))
    {
      if (!is.list(paramk))
        paramk <- list(name=names(paramk),value=as.vector(paramk))
      pk <- rep(paramk$value,N)
      pk <- t(matrix(pk,ncol=N))
      pk <- cbind(uN,pk)
      pk <- data.frame(pk)
      names(pk) <- c('id',paramk$name)
      paramk <- pk
    }
    else
    {
      if (!isfield(paramk,"id"))
      {
        paramk <- paramk[rep(1,each=N),] 
        paramk$id <- uN
      }
      else
      {
        idk <- sort(unique(paramk$id))
        if (any(idk != uN))
        {
          #id.ori[[length(id.ori)+1]]<-idk    
          id.ori<-c(id.ori,idk) 
          paramk$id <- match(paramk$id, idk)
        }
      }
    }
    for (i in uN)
    {
      pki <- paramk[paramk$id==i,]
      pki$id <- i
      if (is.null(parameter[[i]]))
        parameter[[i]] <- pki
      else
        parameter[[i]] <- merge(parameter[[i]],pki)
    }
  }
  # r <- list(parameter=parameter, id=id.ori)
  r <- list(parameter=parameter, id=unique(id.ori))
  return(r)
}
#-----------------------------
merge.parameter <- function(parameter) 
{
  if (!is.null(parameter[[1]]))
  {
    N <- length(parameter)
    pp <- parameter[[1]]
    if (N>1)
    {
      for (i in seq(2,N))
        pp <- rbind(pp, parameter[[i]])
    }
    parameter <- pp[with(pp,order(id)),]
    #    parameter$id <- NULL  
  }else{
    parameter=NULL
  }
  return(parameter)
}

#-----------------------------
format.output <- function(output, out,uN) 
{
  N <- length(uN)
  if (!is.null(names(out)))
    out=list(out) 
  
  ioutput <- output$individual
  goutput <- output$group
  id.ori <- NULL
  for (k in seq(1,length(out)))
  {
    outk <- out[[k]]
    if (is.data.frame(outk))
    {
      n.outk <- names(outk)
      outk <- list(name=n.outk[!(n.outk %in% c("id","time"))], 
                   time = outk[,c("id","time")])
    }
    if (!isfield(outk,"name"))
      outk <- list(name=outk)
    okname <- unlist(outk$name)
    nok <- length(okname)
    if (isfield(outk,'colNames'))
    {
      pp=as.data.frame(outk$value)
      names(pp)=outk$colNames
      outk$design=pp
    }
    if (!is.null(outk$time) && ("id" %in% names(outk$time)))
    {      
      # id <- sort(unique(as.factor(outk$time$id)))
      # id <- sort(unique(outk$time$id))
      # id <- output$id
      # if (length(id) !=length(uN) ||  any(id != uN))
      #   id.ori<-sort(unique(c(id.ori,id)))
      # idk <- match(id,levels(id))
      # idk <- match(sort(unique(outk$time$id)),id)
      idk <- match(sort(unique(outk$time$id)),uN)
      for (i in seq(1:length(idk)))
      {
        ji <- which(outk$time$id == idk[i])
        ti <- sort(outk$time$time[ji])
        oitime <- vector("list" , nok)
        for (j in seq(1,nok))
          oitime[[j]] <- ti
        ik <- idk[i]
        goutput[[ik]]$name=c(goutput[[ik]]$name,okname)
        ioutput[[ik]]$name=c(ioutput[[ik]]$name,okname)
        ioutput[[ik]]$time=c(ioutput[[ik]]$time,oitime)
      }
    }
    else
    { 
      for (i in uN)
      {
        if (isfield(outk,"time"))
        {
          oitime <- vector("list" , nok)
          for (j in seq(1,nok))
            oitime[[j]] <- sort(outk$time)
          ioutput[[i]]$time=c(ioutput[[i]]$time,oitime)
          ioutput[[i]]$name=c(ioutput[[i]]$name,okname)
        }
        goutput[[i]]$name=c(goutput[[i]]$name,okname)
      }
    }
  }
  output$individual <- ioutput 
  output$group <- goutput
  #r <- list(output=output, id=id.ori)
  r <- list(output=output, id=unique(id.ori))
  return(r)
}

#-----------------------------
format.regressor <- function(reg, uN) 
{
  N <- length(uN)
  if (!is.null(names(reg))) 
    reg=list(reg) 
  
  regressor <- vector("list",length(reg))
  #id.ori <- list()
  id.ori <- NULL
  for (k in seq(1,length(reg)))
  {
    regk <- reg[[k]]
    if (!is.data.frame(regk))
    {
      if (isfield(regk,'colNames'))
      {
        mk <- regk$value
        colNames=regk$colNames
      }
      else
      {
        nk <- length(regk$time)
        idk <- rep(uN,each=nk)
        mk <- cbind(regk$time,regk$value)
        mk <- mk[rep(1:nk,N),]
        mk <- cbind(idk, mk)
        colNames <- c("id","time",regk$name)
      }
      regk <- data.frame(mk)
      names(regk) <- colNames
    }
    else
    {
      #       colNames <- names(regk)
      #       #      mk <- data.matrix(regk)
      #       mk <- matrix(as.numeric(unlist(regk)),nrow=nrow(regk))
      #       mk <- as.data.frame(mk)
      #       names(mk) <- colNames
      mk <- regk
      if (is.null(mk$id))
        mk <- cbind(list(id=1),mk)
      # idk <- sort(unique(mk$id))
      # if (any(idk != uN)) 
      # {
      #id.ori[[length(id.ori)+1]]<-idk 
      # id.ori<-c(id.ori,idk) 
      # mk$id <- match(mk$id, idk)
      # }
      regk <- mk
    }
    # idu <- unique(regk$id)
    # regk$id <- match( regk$id , idu )
    regk <- regk[order(regk$id,regk$time),]
    regressor[[k]] <- regk
  }
  # r <- list(regressor=regressor, id=id.ori)
  r <- list(regressor=regressor, id=unique(id.ori))
  return(r)
}

#-----------------------------
merge.regressor <- function(reg, N) 
{
  if (length(reg)>1)
  {
    r <- NULL
    for (i in seq(1,N))
    {
      ri <- NULL
      for (k in seq(1,length(reg)))
      {
        rik <- reg[[k]][reg[[k]]$id==i,]
        if (nrow(rik)>0)
          if (is.null(ri))
            ri <- rik
          else
          {
            ri <- merge(ri,rik,
                        by.x=c("id","time"), by.y=c("id","time"),
                        all.x=TRUE,all.y=TRUE) 
            ri <- ri[order(ri$time),]
          }
      }
      r <- rbind(r,ri)
    }
  }
  else
    r <- reg[[1]]
  m <- matrix(as.numeric(unlist(r)),nrow=nrow(r))
  colNames <- names(r)
  colTypes <- rep("x",length(colNames))
  colTypes[colNames=="id"] <- "id"
  colTypes[colNames=="time"] <- "time"
  regressor <- list(colNames=colNames, colTypes=colTypes, value=m)
  return(regressor)
}

#-----------------------------
format.varlevel <- function(var, uN) 
{
  N <- length(uN)
  if (!is.null(names(var))) 
    var=list(var) 
  
  varlevel <- vector("list",length(var))
  id.ori <- NULL
  for (k in seq(1,length(var))){
    vark <- var[[k]]
    if (!is.data.frame(vark)){
      if (isfield(vark,'colNames')){
        mk <- vark$value
        colNames=vark$colNames
      }else{
        nk <- length(vark$time)
        idk <- rep(uN,each=nk)
        occk <- rep(1:nk,N)
        mk <- vark$time[occk]
        mk <- cbind(idk, mk, occk)
        colNames <- c("id","time",vark$name)
      }
      vark <- data.frame(mk)
      names(vark) <- colNames
    }
    # else
    #   {
    #   idk <- sort(unique(vark$id))
    #   if (any(idk != uN)) {
    #     id.ori[[length(id.ori)+1]]<-idk 
    #     vark$id <- match(vark$id, idk)
    #   }
    # }
    varlevel[[k]] <- vark
  }
  r <- list(varlevel=varlevel, id=id.ori)
  return(r)
}

#-----------------------------
merge.varlevel <- function(var, N) 
{
  if (length(var)>1)
  {
    r <- NULL
    for (i in seq(1,N))
    {
      ri <- NULL
      for (k in seq(1,length(var)))
      {
        rik <- var[[k]][var[[k]]$id==i,]
        if (nrow(rik)>0)
          if (is.null(ri))
          {
            ri <- rik
          }
        else
        {
          ri <- merge(ri,rik,
                      by.x=c("id","time"), by.y=c("id","time"),
                      all.x=TRUE,all.y=TRUE) 
          for (j in seq(1,ncol(ri)))
          {
            ij <- c(which(!is.na(ri[,j])),nrow(ri)+1)
            
            for (l in seq(1,length(ij)-1))
            {
              ijl <- ij[l]
              if (ij[l+1]-ijl>1)
                ri[seq(ijl+1,ij[l+1]-1),j] <- ri[ijl,j]
            }
          }
        }
      }
      r <- rbind(r,ri)
    }
  }else{
    r <- var[[1]]
  }
  
  m <- matrix(as.numeric(unlist(r)),nrow=nrow(r))
  #m <- data.matrix(r)
  colNames <- names(r)
  colTypes <- rep("OCC",length(colNames))
  colTypes[colNames=="id"] <- "id"
  colTypes[colNames=="time"] <- "time"
  varlevel <- list(colNames=colNames, colTypes=colTypes, value=m)
  return(varlevel)
}

#-----------------------------
parse.model <- function(model) 
{
  con        = file(model, open = "r")
  lines      = readLines(con, warn=FALSE)
  close(con)
  ip <- grep(";",lines, fixed=TRUE)
  if (length(ip)>0)
  {
    for (k in (1:length(ip)))
    {
      ll <- lines[ip[k]]
      il <- regexpr(";",ll)
      il <- il[1]
      if (il>1)
        lines[ip[k]] <- substr(ll,start=1,stop=(il-1))
      else
        lines[ip[k]] <- ''
    }
  }
  
  i1 <- grep("[POPULATION]",lines, fixed=TRUE)
  i2 <- grep("[COVARIATE]",lines, fixed=TRUE)
  i3 <- grep("[INDIVIDUAL]",lines, fixed=TRUE)
  i4 <- grep("[LONGITUDINAL]",lines, fixed=TRUE)
  
  level <- NULL
  if (length(i1)>0)
    level <- c(level,'population')
  if (length(i2)>0)
    level <- c(level,'covariate')
  if (length(i3)>0)
    level <- c(level,'individual')
  if (length(i4)>0)
    level <- c(level,'longitudinal')
  
  fs <- fsort(c(i1, i2, i3, i4))
  model.info <- list(level=level[fs[[2]]])  
  return(model.info)
  
}  
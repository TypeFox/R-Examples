simulR <- function(argList)
{
  model <- argList$DATA$model
  if (!exists(model, mode="function"))
    source(model)
  model <- basename(file_path_sans_ext(model))
  if (!is.null(argList$SETTINGS) && !is.null(argList$SETTINGS$seed))
    set.seed(argList$SETTINGS$seed)
  group <- argList$DATA$group
  parameter <- argList$DATA$individual_parameters
  
  if (!is.null(argList$DATA$depot)){
    iop.depot <- TRUE
    depot <- argList$DATA$depot
    n.depot <- length(depot)
    target <- vector(length=n.depot)
    for (k in (1:n.depot))
      target[depot[[k]]$type] <- depot[[k]]$target
  }else
    iop.depot <- FALSE
  
  K <- length(group)
  i <- 0
  r <- list()
  for (k in (1:K)){
    pk <- parameter$value[k,]
    names(pk)=parameter$name
#      pk <- as.list(pk)
    group[[k]]$level=NULL
    trk <- argList$DATA$individual[[k]]$dose
    if (iop.depot){
      trk <- list(data=data.frame(time=trk$time, value=trk$amount, 
                             method="add", var=target[trk$type]))
    }
    
   ok <- argList$DATA$group[[k]]$output
    respk <- argList$DATA$individual[[k]]$response$name
    rtk <- argList$DATA$individual[[k]]$response$time
    outk <- list()
    jparam <- NULL
    for (j in 1:length(ok)){
      outk[[j]] <- list(name=ok[j])
      jj <- which(ok[j]==respk)
      if (length(jj)>0)
        outk[[j]]$time <- rtk[[jj]]
      else
        jparam <- c(jparam,j)
    }
    if (!is.null(argList$DATA$regressor)){
      reg <- argList$DATA$regressor
      jid <- which(reg$colTypes=="id")
      ik <- which(reg$value[,jid]==k)
      regk <- data.frame(reg$value[ik,])
      names(regk) <- reg$colNames
      regk$id <- NULL
      reg$colTypes <- reg$colTypes[-jid]
      reg$colNames <- reg$colNames[-jid]
      jx <- which(reg$colTypes=="x")
      jt <- which(reg$colTypes=="time")
      rk <- list()
      for (j in (1:length(jx))){
        ix <- which(!is.nan(regk[,jx[j]]))
        rk[[j]] <- list(name=reg$colNames[jx[j]],
                        time=regk[ix,jt],
                        value=regk[ix,jx[j]])
      }
    }else{
      rk <- NULL
    }
    Nk <- group[[k]]$size
    dki=data.frame(group=k)
    ri <- NULL
    if (!is.null(trk)){
      if (!is.null(rk)){
        ttt<- paste0('ri=',model,'(parameter=pk, dose=trk, time=rtk, regressor=rk)')
      }else{
        ttt<- paste0('ri=',model,'(parameter=pk, dose=trk, time=rtk)')
      }
    }else{
      if (!is.null(rk)){
        ttt<- paste0('ri=',model,'(parameter=pk, time=rtk, regressor=rk)')
      }else{
        ttt<- paste0('ri=',model,'(parameter=pk, time=rtk)')
      }        
    } 
    for (ik in (1:Nk)){
      i <- i+1
      eval(parse(text=ttt))
      dki$id=i
      for (j in (1:length(ok))){
        rij <- merge(dki,ri[[j]])
        if (i==1){
          r[[j]] <- rij
        }else{
          r[[j]] <- rbind(r[[j]],rij)
        }
      }
    }
  }
  names(r) <- ok
  if (length(jparam)>0){
    for (j in (1:length(jparam)))
      names(r[[jparam[j]]])[3] <- ok[jparam[j]]
  }
  if (i==1){
    for (j in (1:length(ok)))
      r[[j]]$id=NULL
  }else{
    r[[j]]$id <- as.factor(r[[j]]$id)
  }
    
  if (K==1){
    for (j in (1:length(ok)))
      r[[j]]$group=NULL
  }else{
    r[[j]]$group <- as.factor(r[[j]]$group)
  }
if (length(jparam)>1){
    p <- r[[jparam[1]]]
    for (j in (2:length(jparam[j])))
      p <- merge(p,r[[jparam[j]]])
    r$parameter <- p
  }
  return(r)
}
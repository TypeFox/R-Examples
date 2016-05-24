#' @export
arfimaPrep <-
  function(data, timevar, varlist.mean, varlist.fd, varlist.xdif, varlist.ydif
           , d="Hurst", arma=NULL, ecmformula=NULL, decm="Hurst", drop=5, ...){
    # warnings / error messages
    if(missing(timevar)) stop("timevar must be specified!")
    if(length(unique(diff(sort(unique(data[,timevar])))))>1) warning("Time series is not evenly spaced!")
    
    # step 1: aggregate .dif+.fd variables over timevar
    data.mean <- NULL
    if (length(varlist.mean)>0){
      data.mean <- aggregate(data[,varlist.mean],list(data[,timevar]), mean, na.rm=T)
      if (length(varlist.mean)==1) names(data.mean)[2] <- varlist.mean
    }
    
    # step 2: fractionally difference .fd variables (see fd function)
    data.fd=NULL; res.d=NULL
    if (length(varlist.fd)>0){
      data.fd <- data.frame(cbind(data.mean[,1],data.mean[,varlist.fd]))
      if (length(varlist.fd)==1) names(data.fd)[2] <- varlist.fd
      
      if(is.character(d)) {
        for(i in 1:length(varlist.fd)){
          tmp <- fd(data.mean[,varlist.fd[i]],dval=d)
          data.fd[,varlist.fd[i]] <- tmp$series
          res.d <- rbind(res.d,c(varlist.fd[i],tmp$estimator,tmp$value))
          rm(tmp)
        }
      } else if(length(d)==length(varlist.fd)){
        if (length(setdiff(varlist.fd,names(d)))!=0) {
          stop("List of d-parameters does not fit to model specification!")
        } else for(i in 1:length(varlist.fd)){
          tmp <- fd(data.mean[,varlist.fd[i]],dval=d[[grep(varlist.fd[i], names(d),fixed=T)]])
          data.fd[,varlist.fd[i]] <- tmp$series
          res.d <- rbind(res.d,c(varlist.fd[i],tmp$estimator,tmp$value))
          rm(tmp)
        }
      } else stop("List of d-parameters does not have correct length!")
    }
    
    # step 3: estimate AR/MA parameter for selected variables
    res.arma <- NULL
    if(!is.null(arma)){
      res.arma <- vector("list",length=length(arma))
      if(length(intersect(varlist.fd,names(arma)))!=length(arma)){
        stop("List for AR/MA estimation does not fit to model specification!")
      }
      for(i in 1:length(arma)){
          if(class(arma[[i]])=="numeric"){
              tmp <- arima(data.fd[,names(arma)[i]]
                           , order=c(arma[[names(arma)[i]]][1],0,arma[[names(arma)[i]]][2])
                           , include.mean = FALSE)
          }
          if(class(arma[[i]])=="list"){
            if (max(arma[[i]][[1]]) == length(arma[[i]][[1]]) & 
	    max(arma[[i]][[2]]) == length(arma[[i]][[2]])) {
                tmp <- arima(data.fd[,names(arma)[i]]
                           , order=c(max(arma[[i]][[1]]),0,max(arma[[i]][[2]]))
                           , include.mean = FALSE)
	    } else {
                fixed <- rep(0, max(arma[[i]][[1]]) + max(arma[[i]][[2]]))
                fixed[c(arma[[i]][[1]],max(arma[[i]][[1]]) + arma[[i]][[2]])] <- NA
                tmp <- arima(data.fd[,names(arma)[i]]
                           , order=c(max(arma[[i]][[1]]),0,max(arma[[i]][[2]]))
                           , include.mean = FALSE
                           , fixed = fixed, transform.pars = FALSE)	
            }
          }
	data.fd[,names(arma[i])] <- tmp$residuals
        res.arma[[i]] <- tmp
        names(res.arma)[i] <- names(arma[i])
        rm(tmp)
      }
    }
    
    # step 4: adjust varnames
    if (!is.null(data.mean)){
      names(data.mean) <- paste0(names(data.mean),".mean")
      names(data.mean)[1] <- timevar
    }
    if (!is.null(data.fd)){
      names(data.fd) <- paste0(names(data.fd),".fd")
      names(data.fd)[1] <- timevar
    }
    
    # step 5: calculate ecm
    res.ecm=NULL
    if (!is.null(ecmformula)){
      res.ecm <- lm(ecmformula,data=data.mean)
      tmp  <- fd(residuals(res.ecm), decm)
      ecm <- tmp$series
      ecm <- c(NA,ecm[1:(length(ecm)-1)])
      data.fd <- cbind(data.fd,ecm)
      res.d <- rbind(res.d,c("ecm",tmp$estimator,tmp$value))
      rm(tmp)
    }
    
    # step 6: merge datasets
    if (!is.null(data.fd)){
      data.merged <- merge(data.mean, data.fd, by=timevar)
      data.merged <- merge(data, data.merged, by=timevar)
    }
    else if (!is.null(data.mean)){
      data.merged <- merge(data, data.mean, by=timevar)
    }
    else data.merged <- data
    
    # step 7: calculate dif value for dv
    if (length(varlist.ydif)>0){
      for (i in 1:length(varlist.ydif)){
        ydif <- data.merged[,varlist.ydif[i]]-(data.merged[,grep(paste0(varlist.ydif[i],".mean")
                                                               ,names(data.merged))]
                                             -data.merged[,grep(paste0(varlist.ydif[i],".fd")
                                                                ,names(data.merged))])
        data.merged <- cbind(data.merged, ydif)
      }
      names(data.merged)[(ncol(data.merged)
                         -length(varlist.ydif)+1):ncol(data.merged)] <- paste0(varlist.ydif,".ydif")
    }
    
    # step 8: calculate dif values for ivs
    if (length(varlist.xdif)>0){
      for (i in 1:length(varlist.xdif)){
        xdif <- data.merged[,varlist.xdif[i]]-data.merged[,grep(paste0(varlist.xdif[i],".mean")
                                                              ,names(data.merged))]
        data.merged <- cbind(data.merged, xdif)
     }
      names(data.merged)[(ncol(data.merged)
                          -length(varlist.xdif)+1):ncol(data.merged)] <- paste0(varlist.xdif,".xdif")
    }
    
    # step 9: drop certain number of initial observations
    data.merged <- subset(data.merged
                          , data.merged[,timevar]>(min(data.merged[,timevar])+drop-1))
    
    # return data
    if (!is.null(res.d)){
      rownames(res.d) <- res.d[,1]
      res.d <- data.frame(res.d[,2],as.numeric(res.d[,3]))
      colnames(res.d) <- c("Method","Estimate")
    }
    out <- list(data.mean = data.mean, data.fd = data.fd, data.merged = data.merged
                , d = res.d, arma = res.arma, ecm = res.ecm
                )
  }

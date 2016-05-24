getOutputTable = function(object) {
	what = object$outputWhat
	pred = object$predictions
	nCols = length(what)
	blockWhat = object$blockWhat
  bCols = length(names(blockWhat))
  nsim = ifelse("nsim" %in% names(what),what$nsim,0)
  out = matrix(NA, length(pred[[1]]), nCols+bCols+max(nsim-1,0))
	ns = c("x", "y")
	if (dim(coordinates(pred))[2] == 3) ns = c("x", "y", "z")
  ntype = NULL
  if (all(c("var1.pred","var1.var") %in% names(pred))) ntype = "gstat" 
  if (length(grep("sim",names(pred)))>0) {
    ntype = c(ntype,"simul")
    sims = pred@data[,grep("sim",names(pred@data))]
  }
	for (i in 1:nCols) {
		if (names(what)[i] == "mean") {
		  if ("mean" %in% names(pred)) {
        out[,i] = pred$mean
  		} else if ("var1.pred" %in% names(pred)) {
        out[,i] = pred$var1.pred
  		} else if ("simul" %in% ntype) {
        out[,i] = apply(sims,MARGIN=1,FUN=function(arr) mean(arr))
  		} else {
  		  stop(paste("Warning: Could not create variable mean, returning NA"))
  		  out[,i] = NA
      }
      ns = c(ns, "mean")
  	} else if (names(what)[i] == "variance") {
		  if ("variance" %in% names(pred)) {
        out[,i] = pred$variance
  		} else if ("var1.var" %in% names(pred)) {
        out[,i] = pred$var1.var
  		} else if ("simul" %in% ntype) {
        out[,i] = apply(sims,MARGIN=1,FUN=function(arr) var(arr))
  		} else {
  		  stop(paste("Warning: Could not create variable variance, returning NA"))
  		  out[,i] = NA
      }
  		ns = c(ns, "variance")
  	} else if (names(what)[i] == "quantile") {
      qname = paste("quantile",what[[i]],sep="")
		  if (length(grep("quantile",names(pred)))>0) {
        out[,i] = pred[[qname]]
  		} else if ("gstat" %in% ntype) {
        out[,i] = qnorm(what[[i]], pred$var1.pred, sqrt(pred$var1.var))
  		} else if ("simul" %in% ntype) {
        out[,i] = apply(sims,MARGIN=1,FUN=function(arr) quantile(arr,what[[i]]))
  		} else {
  		  stop(paste("Warning: Could not create variable quantile = ", what[i],", returning NA"))
  		  out[,i] = NA
      }
  		ns = c(ns, qname)
  	} else if (names(what)[i] == "cumdistr") {
      cname = paste("cumdistr",what[[i]],sep="")
		  if (length(grep("cumdistr",names(pred)))>0) {
        out[,i] = pred[[cname]]
  		} else if ("gstat" %in% ntype) {
        out[,i] = pnorm(what[[i]], pred$var1.pred, sqrt(pred$var1.var))
  		} else if ("simul" %in% ntype) {
        out[,i] = apply(sims,MARGIN=1,FUN=function(arr) sum(I(arr) < what[[i]])/length(arr))
  		} else {
  		  stop(paste("Warning: Could not create variable cumdistr = ", what[i],", returning NA"))
  		  out[,i] = NA
      }
  		ns = c(ns, cname)
  	} else if (names(what)[i] == "excprob") {
      ename = paste("excprob",what[[i]],sep="")
		  if (length(grep("excprob",names(pred)))>0) {
        out[,i] = pred[[ename]]
  		} else if ("gstat" %in% ntype) {
        out[,i] = 1 - pnorm(what[[i]], pred$var1.pred, sqrt(pred$var1.var))
  		} else if ("simul" %in% ntype) {
        out[,i] = apply(sims,MARGIN=1,FUN=function(arr) sum(I(arr) > what[[i]])/length(arr))
  		} else {
  		  stop(paste("Warning: Could not create variable excprob = ", what[i],", returning NA"))
  		  out[,i] = NA
      }
 	 		ns = c(ns, ename)
  	} else if (names(what)[i] == "MOK") {
      mname = paste("MOK",what[[i]],sep="")
		  if (length(grep("MOK",names(pred)))>0) {
        out[,i] = pred[[mname]]
  		} else {
  		  stop(paste("Warning: Could not create variable MOK = ", what[i],", returning NA"))
  		  out[,i] = NA
      }
 	 		ns = c(ns, mname) 	 		
  	} else if (names(what)[i] == "IWQSEL") {
      iname = paste("IWQSEL",what[[i]],sep="")
		  if (length(grep("IWQSEL",names(pred)))>0) {
        out[,i] = pred[[iname]]
  		} else {
  		  stop(paste("Warning: Could not create variable IWQSEL = ", what[i],", returning NA"))
  		  out[,i] = NA
      }
 	 		ns = c(ns, iname) 	
  	} else if (names(what)[i] == "nsim") {
      if ("simul" %in% ntype && dim(sims)[2] >= nsim) {
        out[,i:(i+nsim-1)] = as.matrix(sims[,1:nsim])         		
      } else {
        stop(paste("Warning: Number of simulations is zero or smaller than nsim = ",nsim))
      }
      ns = c(ns,names(sims))
  	} else
	 		stop(paste("unknown request: ", names(what)[i]))
	}
# Block Predictions:
	i = 0
  for (j in numeric(length(names(blockWhat)))) {  
	  i = i+1
    ib = nCols+i
  	if (names(blockWhat)[i] == "fat") {
      vname = paste("fat",blockWhat[[i]],sep="")
		  if (length(grep("fat",names(pred)))>0) {
		    out[,ib] = pred[[vname]]
  		} else {
  		  stop(paste("Warning: Could not create variable fat, returning NA"))
  		  out[,ib] = NA
      }
      ns = c(ns, vname)
    } else if (names(blockWhat)[i] == "fatVar") {
      vname = paste("fatVar",blockWhat[[i]],sep="")
		  if (length(grep("fatVar",names(pred)))>0) {
		    out[,ib] = pred[[vname]]
  		} else {
  		  stop(paste("Warning: Could not create variable fatVar, returning NA"))
  		  out[,ib] = NA
      }
      ns = c(ns, vname)
    } else if (names(blockWhat)[i] == "blockMax") {
      if ("blockMax" %in% names(pred)) {
		    out[,ib] = pred$blockMax
  		} else {
  		  stop(paste("Warning: Could not create variable blockMax, returning NA"))
  		  out[,ib] = NA
      }
      ns = c(ns, "blockMax")
    } else if (names(blockWhat)[i] == "blockMaxVar") {
      if ("blockMaxVar" %in% names(pred)) {
		    out[,ib] = pred$blockMaxVar
  		} else {
  		  stop(paste("Warning: Could not create variable blockMaxVar, returning NA"))
  		  out[,ib] = NA
      }
      ns = c(ns, "blockMaxVar")
    } else if (names(blockWhat)[i] == "blockMin") {
      if ("blockMin" %in% names(pred)) {
        out[,ib] = pred$blockMin
  		} else {
  		  stop(paste("Warning: Could not create variable blockMin, returning NA"))
  		  out[,ib] = NA
      }
      ns = c(ns, "blockMin")
    } else if (names(blockWhat)[i] == "blockMinVar") {
      vname = paste("blockMinVar",blockWhat[[i]],sep="")
		  if (length(grep("blockMinVar",names(pred)))>0) {
		    out[,ib] = pred$blockMinVar
  		} else {
  		  stop(paste("Warning: Could not create variable blockMinVar, returning NA"))
  		  out[,ib] = NA
      }
      ns = c(ns, "blockMinVar")
  	} else if (names(blockWhat)[i] != "none") {
	 		stop(paste("unknown request: ", names(blockWhat)[i]))
    }
	}

  out = signif(out,digits = 6)
 	ret = cbind(coordinates(object$predictions), out)
  colnames(ret) = ns
	ret
}

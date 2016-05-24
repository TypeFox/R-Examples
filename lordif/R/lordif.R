lordif <-
function(resp.data,group,selection=NULL,criterion=c("Chisqr","R2","Beta"),pseudo.R2=c("McFadden","Nagelkerke","CoxSnell"),alpha=0.01,beta.change=0.1,R2.change=0.02,
           maxIter=10,minCell=5,minTheta=-4.0,maxTheta=4.0,inc=0.1,control=list(),model="GRM",anchor=NULL,MonteCarlo=FALSE,nr=100,weights=NULL,normwt=TRUE) {
    call<-match.call()
    criterion<-match.arg(criterion)
    pseudo.R2<-match.arg(pseudo.R2)
    tni<-ncol(resp.data)
    if (!(criterion %in% c("Chisqr","R2","Beta"))) {
      warning("criterion must be one of the following: \"Chisqr\", \"R2\", or \"Beta\"; will be reset to \"Chisqr\"")
      criterion<-"Chisqr"
    }
    if (!(pseudo.R2 %in% c("McFadden","Nagelkerke","CoxSnell"))) {
      warning("pseudo.R2 must be one of the following \"McFadden\", \"Nagelkerke\", or \"CoxSnell\"; will be reset to \"McFadden\"")
      pseudo.R2<-"McFadden"
    }
    if (alpha<=0 || alpha>1) {
      warning ("alpha must be > 0 & < 1; will be reset to .01")
      alpha<-.01
    }
    if (beta.change<=0 || beta.change>=1) {
      warning ("beta.change must be > 0 & < 1; will be reset to .1")
      beta.change<-.1
    }
    if (R2.change<=0 || R2.change>=1) {
      warning("R2.change must be > 0 & < 1; will be reset to .02")
      R2.change<-.02
    }
    if (maxIter<1) {
      warning("maxIter must be >= 1; will be reset to 10")
      maxIter<-10
    }
    if (minCell<1) {
      warning("minCell must be >= 1; will be reset to 5")
      minCell<-5
    }
    if (minTheta>=maxTheta) {
      warning("minTheta must be < maxTheta; will be reset to default")
      minTheta<--4;maxTheta<-4
    }
    if (inc<=0) {
      warning("inc must be > 0; will be reset to .1")
      inc<-.1
    }
    if (nrow(resp.data) != length(group)) stop("nrow of resp.data and length of group are non-conformable")
    if (length(selection)==0) selection<-1:tni
    else {
      selection <- unique(selection)
      if (!all(selection %in% 1:tni)) {
        warning("selection is not in the full set; all items will be selected")
        selection <- 1:tni
      }
    }
    if (!(model %in% c("GRM","GPCM"))) {
      warning("model must be either \"GRM\", or \"GPCM\"; will be reset to \"GRM\"")
      model<-"GRM"
    }
    if(length(anchor) > 0) {
      anchor<-unique(anchor)
      if (!all(anchor %in% selection)) {
        warning("bad anchor items; no anchor items will be used")
        anchor<-NULL
      } 
    }
    ni<-length(selection)
    good<-!(is.na(group) | rowSums(!is.na(resp.data[,selection]))==0)
    resp.data<-resp.data[good,]
    group<-group[good]
    nobs<-nrow(resp.data)
    if (!is.null(weights)) {
      weights<-weights[good]
      if (nobs != length(weights)) stop("nrow of resp.data and length of weights are non-conformable")
      if (any(weights<0)) {
        warning("some/all elements of weights are negative; will be reset to 0")
        weights[weights<0]<-0
        if (all(weights==0)) {
          warning("all elements of weights are 0; will be reset to NULL")
          weights<-NULL
        }
      }
      if (sum(weights)>0 && normwt) weights <- weights * nobs / sum(weights)
    }
    if (nobs != length(group)) stop("nrow of resp.data and length of group are non-conformable") 
    if (ni<4) stop("number of items must be at least 4") 
    compare<-function(x,table) {
      for (i in 1:nrow(table)) {
        if (all(table[i,]==x)) return(TRUE)
      }
      return(FALSE)
    }
    resp.recoded<-data.frame(matrix(NA,nobs,ni)) 
    names(resp.recoded)<-paste0("I",selection) 
    for (i in 1:ni) {
      resp.recoded[,i]<-collapse(resp.data[,selection[i]],group,minCell) 
    }
    ncat<-as.numeric(apply(resp.recoded,2,max,na.rm=T))
    ng<-length(table(group))
    theta.grid<-seq(minTheta,maxTheta,inc)
    item.labels<-names(resp.recoded)
    meanraw<-apply(resp.recoded,1,mean,na.rm=T)
    outraw<-rundif(selection,resp.recoded,meanraw,group,criterion,alpha,beta.change,pseudo.R2,R2.change,weights)
    calib<-mirt(resp.recoded,1,itemtype=ifelse(model=="GPCM","gpcm","graded"),technical=control)
    ipar<-extract(calib)
    row.names(ipar)<-item.labels
    theta<-calctheta(ipar,resp.recoded,theta.grid,model=model)
    if (length(anchor)==0) {
      out<-rundif(selection,resp.recoded,theta$theta,group,criterion,alpha,beta.change,pseudo.R2,R2.change,weights)
      flag.matrix<-rbind(logical(ni),!logical(ni))
      flags<-out$flag 
      iter<-1
      ndif<-sum(flags)
      cat(paste0(" (mirt) | Iteration: ",iter,", ",ndif," items flagged for DIF (",paste0(selection[flags],collapse=","),")\n"))
      if (ndif==ni) {
        warning("all items got flagged for DIF - stopping\n")
        theta.sparse = NULL
        ipar.sparse = NULL
      } else if (ndif==0) {
        warning("no items got flagged for DIF - stopping\n")
        theta.sparse = NULL
        ipar.sparse = NULL
      }
      if (!compare(flags,flag.matrix)) {
        repeat {
          iter<-iter+1
          flag.matrix<-rbind(flag.matrix,flags)
          sparse.matrix<-separate(resp.recoded,flags,group)
          calib.sparse<-mirt(sparse.matrix,1,itemtype=ifelse(model=="GPCM","gpcm","graded"),technical=control)
          ipar.sparse<-extract(calib.sparse) 
          eqconst<-equate(ipar[!flags,],ipar.sparse[1:sum(!flags),],theta.grid,model=model)
          ipar.sparse[,1]<-ipar.sparse[,1]/eqconst[1]
          ipar.sparse[,2:ncol(ipar.sparse)]<-ipar.sparse[,2:ncol(ipar.sparse)]*eqconst[1]+eqconst[2]
          theta.sparse<-calctheta(ipar.sparse,sparse.matrix,theta.grid,model=model)
          pre.flags<-flags
          out<-rundif(selection,resp.recoded,theta.sparse$theta,group,criterion,alpha,beta.change,pseudo.R2,R2.change,weights)
          flags<-out$flag
          ndif<-sum(flags)
          cat(paste0(" (mirt) | Iteration: ",iter,", ",ndif," items flagged for DIF (",paste0(selection[flags],collapse=","),")\n"))
          if (ndif==ni) {
            warning("all items got flagged for DIF - stopping\n")
            break
          }
          if (compare(flags,flag.matrix) || iter==maxIter) {
            if (!all(pre.flags==flags)) {
              sparse.matrix<-separate(resp.recoded,flags,group)
              calib.sparse<-mirt(sparse.matrix,1,itemtype=ifelse(model=="GPCM","gpcm","graded"),technical=control)
              ipar.sparse<-extract(calib.sparse) 
              eqconst<-equate(ipar[!flags,],ipar.sparse[1:sum(!flags),],theta.grid,model=model)
              ipar.sparse[,1]<-ipar.sparse[,1]/eqconst[1]
              ipar.sparse[,2:ncol(ipar.sparse)]<-ipar.sparse[,2:ncol(ipar.sparse)]*eqconst[1]+eqconst[2]
              theta.sparse<-calctheta(ipar.sparse,sparse.matrix,theta.grid,model=model)
            }
            break
          }
        }
        if (!compare(flags,flag.matrix) && iter==maxIter) {
          sparse.matrix<-separate(resp.recoded,flags,group)
          calib.sparse<-mirt(sparse.matrix,1,itemtype=ifelse(model=="GPCM","gpcm","graded"),technical=control)
          ipar.sparse<-extract(calib.sparse)
          eqconst<-equate(ipar[!flags,],ipar.sparse[1:sum(!flags),],theta.grid,model=model) 
          ipar.sparse[,1]<-ipar.sparse[,1]/eqconst[1]
          ipar.sparse[,2:ncol(ipar.sparse)]<-ipar.sparse[,2:ncol(ipar.sparse)]*eqconst[1]+eqconst[2]
          theta.sparse<-calctheta(ipar.sparse,sparse.matrix,theta.grid,model=model) 
        }
        row.names(ipar.sparse)<-names(sparse.matrix)
      }
    } else {
      iter<-0
      flags<-rep(TRUE,ni)
      flags[anchor]<-FALSE
      sparse.matrix<-separate(resp.recoded,flags,group)
      cat(" (mirt)\n")
      calib.sparse<-mirt(sparse.matrix,1,itemtype=ifelse(model=="GPCM","gpcm","graded"),technical=control)
      ipar.sparse<-extract(calib.sparse)
      eqconst<-equate(ipar[!flags,],ipar.sparse[1:sum(!flags),],theta.grid,model=model) 
      ipar.sparse[,1]<-ipar.sparse[,1]/eqconst[1]
      ipar.sparse[,2:ncol(ipar.sparse)]<-ipar.sparse[,2:ncol(ipar.sparse)]*eqconst[1]+eqconst[2]
      theta.sparse<-calctheta(ipar.sparse,sparse.matrix,theta.grid,model=model)
      out<-rundif(selection,resp.recoded,theta.sparse$theta,group,criterion,alpha,beta.change,pseudo.R2,R2.change,weights)
      flags<-out$flag
      ndif<-sum(flags)
      cat(" (mirt)\n")
      cat("anchor items: ",paste0(anchor,collapse=","),"\n")
      cat(paste0(ndif," items flagged for DIF (",paste0(selection[flags],collapse=","),")\n"))
      row.names(ipar.sparse)<-names(sparse.matrix)
    }
    output<-list(call=call,options=list(model=model,criterion=criterion,pseudo.R2=pseudo.R2,alpha=alpha,beta.change=beta.change,R2.change=R2.change,
                                        maxIter=maxIter,minCell=minCell,minTheta=minTheta,maxTheta=maxTheta,inc=inc,control=control),selection=selection,stats=out$stats,
                 flag=out$flag,recoded=resp.recoded,group=group,ng=ng,ni=ni,ncat=ncat,calib=theta,calib.sparse=theta.sparse,weights=weights,
                 iteration=iter,ipar=ipar,ipar.sparse=ipar.sparse,stats.raw=outraw$stats,meanraw=meanraw,flag.raw=outraw$flag,DFIT=NULL,anchor=anchor,MonteCarlo=NULL)
    class(output)<-"lordif"
    if (MonteCarlo) {
      cat("Monte Carlo simulation\n")
      MC<-montecarlo(output,alpha=alpha,nr=nr)
      cat("\n")
      if (toupper(output$options$criterion)=="CHISQR") {
        flag<-output$stats$chi12 < MC$cutoff$chi12 | output$stats$chi13 < MC$cutoff$chi13 | output$stats$chi23 < MC$cutoff$chi23
      } else if (toupper(output$options$criterion)=="BETA") {
        flag<-output$stats$beta12 > MC$cutoff$beta12
      } else if (toupper(output$options$criterion)=="R2") {
        if (toupper(output$options$pseudo.R2)=="MCFADDEN") {
          flag<-output$stats$pseudo13.McFadden > MC$cutoff$pseudo13.McFadden 
        } else if (toupper(output$options$pseudo.R2)=="NAGELKERKE") {
          flag<-output$stats$pseudo13.Nagelkerke > MC$cutoff$pseudo13.Nagelkerke
        } else if (toupper(output$options$pseudo.R2)=="COXSNELL") {
          flag<-output$stats$pseudo13.CoxSnell > MC$cutoff$pseudo13.CoxSnell
        }
      }
      cat(paste0(sum(flag)," items flagged for DIF by Monte Carlo thresholds (",paste0(selection[flag],collapse=","),")\n"))
      output[["MonteCarlo"]]<-MC
    }
    return(output)
  }

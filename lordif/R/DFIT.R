DFIT <-
function(obj) {
    if (class(obj)!="lordif") stop(paste(deparse(substitute(obj))," must be of class lordif"))
    options<-obj$options
    model<-options$model
    control<-options$control
    ncat<-obj$ncat
    flags<-obj$flag
    group<-obj$group
    resp.recoded<-obj$recoded
    selection<-obj$selection
    ni<-obj$ni
    group.names<-names(table(group))
    ng<-obj$ng 
    group.size<-table(group)
    theta.grid<-seq(obj$options$minTheta,obj$options$maxTheta,obj$options$inc)
    item.labels<-names(resp.recoded)
    NCDIF<-data.frame(matrix(NA,ni,ng-1))
    names(NCDIF)<-paste("NCDIF",group.names[-1],sep=".")
    row.names(NCDIF)<-item.labels
    CDIF<-data.frame(matrix(NA,ni,ng-1))
    names(CDIF)<-paste("CDIF",group.names[-1],sep=".")
    row.names(CDIF)<-item.labels
    DTF<-numeric(ng-1)
    cat("DFIT Analysis\n")
    resp.recoded.group<-split(resp.recoded,group)
    calib.group<-list()
    tcc.group<-list()
    tcc.group[["theta"]]<-theta.grid
    for (g in 1:ng) {
      cat(paste("Group:",group.names[g]),"\n")
      calib.group[group.names[g]]<-mirt(resp.recoded.group[[group.names[g]]],1,itemtype=ifelse(model=="GPCM","gpcm","graded"),technical=control)
      cat(" (mirt)\n")
    }
    ipar.group<-lapply(calib.group,extract)  
    for (g in 2:ng) {
      eqconst<-equate(ipar.group[[1]][!flags,],ipar.group[[g]][!flags,],theta.grid,model=model) 
      ipar.group[[g]][,1]<-ipar.group[[g]][,1]/eqconst[1]
      ipar.group[[g]][,2:ncol(ipar.group[[g]])]<-ipar.group[[g]][,2:ncol(ipar.group[[g]])]*eqconst[1]+eqconst[2]
    }
    for (g in 1:ng) {
      tcc.group[[group.names[g]]]<-as.vector(tcc(ipar.group[[g]][,1],ipar.group[[g]][,2:ncol(ipar.group[[g]])],theta.grid,model=model))
    }
    for (g in 2:ng) {
      theta.group<-calctheta(ipar.group[[g]],resp.recoded.group[[g]],theta.grid,model=model)$theta
      pp.R<-calcprob(ipar.group[[1]],theta.group,model=model)
      pp.F<-calcprob(ipar.group[[g]],theta.group,model=model)
      di<-matrix(NA,group.size[g],ni)
      for (i in 1:ni) {
        di[,i]<-rowSums(pp.F[,i,1:ncat[i]]*(col(pp.F[,i,1:ncat[i]])-1))-rowSums(pp.R[,i,1:ncat[i]]*(col(pp.R[,i,1:ncat[i]])-1))
      }
      Di<-rowSums(di)
      mean.D<-mean(Di)
      mean.di<-colMeans(di)
      for (i in 1:ni) {
        CDIF[i,g-1]<-cov(di[,i],Di)+mean.di[i]*mean.D
        NCDIF[i,g-1]<-var(di[,i])+mean.di[i]^2
      }
      DTF[g-1]<-sum(CDIF[,g-1])
      cat(paste0("DTF (",group.names[g],") = ", round(DTF[g-1],digits=4),"\n"))
    }
    out<-obj
    out[["DFIT"]]<-list(CDIF=CDIF,NCDIF=NCDIF,DTF=DTF,ipar=ipar.group,TCC=tcc.group)
    class(out)<-"lordif"
    return(out)
  }

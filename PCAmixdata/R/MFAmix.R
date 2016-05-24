MFAmix<-function(data, groups, name.groups, ndim=5, rename.level=FALSE, 
                 graph=TRUE, axes=c(1,2))
{
  
  cl <- match.call()
  
  if(length(groups)!=ncol(data))
    stop("\"groups\" must be a vector of size the number of variables in \"data\"")
  
  #test that name.groups doesn't contain special characters
  ch<-paste(name.groups,collapse="")
  ch<-strsplit(ch, split="")
  ch<-tolower(unique(unlist(ch)))
  res.ch<-unique(is.element(ch,c(letters,0,1,2,3,4,5,6,7,8,9,"_")))
  if(length(res.ch)==2)
    stop("In \"name.groups\" spaces and special characters are not allowed.")
  
  n<-nrow(data)
  
  nbr.groups<-length(unique(groups))
  
  if (length(name.groups)!=nbr.groups)
    stop("invalid length of \"name.groups\"")
  
  Lst.groups<-splitgroups(base=data,groups=groups,name.groups=name.groups) 
  long.groups<-sapply(Lst.groups,ncol)
  typ.groups<-unlist(sapply(Lst.groups,splitmix)[3,])
  
  DATA.ord<-data.frame(matrix(NA,ncol=ncol(data),nrow=nrow(data)))
  init<-0
  for(g in 1:nbr.groups)
  {
    DATA.ord[,c((1+init):(init+ncol(Lst.groups[[g]])))]<-Lst.groups[[g]]
    init<-init+ncol(Lst.groups[[g]])
  }
  colnames(DATA.ord)<-unlist(sapply(Lst.groups,colnames))
  rownames(DATA.ord)<-rownames(data)
  groups.ord<-NULL
  for(g in 1:nbr.groups){
    groups.ord<-c(groups.ord,rep(g,ncol(Lst.groups[[g]])))
  }
  groups<-groups.ord
  data<-DATA.ord
  
  ordre.var.data<-data.frame(colnames(data),seq(from=1,by=1,to=ncol(data)))
  colnames(ordre.var.data)<-c("nom.var","ordre.var")
  
  tab.indic.names<-cbind(colnames(data),rep(name.groups,long.groups),rep(typ.groups,long.groups),rep(long.groups,long.groups))
  rownames(tab.indic.names)<-NULL
  tab.indic.names<-data.frame(tab.indic.names)
  colnames(tab.indic.names)<-c("var","groups","type","nb.var")
  
  Res.separe.pcamix<-list()
  Lst.groups.stand<-list()
  
  for(i in 1:nbr.groups){
    base.qt<-splitmix(Lst.groups[[i]])$X.quanti
    base.ql<-splitmix(Lst.groups[[i]])$X.quali
    Res.separe.pcamix[[i]]<-PCAmix(X.quanti=base.qt, X.quali=base.ql, ndim=ndim, rename.level=rename.level, graph=F)
    Lst.groups.stand[[i]]<-Res.separe.pcamix[[i]]$Z
  }
  names(Res.separe.pcamix)<-names(Lst.groups.stand)<-name.groups
  
  eig.groups<-list()
  for(i in 1:nbr.groups){
    eig.groups[[i]]<-Res.separe.pcamix[[i]]$eig[1]
  }
  eig.groups<-unlist(eig.groups)
  names(eig.groups)<-paste(name.groups)
  
  ponde.qt<-NULL
  base.qt<-splitmix(data)$X.quanti
  indic.qt.groups<-tab.indic.names[match(colnames(base.qt),tab.indic.names$var),2]
  ponde.qt<-eig.groups[match(indic.qt.groups,names(eig.groups))]
  
  ponde.ql<-NULL
  base.ql<-splitmix(data)$X.quali
  if (is.null(base.ql)==FALSE){
    indic.ql.groups<-tab.indic.names[match(colnames(base.ql),tab.indic.names$var),2]
    ponde.ql<-eig.groups[match(indic.ql.groups,names(eig.groups))]
    ponde.ql<-rep(ponde.ql,apply(base.ql,2,nb.level))
  }
  
  ponderation<-c(ponde.qt,ponde.ql)
  ponderation<-1/ponderation
  Res.total<-PCAmix(X.quanti=base.qt,X.quali=base.ql,ndim=ndim, rename.level=rename.level, graph=FALSE,weight.col=ponderation)
  
  sqload.order<-data.frame(matrix(NA,ncol=ncol(Res.total$sqload),nrow=nrow(Res.total$sqload)))
  for (i in 1:nrow(Res.total$sqload)){
    index<-as.character(ordre.var.data[i,1])
    sqload.order[i,]<-Res.total$sqload[index,]
  }
  rownames(sqload.order)<-ordre.var.data[,1]
  colnames(sqload.order)<-colnames(Res.total$sqload)
  Res.total$sqload<-sqload.order
  
  ponderation.groups<-ponderation[match(name.groups,names(ponderation))]
  ndim<-ncol(Res.total$ind$coord)
  
  #partial individuals
  data.partiel <- vector(mode = "list", length = nbr.groups)
  names(data.partiel) <- name.groups
  
  for (g in 1:nbr.groups){
    col.interet<-rownames(Res.separe.pcamix[[g]]$V)
    data.partiel[[g]]<-data.frame(Res.total$W[,col.interet])
    colnames(data.partiel[[g]])<-col.interet
  }
  
  V<-Res.total$V
  M<-Res.total$M
  
  res.ind.partiel <- vector(mode = "list", length = nbr.groups)
  names(res.ind.partiel) <- name.groups
  for (g in 1:nbr.groups) {
    coord.ind.sup<-nbr.groups*data.partiel[[g]]
    coord.ind.sup<-sweep(coord.ind.sup,2,STATS=M[colnames(coord.ind.sup)],FUN="*")
    coord.ind.sup<-as.matrix(coord.ind.sup)%*%V[colnames(coord.ind.sup),]
    res.ind.partiel[[g]]$coord.sup <- coord.ind.sup
  }
  
  ndim.max.groups<-NULL
  for (g in 1:nbr.groups){
    ndim.max.groups<-c(ndim.max.groups,ncol(res.ind.partiel[[g]]$coord.sup))
  }
  
  nom.ligne <- NULL
  for (i in 1:n) {
    ind.tmp <- rownames(data)[i]
    nom.ligne <- c(nom.ligne, paste(ind.tmp, name.groups, sep = "."))
  }
  coord.ind.partiel <- as.data.frame(matrix(NA, (n *nbr.groups), ndim))
  rownames(coord.ind.partiel) <- nom.ligne
  colnames(coord.ind.partiel) <- paste("Dim", c(1:ndim), sep = ".")  
  liste.ligne <- seq(1, n* nbr.groups, by = nbr.groups)  
  for (g in 1:nbr.groups){
    coord.ind.partiel[liste.ligne +  g - 1, ] <- res.ind.partiel[[g]]$coord.sup[,1:ndim.max.groups[g]]
  }
  
  #inertia
  Inertie.tot <- vector(length = ndim)
  for (g in 1:nbr.groups){
    Inertie.tot <- Inertie.tot + apply(res.ind.partiel[[g]]$coord.sup^2 * n, 2, sum)}
  
  rap.inertie <- apply(Res.total$ind$coord^2 * n, 2, sum) * nbr.groups/Inertie.tot
  
  #partial axes
  Res.partial.axes.groups<-function(nom.groups){
    score.sep<-eval(parse(text=paste("Res.separe.pcamix$",nom.groups,"$ind$coord",sep="")))
    colnames(score.sep)<-paste(colnames(score.sep),nom.groups,sep=".")
    score.global<-Res.total$ind$coord
    cor(score.sep,score.global)  
  }
  
  partial.axes.coord<-NULL
  for (g in 1:nbr.groups) {
    a<-Res.partial.axes.groups(name.groups[g])
    partial.axes.coord<-rbind(partial.axes.coord,a)
  }
  partial.axes.coord
  
  #variables contributions                                  
  contrib.total<-rbind(Res.total$quanti$contrib,Res.total$quali$contrib)
  contrib.total<-contrib.total[match(as.vector(tab.indic.names$var),rownames(contrib.total)),]
  
  #groups
  contrib.groups<-data.frame(matrix(NA,nrow=nbr.groups,ncol=ndim))
  a<-0
  for (i in 1:nbr.groups){
    if (is.vector(contrib.total[(a+1):(a+long.groups[i]),]))
    {contrib.groups[i,]<-contrib.total[(a+1):(a+long.groups[i]),]} else
    { contrib.groups[i,]<-apply(contrib.total[(a+1):(a+long.groups[i]),],2,sum)
    }
    a<-a+long.groups[i]
  }
  rownames(contrib.groups)<-name.groups
  colnames(contrib.groups)<-colnames(partial.axes.coord)
  
  
  coord.groups<-contrib.groups
  contrib.groups.pct<-sweep(coord.groups,2,STATS=Res.total$eig[1:ndim,1],FUN="/")*100
  
  dist.groups<-NULL
  for (i in 1:nbr.groups){
    valP<-Res.separe.pcamix[[i]]$eig[,1]
    dis<-valP/valP[1]
    dis<-sum(dis^2)
    dist.groups<-c(dist.groups,dis)
  }
  names(dist.groups)<-name.groups
  
  cos2.groups<-sweep(coord.groups^2,1,STATS=dist.groups,FUN="/")
  
  
  Lg.groups<-Lg.pond(Lst.groups.stand,ponderation.groups)
  RV.groups<-RV.pond(Lst.groups.stand,ponderation.groups)
  
  res.groups<-list(Lg=Lg.groups,RV=RV.groups,coord=coord.groups,contrib.pct=contrib.groups.pct,dist2=dist.groups,cos2=cos2.groups)
  res.partial.axes<-list(coord=partial.axes.coord)  
  Res.total$ind$coord.partial<-coord.ind.partiel
  
  #separated analysis
  Recap.eig<-list()
  for (i in 1:nbr.groups){
    if (nrow(Res.separe.pcamix[[i]]$eig)<ndim){
      Recap.eig[[i]]<-c(Res.separe.pcamix[[i]]$eig[,1],rep(NA,ndim-nrow(Res.separe.pcamix[[i]]$eig)))
    } else       {
      Recap.eig[[i]]<-c(Res.separe.pcamix[[i]]$eig[1:ndim,1])
    }
  }
  
  Recap.eig.frame<-matrix(NA,ncol=nbr.groups,nrow=ndim)
  for (i in 1:nbr.groups){
    Recap.eig.frame[,i]<-Recap.eig[[i]]
  }
  colnames(Recap.eig.frame)<-name.groups
  rownames(Recap.eig.frame)<- paste("dim", 1:ndim, sep =" ")
  
  res<-list(call=cl,
            eig=Res.total$eig,
            eig.separate=Recap.eig.frame,
            separate.analyses=Res.separe.pcamix,
            groups=res.groups,
            partial.axes=res.partial.axes,
            inertia.ratio=rap.inertie,
            ind=Res.total$ind,
            quanti=Res.total$quanti,        
            levels=Res.total$levels,
            quali=Res.total$quali,
            coef=Res.total$coef,
            sqload=Res.total$sqload,
            global.pca=Res.total,
            ind.partial=res.ind.partiel,
            lst.groups=groups
  )
  
  
  class(res)<-c("MFAmix","list")
  
  if (graph) {  
    plot.MFAmix(res,axes=axes,choice="axes",coloring.var="groups")
    plot.MFAmix(res,axes=axes,choice="groups",coloring.var="groups")
    plot.MFAmix(res,axes=axes,choice="ind",cex=0.8)
    plot.MFAmix(res,axes=axes,choice="sqload", coloring.var="groups")
    
    if (!is.null(Res.total$quanti$coord)){
      plot.MFAmix(res,axes=axes,choice="cor", coloring.var="groups")
    }
    if (!is.null(Res.total$levels$coord)){
      plot.MFAmix(res,axes=axes,choice="levels",coloring.var="groups")
    }
  }
  
  return(res)
}

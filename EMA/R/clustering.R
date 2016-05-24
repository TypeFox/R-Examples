#########
##
## Clustering Plots
##
########

clustering.plot <- function(tree, tree.sup, data=NULL, lab, lab.sup, dendro=TRUE, dendro.sup=TRUE, title="", scale="row", heatcol, names=TRUE, names.sup=TRUE, names.dist=TRUE, trim.heatmap=1, palette="rainbow", legend=TRUE, legend.pos="topright", ...){

  if (missing(tree)){
    stop("At least, one tree is required.")
  }

  if ((!missing(tree.sup) && missing(data))){
    stop("Two-ways clustering need two trees and data matrix.")
  }
  ##check that label exist
  if(! "order.lab" %in% names(tree))
      stop("order.lab missing in tree object. Check you data colnames.") 

  if (!missing(lab) && is.vector(lab)){
    lab <- matrix(lab, ncol=1)
  }
  if (!missing(lab) && is.data.frame(lab)){
      for (i in 1:ncol(lab)) lab[, i] <- as.character(lab[, i])
      lab <- as.matrix(lab)
  }
  
  if (!missing(lab.sup) && is.vector(lab.sup)){
    lab.sup <- matrix(lab.sup, ncol=1)
  }
  if (!missing(lab.sup) && is.data.frame(lab.sup)){
      for (i in 1:ncol(lab.sup)) lab.sup[, i] <- as.character(lab.sup[, i])
      lab.sup <- as.matrix(lab.sup)
  }
  
  op <- par(no.readonly = TRUE)
  #op<-par(mar=c(4,4,3,3)) 
  on.exit(par(op))

  ## Heatmap display
  if (!is.null(data)){
    if (!is.matrix(data)){
      data <- as.matrix(data)
    }

    #par(mar=c(5,4,3,9))
    mean.ori<-mean(data, na.rm=TRUE)
    
    ##Always center (to -1/1 scaling)
    data <- data-mean(data, na.rm=TRUE)
    
    ##Heatmap
    if (scale=="row"){
      data=t(scale(t(data)))
      mean.ori<-0
      if (tree$method=="euclidean"){
        warnings("We do not advice to scale the data, when you use an euclidean distance ...")
      }
    }
    else if (scale=="column"){
      data=scale(data)
      mean.ori<-0
    }
    
    ##Trimming
    q <- quantile(data, c((1-trim.heatmap), trim.heatmap), na.rm=TRUE)
    data[data < q[1]] = q[1]
    data[data > q[2]] = q[2]
    
    ##-1/1
    maxi <- max(data, na.rm=TRUE)
    mini <- min(data, na.rm=TRUE)
    data[!is.na(data) & data > 0] <- data[!is.na(data) & data > 0]/maxi
    data[!is.na(data) & data < 0] <- -data[!is.na(data) & data < 0]/mini
    maxi.real <- q[2]+mean.ori
    mini.real <- q[1]+mean.ori

    ##RowSideColors
    if (!missing(lab.sup)){
      rowcolorlab<-t(as.colors(t(lab.sup), palette=palette))
      if (ncol(lab.sup)==1){
        rowcolorlab<-cbind(rowcolorlab,rowcolorlab)
      }
      rownames(rowcolorlab)<-NULL
    }
    ##ColSideColors
    if (!missing(lab)){
      colcolorlab<-as.colors(lab, palette=palette)
      if (ncol(lab)==1){
        colcolorlab<-cbind(colcolorlab,colcolorlab)
      }
      colnames(colcolorlab)<-NULL
    }

    ##Dendrograms
    d.t <- d.ts <-NA
    if (dendro){
      d.t <- as.dendrogram(as.hclust(tree))
    }
    if (dendro.sup & !missing(tree.sup)){
      d.ts <- as.dendrogram(as.hclust(tree.sup))
    }

    ##Margins
    bmar=floor(max(nchar(colnames(data))/2))
    if (names.sup)
        rmar=floor(max(nchar(rownames(data))/2))
    else
        rmar=4

    ##Text & legend
    n <- n.sup <-  NULL
    if (!names){
      n <- rep("", length(tree$order))
    }
    if (!names.sup){
        if (!missing(tree.sup)){
            n.sup <- rep("", length(tree.sup$order))
        }
        else{
            n.sup <- rep("", nrow(data)) 
        }
    }
   
    if (names.dist){
      xtitle <- paste("Hierarchical Clustering :",tree$method,attr(tree$diss,"Metric"))
      bmar <- bmar+3
      if (dendro.sup & !missing(tree.sup)){
        ytitle <- paste("Hierarchical Clustering :",tree.sup$method,attr(tree.sup$diss,"Metric"), sep=" ")
        rmar <- rmar+3
      }else{
        ytitle <- NA
      }
    }
    else{
      xtitle <- NA
      ytitle <- NA
    }
    if (missing(heatcol)) 
      heatcol <- myPalette(low = "green", high = "red", mid ="black", k=50)
    
    ##Heatmap
    if (!missing(lab) && !missing(lab.sup)){
        heatmap.plus(data, Rowv=d.ts, Colv=d.t, margins=c(bmar,rmar), cexCol=0.9, cexRow=0.9, scale="none", cex.main=0.7, main=title, ylab=ytitle, xlab=xtitle, RowSideColors=rowcolorlab,ColSideColors=colcolorlab, labRow=n.sup, labCol=n, col=heatcol,...)
    }
    else if (!missing(lab.sup)){
        heatmap.plus(data, Rowv=d.ts, Colv=d.t, margins=c(bmar,rmar), cexCol=0.9, cexRow=0.9, scale="none",  cex.main=0.7, main=title, ylab=ytitle, xlab=xtitle, RowSideColors=rowcolorlab, labRow=n.sup, labCol=n, col=heatcol,...)
    }
    else if (!missing(lab)){
      heatmap.plus(data, Rowv=d.ts, Colv=d.t, margins=c(bmar,rmar), cexCol=0.9, cexRow=0.9, scale="none", cex.main=0.7,  main=title, ylab=ytitle, xlab=xtitle, ColSideColors=colcolorlab, labRow=n.sup, labCol=n, col=heatcol,...)
    }
    else{
        heatmap.plus(data, Rowv=d.ts, Colv=d.t,  scale="none", cexCol=0.9, cexRow=0.9, cex.main=0.7, main=title,  ylab=ytitle, xlab=xtitle, labRow=n.sup, labCol=n, col=heatcol, margins=c(bmar,rmar),...)
    }

    ##Label legend
    if (legend && (!missing(lab) || !missing(lab.sup))){
      op1 <- par(fig = c(0.05,0.2,0.78,0.99), new=TRUE)
      op2 <- par(mar=c(0,0,0,0))
      lall<-c()
      if (!missing(lab)){
         lall<-c(lall, as.vector(lab))
      }
      if (!missing(lab.sup)){
         lall<-c(lall, as.vector(lab.sup))
      }
      leg <- lall
      if (is.element(NA,leg)){
          leg[which(is.na(leg))] <- "NA"
      }
         
      legend(legend.pos,legend=unique(leg), fill=unique(as.colors(lall, palette=palette)), bty="n", cex=0.7, text.col="gray50")
      par(op1)
      par(op2)
    }
    ##Heatmap legend
    if (legend){
      op1 <- par(fig = c(0.93,0.95,0.8,0.95), new=TRUE)
      op2 <- par(mar=c(0,0,0,0))
      iml<-matrix(seq(from=mini.real, to=maxi.real, length.out=10), nrow=1)
      image(x=1, y=1:10, z=iml, axes=FALSE, ylab="", xlab="", col=heatcol, ...)
      if (trim.heatmap != 1){
	axis(side=4, labels=c(paste("<",round(mini.real,1)),round((mini.real + maxi.real)/2,1),paste(">",round(maxi.real,1))), at=c(1,5.5,10), lwd=0.5, las=1, cex.axis=0.7, tick=FALSE,line=-0.7, ...)
      }else{
	axis(side=4, labels=c(round(mini.real,1),round((mini.real + maxi.real)/2,1),round(maxi.real,1)), at=c(1,5.5,10), lwd=0.5, las=1, cex.axis=0.7, tick=FALSE,line=-0.7, ...)
      }
      par(op1)
      par(op2)
    }
  }
  
  ##One tree
  else{

      if (names){
          bmar=floor(max(nchar(tree$order.lab))/2)
      }else{
          bmar=0
      }
    par(mar=c(bmar,4,2,1), font.lab=2) ##3,2,1
    if(!missing(lab)){
        if (names){
            pl=0.1+0.1*ceiling(ncol(lab)/3)
        }else{
            pl=0.08
        }
        layout(matrix(c(1:2), 2, 1, byrow=TRUE), heights=c(1-pl,pl))#, respect=TRUE)
    }
    else{
	layout(matrix(c(1:2), 2, 1, byrow=TRUE), heights=c(0.9,0.1))#, respect=TRUE)
    }
    plot(as.dendrogram(as.hclust(tree)), edgePar=list(col=c('blue', 'black')), leaflab="none", dLeaf=NULL,horiz=FALSE, frame.plot=FALSE, ylab=paste(tree$method,"linkage"), cex.lab=0.8, main=title, cex.main=0.8, ...)
    parmar=par("mar")
    paru=par("usr")

  if (legend && !missing(lab)){
      leg <- lab
      if (is.element(NA,leg)){
          leg[which(is.na(leg))] <- "NA"
      }
      
      legend(legend.pos,legend=unique(as.vector(leg)), fill=unique(as.colors(as.vector(lab), palette=palette)), bty="n", cex=0.7, text.col="gray50")
    }

    if (names){
      mtext(text=tree$order.lab, at=1:length(tree$order), side=1, line=1, col="black", las=2, cex=0.65, font=1, ...)
    }

    if (names.dist)
        mtext(paste("Hierarchical Clustering :",tree$method,attr(tree$diss,"Metric"), sep=" "), at=length(tree$order)/2, side=3, las=1, col="gray50", cex=0.8)
    
    if (!missing(lab)){
      lab.num <- lab
      ind <- 0
      for (i in unique(as.vector(lab))){
        if (is.na(i)){
            lab.num[which(is.na(lab))]<-ind
        }
        else{
            lab.num[which(lab==i)]<-ind
        }
        ind <- ind+1
      }
      im<-matrix(as.numeric(lab.num),ncol=ncol(lab), nrow=nrow(lab))
      im <- as.matrix(im[tree$order,])
      
      par(mar=c(2,parmar[2],parmar[3],parmar[4]))
      #par(mar=c(2,parmar[2],0,parmar[4]))

      image(x=1:length(tree$order), y=1:ncol(im),xlim=c(paru[1],paru[2]), z=im,axes=FALSE, ylab="", xlab="", col=unique(as.colors(as.vector(lab),palette=palette)),...)
  
      if (ncol(lab)>1){
        axis(side=2, labels=colnames(lab), at=1:ncol(im), lwd=0.5, las=1, cex.axis=0.7, tick=FALSE, ...)
        #text(x=rep(-1,ncol(im)),y=1:ncol(im), labels=colnames(lab), cex=.7, pos=4)
      }
      else{
        axis(side=4, labels=FALSE, at=1:ncol(im), lwd=0.5, las=1, cex.axis=0.7, tick=FALSE, ...)
      }

    }
  }
 }

#########
##
## Clustering Functions
##
########

clustering <-  function(data, metric="euclidean", method="ward", nb) {

  METHOD <- c("average", "single", "complete", "ward", "weighted","diana","kcentroids")
  mea <- pmatch(method, METHOD)
  
  if (is.na(mea)) stop("Error : Unknown Linkage.")
  if (!is.na(mea) && mea != 7)  DIS <- clust.dist(data, metric)

  if (mea<=5) HIERA <- agnes(DIS, method=method, diss=TRUE, keep.diss=TRUE)
  else if (mea==6) HIERA <- diana(DIS,diss=TRUE)
  else if (mea==7) {
    if (missing(nb)) stop("Number of clusters 'k' have to be selected for kcentroids algorithm")
    HIERA<-list()
    if (metric=="euclidean") {
      for (i in nb) HIERA <- c(HIERA, km=kmeans(t(data), centers=i)$cluster )
    }
    else{
      DIS <- clust.dist(data, metric)
      for (i in nb) HIERA <- c(HIERA, pm=pam(DIS$DIS, k=i, diss=TRUE)$cluster )
    }
  }
 
  attr(HIERA$diss,"Metric")<-attr(DIS,"Metric")
  
  return(HIERA)             
}

clustering.kmeans <-function(data, N=100, iter.max=20, title="Kmeans - Hierarchical Clustering", dist.s="pearson", dist.g="pearsonabs", method="ward") {

  print ("Running kmeans ...")
  km <- kmeans(data, centers=N, iter.max=iter.max)
  data2cluster<-km$centers

  print ("Running clustering ...")
  print (paste(dist.s,method))
  c.sample <- clustering(data2cluster, metric=dist.s, method=method)
  print (paste(dist.g, method))
  c.gene <- clustering(t(data2cluster), metric=dist.g, method=method)
  clustering.plot(tree=c.sample, tree.sup=c.gene, data=data2cluster, title=title)
  
  return(list(c.km=km, c.sample=c.sample, c.kcenters=c.gene))
}



#########
##
## Clustering Distances
##
########

clust.dist <- function(mat, meth.dis="euclidean") { 
  
  MEASURE <- c("euclidean", "manhattan", "pearson", "pearsonabs", "spearman", "spearmanabs", "jaccard", "dice")
  mea <- pmatch(meth.dis, MEASURE)
  
  if (is.na(mea)) stop("Error :Unknown Metric.")
  if (mea==1) DIS <- as.matrix(daisy(t(mat), metric="euclidean"))      
  if (mea==2) DIS <- as.matrix(daisy(t(mat), metric="manhattan"))    
  if (mea==3) DIS <- (1-cor(mat, method="pearson"))/2
  if (mea==4) DIS <- 1-abs(cor(mat, method="pearson"))
  if (mea==5) DIS <- (1-cor(mat, method="spearman"))/2
  if (mea==6) DIS <- 1-abs(cor(mat, method="spearman"))
  if (mea==7) DIS <- jaccard(mat)                                    
  if (mea==8) DIS <- dice(mat)                                        
  attr(DIS,"Metric")<-meth.dis

  return(DIS)
}

jaccard <- function(mat){
  out <- matrix(NA, ncol=ncol(mat), nrow=ncol(mat))
  colnames(out) <- colnames(mat)
  rownames(out) <- colnames(mat)

  for (i in 1:(ncol(mat))){
    for (j in 1:(ncol(mat))){
      a <- mat[,i]
      b <- mat[,j]
      out[i,j] <- (sum(b, na.rm=TRUE) + sum(a, na.rm=TRUE) -2*sum(a & b, na.rm=TRUE))/(sum(a, na.rm=TRUE) + sum(b, na.rm=TRUE) - sum(a & b, na.rm=TRUE))
    }
  }
  return (out)
}

dice <- function(mat) {
  out <- matrix(NA, ncol=ncol(mat), nrow=ncol(mat))
  colnames(out) <- colnames(mat)
  rownames(out) <- colnames(mat)

  for (i in 1:(ncol(mat))){
    for (j in 1:(ncol(mat))){
      a <- mat[,i]
      b <- mat[,j]
      out[i,j] <- (sum(abs(b), na.rm=TRUE) + sum(abs(a), na.rm=TRUE) -2*sum(a & b, na.rm=TRUE))/(sum(abs(a), na.rm=TRUE) + sum(abs(b), na.rm=TRUE))
    }
  }
  return (out)
}



#########
##
## Clustering Stability
##
########

eval.stability.clustering<-  function(X,nb=c(2:4),f=0.8,nsub=10,s0=0.98, list_DIS=c("euclidean","pearson"), list_ALGO=c("average","complete","ward"), pdfname = NULL, verbose = TRUE){
 
  Transform.vector.to.list <-function (v) {
    if (is.integer(v) == FALSE)
      stop("Transform.vector.to.list: the elements of the input vector must be integers", call.=FALSE);
    ## relabeling of the classes in order to avoid (that is the labels of the classes will be consecutive integers)
    n.examples <- length(v);
    new.v <- integer(n.examples); ## vector with relabeled elements
    max.class <- max(v);
    v.index <- integer(max.class); ## vector of the class indices for label translation
    label.class <- 0;
    
    for (i in 1:n.examples) {
      old.label <- v[i]; # old label of the ith example
      if (v.index[old.label] == 0) {
        label.class <- label.class + 1;
        v.index[old.label] <- label.class;
      }
      new.v[i] <- v.index[old.label];	
    }
    ## building of the list of clusters
    cl =list();
    for (i in 1:n.examples) {
      if (length(cl) < new.v[i]) 
        cl[[new.v[i]]] <- i
      else  
        cl[[new.v[i]]][length(cl[[new.v[i]]]) + 1] <- i;
    }
    return(cl);
  }
  
  Do.boolean.membership.matrix <-function(cl, dim.M, examplelabels) {
    M <- matrix(integer(dim.M*dim.M), nrow=dim.M);
    colnames(M) <- rownames(M) <- examplelabels;
    singletons <- integer(dim.M);  
    c <- length(cl); # number of clusters 
    for (j in 1:c) {
      n.ex <- length(cl[[j]]);
      if (n.ex == 1)
        singletons[cl[[j]][1]] <- 1
      else {
        for (x1 in 1:(n.ex-1)) {
          for (x2 in (x1+1):n.ex) {
            x <- cl[[j]][x1];
            y <- cl[[j]][x2];
            M[x,y] <- 1;
          }
        }
      }
    }
    for (x1 in 1:(dim.M-1)) 
      for (x2 in (x1+1):dim.M) 
        M[x2,x1] <- M[x1,x2];
    for (x in 1:(dim.M)) 
      M[x,x] <- singletons[x];
    return(M);
  }
  
  
  Compute.Chi.sq<-function (M, s0){
    n <- ncol(M)
    K <- nrow(M)
    x <- numeric(K)
    for (k in 1:K) {
      for (j in 1:n) {
        if (is.na(M[k, j])==FALSE & M[k, j] > s0) x[k] <- x[k] + 1
      }
    }
    theta <- sum(x)/(n * K)
    
    if (theta == 0) {
      p.value = NA
    } else if (theta == 1) {
      p.value = 1
    } else {
      chi.statistic <- sum((x - n * theta)^2)/(n * theta * (1 - theta))
      p.value <- 1 - pchisq(chi.statistic, K - 1)
    }
    return(p.value)
  }
  
  
  teste.stab<-function (sim.matrix, s0 ){
    n.clusterings <- nrow(sim.matrix)
    n.measures <- ncol(sim.matrix)
    ordered.clusterings <- integer(n.clusterings)
    p.value <- numeric(n.clusterings)
    means <- numeric(n.clusterings)
    variance <- numeric(n.clusterings)
    
    means <- mean(as.data.frame(t(sim.matrix)),na.rm=TRUE)
    for (i in 1:n.clusterings) variance[i] <- var(sim.matrix[i,],na.rm=TRUE)
    
    sorted.means <- sort(means, decreasing = TRUE)
    sorted.indices <- order(means, decreasing = TRUE)
    means <- sorted.means
    variance <- variance[sorted.indices]
    ordered.sim.matrix <- sim.matrix[sorted.indices, ]
    
    ordered.clusterings <- sorted.indices 
    p.value[1] <- 1
    for (k in n.clusterings:2)
      p.value[k] <- Compute.Chi.sq(ordered.sim.matrix[1:k,], s0)
    
    d <- data.frame(ordered.clusterings = ordered.clusterings, p.value = p.value, means = means, variance = variance)
    rownames(d) <- 1:n.clusterings
    
    return(d)
  }
  
  t0=Sys.time()
    
  VALID.STAB<-matrix(0,0,4)
  
  for (j in list_DIS){
    for (q in list_ALGO){
      
      ## STABILITY
      ## with  partition   
      n <- ncol(X)
      n.sub.ex <- ceiling(n * f)
      
      for (t in 1:nsub) {
        Jsim.vector <- rep(NA,length(nb))         
        sub1 <- sample(n, n.sub.ex)
        sub2 <- sample(n, n.sub.ex)
        Xsub1 <- X[, sub1]
        colnames(Xsub1) <- sub1
        Xsub2 <- X[, sub2]
        colnames(Xsub2) <- sub2
        
        if (q=="kcentroids"){
          S1 <- clustering(Xsub1,j,q, nb=nb)
          S2 <- clustering(Xsub2,j,q, nb=nb)
        }
        else{
          S1 <- clustering(Xsub1,j,q)
          S2 <- clustering(Xsub2,j,q)
        }
     
        ##Works on 2 partitions 80%
        for (l in nb){
          if (q!="kcentroids") {
            t1 <- cutree(as.hclust(S1), k=l)
            cl1<-Transform.vector.to.list(t1)
            t2 <- cutree(as.hclust(S2), k=l)
            cl2<-Transform.vector.to.list(t2)
            
            M1<-Do.boolean.membership.matrix(cl1, n.sub.ex, sub1)
            M2<-Do.boolean.membership.matrix(cl2, n.sub.ex, sub2)
          }  
          else{
            cl1<-Transform.vector.to.list(S1[[l-1]])   
            M1<-Do.boolean.membership.matrix(cl1, n.sub.ex, sub1)
            cl2<-Transform.vector.to.list(S2[[l-1]])   
            M2<-Do.boolean.membership.matrix(cl2, n.sub.ex, sub2)
          }
          
          sub.common <- intersect(sub1, sub2)
          label.examples <- as.character(sub.common)
          M1 <- M1[label.examples, label.examples]
          M2 <- M2[label.examples, label.examples]
          
          ## calculate Jaccard coefficient
          Jsim.vector[l-1] <-sum(M1 * M2)/(sum(M1 * M1) + sum(M2 * M2) - sum(M1 * M2))
        }    
        
        STAB<-data.frame(dist=rep(j,length(nb)), algo=rep(q,length(nb)), nb, res=rep(t,length(nb)), Jsim.vector=Jsim.vector)
        VALID.STAB<-rbind(VALID.STAB,STAB)
      } # t resampling
    } # q algo
  } # j linkage
  
  t1=Sys.time()
  tdiff=t1-t0

  ##================================================
  ##                CHI-2 TEST - Which one is the most stable ?
  ##================================================
  
  NAMES   <-list()
  RESULTS <-list()
  NBSET   <-c()
  
  for (i in nb) {
    VALID.STABsub<-VALID.STAB[which(VALID.STAB$nb==i),]
    
    sim.matrix   <- matrix(VALID.STABsub$Jsim.vector, ncol = nsub, byrow = TRUE )
    names.matrix <- matrix(paste(paste(VALID.STABsub$nom,VALID.STABsub$dist),VALID.STABsub$algo),ncol = nsub, byrow = TRUE )
    
    D<-teste.stab(sim.matrix, s0 ) 
    NAMES   <-c(NAMES,list(names.matrix[,1]))
    RESULTS <-c(RESULTS,list(D))
    NBSET   <-c(NBSET,paste("nb.clust=",i))
  }
  
  for( i in 1:length(RESULTS)) {
    corr.p.value<-p.adjust(RESULTS[[i]]$p.value,method ="holm")
    RESULTS[[i]]<-cbind(RESULTS[[i]], corr.p.value)
  }
  
  RES <-list()
  NAM <-list()
  for( i in 1:length(RESULTS)) {
    res<- RESULTS[[i]]
    nam<- NAMES[[i]]
    res1<-res[which(res$corr.p.value>0.05),]
    RES <-c(RES,list(res1))
    nam.ord<-nam[res1$ordered.clusterings]
    NAM <-c(NAM,list(nam.ord))
  }
  
  stab.methods<-list()
  for (i in 1:(max(nb)-1))  stab.methods<-c(stab.methods, c(nclass=i+1,data.frame(methods=NAMES[[i]][RES[[i]][,1]],  p.value=RES[[i]][,2] ))) 
  if ( verbose == TRUE ) print(stab.methods) 
  
  ##================================================
  ##                GRAPHICS
  ##================================================
  
  spp<-NULL
  for (w in 1:length(NAM)) spp<-c(spp,NAM[[w]])
  PERC<-sort(   table(spp)*100/(max(nb)-1)  ,decreasing=TRUE)
  
  if (!is.null(pdfname)) {
    pdf(paste(pdfname,".pdf",sep=""))
  }
  par(mfrow=c(1,1))
  par(mar=c(12,4,4,2)+0.1)
  
  barplot(PERC, las=2, ylab ="%",ylim=c(0,100), cex.names = 0.75, axes=FALSE)   
  axis(2)
  mtext("Frequency of stable methods")
  
  if (!is.null(pdfname)) {
    dev.off()
  }
  

  t2=Sys.time()
  tdiff=t2-t1
  print(tdiff)
  
  stab.methods
}



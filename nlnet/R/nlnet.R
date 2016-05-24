# source this code to conduct the non-linear relations network reconstruction.
# input: the data matrix with no missing values
# min.fdr.cutoff: the minimun value of the local false discovery cutoff in establishing links between genes
# max.fdr.cutoff: the maximun value of the local false discovery cutoff in establishing links between genes
# conn.proportion: determine how much connections between genes.
# gene.fdr.plot: whether plot a figure with estimated densities, distribution functions, and (local) false discovery rates
# min.module.size: the min number of genes together as a module.
# gene.community.method: it provides three kinds of community detection method:
# "mutilevel", "label.propagation" and "leading.eigenvector".
# use.normal.approx: whether to use the normal approximation to for the null hypothesis. If TRUE, normal approximation is used for every feature, AND all covariances are assumed to be zero. If FALSE, generates permutation based null distribution - mean vector and a variance-covariance matrix.
# normalization: the normalization method for the array. There are three choices - "standardize" means removing the mean of each row and make the standard deviation one; "normal_score" means normal score transformation; "none" means do nothing. In that case we still assume some normalization has been done by the user such that each row has approximately mean 0 and sd 1. 
# If TRUE, normal approximation is used for every feature, AND all covariances are assumed to be zero. If FALSE, generates permutation based null distribution - mean vector and a variance-covariance matrix.
# plot.method: it provides three kinds of ploting method:
# "none" means no graph "communitygraph" means ploting community with graph, "graph" means ploting graph, "membership" means ploting membership of the community
# return the community and its membership
################################################################################

nlnet<-function(input, min.fdr.cutoff=0.05,max.fdr.cutoff=0.2,conn.proportion=0.007,gene.fdr.plot=FALSE,min.module.size=0,gene.community.method="multilevel",use.normal.approx=FALSE,normalization="standardize",plot.method="communitygraph")
{
    
    normrow<-function(array)
    {
        m<-apply(array,1,mean,na.rm=T)
        s<-apply(array,1,sd,na.rm=T)
        array<-(array-m)/s
        return(array)
    }
    
    gene.specific.null<-function(array, B=500)
    {
        null.mat<-matrix(0, nrow=nrow(array), ncol=B)
        l<-ncol(array)
        d.array<-array[,1:(l-1)]
        for(i in 1:B)
        {
            this.order<-sample(l, l, replace=FALSE)
            for(j in 1:(l-1)) d.array[,j]<-abs(array[,this.order[j+1]]-array[,this.order[j]])
            null.mat[,i]<-apply(d.array, 1, sum)
        }
        r<-cbind(apply(null.mat, 1, mean), apply(null.mat, 1, sd))
        return(r)
    }
    
    scol.matrix.order<-function(array,x) # x is the vector, a is the matrix, find ordered distance of rows.of.a|x
    {
        if(is.null(nrow(array)) | nrow(array) == 1)
        {
            array<-as.vector(array)
            array<-array[order(x)]
            d<-array[2:length(array)]-array[1:(length(array)-1)]
            dd<-sum(abs(d),na.rm=T)
        }else{
            array<-array[,order(x)]
            d<-array[,2:ncol(array)]-array[,1:(ncol(array)-1)]
            dd<-apply(abs(d),1,sum,na.rm=T)
        }
        return(dd)
    }
    
    scol.matrix<-function(a, direction=2)  # when direction is 1, scol.matrix[i,j] = SCOL(a[i,], a[j,]), j|i
    {
        
        rdmat<-matrix(0, ncol=nrow(a), nrow=nrow(a))
        for(j in 1:nrow(a))
        {
            rdmat[j,]<-scol.matrix.order(a, a[j,])
        }
        
        if(direction == 2)
        {
            rdmat.diff<-rdmat-t(rdmat)
            sel<-which(rdmat.diff > 0)
            rdmat[sel]<-t(rdmat)[sel]
        }
        return(rdmat)
    }
    gene.specific.p<-function(null.distr, new.d)
    {
        for(i in 1:length(new.d))
        {
            new.d[i]<-pnorm(new.d[i], mean=null.distr[i,1], sd=null.distr[i,2], lower.tail=TRUE)
        }
        return(new.d)
    }
    gene.specific.q<-function(new.d)
    {
        for(i in 1:length(new.d))
        {
            new.d[i]<-qnorm(new.d[i], mean=0, sd=1, lower.tail=TRUE)
        }
        return(new.d)
    }
    normscore.row<-function(a)
    {
        #library(coin)
        b<-t(apply(a, 1, normal_trafo))
        return(b)
    }
    
    #library(fdrtool)
    #library(igraph)
    
    n<-nrow(input)
    m<-ncol(input)
    if(use.normal.approx & normalization=="none") normalization<-"standardize"
    
    if(normalization=="normal_score")
    {
        array<-normscore.row(input)
    }else{
        if(normalization=="standardize")
        {
            array<-normrow(input)
        }else{
            if(normalization=="none")
            {
                array<-input
            }else{
                return("wrong normalization method.")
            }
        }
    }
    orig.array <-array
    if(!use.normal.approx)
    {
        null.distr<-gene.specific.null(array)
    }else{
        m.emp<-2/sqrt(pi)*(m-1)
        s.emp<-sqrt(m-1)*4*(2*pi+3*sqrt(3)-9)/3/pi
        null.distr<-cbind(rep(m.emp, nrow(array)), rep(s.emp, nrow(array)))
    }
    sim.mat<-scol.matrix(array,direction=1)  ## similarity matrix by SCOL, asymmetric, column given row
    d.mat<-sim.mat
    d.tmp.mat<-d.mat
    for(i in 1:nrow(sim.mat)) d.tmp.mat[i,]<-gene.specific.p(null.distr, sim.mat[i,])
    for(i in 1:nrow(sim.mat)) d.mat[i,]<-gene.specific.q(d.tmp.mat[i,])
    diag(d.mat)<-0
    gene.rel.mat<-matrix(0,nrow=n,ncol=n)#the matrix store the relationship between two genes, 1 means having relationship while 0 means no relationship
    gene.fdr.mat<-matrix(0,nrow=n,ncol=n)
    for(i in 1:n) {
      sim.vec<-d.mat[i,]
      suppressWarnings(t.locfdr<-fdrtool(sim.vec, statistic="normal", plot=gene.fdr.plot,color.figure=TRUE,verbose=FALSE,cutoff.method="pct0",pct0=0.75))
      t.row.lfdr<-as.vector(t.locfdr$lfdr)
      gene.fdr.mat[i,]<-t.row.lfdr
    }
    #now dynamic adjust the fdr cutoff
    last.fdr.cutoff<-quantile(as.vector(gene.fdr.mat),conn.proportion)
    cat('------fdr cutoff aimed ',last.fdr.cutoff,'-------\n')
    if(last.fdr.cutoff > max.fdr.cutoff){
        gene.fdr.cutoff<- max.fdr.cutoff
    }else{
     if(last.fdr.cutoff < min.fdr.cutoff){
        gene.fdr.cutoff<- min.fdr.cutoff
      }else{
        gene.fdr.cutoff<-last.fdr.cutoff
      }
    }
    cat('------fdr cutoff real ',gene.fdr.cutoff,'-------\n')
    for(i in 1:n){
      for(j in 1:n){
        if(gene.fdr.mat[i,j] < gene.fdr.cutoff){
          gene.rel.mat[i,j]<-1 
        } 
      }
    }
    gene.rel.mat<-t(gene.rel.mat)+gene.rel.mat ##symmetric
    
    gene.graph<-graph.adjacency(gene.rel.mat, mode="undirected", weighted=NULL)
    if(gene.community.method=="multilevel"){
        commu<-multilevel.community(gene.graph, weights=NA)
    }else if(gene.community.method=="label.propagation"){
        commu<-label.propagation.community(gene.graph)
    }else if(gene.community.method=="leading.eigenvector"){
      errormsg = "yes"
      graph.tmp<-gene.graph
      commu<-tryCatch({
        commu<-leading.eigenvector.community(graph.tmp)
        return(commu)
      },error = function(e) {
        cat('---max interations reached, so restrict steps to 100')
        errormsg2 = tryCatch({
          commu<-leading.eigenvector.community(graph.tmp,steps=100)
          return(commu)
        },error = function(e){
          cat('---max interations reached, so restrict steps to 10')
          commu<-leading.eigenvector.community(graph.tmp,steps=10)
          return(commu)
        })
      })
    }else{
      print("can not find the method, use multilevel as the default one")
      commu<-multilevel.community(gene.graph, weights=NA)
    }
    
    mem<-commu$membership
    ta<-table(mem)
    for(i in 1 : length(ta)){
      if(ta[i] < min.module.size){
        ##now replace the value in the mem that equals i to zero
        for(j in 1 : length(mem)){
          if(mem[j] == i){
            mem[j]<-0  
          }
        }
      }
    }
    ## then we adjust the index of community after we have set some value into zero
    ta.new<-table(mem)
    data.frame<-as.data.frame(ta.new)$mem
    if(data.frame[1] == 0){
      for(i in 1 : length(ta.new)){
        if(i > 1){
          index.origin <- data.frame[i]
          index.new <- i-1
          for(j in 1 : length(mem)){
            if(mem[j] == index.origin){
              mem[j]<-index.new
            }
          }
        }
      }
    }
    if(plot.method=="none"){
        
    }else if(plot.method=="communitygraph"){
      plot(commu,gene.graph)
    }else if(plot.method=="graph"){
      plot.igraph(gene.graph)
    }else if(plot.method=="membership"){ 
        plot(mem) 
    }
    graph.community<-new("list")
    graph.community$algorithm<-gene.community.method
    graph.community$graph<-gene.graph
    graph.community$community<-mem
    return (graph.community)
}


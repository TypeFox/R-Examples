#implementation of K-profiles clustering
#contact: tianwei.yu@emory.edu
#parameters:
#dataset: the data matrix with genes in the row and samples in the column
#nCluster: the number of clusters K
#maxIter: the maximum number of iterations
#p.max: the starting p-value cutoff to exclude noise genes
#p.min: the final p-value cutoff to exclude noise genes

KPC<-function(dataset, nCluster, maxIter =100, p.max=0.2, p.min=0.05){
    
    #library(TSP)
    
    gene.specific.null<-function(array, B=200)
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
    
    
    null.mat<- gene.specific.null(dataset)  #mean and sd for each gene
    
    dims= dim(dataset)
    cat('dim:',dims,'\n')
    ngene= dims[1]
    nsample= dims[2]
    
    if (nCluster > ngene){
        sprintf('datasize %d smaller than the amount of clusters %d\n',ngene,nCluster)
        q()
    }
    
    # initialization randomly
    indCluster= sample(1:nCluster, ngene, replace=TRUE)
    path= list()
    
    for(i in 1:nCluster){
        d= dist(t(dataset[indCluster==i,]))
        tsp= TSP(d)
        set.seed(999999)
        tour= solve_TSP(tsp, method="nn")
        path[[i]]= as.integer(tour)
    }
    
    p_max= p.max           # initial p-value threshold
    p_min= p.min           # more strict p-value
    p_val= p_max
    
    iloop= 0
    unclustered=c()
    while (TRUE){
        preindCluster=indCluster
        isEmpty= rep(TRUE, nCluster)
        unclustered=c()
        
        iloop= iloop+1
        if (iloop%%5 == 1){
            cat('------iteration ',iloop,'-------\n')
        }
        
        dist.mat<-matrix(0, nrow=nrow(dataset), ncol=nCluster)
        
        for(i in 1:nCluster){
            this.data<-dataset[,path[[i]]]
            this.d<-this.data[,2:ncol(this.data)]-this.data[,1:(ncol(this.data)-1)]
            this.d<-abs(this.d)
            dist.mat[,i]<-apply(this.d,1,sum)
        }
        
        indCluster<-max.col(-dist.mat, "last")
        p.rec<-pnorm((dist.mat[cbind(1:nrow(dist.mat), indCluster)]-null.mat[,1])/null.mat[,2], lower.tail=TRUE)
        
        indCluster[p.rec>p_val]<- -1
        if(p_val>p_min) p_val <- p_val- (p_max-p_min)/min(50, maxIter)
        
        for(i in 1:nCluster){
            if(sum(indCluster==i)>0){
                d= dist(t(dataset[indCluster==i,]))
                tsp= TSP(d)
                set.seed(999999)
                tour= solve_TSP(tsp, method="nn")
                path[[i]]= as.integer(tour)
            }
        }
        
        if(sum(preindCluster == indCluster)/nrow(dataset) > 0.99 || iloop> maxIter) break
    }
    
    r<-new("list")
    r$cluster= indCluster
    r$p.list=p.rec
    return (r)
}


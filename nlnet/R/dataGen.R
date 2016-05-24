############# the function to generate simulated dataset
# n.gene: the number of rows of the matrix 
# n.samples: the number of columns of the matrix
# n.grps: the number of hidden clusters
# aver.grp.size: averge number of genes in a cluster
# n.fun.types: number of function types to use
# epsilon: noise level
# n.depend: data generation dependence structure. can be 0, 1, 2
data.gen<-function(n.genes=100, n.samples=10, n.grps=10, aver.grp.size=10, n.fun.types=3, epsilon=1, n.depend=0)
{
  set.seed(Sys.time())
  link<-function(x, type)
  {
    x<-(x-mean(x))/sd(x)
    if(type == 1) return(x)
    if(type == 2) return(sin(2*x))
    if(type == 3) return(x^2)
    if(type == 4) return(abs(x))
    if(type == 5) return(x^3)
    if(type == 6) return(atan(4*x))
  }
  
  a<-matrix(rnorm(n.genes*n.samples),ncol=n.samples)
  curr.count<-0
  g<-new("list")
  for(i in 1:n.grps)
  {
    #		this.size<-rpois(1, aver.grp.size)
    this.size<-aver.grp.size
    if(this.size < 2) this.size<-2
    
    this.mat<-matrix(0, nrow=this.size, ncol=n.samples)
    this.mat[1,]<-rnorm(n.samples)
    for(j in 2:this.size)
    {
      if(n.depend==0)
      {
        this.basis<-c(1, rep(0,j-2))
      }else{
        #				this.basis<-sample(c(1,0), j-1, replace=T, prob=c(min(1, n.depend/(j-1)), 1-min(1, n.depend/(j-1))))
        if(j-1 <= n.depend) 
        {
          this.basis<-rep(1, j-1)
        }else{
          this.basis<-sample(c(rep(1, n.depend), rep(0,j-1-n.depend)), j-1, replace=F)
        }
        
      }
      if(sum(this.basis) > 0)
      {
        x<-rep(0,n.samples)
        for(k in which(this.basis == 1))
        {
          x<-x+link(this.mat[k,], sample(n.fun.types,1))*runif(1,min=-1,max=1)
        }
        #				x[x>quantile(x, 0.95)]<-quantile(x, 0.95)
        #				x[x<quantile(x, 0.05)]<-quantile(x, 0.05)
        this.mat[j,]<-x
        this.mat[j,]<-(this.mat[j,]-mean(this.mat[j,]))/sd(this.mat[j,])
      }else{
        this.mat[j,]<-rnorm(n.samples)
      }
    }
    if(n.depend == 0)
    {
      this.mat[1,]<-link(this.mat[1,], sample(n.fun.types,1))
      this.mat[1,]<-(this.mat[1,]-mean(this.mat[1,]))/sd(this.mat[1,])
    }
    
    if(curr.count+this.size <= n.genes)
    {
      a[(curr.count+1):(curr.count+this.size),]<-this.mat
      g[[length(g)+1]]<-(curr.count+1):(curr.count+this.size)
    }
    curr.count<-curr.count+this.size		
  }
  a<-a+matrix(rnorm(n.genes*n.samples, sd=epsilon),ncol=n.samples)
  
  g2<-rep(0, nrow(a))
  for(i in 1:length(g)) g2[g[[i]]]<-i
  
  r<-new("list")
  r$data<-a
  r$grps<-g2
  return(r)
}
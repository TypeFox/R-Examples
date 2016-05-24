"mkle" <- function(data,bw=2*sd(data),kernel=c("gaussian", "epanechnikov", "rectangular", 
                   "triangular", "biweight", "cosine", "optcosine"),gridsize=2^14) {

if(bw<0)stop('Bandwidth must be positive')
if(bw==0)return(mean(data))
if(gridsize<4)stop('gridsize must be at least 4')
if(gridsize<2^11)warning('gridsize might be too small for reliable results')

kernel <- match.arg(kernel)
kde<-density(data,bw=bw,kernel=kernel,from=min(data)-2*bw,to=max(data)+2*bw,n=gridsize)
min<-kde$x[1]
grid<-kde$x[2]-min

return(mean(data)+optimize(klik,maximum=TRUE,lower=-2*bw,upper=2*bw,data=data,kde=kde,grid=grid,min=min)$maximum)

}


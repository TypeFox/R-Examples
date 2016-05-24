Silh <-
function(p){
# %%%%Silhouette for probabilistic clustering methods
# %%%%%INPUT
# %p probability matrix
# %%%%%OUTPUT
# % Silhouettes
  n=nrow(p)
nc=ncol(p)

la=max.col(p)
ds=disS(p)

m=cbind(la,ds)
m2=m[sort.list(m[,2]), ]
m=m2[sort.list(m2[,1]), ]

ss=m[,2]

ll=table(la)
pl=cumsum(ll)
#pl=[0,pl(1:nc-1)];
ll=round(ll/2)
  pl=c(0,pl[1:(nc-1)])
tic=ll+pl
#tcks=ll+pl;
  barplot(ss,space=0, main="Shiluette plot", horiz=TRUE,xlab='Silhouette Value', ylab='Cluster',xpd=F)
axis(2,at=tic, labels=(1:nc))

}

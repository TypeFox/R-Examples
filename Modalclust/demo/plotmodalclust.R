data(disc2d.hmac)

### hard plot
hard.hmac(disc2d.hmac,n.cluster=2)
#hard.hmac(disc2d.hmac) press "enter" to give the hard plot of next level

### soft plot 
soft.hmac(disc2d.hmac,level=3)

### contour plot
contour(disc2d.hmac,n.cluster=3)

### plot cluster tree 
plot(disc2d.hmac)

### choose the cluster near the origin
choose.cluster(disc2d.hmac,x=c(0,0),n.cluster=3)


require(frontier)
data(front41Data)


summary(sfa(log( output ) ~ log( capital ) + log( labour ), data=front41Data))
summary(spfrontier(log( output ) ~ log( capital ) + log( labour ), data=front41Data))


sfa <- sfa(log( output ) ~ log( capital ) + log( labour ), data=front41Data, truncNorm=T)
summary(sfa)
ssf <- spfrontier(log( output ) ~ log( capital ) + log( labour ), data=front41Data, inefficiency="truncated")
summary(ssf)

logLikelihood(log( output ) ~ log( capital ) + log( labour ), data=front41Data, inefficiency="truncated",values=c(0.46453,0.28327,0.54098,0.22762,0.91568,-2.84150))

mu <- seq(-3,0, by=0.01)
sigma<- seq(0.3,1, by=0.01)
z<-matrix(nrow=length(mu),ncol=length(sigma))
i<-1
j<-1
for(m in mu){
    for(s in sigma){
        z[i,j]<-logLikelihood(log( output ) ~ log( capital ) + log( labour ), data=front41Data, inefficiency="truncated",values=c(0.46453,0.28327,0.54098,0.22762,s,m))
        j <- j + 1
    }
    i <- i + 1
    j <- 1
}
z[z < -1000] <- Inf
z[z == Inf] <- min(z)

contour(mu, sigma, z, nlevels=30, col=rev(terrain.colors(50)), xlab=expression(italic(mu)), ylab=expression(italic(sigma[u])), labcex=0.8)

x<-c(-2.6491,-1.4106)
y<-c(0.8903,0.7053)
points(x,y)
text(x, y,labels=c("spfrontier", "frontier"),pos=4)
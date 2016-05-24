# Example 2:  bivariate beta with repeated measures
library(zoib)
data("BiRepeated", package = "zoib")

eg2 <- zoib(y1|y2 ~ x|1|x, data= BiRepeated, n.response=2,
            random=1, EUID= BiRepeated$id,
            zero.inflation = FALSE, one.inflation = FALSE,				
            prior.Sigma = "VC.unif", n.iter=5, n.thin=1, n.burn=1)
coeff<- eg2$coeff
ypred<- eg2$ypred
Xb<- eg2$Xb
Xd<- eg2$Xd
Xb0<- eg2$Xb0

if(0){
eg2 <- zoib(y1 ~x|1|x, data= BiRepeated, n.response=1,
            random=1, EUID= BiRepeated$id, joint=FALSE,
            zero.inflation = FALSE, one.inflation = FALSE,  			
            prior.Sigma = "UN.halfcauchy", n.iter=600, n.thin=5, n.burn=10,            
            inits=list(list(b0=NULL,b1=NULL,b=matrix(c(-1.3,-2.6),2,1),
                            d=matrix(1.75,1,1),sigma=c(0.16,0.25),R=c(1,0.12,1)),
                       list(b0=NULL,b1=NULL,b=matrix(c(-0.7,-1.4),2,1),
                            d=matrix(3.25,1,1),sigma=c(0.25,0.16),R=c(1,0.08,1)))
)
coeff<- eg2$coeff
traceplot(coeff)


eg2 <- zoib(y1|y2 ~ x|1|x, data= BiRepeated, n.response=2,
            random=1, EUID= BiRepeated$id,
            zero.inflation = FALSE, one.inflation = FALSE,  			
            prior.Sigma = "UN.halfcauchy", n.iter=100, n.thin=2, n.burn=1,
            inits=list(list(b0=NULL,b1=NULL,b=matrix(c(-1.3,-2.6,0.5,1.0),2,2),
                            d=matrix(c(1.75,2),1,2),sigma=c(0.16,0.25),R=c(1,0.12,1)),
                       list(b0=NULL,b1=NULL,b=matrix(c(-0.7,-1.4,0.5,1.0),2,2),
                            d=matrix(c(3.25,2),1,2),sigma=c(0.25,0.16),R=c(1,0.08,1)))
)
coeff<- eg2$coeff
traceplot(coeff)            
}
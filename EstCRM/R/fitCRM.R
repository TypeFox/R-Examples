fitCRM <-
function(data,ipar,est.thetas,max.item,group=20) {

n=ncol(data)
N=nrow(data)

if(is.data.frame(data)==FALSE) stop("The input response data is not a data frame.
Please use as.data.frame() and convert your response data to a data frame object before the analysis")

if(dim(data)[2]!=length(max.item)) stop("The number of columns in the data is not equal to the length of max.item vector")

if(dim(data)[2]!=nrow(ipar)) stop("The number of columns in the data is not equal to the number of rows in the item parameter matrix")

for(i in 1:n) {
if(max(na.omit(data[,i]))> max.item[i]) stop("The column ",i," has values higher than the maximum available score in the
user specified max.item vector. Please check and clean your data.")
}

if(class(est.thetas)!="CRMtheta") stop("The estimated thetas is not an object created by EstCRMperson()")


data$theta <- est.thetas$thetas[,2]
data <- data[order(data$theta),]
data$groups <- cut2(data$theta,g=group)
labels <- unique(data$groups)


expected <- function(t,item) { 
prob <- function(t,x) {
v1=ipar[item,1]*(t-ipar[item,2]-((1/ipar[item,3])*log((x-.5)/(max.item[item]-x+.5))))
v2=ipar[item,1]*(t-ipar[item,2]-((1/ipar[item,3])*log((x+.5)/(max.item[item]-x-.5))))
integrate(dnorm,v2,v1)$value
}
all.scores <- 1:(max.item[item]-1)
all.probs <- c()
for(i in 1:length(all.scores)) { all.probs[i]=prob(t,all.scores[i]) }
sum(all.scores*all.probs)
} 

fitindex <- vector("list",n)                
median.thetas.list <- vector("list",n)     
mean.obs.score.list <- vector("list",n)    
sd.obs.score.list <- vector("list",n)      
exp.score.list <- vector("list",n)         

for(item in 1:n) { 

median.thetas <- c()
exp.score <- c()
mean.obs.score <- c()
sd.obs.score <- c()

for(i in 1:length(labels)) { 

median.thetas[i]= median(data[which(data$groups==labels[i]),]$theta,na.rm=TRUE)
exp.score[i]= expected(median.thetas[i],item)
mean.obs.score[i]= mean(data[which(data$groups==labels[i]),item],na.rm=TRUE)
sd.obs.score[i]= sqrt(var(data[which(data$groups==labels[i]),item],na.rm=TRUE)/50)
}

median.thetas.list[[item]] <- median.thetas
mean.obs.score.list[[item]] <- mean.obs.score
sd.obs.score.list[[item]] <- sd.obs.score
exp.score.list[[item]] <- exp.score 
fitindex[[item]]<-(mean.obs.score - exp.score)/sd.obs.score
}

fits <- as.data.frame(matrix(nrow=length(labels),ncol=n+1))
colnames(fits) <- c("Interval",colnames(data)[1:n])
fits[,1]=labels
for(i in 2:(n+1)) { fits[,i]=fitindex[[i-1]] }

emp.plots <- vector("list",nrow(ipar))
for(u in 1:nrow(ipar)) {
emp.prop <- cbind(expand.grid(median.thetas.list[[u]],1:max.item[u]),labels)
emp.prop$prob <- NA
for(i in 1:dim(emp.prop)[1]) {
emp.prop[i,4]=length(which(round(data[which(data$groups==emp.prop[i,3]),u],0)==emp.prop[i,2]))/length(which(data$groups==emp.prop[i,3]))
}
emp.plots [[u]] <- wireframe(emp.prop[,4]~emp.prop[,1]*emp.prop[,2],xlab="Ability Scale",ylab="Response Scale",
zlab="Prob",zlim=c(0,max(emp.prop[,4])+.05),screen=list(z =-50, x = -70),
scales =list(arrows = FALSE,tck=.5),main=paste("Category Response Curves - Item ",u,sep=""))
}
CRM.fit <- list(fit.stat=fits,emp.irf=emp.plots)
return(CRM.fit)
}


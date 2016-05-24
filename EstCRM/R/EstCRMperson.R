EstCRMperson <-
function(data,ipar,min.item,max.item) { 

N=nrow(data)
n=ncol(data)

if(is.data.frame(data)==FALSE) stop("The input response data is not a data frame.
Please use as.data.frame() and convert your response data to a data frame object before the analysis")
if(length(max.item)!=length(min.item)) stop("The length of max.item vector is not equal to the length of 
min.item vector. Please check your inputs")
if(dim(data)[2]!=length(max.item)) stop("The number of columns in the data is not equal to the length of max.item vector")
if(dim(data)[2]!=length(min.item)) stop("The number of columns in the data is not equal to the length of min.item vector")
for(i in 1:n) {
if(max(na.omit(data[,i]))> max.item[i]) stop("The column ",i," has values higher than the maximum available score in the
user specified max.item vector. Please check and clean your data.")
}
for(i in 1:n) {
if(min(na.omit(data[,i]))< min.item[i]) stop("The column ",i," has values smaller than the minimum available score in the
user specified min.item vector. Please check and clean your data.")
}

#############

for(i in 1:n){data[,i]= data[,i]-min.item[i]}
max.item <- max.item-min.item
min.item <- c(0,0,0,0)
for(i in 1:n) {
if(length(which(data[,i]==max.item[i]))!=0) {
data[which(data[,i]==max.item[i]),i]=max.item[i]-.01
}
if(length(which(data[,i]==0))!=0) {
data[which(data[,i]==0),i]=.01
}
}
data.original <- data

for(i in 1:n){data[,i]= log(data[,i]/(max.item[i]-data[,i]))}

###############

theta.vector <- c()
for(k in 1:nrow(data)){
alll <- c()
itemss <- which(data[k,]!="NA")
for(i in 1:length(itemss)){
alll[i]=ipar[itemss[i],1]^2*(ipar[itemss[i],2]+data[k,itemss[i]]/ipar[itemss[i],3])
}
theta.vector[k]=sum(alll)/sum(ipar[itemss,1]^2)
}
theta.vector <- (theta.vector-mean(theta.vector,na.rm=TRUE))/sd(theta.vector)
se=1/sqrt(sum(ipar[itemss,1]^2))
thetas <- cbind(c(1:N),theta.vector,se)
colnames(thetas) <- c("ID","Theta Est.","SE")
thetas <- list(thetas=thetas)
class(thetas) <- "CRMtheta"
return(thetas)
}


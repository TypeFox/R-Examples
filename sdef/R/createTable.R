createTable <-
function(output.ratio,output.bay,dir=getwd(),h=NULL){

if(output.ratio$pvalue==TRUE){
matrix.results =  cbind(output.ratio$h,output.ratio$ratio,round(output.bay,3),output.ratio$Common,output.ratio$DE)
}
if(output.ratio$pvalue==FALSE){
matrix.results =  cbind(1-output.ratio$h,output.ratio$ratio,round(output.bay,3),output.ratio$Common,output.ratio$DE)
}

lists = dim(output.ratio$DE)[2]
namesDE = paste("O",seq(1,lists),rep("+",lists),sep="")
nameCommon="O"
for(i in 1:lists){
nameCommon=paste(nameCommon,"1",sep="")
}

names.matrix = c("h","T",colnames(output.bay)[1],colnames(output.bay)[2],colnames(output.bay)[3],nameCommon,namesDE)
dimnames(matrix.results)[[2]]<-names.matrix

setwd(dir)
write.csv(matrix.results,"Output.csv",row.names=FALSE)

#Decision rules:
#1) Maximum for CI not including 1
if(length(output.bay[round(output.bay[,1],3)>1,2])==0){
cat("WARNING: the requested contrast is under-represented in the data (Rmax<1)\n")
}
if(length(output.bay[round(output.bay[,1],3)>1,2])>0){
max.R = max(output.bay[round(output.bay[,1],3)>1,2])
max.T = max(output.ratio$ratios)
h.Tmax=output.ratio$h[output.ratio$ratios==max.T]
h.Rmax=output.ratio$h[output.bay[,2]==max.R]
maximum1 = matrix.results[matrix.results[,4]==round(max.R,3),]
maximum=c(maximum1[1],maximum1[2],maximum1[3],maximum1[4],maximum1[5],maximum1[6])
for(i in 1:lists){
maximum=c(maximum,maximum1[6+i])
}
maximum=matrix(maximum,nrow=1,ncol=6+lists)
dimnames(maximum)[[2]]<-names.matrix

if(h.Tmax != h.Rmax){
maximum1.T = matrix.results[matrix.results[,2]==max.T,]
maximum.T=c(maximum1.T[1],maximum1.T[2],maximum1.T[3],maximum1.T[4],maximum1.T[5],maximum1.T[6])
for(i in 1:lists){
maximum.T=c(maximum.T,maximum1.T[6+i])
}
maximum.T=matrix(maximum.T,nrow=1,ncol=6+lists)
dimnames(maximum.T)[[2]]<-names.matrix
}

if(length(output.bay[output.bay[round(output.bay[,1],3)>1,2]>=2,1])>0){
#2) Rule 2
R2 = max(matrix.results[round(output.bay[,2],3)>=2 & round(output.bay[,1],3)>1 ,1])
rule2_temp = matrix.results[matrix.results[,1]==R2,]
rule2=c(rule2_temp[1],rule2_temp[2],rule2_temp[3],rule2_temp[4],rule2_temp[5],rule2_temp[6])
for(i in 1:lists){
rule2=c(rule2,rule2_temp[6+i])
}
rule2=matrix(rule2,nrow=1,ncol=6+lists)
dimnames(rule2)[[2]]<-names.matrix
if(is.null(h)){
if(h.Tmax == h.Rmax){
return(list(maximum.R=maximum,rule2=rule2))
}
if(h.Tmax != h.Rmax)
return(list(maximum.R=maximum, maximum.T=maximum.T,rule2=rule2))
}
####ruleh############################################
if(!is.null(h)){
ruleh=list()
for(j in 1:length(h)){
ruleh_temp = matrix.results[matrix.results[,1]==h[j],]
ruleh_temp2=c(ruleh_temp[1],ruleh_temp[2],ruleh_temp[3],ruleh_temp[4],ruleh_temp[5],ruleh_temp[6])
for(i in 1:lists){
ruleh_temp2=c(ruleh_temp2,ruleh_temp[6+i])
}
ruleh_temp2=matrix(ruleh_temp2,nrow=1,ncol=6+lists)
dimnames(ruleh_temp2)[[2]]<-names.matrix
ruleh[[j]]<-ruleh_temp2
}

if(h.Tmax == h.Rmax){
return(list(maximum.R=maximum,rule2=rule2,User=ruleh))
}
if(h.Tmax != h.Rmax){
return(list(maximum.R=maximum,maximum.T=maximum.T,rule2=rule2,User=ruleh))
}

}
#####################################################
}

if(length(output.bay[output.bay[round(output.bay[,1],3)>1,2]>=2,1])==0){
if(is.null(h)){
if(h.Tmax == h.Rmax){
return(list(maximum.R=maximum))
}
if(h.Tmax != h.Rmax){
return(list(maximum.R=maximum,maximum.T=maximum.T))
}
}
####ruleh############################################
if(!is.null(h)){
ruleh=list()
for(j in 1:length(h)){
ruleh_temp = matrix.results[matrix.results[,1]==h[j],]
ruleh_temp2=c(ruleh_temp[1],ruleh_temp[2],ruleh_temp[3],ruleh_temp[4],ruleh_temp[5],ruleh_temp[6])
for(i in 1:lists){
ruleh_temp2=c(ruleh_temp2,ruleh_temp[6+i])
}
ruleh_temp2=matrix(ruleh_temp2,nrow=1,ncol=6+lists)
dimnames(ruleh_temp2)[[2]]<-names.matrix
ruleh[[j]]<-ruleh_temp2
}

if(h.Tmax == h.Rmax){
return(list(maximum.R=maximum,User=ruleh))
}
if(h.Tmax != h.Rmax){
return(list(maximum.R=maximum,maximum.T,User=ruleh))
}

}
#####################################################
}

}
}


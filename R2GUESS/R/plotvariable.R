plotvariable <-
function(x,threshold.model=0.05,file.annotation=NULL){
if(is.null(x$label.Y)) Pheno <- paste("Y",1:x$q,sep="",collapse="_") else Pheno <- paste("Y: ",paste(x$label.Y,collapse="/"),sep="") 

label.X <- as.character(1:x$p)
if(!is.null(file.annotation)){
#annot <- read.table(paste(x$path.input,file.annotation,sep=""),header=TRUE)
annot <- read.table(file.annotation,header=TRUE)
  
label.X <- as.character(annot[,1])
} 

if(threshold.model<1) {
ListSelect <- which(x$BestModels$postProb>threshold.model) 
title1 <- paste(Pheno,"\n Variables selected in each model (posterior > ",threshold.model,")",sep="")}else{
title1 <- paste(Pheno,"\n Variables selected in each model (",threshold.model," best models)",sep="")
ListSelect <-1:threshold.model
}


vect.variable <- unique(unlist(strsplit(x$BestModels$modelName[ListSelect]," ")))

if(is.null(vect.variable)){
cat("No model has been selected with a threshold greater than ", threshold.model, "\n")
res.var=NULL}
else{
if(any(is.na(vect.variable))&length(vect.variable)==1){
cat("The model selected with a threshold: ", threshold.model,"contains the NULL model \n ")
res.var=NULL
}else{
list.best.model <- strsplit(x$BestModels$modelName[ListSelect]," ")
res <- t(sapply(list.best.model,FUN=function(x){vect.variable%in%x}))

colnames(res) <- vect.variable

if(is.null(x$MAP.file)) labx <- label.X[as.numeric(vect.variable)] else labx <- vect.variable


par(mar=c(10,4,4,2))
image(x=1:length(vect.variable),y=1:length(ListSelect),t(res),axes=FALSE,xlab="",ylab="Model",main=title1,col=rev(heat.colors(10)))
par(las=2)
par(cex.axis=1)
axis(1, at = 1:length(vect.variable),labx)
par(cex.axis=0.5)
axis(2, at = 1:length(ListSelect),paste("M",1:length(ListSelect),sep=""))
res.var <- na.omit(vect.variable)
grid(nx=length(vect.variable),ny=length(ListSelect),lty=1,col="gray")
}}

return(output<-list(var.best.models=labx))
}

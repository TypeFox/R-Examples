plotcim <- function(x,select.variable,labelX=NULL){

if(is.null(x$label.Y)) labelY <- paste("Y",1:x$q,sep="") else labelY <- x$label.Y

var.gene <- select.variable
if(length(var.gene)<2){cat("more than on variable need to be selected")
mat.sim <- NULL
}else{
datax <- data.frame(matrix(scan(file=file.path(x$path.input,x$dataX),skip=2,sep=""),ncol=x$p,byrow=TRUE))
if(is.null(labelX)) colnames(datax) <- 1:x$p else colnames(datax) <- labelX
if(!is.null(x$MAP.file)){
if(is.data.frame(x$MAP.file)){annot<- x$MAP.file}else{
annot <- read.table(file.path(x$path.input,x$MAP.file),header=TRUE)}
colnames(datax) <- annot[,1]
}

matY <- data.frame(matrix(scan(file=file.path(x$path.input,x$dataY),skip=2,sep=""),ncol=x$q,byrow=TRUE))

matX <- datax[,var.gene]


mat.sim <- matrix(NA,ncol=length(var.gene),nrow=length(labelY))
for (i in 1:length(var.gene)){
for (j in 1:length(labelY)){
mat.sim[j,i] <- cor(matX[,i],matY[,j])
}}
cim(mat.sim,row.names = labelY,col.names = var.gene)}
return(result=list(mat.sim=mat.sim,labelY=labelY,labelX=var.gene))
}

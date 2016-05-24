inline <- function(source,layer='',attributes,type='l',label='',col='')
{
#read data
#read data
if(layer=="gisobject"||is.null(layer)||layer=="gisdata"||layer==""){
object <- source
}
else{

file=c("[.]csv$","[.]txt$","[.]tab$","[.]dat$")
i=1
nonspatial=0
while (i <= length(file)) {
log=grep(file[i],source)
if(length(log)>0)
{
file.type=file[i]
thisfile=source
nonspatial=nonspatial+1
}
else{
log=grep(file[i],layer)
if(length(log)>0){
file.type=file[i]
thisfile=paste(source,layer,sep = "/")
nonspatial=nonspatial+1
}
}
i=i+1
}
if(nonspatial<1){
if(layer==""||layer=="nofile"||layer=="notfile"||layer=="nonfile"){
object=source
}
else{
#check input data
if(source==""){
print("Full Folder path is needed:source")
return
}
if(layer==""){
print("Map name is needed: layer")
return
}

object <- readOGR(source, layer)
}
}
else{

log=grep('[.]csv$',file.type)
if(length(log)>0)
{
file.type=file[i]
object=read.csv(thisfile,header=TRUE)

}
else{
object=read.table(thisfile,header=TRUE)
}

}
}
#line chart Region

#split
attribute=strsplit(attributes,",")
attribute1=attribute[[1]][1]
attribute2=attribute[[1]][2]

slice1=object[[attribute1]]
slice2=object[[attribute2]]
slices=object[[attribute2]]

main=label
xlab=toupper(attributes)
#default label
if(label==""){
main=paste("Line plot of", attributes)
}

#label only names without %

#define colours
colour=col
if(col=="")
{
colour="black"
}
ylab=attribute2
#construct line
if(type=="map"||type=="spatial"){
if(attribute2!=""){
#par(mfrow=c(1,2)) 
lists=c(attribute1,attribute2)
plot(spplot(object,lists,scales = list(draw = TRUE)))
plot(spplot(object,attribute1,scales = list(draw = TRUE)))
}
plot(spplot(object,attribute2,scales = list(draw = TRUE)))
}
else{

plot(slice1,slice2, xlab=attribute1,ylab=attribute2,col=colour,  
main=main,type="n")
}
lines(slice1,slice2,type=type)
data=c(object[[attribute1]],object[[attribute2]])
list(data=data,table=slices,source=source,layer=layer,
attributes=attributes,type=type,label=label,colour=col,attribute1=attribute1,
attribute2=attribute2)
}

#generic function
linemap<-function(source,layer='',attributes,type='l',label="",col='') UseMethod ("linemap")

#default function
linemap.default<-function(source,layer='',attributes,type='l',label="",col='')
{

factor<-inline(source,layer,attributes,type,label,col)

factor$call<-match.call()

class(factor)<-"linemap"
factor
}

#summary function
summary.linemap<-function(object,...)
{
x<-object
cat("Call:\n")
print(x$call)
#cat("\n Summary  \n")
#print(summary(x$data))

}

#print function
print.linemap<-function(x,...)
{
cat("Call:\n")
print(x$call)
#cat("\nCall: Data\n")
print(x$data)

#cat("\n Summary: attribute1  \n")
#print(summary(x$attribute1))

#cat("\n Summary: attribute2  \n")
#print(summary(x$attribute2))


}

#plot function
plot.linemap<-function(x,...)
{
linemap(source=x$source,layer=x$layer,attributes=x$attributes,type=x$type,col=x$colour,label=x$label)
}

#graph=linemap("E:/geokings/advancegis/R-process", type="map","KwabibiremAtiwa",attributes='Age,FArea',col='black',label="")


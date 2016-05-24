inbox <- function(source,layer='',attribute,type='',label='',col='',factor="")
{
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
#box chart Region
slices=object[[attribute]]
pct <- round(slices/sum(slices)*100)
lbl=names(slices)
lbls <- paste(lbl, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="")
main=label
ylab=toupper(attribute)
#default label
if(label==""){
main=paste("Box plot of", attribute)
}

#label only names without %

#define colours
colour=col
if(col==""){
colour="light blue"
}
if(col=="gray"||col=="black"||col=="white"){
colour=gray(seq(0.4,1.0,length=6))
}

#set type
notch=FALSE
if(type=="notch"){
notch=TRUE
}

#insert factors
#construct bar
if(type=="map"||type=="spatial"){
if(factor!=""){
#par(mfrow=c(1,2)) 
lists=c(attribute,factor)
plot(spplot(object,lists,scales = list(draw = TRUE)))
plot(spplot(object,factor,scales = list(draw = TRUE)))
}
plot(spplot(object,attribute,scales = list(draw = TRUE)))
}
else{

if(factor!=""){
factors=toupper(object[[factor]])
xlab=toupper(factor)
ylab=toupper(attribute)
boxplot(slices~factors, col=colour,  
main=main,xlab=xlab,ylab=ylab,notch=notch)
}
else{
#construct box
boxplot(slices, col=colour,  
main=main,ylab=ylab,notch=notch)
}
}
data=object[[attribute]]
list(data=data,table=slices,source=source,layer=layer,
attribute=attribute,type=type,label=label,factor=factor,colour=col)
}

#generic function
boxmap<-function(source,layer='',attribute,type='',label="",col='',factor="") UseMethod ("boxmap")

#default function
boxmap.default<-function(source,layer,attribute,type='',label="",col='',factor="")
{

factor<-inbox(source,layer,attribute,type,label,col,factor)

factor$call<-match.call()

class(factor)<-"boxmap"
factor
}

#summary function
summary.boxmap<-function(object,...)
{
x<-object
cat("Call:\n")
print(x$call)
cat("\n Summary statistics ")
print(x$layer)

summary(x$data)
}

#print function
print.boxmap<-function(x,...)
{
cat("Call:\n")
print(x$call)

cat("\n Summary statistics ")
print(x$layer)
summary(x$data)

cat("\nCall: Data\n")
print(x$data)
}

#plot function
plot.boxmap<-function(x,...)
{
boxmap(source=x$source,layer=x$layer,attribute=x$attribute,type=x$type,col=x$colour,label=x$label,factor=x$factor)
}
#graph=boxmap("E:/geokings/advancegis/R-process", "KwabibiremAtiwa",attribute='Age',col='black',label="",factor="Sex",type="group")
#boxmap(source=meuse,layer="nofile",attribute='zinc',factor="landuse",type="notch",col="red")



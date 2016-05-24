inbar <- function(source,layer='',attribute,type='',label='',col='',factor='')
{
#read data

if(layer=="gisobject"||is.null(layer)||layer=="gisdata"||layer==""){
object <- source
}
else
{
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
#bar chart Region
slices=table(toupper(object[[attribute]]))
if(factor!=""){
slices=table(toupper(object[[attribute]]),toupper(object[[factor]]))
}

pct <- round(slices/sum(slices)*100)
lbl=names(slices)
lbls <- paste(lbl, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="")
main=label
xlab=toupper(attribute)
#default label
if(label==""){
main=paste("Bar Chart of", attribute)
}
legend.text=TRUE
beside <- FALSE
#label only names without %
if(type=="simple"){
lbls <- names(slices)
legend.text=FALSE
}
if(type=="group"){
beside <- TRUE
}
if(type=="stack"){
beside <- FALSE
}



#define colours
colour=col

if(col=="gray"||col=="black"||col=="white"){
colour=gray(seq(0.4,1.0,length=6))
}

if(col==""&&type=="simple"){
colour="light blue"
}
if(col==""){
colour=rainbow(length(lbls))

}


#construct bar
if(type=="map"||type=="spatial"){
plot(spplot(object[attribute],scales = list(draw = TRUE)))
}
else{

if(type=="simple"||type==""){
barplot(slices, main=main,xlab=xlab,cex.names=0.6)
}
else{
barplot(slices,  
main=main,xlab=xlab,legend = rownames(slices), beside=beside,cex.names=0.6)
}
}

data=object[[attribute]]
list(data=data,table=slices,source=source,layer=layer,
attribute=attribute,type=type,label=label,colour=col,factor=factor)
}


#generic function
#barmap<-function(source,layer='',attribute,type='',label='',col='',factor='') UseMethod ("barmap")
barmap<-function(source,layer='',attribute,type='',label="",col='',factor='') UseMethod ("barmap")

#default function
barmap.default<-function(source,layer='',attribute,type='',label="",col='',factor='')
{

factor<-inbar(source,layer,attribute,type,label,col,factor)

factor$call<-match.call()

class(factor)<-"barmap"
factor
}

#summary function
summary.barmap<-function(object,...)
{
x<-object
cat("Call:\n")
print(x$call)
cat("\n Table  \n")
print(x$table)

}

#print function
print.barmap<-function(x,...)
{
cat("Call:\n")
print(x$call)
cat("\nCall: Data\n")
print(x$data)

cat("\n Table  \n")
print(x$table)

}

#plot function
plot.barmap<-function(x,...)
{
barmap(source=x$source,layer=x$layer,attribute=x$attribute,type=x$type,col=x$colour,label=x$label,factor=x$factor)
}
#graph=barmap("E:/geokings/advancegis/R-process","KwabibiremAtiwa",attribute='Age',factor="")
#barmap(source=meuse,attribute='lime',factor="landuse",type="group",col="red")

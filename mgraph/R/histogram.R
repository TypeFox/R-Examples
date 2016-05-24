inhist <- function(source,layer='',attribute,type='',label='',col='',trans='')
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
#hist chart Region
slices=object[[attribute]]
pct <- round(slices/sum(slices)*100)
lbl=names(slices)
lbls <- paste(lbl, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="")
main=label
xlab=toupper(attribute)
#default label
if(label==""){
main=paste("Histogram of", attribute)
}


#define colours
colour=col
if(col==""){
colour="light blue"
}
#transformation
if(trans=="log"){
slices=log(slices)
}

if(trans=="sqrt"){
slices=sqrt(slices)
}

if(trans=="scale"){
slices=scale(slices)
}


#construct hist
if(type=="map"){
plot(spplot(object[attribute],scales = list(draw = TRUE)))
}
else{

if(type!="density"){

h=hist(slices, col=colour,  
main=main,xlab=xlab)
}
#label only names without %
if(type=="normal"||type=="Normal"){
xfit<-seq(min(slices),max(slices),length=40)
yfit<-dnorm(xfit,mean=mean(slices),sd=sd(slices))
yfit <- yfit*diff(h$mids[1:2])*length(slices)
lines(xfit, yfit, col="black", lwd=2) }

if(type=="density"||type=="Density"){
main=paste("Density Plot of", attribute)

d <- density(slices)
plot(d, main=main)
polygon(d, col=col, border="black")  }

}
data=object[[attribute]]
list(data=data,table=slices,source=source,layer=layer,
attribute=attribute,type=type,label=label,colour=col)
}

#generic function
histmap<-function(source,layer='',attribute,type='',label="",col='',trans='') UseMethod ("histmap")

#default function
histmap.default<-function(source,layer='',attribute,type='',label="",col='',trans='')
{

factor<-inhist(source,layer,attribute,type,label,col,trans)

factor$call<-match.call()

class(factor)<-"histmap"
factor
}

#summary function
summary.histmap<-function(object,...)
{
x<-object
cat("Call:\n")
print(x$call)
cat("\n Summary  \n")
print(summary(x$call))


}

#print function
print.histmap<-function(x,...)
{
cat("Call:\n")
print(x$call)
cat("\nCall: Data\n")
print(x$data)

cat("\n Table  \n")
print(table(x$table))

}

#plot function
plot.histmap<-function(x,...)
{
histmap(source=x$source,layer=x$layer,attribute=x$attribute,type=x$type,col=x$colour,label=x$label)
}

#graph=histmap("E:/geokings/advancegis/R-process","KwabibiremAtiwa",attribute='Age',trans='scale',type="map")


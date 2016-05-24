indot <- function(source,layer='',attribute,type='',label='',col='',symbol='')
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
#dot chart Region
slices=object[[attribute]]
pct <- round(slices/sum(slices)*100)
lbl=names(slices)
lbls <- paste(lbl, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="")
main=label
xlab=toupper(attribute)
#default label
if(label==""){
main=paste("Dot Plot of", attribute)
}


#define colors
color=col
if(col==""){
color="black"
}

#construct dot
if(type=="map"){
#plot(spplot(object[attribute],main=main,scales = list(draw = TRUE)))
#convert to points
if(symbol==""){
symbol="o"
}
scot_SP_LL <- SpatialPointsDataFrame(coordinates(object), proj4string = CRS(as.character(NA)), data = as(object, "data.frame")[c(attribute)])
 plot(bubble(scot_SP_LL,col=color, main=main,attribute,pch = symbol,scales = list(draw = TRUE)))

}
else{


dotchart(slices, labels=row.names(object),cex=.7,color=color,main=main,xlab=xlab)

}
data=object[[attribute]]
list(data=data,table=slices,source=source,layer=layer,
attribute=attribute,type=type,label=label,color=col)
}

#generic function
dotmap<-function(source,layer='',attribute,type='',label="",col='',symbol='') UseMethod ("dotmap")

#default function
dotmap.default<-function(source,layer='',attribute,type='',label="",col='',symbol='')
{

factor<-indot(source,layer,attribute,type,label,col,symbol)

factor$call<-match.call()

class(factor)<-"dotmap"
factor
}

#summary function
summary.dotmap<-function(object,...)
{
x<-object
cat("Call:\n")
print(x$call)
cat("\n Summary  \n")
print(summary(x$call))


}

#print function
print.dotmap<-function(x,...)
{
cat("Call:\n")
print(x$call)
cat("\nCall: Data\n")
print(x$data)

cat("\n Table  \n")
print(table(x$table))

}

#plot function
plot.dotmap<-function(x,...)
{
dotmap(source=x$source,layer=x$layer,attribute=x$attribute,type=x$type,col=x$color,label=x$label)
}

#graph=dotmap("E:/geokings/advancegis/R-process","KwabibiremAtiwa",attribute='PlantPop',type="map",label="Plant popu")


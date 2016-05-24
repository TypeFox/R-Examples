
#################################################################################################
# function to draw a map
draw.map<-function(a=-10,b=30,c=23,d=65, bbox=NULL, fill=TRUE, col.land=grey(0.5),col.water="white", border=NA,
      detail=FALSE, box=TRUE, axes=FALSE, las=1, cex.axis=0.8, dist.axislab=0.2, whichaxis=c(1:4),
      tck=-0.005, mercator=FALSE, mar=rep(0.5, 4), asp=NA){
# a = longitude of left margin
# b = longitude of right margin
# c = latitude of lower margin
# d = latitude of upper margin
# bbox = an alternative way of specifying the borders, to be consistent with draw.recmap
# fill = if TRUE areas are filled with colors
# col.land = color of land
# col.water = color of water
# border = if NA no borders are painted
# detail = if true, the coast lines are drawn
# box = should the map be boxed
# mercator = if true, mercatorprojection is used
# asp = if set to 1 x and y axes are scaled equally (not recommended if mercator=TRUE)
#-----------------------------------------------------------------------------------
# log of changes
# 07. 3. 2014 fk: data files provided as rda instead of txt files
# 22. 1. 2014 fk: argument asp added so that it is possible to choose this 
#                 option (formerly it was 1 if mercator was TRUE and NA otherwise)
#-----------------------------------------------------------------------------------
if(!is.null(bbox)) {
a <- bbox[1]; b <- bbox[2]; c <- bbox[3]; d <- bbox[4]
}
#data(coastEU, envir = environment())
#data(coastpaleo, envir = environment())
#filename1 <- system.file("data", "coastEU.rda", package = "birdring")
#filename2 <- system.file("data", "coastpaleo.rda", package = "birdring")
#coastpaleo <- read.table(filename2, header=TRUE)
coast<- coastEU #read.table(filename1, header=TRUE)
paleo<-(a<min(coast$x, na.rm=TRUE)|b>max(coast$x, na.rm=TRUE)|c<min(coast$y, na.rm=TRUE)|d>max(coast$y, na.rm=TRUE))
if(paleo) coast<-coastpaleo
if(mercator) coast$y <- mercatorlat(coast$y)
if(mercator) c <- mercatorlat(c)
if(mercator) d <- mercatorlat(d)
par(mar=mar)
plot(coast$x, coast$y, type="n", xlab="", ylab="", axes=FALSE, xlim=c(a, b), ylim=c(c, d), xaxs="i", yaxs="i",asp=asp)
rect(-180,-90,180,90, col=col.water, border=NA)
if(axes){
par(mgp=c(3,dist.axislab, 0))
if(is.element(1, whichaxis)){
  axis1 <- axis(1, tck=0, labels=NA, line=-0.01)
  axis(1, at=axis1, labels=axis1, line=-0.01, las=las, cex.axis=cex.axis, tck=tck)
  }
  
if(is.element(2, whichaxis)){
if(!mercator){ 
  axis2 <- axis(2, tck=0, labels=NA, line=-0.01)
  axis(2, at=axis2, labels=axis2, las=las, cex.axis=cex.axis, tck=tck,  line=-0.01)
  }
if(mercator) axis(2, at=mercatorlat(seq(-60, 70, by=10)),  line=-0.01, labels=seq(-60, 70, by=10), cex.axis=cex.axis, las=las, tck=tck)
}
if(is.element(3, whichaxis)){
  axis3 <- axis(3, tck=0, labels=NA, line=-0.01) 
  axis(3, at=axis3, labels=axis3, las=las, cex.axis=cex.axis, tck=tck, line=-0.01)
  }
if(is.element(4, whichaxis)){ 
if(!mercator){ 
  axis4 <- axis(4, tck=0, labels=NA, line=-0.01)
  axis(4, at=axis4, labels=axis4, las=las, cex.axis=cex.axis, tck=tck, line=-0.01)
  }
if(mercator) axis(4, at=mercatorlat(seq(-60, 70, by=10)), line=-0.01,  labels=seq(-60, 70, by=10), cex.axis=cex.axis, las=las, tck=tck)
}
}

if(!paleo){
  coast$name<-as.character(coast$name)
  coast$name[coast$name=="mainland"&!is.na(coast$name)]<-"aaaa"
  if(fill) {for(ind in levels(factor(coast$name))){
  index<-coast$name==ind&!is.na(coast$name)
  ifelse(is.element(ind, c("caspian", "aral", "baikal", "lagoda","peipus", "vanarn")), col<-col.water, col<-col.land)
  polygon(coast$x[index], coast$y[index], col=col, border=border)
  }}
  }
if(paleo){
  index<-coast$kategorie=="m"
  if(fill) polygon(coast$x[index], coast$y[index], col=col.land, border=border)
  index<-is.element(coast$kategorie,c("s", "i"))
  if(fill){ for(ind in levels(factor(coast$name[index]))){
  index1<-coast$name==ind&!is.na(coast$name)
  polygon(coast$x[index1], coast$y[index1], col=ifelse(coast$kategorie[index1][1]=="s", col.water, col.land), border=border)
  }}
  }

if(detail)lines(coast$x, coast$y)
if(box) box()
}
##################################################################################################


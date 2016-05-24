choropleth.plot <-
function(sp,dem="P0010001",cuts=list("quantile",seq(0, 1, 0.25)),color=list(fun="hsv",attr=list(h = c(.4,.5,.6,.7), s = .6, v = .6, alpha=1)),main=NULL,sub="Quantiles (equal frequency)",legend=list(pos="bottomleft",title="Population Count"),type=NULL,...){
	
color.map<- function(x,dem,y=NULL){
	l.poly<-length(x@polygons)
	dem.num<- cut(dem, breaks=ceiling(do.call(cuts[[1]],list(x=dem,probs=cuts[[2]]))),dig.lab = 6)
	dem.num[which(is.na(dem.num)==TRUE)]<-levels(dem.num)[1]
	l.uc<-length(table(dem.num))
if(is.null(y)){
	col.heat<-do.call(color$fun,color$attr)
}else{
	col.heat<-y
	}
dem.col<-cbind(col.heat,names(table(dem.num)))
colors.dem<-vector(length=l.poly)
for(i in 1:l.uc){
	colors.dem[which(dem.num==dem.col[i,2])]<-dem.col[i,1]
	}
out<-list(colors=colors.dem,dem.cut=dem.col[,2],table.colors=dem.col[,1])
out
}	
	
colors.use<-color.map(sp,sp[[dem]])
col<-colors.use$color

args <- list(x=sp,...,col=colors.use$color)
do.call("plot", args)
title(main=main,sub=sub)
legend(legend$pos,legend=colors.use$dem.cut,fill=colors.use$table.colors,bty="o",title=legend$title,bg="white")
}


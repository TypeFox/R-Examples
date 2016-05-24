#' @rdname pcps
#' @encoding UTF-8
#' @export
plot.pcps<-function(x,display=c("text","points"),groups,showlabel=TRUE,choices=c(1,2),...){
	sco<-scores.pcps(x,choices = choices)
	plot(sco$scores.sites,type="n",ylim=c(min(sco$scores.sites[,2],sco$scores.species[,2],na.rm=TRUE)-0.05, max(sco$scores.sites[,2],sco$scores.species[,2],na.rm=TRUE)+0.05),xlim=c(min(sco$scores.sites[,1],sco$scores.species[,1],na.rm=TRUE)-0.05,max(sco$scores.sites[,1],sco$sco[,1],na.rm=TRUE)+0.05),...)
	if(display=="text"){
		text(sco$scores.sites,labels=rownames(sco$scores.sites),...) 
	}
	if(display=="points"){
		points(sco$scores.sites,...)
	}
	vegan::ordispider(sco$scores.species,groups=groups,label=showlabel,...)
	if(showlabel){
		g1<-ifelse(table(groups)==1,1,0)
		g1_groups<-names(g1)[g1==1]
		if(sum(g1)>0){
			for(i in 1:sum(g1)){
				position<-which(groups==g1_groups[i])
				vegan::ordilabel(sco$scores.species[position,],labels=groups[position],...)
			}	
		}
	}
}
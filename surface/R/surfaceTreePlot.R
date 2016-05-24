surfaceTreePlot <-
function(tree,hansenfit,cols=NULL,convcol=TRUE,labelshifts=FALSE,...){

	fit<-hansenfit$fit[[1]]
	otree<-as(fit,"data.frame")
	otree<-data.frame(otree,shifts=rep(NA,length(otree$nodes)))
	otree$shifts[match(names(hansenfit$savedshifts),otree$nodes)]<-1:length(hansenfit$savedshifts)
	ntip<-(dim(otree)[1]+1)/2;nnode<-ntip-1
	otree2<-otree[match(c(tree$tip.label,tree$node.label),otree$labels),]
	otree2<-otree2[tree$edge[,2],]

if(length(cols)==1)cols<-rep(cols,length(unique(hansenfit$savedshifts)))
if(is.null(cols)){
		xx<-summary(factor(hansenfit$savedshifts))
	if(convcol){
		cols<-character(length(xx))
		cols[xx>1]<-rainbow(sum(xx>1))
		if(any(xx==1))	
			cols[xx==1]<-c("black",grey(seq(0.7,0.3,length.out=sum(xx==1)-1)))
	}else{
		cols<-c("black",rainbow(length(xx)-1))
	}	}
	edgecols<-cols[as.numeric(factor(otree2[,5]))]

	plot(tree,edge.color=edgecols,...)
	if(labelshifts){
	nodelabels(node=tree$edge[,2][which(!is.na(otree2$shifts))],bg="white", text=otree2$shifts[!is.na(otree2$shifts)],cex=0.6,frame="circle")
		}
	}

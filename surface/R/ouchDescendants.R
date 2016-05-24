ouchDescendants <-
function(node,otree){
	newdesc<-desc<-which(otree@ancestors==node)
	if(length(desc)==0)return(NA)
	while(length(newdesc)>0){
		newdesc<-which(otree@ancestors%in%newdesc)
		desc<-c(desc,newdesc)	
	}
	return(desc)
	}

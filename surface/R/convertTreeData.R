convertTreeData <-
function(tree,dat){
	otree<-ape2ouch(tree,scale=F)
	odata<-dat[as.character(as(otree,"data.frame")$labels),,drop=F]
	rownames(odata)<-otree@nodes
	return(list(otree,odata))
	}
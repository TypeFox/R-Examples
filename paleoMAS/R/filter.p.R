filter.p <-
function(x,rare=TRUE,presen=1,persist=0.05)
{
{
	if(rare==TRUE){
		filter<-matrix(nrow=ncol(x),ncol=3,
			dimnames=list(colnames(x),
			c("n","n over minimum presence","quality")))
		ifelse(x>0,1,0)->pres.abs
		apply(pres.abs,2,sum)->filter[,1]
		ifelse(x>presen,1,0)->presence
		apply(presence,2,sum)->persistence
		persistence->filter[,2]
		a<-as.integer(persist*nrow(x))
		ifelse(filter[,2]>=a,1,0)->filter[,3]
		which(as.matrix(persistence)>=a)->filter1
		x[,filter1]->filtered
		matrix(nrow=3,ncol=1,dimnames=list(c("%","minimum n",
			"# taxa"),"value"))->result
		result[1,]<-presen
		result[2,]<-a
		result[3,]<-ncol(filtered)
		}
	else{
		ifelse(x>0,1,0)->pres.abs
		apply(pres.abs,2,sum)->filter1
		a<-as.integer(persist*nrow(x))
		names(which(filter1<=a))->filter2
		x[,filter2]->filtered
		matrix(nrow=2,ncol=1,dimnames=list(c("maximum n",
			"#taxa"),"value"))->result
		result[1,]<-a
		result[2,]<-ncol(filtered)
		}
}
if(rare==TRUE){
	results<-list(filtered,filter,result)
	names(results)<-c("filtered","filter","result")
	return(results)
	}
else{
	results<-list(filtered,result)
	names(results)<-c("filtered","result")
	return(results)
	}

}


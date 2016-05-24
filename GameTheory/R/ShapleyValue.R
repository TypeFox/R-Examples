ShapleyValue <-
function(x,Names=NULL){
	
	#####
	as.matrix(x$Lex)->z
	as.vector(z)->z
	#####
	coalitions<-set.func(c(0,z))
	SV<-Shapley.value(coalitions)
	SV<-as.matrix(SV)
	rownames(SV)<-Names
	barplot(SV,col=c(
		"#63A5DB",
		"#005B9A",
		"#734A75",
	    "#B80606",
	    "#E38030"),axes=TRUE,family="Times",main="Shapley Value",beside=TRUE,names.arg=Names,space=0.25)
    colnames(SV)<-"Shapley Value"
    axis(2,family="Times")
    Output<-list(SV=SV)
    class(Output)<-"ShapleyValue"
	return(Output)
}

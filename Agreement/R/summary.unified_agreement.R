summary.unified_agreement <-
function(object,...){
	if (object$Input[2]>1){
		print=rbind(object$Intra,object$Inter,object$Total)
		print=cbind(rbind("Intra"," "," ","Inter"," "," ","Total"," "," "),print)
		colnames(print)=c("Type","Statistics","CCC","Precision","Accuracy","TDI","CP","RBS");
		rownames(print)=c("","","","","","","","","");
		if (object$Inputs[3]==0)	print=print[,1:5]
		print(print,quote=FALSE, na.print=".",...);
	}else{
		print=object$Stat;
		colnames(print)=c("Statistics","CCC","Precision","Accuracy","TDI","CP","RBS");
		rownames(print)=c("","","");
		if (object$Inputs[3]==0)	print=print[,1:4]
		print(print,quote=FALSE, na.print=".",...);
	}
}


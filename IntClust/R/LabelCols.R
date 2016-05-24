LabelCols <- function(x,Sel1,Sel2=NULL,col1,col2=NULL) {
	colfunc=function(x,Sel1,Sel2,col1,col2){
		if (x %in% Sel1){
			return(col1)
		}
		else if(x %in% Sel2){
			return(col2)
		}
		else{
			return("black")
		}
	}
	
	if (stats::is.leaf(x)) {
		## fetch label
		label <- attr(x, "label") 
		## set label color to red for SelF, to black otherwise
		attr(x, "nodePar") <- list(pch=NA,lab.col=colfunc(label,Sel1,Sel2,col1,col2),lab.cex=0.9,font=2)
		attr(x, "edgePar") <- list(lwd=2,col=colfunc(label,Sel1,Sel2,col1,col2))
	}
	return(x)
}

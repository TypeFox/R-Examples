summary.tir_iir <-
function(object,...){
	cat("TIR:\n");
	cat("------------------------------------------\n");
	for (i in 1:length(object$Model$TIR)){
		cat("  ");cat(object$Model$TIR[i]);cat(":\n");
		print <- cbind(round(object$TIR$TIR[i],3), round(object$TIR$TIR_upper[i],3), object$TIR$TIR_a);
		colnames(print) <- c("Estimate:","95% Conf. Limit:","Compared to:");
		rownames(print) <- "";
		print(print,quote=FALSE, na.print=".",...);
		cat("------------------------------------------\n");
	}
	cat("IIR:\n");
	cat("------------------------------------------\n");
	for (i in 1:length(object$Model$IIR)){
		cat("  ");cat(object$Model$IIR[i]);cat(":\n");
		IIR_conf <- paste("(",round(object$IIR$IIR_lower[i],3),",",round(object$IIR$IIR_upper[i],3),")",sep="")
		print <- cbind(round(object$IIR$IIR[i],3), IIR_conf, 1.00);
		colnames(print) <- c("Estimate:","95% Conf. Limit:","Compared to:");
		rownames(print) <- "";
		print(print,quote=FALSE, na.print=".",...);
		cat("------------------------------------------\n");
	}
}


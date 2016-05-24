print.SJ.dag <-
function(x,...){	
	cat("Object of class \"SJ.dag\"\n")
	cat("Contains the graphs:\n")
	if(class(x$graph) == "list"){ 
		print(data.frame("lambda" = signif(x$lambda,2), "Edges" = unlist(lapply(x$graph,ecount)), "bic" = x$bic))
		}	
	else print(data.frame("lambda" = signif(x$lambda,2), "Edges" = ecount(x$graph), "bic" = x$bic))
	cat("Call:\n\t");print(x$cl)
	invisible(x)
}

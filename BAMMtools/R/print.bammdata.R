print.bammdata = function(x, ...)
{
	print.phylo(as.phylo.bammdata(x));
	nsamples <- length(x$eventData);
	cat(paste("\nPosterior samples:", nsamples,"\n\n"));
	cat("List elements:\n");
	cat("\t",names(x)[1:10]);
	cat("\t",names(x)[11:length(x)]);
	cat('\n');
}

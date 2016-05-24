slice_show <- function(slices_obj, main="Counts in each slice", 
                       xlab="Slices", ylab="Percentage")
{
	s_n <- nrow(slices_obj)
	l_n <- ncol(slices_obj)
	s_v <- rownames(slices_obj)
	l_v <- colnames(slices_obj)[-l_n]
	if(s_n == 1){
		smat <- matrix(slices_obj[, -l_n], nrow = 1)
	}else{
		smat <- slices_obj[, -l_n]
	}
	ns <- slices_obj[, l_n]
	l_n <- l_n - 1
	s_h <- paste(rep("Slice", s_n), 1:s_n, sep = " ")
	slice <- paste(rep(s_h, l_n), " (n_", rep(1:s_n, l_n), " = ", rep(ns, l_n), ")", sep="")
	level <- rep(l_v, each = s_n)
	count <- as.vector(smat)
	percentage <- as.vector(smat) / rep(ns, l_n)
	pmat <- data.frame(slice, level, count, percentage)
	gsp <- ggplot(pmat, aes(x = level, y = percentage, ymax = 1, fill = level)) + 
		geom_bar(stat = "identity", position="stack") + 
		scale_y_continuous(labels = percent, limits=c(0, 1)) + 
		scale_fill_hue(name="Types") + 
		ggtitle(main) + xlab(xlab) + ylab(ylab) + 
		geom_text(aes(label=count), position="stack", vjust=-0.5, colour="black", size = 3) + 
		facet_grid(. ~ slice) + 
		theme_bw()
	return(gsp)
}

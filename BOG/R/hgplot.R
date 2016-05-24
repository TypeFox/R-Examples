hgplot <-
function(x){

    stat=x[[1]]

	mgn <- par(oma=c(0.1,0.1,0.1,0.1))

	legend1="# of the observed (o)"
	legend2="# of the expected (e)"
    graphs=internalhgplot(stat)
    plotmat=graphs$plotmat
    label=graphs$label
    margin.y=graphs$margin.y
    e.total=graphs$e.total
    o.total=graphs$o.total
    barplot(plotmat,names.arg=label,beside=TRUE,col=c("gold","brown"),ylim=c(0,margin.y),ylab="Counts",
		xlab="COGs with Adjusted p-values < 0.1",legend.text=c(legend1,legend2),cex.names=1.0)
		par(mgn)
}

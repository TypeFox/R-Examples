davesil.default<- function(ddist,o.hclr,o.relgr,...)  {
	o.davesil<- dsil(ddist,o.hclr,o.relgr)
	o.davesil$call<- match.call()
	cat("Call:\n") 
	class(o.davesil) <- "davesil"
	print(o.davesil$call)
	cat("\nAverage silhouette widths:\n",round(o.davesil$meansilwidth,digits=3),"\n")
	o.davesil
}

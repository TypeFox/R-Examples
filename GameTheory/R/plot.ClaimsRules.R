plot.ClaimsRules<-
function(x,y, ...){
#assignInNamespace("plot.ClaimsRules", plot.ClaimsRules, ns = asNamespace("stats"))
nam<-rownames(x)
N<-paste("Allocations for ",nam[y])
	barplot(x[y,],
	col=c(
		"#FFFDED",
		"#63A5DB",
		"#005B9A",
		"#734A75",
	   "#B80606",
	   "#E38030"),
	family="Times",
	main=N)
}
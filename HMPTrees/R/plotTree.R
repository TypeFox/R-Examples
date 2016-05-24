plotTree <-
function(trees, colors, divisions, main, sub, 
showTipLabel=TRUE, showNodeLabel=FALSE, displayLegend=TRUE){

if(missing(trees))
stop("At least one valid tree of type 'phylo' is required inside a list.")

if(missing(sub))
sub <- ""

genMain <- missing(main)
if(genMain)
main <- NULL
passMain <- length(main) == length(trees)

if(displayLegend)
displayLegend(colors, divisions)

par(cex=.75)

for(i in 1:length(trees)){
tr <- trees[[i]]
bs <- getBranchSizes(tr, colors, divisions)

if(passMain){
main2 <- main[i]
}else if(genMain){
main2 <- names(trees)[i]
}else{
main2 <- main[1]
}

ape::plot.phylo(tr, type="r", root.edge=FALSE, edge.color=bs$edgecol, edge.width=bs$edgewid, 
show.tip.label=showTipLabel, show.node.label=showNodeLabel, main=main2, sub=sub)
}
}

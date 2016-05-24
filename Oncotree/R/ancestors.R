"ancestors" <-
function(otree, vertex){
	v.num <- match(vertex, otree$parent$child)
	curr.anc <- vertex
	curr.vertex <- v.num
	anc.list <- vertex
	while (curr.anc!="Root"){
		curr.vertex <- otree$parent$parent.num[curr.vertex]
		curr.anc <- otree$parent$child[curr.vertex]
		anc.list <- c(anc.list, curr.anc)
	}
	anc.list
}


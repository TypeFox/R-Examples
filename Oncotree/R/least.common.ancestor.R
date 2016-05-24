"least.common.ancestor" <-
function(otree, v1, v2){
	anc1 <- ancestors(otree, v1)
	anc2 <- ancestors(otree, v2)
	overlap <- na.omit(match(rev(anc1),rev(anc2))) #"Root" will always match
	lca <- rev(anc1)[max(overlap)]
	lca
}


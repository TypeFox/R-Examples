`matchIndexes` <-
function(m, string) {
	ans = 1:length(m)
	for (i in 1:length(m)) {
		ans[i] =match(m[i],string)
	}
	return (ans)
}


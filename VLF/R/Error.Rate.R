Error.Rate <-
function(single, shared, spec,seqlength){
	pos2.Sum.single = 0
	pos2.Sum.shared = 0
	pos2.Sum = 0
	i = 2
	while(i <= seqlength){
		pos2.Sum.single = pos2.Sum.single + single[i]
		#pos2.Sum.shared = pos2.Sum.shared + shared[i]
		#pos2.Sum = pos2.Sum + single[i] + shared[i]
		i = i+3
	}
	Rate.single <- pos2.Sum.single/((216*spec) - pos2.Sum.shared)
	Rate.shared <- pos2.Sum.shared/((216*spec) - pos2.Sum.single)
	Rate <- pos2.Sum/(216*spec)
	return(rbind(Rate.single, Rate.shared, Rate))
}

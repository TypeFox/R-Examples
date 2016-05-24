summary.ibd <-
function(segments, merged=TRUE, verbose=TRUE) {
	total = attr(segments, 'total_map_length_Mb')
	disease = all(sapply(segments,function(r) any(r[,'disreg']==1)))

	res <- sapply(segments, function(run) {
		if(merged) r = mergeRuns(run) else r = run
		count = nrow(r)
		all.lengths = as.numeric(r[,'end']-r[,'start'])
		len = sum(all.lengths)
		aver = ifelse(count>0, len/count, 0)
		longest = max(c(0, all.lengths))
		stats = c(count.all = count, fraction.all = len/total*100, average.all = aver, longest.all = longest)
		
		if(disease) {
			r = run[run[,'disreg']==0, , drop=F]
			if(merged) r = mergeRuns(r)
			count = nrow(r)
			rlengths = as.numeric(r[,'end']-r[,'start'])
			len = sum(rlengths)
			aver = ifelse(count>0, len/count, 0)
			longest = max(c(0, rlengths))
			stats = c(stats, count.rand = count, fraction.rand = len/total*100, average.rand = aver, longest.rand = longest)
		
			dis = run[,'disreg']==1
			if(sum(dis)!=1) {print(run); warning("More or less than one disease region")}
			else {
				disrun = run[dis, ]
				dis.length = as.numeric(disrun['end']-disrun['start'])
				stats = c(stats, length.dis = dis.length, rank.dis =rank(-c(dis.length, all.lengths), ties.method="first")[1])
			}
		}
		stats
	})
	no.seg = mean(res['count.all',]==0)
	if (verbose) {
        summary = round(c(rowMeans(res), 'zeroprob' = no.seg*100), 3)
        print(summary)
    }
	invisible(res)
}

	
#merges overlapping and adjacent segments
mergeRuns <- function(r) {
	if (nrow(r)<2) return(r)
	r = r[order(r[,'chrom'], r[,'start'], r[,'end']),]
	mergelist = lapply(unique(r[, 'chrom']), function(i) {
		m = rangeUnion(list(r[ r[, 'chrom']==i, 2:3]))
		m = cbind(i, m)
		colnames(m) = c('chrom', 'start', 'end')
		m
	})
	do.call(rbind, mergelist)
}

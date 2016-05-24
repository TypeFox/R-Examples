# The underlying functions for blending a single landscape

blendOnce = function(
	native, 
	exotic, 
	landscape.name = "",
	warn = FALSE # loess gives annoying warnings; disabled by default
){
	# The main ecoblender workhorse. Returns a blended.landscape object
	
	
	verifyLandscapes(native, exotic)
	
	
	
	# Vegan is fast, but we need to avoid comparing the native species
	# to one another every single time we look at the effect of a new
	# exotic species.  We do this by finding the numerators and 
	# denominators for each pair of native sites, then adding one as 
	# needed for a given exotic species.
	
	# native numerators and denominators
	tops.and.bottoms = findTopsAndBottoms(native) 
	
	# native similarity
	J.Bar = with(tops.and.bottoms, averageRatio(tops, bottoms))
	J.Star = with(tops.and.bottoms, ratioOfAverages(tops, bottoms))
	
	# Change in similarity by exotic species
	similarity.by.exotic = findSimilarityByExotic(
		tops.and.bottoms$tops,
		tops.and.bottoms$bottoms,
		exotic
	)
	delta.J.Bars = similarity.by.exotic$J.Bar - J.Bar
	delta.J.Stars = similarity.by.exotic$J.Star - J.Star
	
	
	
	# threshold and nadir
	crits = findCriticalPoints(exotic, delta.J.Bars, warn)
	
	
	# vegan's value of JBar for full native+exotic
	total.J.Bar = jbar(rbind(native, exotic))
	
	# JStar for full native+exotic
	total.J.Star = jstar(rbind(native, exotic))
	
	
	# x axis for building the "scoop" shaped graphs: 
	scoop.occupancies = seq(0, ncol(exotic), .01)
	
	
	exotic.jstar.tops = choose(round(rowMeans(exotic) * ncol(exotic)), 2)
	exotic.jstar.bottoms = choose(ncol(exotic), 2) - choose(round((1 - rowMeans(exotic)) * ncol(exotic)), 2)
	
	`/`(
		sum(tops.and.bottoms$tops) + exotic.jstar.tops,	
		sum(tops.and.bottoms$bottoms) + exotic.jstar.bottoms
	) - jstar(native)
	
	
	#####
	# Compose output material
	#####
	basic.results = data.frame(
		J.Bar = J.Bar,
		J.Star = J.Star,
		delta.J.Bar = total.J.Bar - J.Bar,
		delta.J.Star = total.J.Star - J.Star,
		R2 = 1 - var(delta.J.Bars - delta.J.Stars)/var(delta.J.Stars),
		threshold = crits$threshold,
		p.Star = `/`(
			(1 + J.Star * (2 * ncol(native) - 1)),
			(ncol(native) * (1 + J.Star))
		),
		nadir = crits$nadir,
		row.names = landscape.name
	)
	
	
	structure(
		c(
			name = landscape.name,
			basic.results,
			list(
				results.table = basic.results,
				species.delta.table = data.frame(
					species = row.names(exotic),
					occupancy = rowSums(exotic)/ncol(exotic),
					delta.J.Bars =  delta.J.Bars,
					delta.J.Stars = delta.J.Stars
				),
				scoop = data.frame(
					x = scoop.occupancies / ncol(exotic),
					y = with(
						tops.and.bottoms, 
						findJ.Stars(
							scoop.occupancies,
							tops,
							bottoms,
							ncol(exotic),
							FALSE
						) - J.Star
					)
				),
				native = native,
				exotic = exotic
			)
		),
		class = "blended.landscape"
	)
}


##################
# helper functions
##################


sitePairs = function(num.sites){
	t(combn(num.sites, 2))
}

findTopsAndBottoms = function(landscape){
	# How many species are shared (tops) and present (bottoms) across
	# all n-choose-2 comparisons?
	
	compare = function(site.pair, landscape){
		left = landscape[ , site.pair[1]]
		right = landscape[ , site.pair[2]] 
		
		
		c(
			sum(left & right), 
			sum(left | right)
		)
	}
	
	out = as.data.frame(
		t(
			apply(
				sitePairs(ncol(landscape)),
				1,
				function(site.pair) compare(site.pair, landscape)
			)
		)
	)
	
	names(out) = c("tops", "bottoms")	
	
	out
}


findSimilarityByExotic = function(tops, bottoms, exotic.landscape){	
	# For each exotic species, put the number of times it's shared on top
	# and the number of times it appears on the bottom.  Then add
	# it to the native totals for each site pair.  Do this for all
	# species (rows).
	
	site.pairs = sitePairs(ncol(exotic.landscape)) 
	
	compareOneSpecies = function(row){
		top = (row[site.pairs[,1]] & row[site.pairs[,2]]) + tops
		bottom = (row[site.pairs[,1]] | row[site.pairs[,2]]) + bottoms
		
		c(
			J.Bar = averageRatio(top, bottom),
			J.Star = ratioOfAverages(top, bottom)
		)
	}
	
	as.data.frame(t(apply(exotic.landscape, 1, compareOneSpecies)))
}



ratioOfAverages = function(tops, bottoms){
	# summing instead of averaging should avoid some small rounding error 
	sum(tops) / sum(bottoms)
}

averageRatio = function(tops, bottoms) mean(tops / bottoms)



findCriticalPoints = function(exotic.landscape, delta.J.Bars, warn = FALSE){
	# use local regression to smooth the data and find where it has a
	# minimum and crosses the x-axis
	
	# vector of exotic occupancy rates
	pe = rowSums(exotic.landscape) / ncol(exotic.landscape)
	
	# Turn off warnings by default; they're confusing and uninformative
	if(warn){}else{
		prior.warn.level = getOption("warn")
		on.exit(options(warn = prior.warn.level))
		options(warn = -1)
	}
	
	# Smooth the data (adding the 0,0 should make loess work better
	# on landscapes with ~4 sites)
	loess.scoop = loess(c(0, delta.J.Bars) ~ c(0, pe))
	
	
	
	if(any(is.nan(loess.scoop$fitted))){
		return(list(threshold = NA, nadir = NA))
	}
	
	nadir = optimize(
		function(x) predict(loess.scoop, x), 
		c(0, max(pe))
	)$minimum
	
	
	if(any(delta.J.Bars > 0) & any(delta.J.Bars < 0)){
		threshold = optimize(
			function(x) abs(predict(loess.scoop, x)), 
			c(nadir, max(pe)) # threshold must be to the right of nadir
		)$minimum
	}else{
		threshold = NA
	}
	
	list(threshold = threshold, nadir = nadir)
}



findJ.Stars = function(occupied, tops, bottoms, n.sites, avg){
	# find J.Star from occupancy and/or native tops/bottoms
	star.tops = function(p, tops, sites){
		sum(tops) + choose(p * sites, 2)
	}
	star.bottoms = function(p, tops, bottoms, sites){
		sum(bottoms) + choose(sites, 2) - choose((1 - p) * sites, 2)
	}
	
	
	p = occupied / n.sites
	
	star.tops = star.tops(p, tops, n.sites)
	star.bottoms = star.bottoms(p, tops, bottoms, n.sites)
	
	if(avg){
		mean(star.tops)/mean(star.bottoms)
	}else{
		star.tops/star.bottoms	
	}
	
	
}





verifyLandscapes = function(native, exotic){
	if(any(is.na(native))){
		stop("the native data table has empty entries.  All entries must be zeros or ones.")
	}
	if(any(is.na(exotic))){
		stop("the exotic data table has empty entries.  All entries must be zeros or ones.")
	}
	
	
	if(!all(unlist(native) %in% c(0, 1))){
		stop("the native landscape contains entries that are not 0 or 1")
	}
	if(!all(unlist(exotic) %in% c(0, 1))){
		stop("the exotic landscape contains entries that are not 0 or 1")
	}
	
	
	if(!all(colnames(native) == colnames(exotic))){
		stop("the sites in the native landscape have different names than in the exotic landscape")
	}
	
	if(any(colSums(native) == 0)){
		stop("one or more sites has no native species.  Similarity is not defined for these sites.")
	}
	
	if(ncol(native) < 2){
		stop("landscapes must have at least two sites for mean similarity to be defined")
	}
	
	
}


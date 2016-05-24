recodeBinary = function(binary.vec, k)
	{
	#Constants and vectors used in the function
	m = length(binary.vec)
	output.vec = rep(0, m)

	#Define the entries of output.vec
	#Move to the right of k
	temp.spot = k
	temp.value = binary.vec[temp.spot]
	while ((temp.value > 0) && (temp.spot < m))
		{
		output.vec[temp.spot] = 1
		temp.spot = temp.spot + 1
		temp.value = binary.vec[temp.spot]
		}
	if ((temp.spot == m) && (temp.value > 0))
		{
		output.vec[temp.spot] = 1
		}

	#Move to the left of k
	temp.spot = k
	temp.value = binary.vec[temp.spot]
	while ((temp.value > 0) && (temp.spot > 1))
		{
		output.vec[temp.spot] = 1
		temp.spot = temp.spot - 1
		temp.value = binary.vec[temp.spot]
		}
	if ((temp.spot == 1) && (temp.value > 0))
		{
		output.vec[temp.spot] = 1
		}
	return(output.vec)
	}

makeCytoband = function(marker.data, annot.file, reformat.annot = FALSE)
	{
	#Reformat annot.file, if necessary
	if (reformat.annot == TRUE)
        {
        autosome.stop = min(which(as.matrix(annot.file[,1]) == "chrX")) - 1
        annot.file = annot.file[1:autosome.stop,]
        new.chr.col = unlist(strsplit(as.matrix(annot.file[,1]), "r", fixed = TRUE))[2*c(1:dim(annot.file)[1])]
        annot.file[,1] = as.data.frame(new.chr.col)
        }
	
	#Constants and vectors used in the function
	chr.vec = intersect(unique(as.matrix(marker.data[,1])), c(1:22))
	output = c()

	#Find the boundaries of the p and q arms for each chromosome, then append the chromosome arm information to output.vec
	for (i in chr.vec)
		{
		temp.annot.rows = which(as.matrix(annot.file[,1]) == i)
		temp.annot = annot.file[temp.annot.rows,]
		temp.pq.vec = substring(as.matrix(temp.annot[,4]), 1, 1)
		temp.q.start = temp.annot[min(which(temp.pq.vec == "q")), 2]
		
		temp.marker.rows = which(marker.data[,1] == i)
		temp.marker = marker.data[temp.marker.rows,]
		temp.num.markers = dim(temp.marker)[1]
		temp.num.p = sum(temp.marker[,2] < temp.q.start)
		temp.num.q = temp.num.markers - temp.num.p
		
		temp.output = c(rep("p", temp.num.p), rep("q", temp.num.q))
		output = c(output, temp.output)
		}

	#Return the output
	return(output)
	}

findNull = function(x, num.perms, random.seed = NULL)
	{
	#Set random seed, if apppropriate
	if (length(random.seed) > 0)
		{
		set.seed(random.seed)
		}

	#Constants and vectors used in the function
	n = dim(x)[1]
	m = dim(x)[2]
	shift.vec = rep(23, num.perms)

	#Perform the cyclic shift procedure, and record the minimum and maximum column sum after each
	#cyclic shift.
	for (j in (1:num.perms))
		{
		cutpoints = sample(c(1:m), size = n, replace = TRUE)
	 	perm.x = matrix(0, n, m)
	 	for (i in (1:n))
			{
			index = cutpoints[i]
			if (index == 1)
				{
				second = c()
				}
			if (index > 1)
				{
				second = x[i, 1:(index - 1)]
				}			
			first = x[i, index:m]
	  		perm.x[i,] = c(first, second)
	  		}
		perm.col.sums = colSums(perm.x)
		shift.vec[j] = max(colSums(perm.x))		
 		}

	#Define and return the output
	return(shift.vec)
	}

peeling = function(x, marker.data, cytoband, k)
	{
	n = dim(x)[1]
 	m = dim(x)[2]
	rowmeans = rowMeans(x)	
	grandmean = mean(rowmeans)
	
	#This section of code creates exceed, a binary n x chr.m matrix, where chr.m is the number of markers
	#in the chromosome containing k, the marker where the peeling procedure begins.  An entry of exceed is 1
	#if the corresponding entry of x belongs to the peak interval
	chr = marker.data[k, 1]
	chr.start = min(which(marker.data[,1] == chr))
	chr.end = max(which(marker.data[,1] == chr))
	chr.m =  chr.end - chr.start + 1
	chr.k = k - chr.start + 1
	chr.data = x[,chr.start:chr.end]
	col.means = colMeans(chr.data)
	col.indicator = as.numeric(col.means > grandmean)
	col.indicator = recodeBinary(col.indicator, chr.k)
	exceed = rep(1, n) %o% col.indicator
	#which.band = substr(cytoband[k], 1, 1)
	which.band = cytoband[k]
	#band.vec = substr(cytoband[chr.start:chr.end], 1, 1)
	band.vec = cytoband[chr.start:chr.end]
	band.ind.vec = as.numeric(band.vec == which.band)
	band.matrix = rep(1, n) %o% band.ind.vec
	exceed = exceed * band.matrix
	
	#This section of code finds the markers that make up the peak interval
	peak.interval = as.numeric(colSums(exceed) > 0)
	chr.left.side = min(which(peak.interval == 1))
	chr.right.side = max(which(peak.interval == 1))
	left.side = chr.left.side + chr.start - 1
	right.side = chr.right.side + chr.start - 1

	#These will be used below
 	columnmean = mean(chr.data[,chr.k])
	
	#This section of code creates exceed.running, a binary n x chr.m matrix.  An entry of exceed.running
	#is 1 if the corresponding entry of x will be peeled, otherwise it is zero.
 	exceed.running = matrix(0, n, chr.m)
 	exceed.running[,chr.k] = as.numeric(chr.data[,chr.k] > rowmeans)
	if (chr.k > 1)
		{
 		for (j in (chr.k - 1):1)
			{
  			exceed.running[,j] = exceed.running[,(j + 1)] * exceed[,j] * as.numeric(chr.data[,j] > rowmeans)
  			}
		}
	if (chr.k < chr.m)
		{
 		for (j in (chr.k + 1):chr.m)
			{
  			exceed.running[,j] = exceed.running[,(j - 1)] * exceed[,j] * as.numeric(chr.data[,j] > rowmeans)
  			}
		}

	#This section of code implements a multiplicative-based shrinkage approach at column
	not.imputed = chr.data * (1 - exceed.running)
	to.be.imputed = chr.data * exceed.running
	num.term1 = n * grandmean
	num.term2 = sum(not.imputed[,chr.k])
	num.term3 = sum(exceed.running[,chr.k] * rowmeans)
	den.term = sum(to.be.imputed[,chr.k]) - num.term3
	tau = (num.term1 - num.term2 - num.term3)/den.term
	imputed = tau * (to.be.imputed - (matrix(rowmeans, n, chr.m) * exceed.running)) + 
		(matrix(rowmeans, n, chr.m) * exceed.running)
	final.matrix = imputed + not.imputed
		
	#Create the output
	x[,chr.start:chr.end] = final.matrix
	output.list = list(x, c(left.side, right.side))	
		
	#Return the output
	return(output.list)
	}

detailedLook = function(x, marker.data, annot.file, num.perms, num.iters, gain.loss = "gain", reformat.annot = FALSE, random.seed = NULL)
	{
	#Constants and vectors used in the function
	n = dim(x)[1]
 	m = dim(x)[2]
	r = dim(marker.data)[2]
	small.marker.data = as.matrix(marker.data[,1:2])
	chrom.vec = small.marker.data[,1]
	cytoband = makeCytoband(small.marker.data, annot.file,reformat.annot)
	marker.matrix = as.matrix(marker.data)
	gain.loss.ind = as.numeric(gain.loss == "gain") - as.numeric(gain.loss == "loss")
	
	#Perform the Detailed Look procedure
	data.matrix = x * gain.loss.ind
	output.matrix = c()
	
	for (i in (1:num.iters))
		{
		null.dist = findNull(data.matrix, num.perms, random.seed)
		col.sums = colSums(data.matrix)
		obs.max = max(col.sums)
		k = which.max(col.sums)
		p.val = min(mean(obs.max < null.dist) + 1/num.perms, 1)	
		peeling.data = peeling(data.matrix, marker.matrix, cytoband, k)
		data.matrix = peeling.data[[1]]
		interval = peeling.data[[2]]
		output.matrix = rbind(output.matrix, c(marker.matrix[k,], k, p.val, marker.matrix[interval[1], 2], marker.matrix[interval[2], 2]))
		}
	colnames(output.matrix) = c(colnames(marker.data), "Marker", "p-Val", "Peak Int. (L)", "Peak Int. (R)")

	#Return output
	return(output.matrix)
	}
	
quickLook = function(x, marker.data, annot.file, num.perms, num.iters, gain.loss = "gain", reformat.annot = FALSE, random.seed = NULL)
	{
	#Constants and vectors used in the function
	n = dim(x)[1]
 	m = dim(x)[2]
	r = dim(marker.data)[2]
	small.marker.data = as.matrix(marker.data[,1:2])
	chrom.vec = small.marker.data[,1]
	cytoband = makeCytoband(small.marker.data, annot.file, reformat.annot)
	marker.matrix = as.matrix(marker.data)
	gain.loss.ind = as.numeric(gain.loss == "gain") - as.numeric(gain.loss == "loss")
	
	#Perform the Detailed Look procedure
	data.matrix = x * gain.loss.ind
	output.matrix = c()
	
	null.dist = findNull(data.matrix, num.perms, random.seed)
	for (i in (1:num.iters))
		{
		col.sums = colSums(data.matrix)
		obs.max = max(col.sums)
		k = which.max(col.sums)
		p.val = min(mean(obs.max < null.dist) + 1/num.perms, 1)	
		peeling.data = peeling(data.matrix, marker.matrix, cytoband, k)
		data.matrix = peeling.data[[1]]
		interval = peeling.data[[2]]
		output.matrix = rbind(output.matrix, c(marker.matrix[k,], k, p.val, marker.matrix[interval[1], 2], marker.matrix[interval[2], 2]))
		}
	colnames(output.matrix) = c(colnames(marker.data), "Marker", "p-Val", "Peak Int. (L)", "Peak Int. (R)")

	#Return output
	return(output.matrix)
	}

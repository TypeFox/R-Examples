make.example.files <- function() {
	file.copy(system.file('data.dat',package='DetSel'),'data.dat')
	file.copy(system.file('data.gen',package='DetSel'),'data.gen')
}

GetData <- function() {
	.C('GetData',PACKAGE = 'DetSel')
}

SimulDiv <- function(p1,p2,n1,n2) {
	.C('SimulDiv',
	as.integer(p1),
	as.integer(p2),
	as.integer(n1),
	as.integer(n2),
	PACKAGE = 'DetSel')
}

read.genepop.file <- function(infile) {
	if(!file.exists(infile)) {
		stop(paste('The file ',infile,' does not exist. Check the input file name...',sep = ''))
	}
	connection <- file(infile, 'r', blocking = FALSE)
	test1 <- readLines(connection)
	close(connection)
	connection <- file(infile, 'r', blocking = TRUE)
	test2 <- readLines(connection,warn = FALSE)
	close(connection)
	if (length(test1) != length(test2)) {
		cat('\n', file = infile, append = TRUE)	# add an empty line to the end of the file
	}
	data <- readLines(infile)					# read the data file (line by line)
	data <- gsub('\t',' ',data)                 # remove the tabulations
	data <- data[-1]							# remove the comments line
	data <- data[!(data == '')]					#? remove any blank line
	lpop <- grep('pop',tolower(data))			# give the positions of 'pop' items
	nbrpop <- length(lpop)						# give the number of populations
	nloci <- (min(lpop) - 1)					# give the number of loci
	last <- length(data) + 1					# give the last line of the file
	data <- data[-lpop]							# remove the 'pop' items
	data <- data[-c(1:nloci)]					# remove the list of loci in the header
	n <- (c(lpop[-1],last) - lpop) - 1			# give the sample size of each population
	o <- c(0,cumsum(n))							# give the cumulative sample sizes over all populations
	m <- matrix(nrow = sum(n),ncol = (nloci + 1)) # this is the data matrix
	cpt <- 0									# this will count the lines in the data matrix
	for (i in 1:nbrpop) {						# loop over populations
    	m[((o[i] + 1):o[i + 1]),1] <- rep(i,n[i]) # give the sub-matrix of data for the jth population 
		for (j in 1:n[i]) {						# loop over individuals within populations
			cpt <- cpt + 1						# increment the line count
			genotypes <- unlist(strsplit(data[cpt],',')) # split a line into the two parts separated by a comma
			genotypes <- genotypes[2]			# the genotypes are in the second part
			genotypes <- unlist(strsplit(genotypes,' ')) # transform the genotypes (encoded as strings) into a list
			genotypes <- genotypes[!(genotypes == '')] # remove any blank element in that list
			genotypes <- as.integer(genotypes)	# transform the genotypes into a vector of integers
			m[cpt,2:(nloci + 1)] <- genotypes	# put the genotypes into the data matrix
		}
	}                                           # in the following, determine the format of the data (number of digits that encode genotypes)
	genotypes <- unlist(strsplit(data[1],',')) 	# split a line into the two parts separated by a comma
	genotypes <- genotypes[2]					# the genotypes are in the second part
	genotypes <- unlist(strsplit(genotypes,' ')) # transform the genotypes (encoded as strings) into a list
	genotypes <- genotypes[!(genotypes == '')]	# remove any blank element in that list
	ploidy <- length(strsplit(genotypes[1],'')[[1]])
	if ((ploidy == 2) | (ploidy == 4)) {digits <- 100} else {
		if ((ploidy == 3) | (ploidy == 6)) {digits <- 1000}
	}
	res <- list(genotypes,format)
	res$genotypes <- m
	res$format <- digits
	return(res)
}

write.detsel.file <- function(data,outfile) {
	pop <- unique(data$genotypes[,1])
	nbrpop <- length(pop)
	nloci <- length(data$genotypes[1,]) - 1
	write(0, file = outfile)
	write(nbrpop, file = outfile,append = TRUE)
	write(nloci, file = outfile,append = TRUE)
	for (i in 1: nloci) {
		genotypes <- data$genotypes[,(i+1)]
		all1 <- floor(genotypes / data$format)
		all2 <- genotypes - (all1 * data$format)
		genes <- c(all1,all2)
		alleles <- unique(genes)
		alleles <- alleles[alleles > 0]
		ntotalleles <- (length(alleles))
		write('', file = outfile,append = TRUE)
		write(ntotalleles, file = outfile,append = TRUE)
		allelescounts <- vector('numeric',ntotalleles)
		for (j in 1:nbrpop) {
			tmp <- genes[data$genotypes[,1] == pop[j]]
			for (k in 1:ntotalleles) {
				allelescounts[k] <- length(tmp[tmp == alleles[k]])
			}
			write(format(allelescounts,width = 5), file = outfile,ncolumns = length(allelescounts),append = TRUE) ###
		}
	}
}

genepop.to.detsel <- function(infile,outfile = 'data.dat') {
	data <- read.genepop.file(infile)
	write.detsel.file(data,outfile)
}

read.data <- function(infile = 'data.dat',dominance = FALSE,maf = 0.99,a = 0.25,b = 0.25) {
	parameterfile <- 'parameters.dat'
	if(!file.exists(infile)) {
		stop(paste('The file ',infile,' does not exist. Check the input file name...',sep = ''))
	}
	write(as.character(infile), file = parameterfile)
	write(as.integer(dominance), file = parameterfile,append = TRUE)
	write(as.double(maf), file = parameterfile,append = TRUE)
	if(dominance) {
		write(as.double(a), file = parameterfile,append = TRUE)
		write(as.double(b), file = parameterfile,append = TRUE)
	}
	if (file.exists('infile.dat')) {
		unlink('infile.dat')
		unlink('plot_*.dat')
		unlink('sample_sizes.dat')
	}
	GetData()
	if(!file.exists('infile.dat')) {
		stop(paste('Problem reading file ',infile,'... Program stopped.',sep = ''))
	}
	data <- get.data.information(infile)
	if (data$biallelic) {
		message(paste('The data file ',infile,' contains ',toString(data$nloci),' biallelic loci, and ',toString(data$npops),' populations',sep = ''))
	} else {
		message(paste('The data file ',infile,' contains ',toString(data$nloci),' loci, with ',min(data$alleles),'-',max(data$alleles),' alleles per locus, and ',toString(data$npops),' populations',sep = ''))
	}
	out <- read.table('infile.dat',skip = 1)
	message('The average values of population-specific measures of differentiation are:')
	message('----------------------------------------------')
	message('Pair\t\tF_1\t\t\tF_2')
	for (i in 1:dim(out)[1]) {
		message(paste(toString(out[i,2]),'-',toString(out[i,3]),'\t\t',format(out[i,4],digits = 3),'\t\t\t',format(out[i,5],digits = 3),sep = ''))
	}
#	message('Pair\t\tF_1\t\t\tF_2\t\t')
#	for (i in 1:dim(out)[1]) {
#		message(paste(toString(out[i,2]),'-',toString(out[i,3]),'\t\t\t',format(out[i,4],digits = 3),'\t\t',format(out[i,5],digits = 3),sep = ''))
#	}
	message('----------------------------------------------')
}

get.data.information <- function(infile) {
	raw <- scan(infile)
	pop_by_rows <- raw[1]
	number_of_alleles <- {}
	index_locus <- 4
	while (index_locus <= length(raw)) {
		tmp <- raw[index_locus]
		number_of_alleles <- append(number_of_alleles,tmp)
		index_locus <- index_locus + (tmp * raw[2]) + 1
	}
	if (length(unique(number_of_alleles)) < 2) {
		if(unique(number_of_alleles) == 2) {
			biallelic = TRUE
		}
	} else {
		biallelic = FALSE
	}
	data <- list(npops = raw[2],nloci = raw[3],alleles = number_of_alleles,biallelic = biallelic)
	return(data)
}

get.simulation.parameters <- function(example) {
	parameterfile <- 'parameters.dat'
	if(!file.exists(parameterfile)) {
		stop(paste('The file ', parameterfile,' does not exist. You must run read.data() first...',sep = ''))
	}
	datafile <- scan(file = parameterfile,what = character(0),skip = 0,n = 1,quiet = TRUE)
	dominance <- scan(file = parameterfile,what = integer(0),skip = 1,n = 1,quiet = TRUE)
	maf <- scan(file = parameterfile,what = double(0),skip = 2,n = 1,quiet = TRUE)
	if(dominance) {
		a <- scan(file = parameterfile,what = double(0),skip = 3,n = 1,quiet = TRUE)
		b <- scan(file = parameterfile,what = double(0),skip = 4,n = 1,quiet = TRUE)
	}
	write(as.character(datafile), file = parameterfile)
	write(as.integer(dominance), file = parameterfile,append = TRUE)
	write(as.double(maf), file = parameterfile,append = TRUE)
	if(dominance) {
		write(as.double(a), file = parameterfile,append = TRUE)
		write(as.double(b), file = parameterfile,append = TRUE)
	}
	if (!example) {
		cat('Total number of simulated points? (default: 500000)')
		n <- scan(file = '',n = 1,what = integer(0),quiet = TRUE)
		if (length(n) == 0) n <- 500000
		write(as.integer(n), file = parameterfile,append = TRUE)
		repeat {
			cat('Average mutation rate?')
			mu <- scan(file = '',n = 1,what = double(0),quiet = TRUE)
			if (length(mu) > 0) break
		}
		write(as.double(mu), file = parameterfile,append = TRUE)
		if(!dominance) {
			repeat {
				cat('Mutation model?')
				cat(' [0 for the infinite allele model; 1 for the stepwise mutation model; any integer k > 1 for a k allele model]')
				mutmod <- scan(file = '',n = 1,what = integer(0),quiet = TRUE)
				if (length(mutmod) > 0) break
			}
			write(as.integer(mutmod), file = parameterfile,append = TRUE)
		}
		cat('Number of sets of parameters? (default: 1)')
		nsets <- scan(file = '',n = 1,what = integer(0),quiet = TRUE)
		if (length(nsets) == 0) nsets <- 1
		parameters <- matrix(nrow = nsets, ncol = 4)
		for (i in 1:nsets) {
			repeat {
				cat(paste('Parameter set ',toString(i),'? (this is a vector with 4 parameters: t,N0,t0,Ne)',sep = ''))
				sets <- scan(file = '',n = 4,what = c(as.matrix(1,4)),quiet = TRUE)
				if (length(sets) == 4) break
			}
			parameters[i,] <- sets
		}
		write(as.matrix(t(parameters)),ncolumns = 4,file = parameterfile,append = TRUE)
	}
	else {
		write(as.integer(10000), file = parameterfile,append = TRUE)
		write(as.double(0.0001), file = parameterfile,append = TRUE)
		write(as.integer(5), file = parameterfile,append = TRUE)
		parameters <- c(100,0,0,20000)
		write(as.matrix(t(parameters)),ncolumns = 4,file = parameterfile,append = TRUE)
	}
}

run.detsel <- function(example = FALSE) {
	get.simulation.parameters(example)
	obs <- read.table('infile.dat',skip = 1)
	all_sample_sizes <- read.table('sample_sizes.dat',skip = 2)
	n <- dim(obs)[1]							# get the number of population pairs
	cpt <- 0									# this is to compute the total number of files that will be created
	for (i in 1:n) {							# loop over population pairs
		pop1 <- obs[i,2]
		pop2 <- obs[i,3]
		if ((obs[i,4] > 0) & (obs[i,5] > 0)) {
			sample_sizes  <- cbind(all_sample_sizes[,(pop1 + 1)],all_sample_sizes[,(pop2 + 1)])
			sample_sizes <- unique(sample_sizes)
			s <- dim(sample_sizes)[1]
			cpt <- cpt + s
		}
	}
	message(paste('The program will now create ',toString(cpt),' simulation files. Please wait, this can take some time...',sep = ''))
	flush.console()
	for (i in 1:n) {							# loop over population pairs
		pop1 <- obs[i,2]
		pop2 <- obs[i,3]
		if ((obs[i,4] > 0) & (obs[i,5] > 0)) {
			sample_sizes  <- cbind(all_sample_sizes[,(pop1 + 1)],all_sample_sizes[,(pop2 + 1)])
			sample_sizes <- unique(sample_sizes)
			s <- dim(sample_sizes)[1]
			for (j in 1:s) {					# loop over sample sizes
				n1 <- sample_sizes[j,1]
				n2 <- sample_sizes[j,2]	
				message(paste('Simulating data in output file: `Pair_',toString(pop1),'_',toString(pop2),'_',toString(n1),'_',toString(n2),'.dat`...',sep = ''))
				flush.console()
				SimulDiv(pop1,pop2,n1,n2)
			}
		}
	}
	message('All the simulations have been completed.');
	target <- read.table('infile.dat',skip = 1)
	realized <- read.table('out.dat',skip = 1)
	message('The difference between observed and simulated values of population-specific measures of differentiation are:')
	message('-------------------------------------------------------------------------')
	message('Pair\t\tF_1 (obs)\tF_1 (sim)\tF_2 (obs)\tF_2 (sim)')
	for (i in 1:n) {
		list <- grep(as.character(target[i,1]),as.character(realized[,1]))
		if (length(list) > 0) {
			message(paste(as.character(target[i,1]),'\t',format(target[i,4],digits = 4),'\t\t',format(mean(realized[list,2]),digits = 4),'\t\t',format(target[i,5],digits = 4),'\t\t',format(mean(realized[list,3]),digits = 4),sep = ''))
		}
	}
	message('-------------------------------------------------------------------------')
}

draw.detsel.graphs <- function(i,j,x.range = c(-1,1),y.range = c(-1,1),n.bins = c(100,100),m = c(2,2),alpha = 0.05,pdf = FALSE,outliers) {
	if (pdf) {
		pdf(file = 'DetSel-outputs.pdf')
	}
	if (missing(i) | missing(j)) {
		infile <- read.table('infile.dat',skip = 1)
		for (n.pairs in 1: dim(infile)[1]) {
			pop.1 <- infile[n.pairs,2]
			pop.2 <- infile[n.pairs,3]
			if ((infile[n.pairs,4] > 0) & (infile[n.pairs,5] > 0)) {
				draw.single.detsel.graph(i = pop.1,j = pop.2,x.range = x.range,y.range = y.range,n.bins = n.bins,m = m,alpha = alpha,pdf = pdf,outliers)
			} else {
			message('multilocus estimates of differentiation in populations ',toString(pop.1),' and ',toString(pop.2),' are negative; cannot draw the graphs...',sep = '')	
		}
		} 
	} else {
		infile <- read.table('infile.dat',skip = 1)
		list <- grep(paste('Pair_',toString(i),'_',toString(j),sep = ''),as.character(infile[,1]))
		if ((infile[list,4] > 0) & (infile[list,5] > 0)) {
			draw.single.detsel.graph(i,j,x.range = x.range,y.range = y.range,n.bins = n.bins,m = m,alpha = alpha,pdf = pdf,outliers)
		} else {
			message('multilocus estimates of differentiation in populations ',toString(i),' and ',toString(j),' are negative; cannot draw the graphs...',sep = '')	
		}
	}
	if (pdf) {
		dev.off()
	}
}

draw.single.detsel.graph <- function(i,j,x.range,y.range,n.bins,m,alpha,pdf,outliers) {
	a <- c(x.range[1],y.range[1])
	b <- c(x.range[2],y.range[2])
	q <- 1 - alpha
	d <- (b - a) / (n.bins - 1)
	lx <- seq(a[1], b[1], by = d[1])
	ly <- seq(a[2], b[2], by = d[2])
	pop1 <- i
	pop2 <- j	
	all_sample_sizes <- read.table('sample_sizes.dat',skip = 2)
	sample_sizes  <- cbind(all_sample_sizes[,(pop1 + 1)],all_sample_sizes[,(pop2 + 1)])
	list_sample_sizes <- unique(sample_sizes)
	s <- dim(list_sample_sizes)[1]
	message(paste('Reading simulation files for populations ',toString(pop1),' and ',toString(pop2),'. Please wait, this can take some time...',sep = ''))
	flush.console()
	data <- {}
	for (j in 1:s) {						# loop over sample sizes
		n1 <- list_sample_sizes[j,1]
		n2 <- list_sample_sizes[j,2]
		tmp <- read.table(paste('Pair_',toString(pop1),'_',toString(pop2),'_',toString(n1),'_',toString(n2),'.dat',sep = ''))
		data <- rbind(data,tmp)
	}
	plotfile <- read.table(paste('plot_',toString(pop1),'_',toString(pop2),'.dat',sep = ''))
	id <- plotfile[,6]
	obs <- cbind(plotfile[,1],plotfile[,2])
	pv <- read.table(paste('P-values_',toString(pop1),'_',toString(pop2),'.dat',sep = ''),skip = 1)
	all <- plotfile[,5]
	nall <- sort(unique(all))
	missing.pvalue <- pv[is.na(pv[,2]),1]
	missing.allele <- plotfile[match(missing.pvalue,plotfile[,6]),5]
	if (length(missing.allele) > 0) {
		for (i in 1:length(missing.allele)) {
			nall <- nall[!nall == missing.allele[i]]
		}
	}
	nbr_pages <- (length(nall) %/% 4) + (length(nall) %% 4)
	message('Plotting graphs...')
	flush.console()
	if (!pdf) {
		dev.new()
	}
	if (length(nall) > 1) {
		par(mfrow = c(2,2))
	}
	for (l in 1:length(nall)) {
		x <- cbind(data[,1][data[,5] == nall[l]],data[,2][data[,5] == nall[l]])
		h <- make.2D.histogram(x,a,b,n.bins)
		f <- ash2(h,m)
		hist <- f$z / sum(f$z)
		freq <- cumulative.distribution.of.probabilities(hist)
		prob <- cbind(freq[,1],(freq[,2] <= q) * q)
		filled.contour3(lx,ly,as.matrix(((hist >=  min(prob[prob[,2] == q,1])))), nlevels = 2,xlab = expression(italic(F)[1]),ylab = expression(italic(F)[2]), main = paste('Marker loci with ',toString(nall[l]),' alleles',sep = ''),col = grey(c(1.0,0.7)))
		nobs <- plotfile[plotfile[,5] == nall[l],1:2]
		nid <- plotfile[plotfile[,5] == nall[l],6]
		pos <- match(nid,pv[,1])
		p <- cbind(nid,pv[pos,2]) #### BUG corrected 14-08-2011
		if (missing(outliers)) {
			selected <- cbind(nobs[p[,2] <= alpha,1],nobs[p[,2] <= alpha,2])
			locus.name <- p[p[,2] <= alpha,1]
			neutral <- cbind(nobs[p[,2] > alpha,1],nobs[p[,2] > alpha,2])
		} else {			
			selected <- cbind(nobs[match(outliers,p[,1]),1],nobs[match(outliers,p[,1]),2])
			locus.name <- p[match(outliers,p[,1]),1]
			neutral <- cbind(nobs[match(setdiff(p[,1],outliers),p[,1]),1],nobs[match(setdiff(p[,1],outliers),p[,1]),2])
		}
		if (length(neutral) > 0) {
			points(neutral,col='black',pch = 16,cex = 0.5)
		}
		if (length(selected) > 0) {
			points(selected,col='black',pch = 8,cex = 0.75)
			text(selected,as.character(locus.name),pos=4,cex=0.6) 
		}
		if (l %/% 4 == l / 4 | length(nall) == 1 | l == length(nall)) {
			title(main = paste('Populations ',toString(pop1),' and ',toString(pop2),' (page ',toString((l - 1) %/% 4 + 1),'/',toString(nbr_pages),')',sep = ''),outer = TRUE,line = -1)
			if (l < length(nall)) {
				if (!pdf) {
					dev.new()
				}
				if (length(nall) > 1) {
					par(mfrow = c(2,2))
				}
			}
		}
	}
	message('Done.')
}
	
compute.p.values <- function(x.range = c(-1,1),y.range = c(-1,1),n.bins = c(100,100),m = c(2,2)) {
	a <- c(x.range[1],y.range[1])
	b <- c(x.range[2],y.range[2])
	message('Computing p-values. Please wait, this can take some time...')
	flush.console()
	pops <- read.table('infile.dat',skip = 1)
	all_sample_sizes <- read.table('sample_sizes.dat',skip = 2)
	n <- dim(pops)[1]
	for (i in 1:n) {						# loop over population pairs
		pop1 <- pops[i,2]
		pop2 <- pops[i,3]
		if ((pops[i,4] > 0) & (pops[i,5] > 0)) {
			pv <- {}
			sample_sizes  <- cbind(all_sample_sizes[,(pop1 + 1)],all_sample_sizes[,(pop2 + 1)])
			list_sample_sizes <- unique(sample_sizes)
			s <- dim(list_sample_sizes)[1]
			plotfile <- read.table(paste('plot_',toString(pop1),'_',toString(pop2),'.dat',sep = ''))
			id <- plotfile[,6]
			obs <- cbind(plotfile[,1],plotfile[,2])
			pv <- matrix(NA,dim(obs)[1],2)
			cpt <- 1
			for (j in 1:s) {					# loop over sample sizes
				n1 <- list_sample_sizes[j,1]
				n2 <- list_sample_sizes[j,2]
				data <- read.table(paste('Pair_',toString(pop1),'_',toString(pop2),'_',toString(n1),'_',toString(n2),'.dat',sep = ''))
				list_loci  <- id[sample_sizes[id,1] == n1 & sample_sizes[id,2] == n2]
				if (length(list_loci) > 0) {
					pos <- match(list_loci,id)
					all <- plotfile[pos,5]
					nall <- unique(all)
					for (l in 1:length(nall)) {
						raw <- cbind(data[,1][data[,5] == nall[l]],data[,2][data[,5] == nall[l]])
						if (nall[l] == 2) {
							hist <- unique(raw)
							list <- array(0,length(hist[,1]))
							for (k in 1:length(hist[,1])) {
								list[k] <- length(raw[,1][raw[,1] == hist[k,1] & raw[,2] == hist[k,2]])
							}
							sub <- list / sum(list)
							list <- sort(unique(sub),decreasing = TRUE)
							pr <- array(0,length(list))
							for (k in 1:length(list)) {
								pr[k] <- list[k] * length(sub[sub == list[k]])
							}
							freq <- cbind(list,pr)
							dist <- cbind(freq[,1],cumsum(freq[,2]))
							nobs <- cbind(plotfile[pos,1][plotfile[pos,5] == nall[l]],plotfile[pos,2][plotfile[pos,5] == nall[l]])
							nid <- plotfile[pos,6][plotfile[pos,5] == nall[l]]
							p <- array(0,length(nobs[,1]))
							for (k in 1:length(nobs[,1])) {
								r1 <- match(nobs[k,1],hist[,1])
								r2 <- match(nobs[k,2],hist[,2])
								if (!(NA %in% r1) & !(NA %in% r2)) {
									if (length(intersect(r1,r2)) > 0) {
										if (nobs[k,1] == hist[r1,1] & nobs[k,2] == hist[r2,2]) {
											lev <- sub[r1]
											x <- seq(1,length(dist[,1]))[dist[,1] == lev]
											p[k] <- dist[x,2]
										} else {
											p[k] <- 1.0
		  								}		
									} else {
										p[k] <- 1.0
									}	
								} else {
									p[k] <- 1.0
								}
 							}
						} else {
							nobs <- plotfile[pos[plotfile[pos,5] == nall[l]],1:2]
							nid <- plotfile[pos[plotfile[pos,5] == nall[l]],6]
							if (dim(raw)[1] > 0) {
								h <- make.2D.histogram(raw,a,b,n.bins)
								f <- ash2(h,m)
								hist <- f$z / sum(f$z)
								freq <- cumulative.distribution.of.probabilities(hist)
								d <- (b - a) / (n.bins - 1)
								v <- trunc((nobs - a) / d) + 1
								lev <- hist[as.matrix(v)]
								p <- (freq[match(as.factor(lev),as.factor(freq[,1])),2]) # need to coerce with 'as.factor'
							} else {
								message(paste('Could not compute the p-value for locus ',toString(nid),' in population pair ',toString(pop1),'-',toString(pop2),sep = ''))
								flush.console()
								p <- NA
							}
						}
	 					q <- cbind(nid,(1 - p))
 	 					pv[cpt:(cpt + length(p) - 1),] <- q
 						cpt <- cpt + length(p)
					}
				}
			}
  			o <- order(pv[,1])
  			loc <- cbind.data.frame('Locus' = pv[o,1])
  			pvalue <- cbind.data.frame('P-value' = as.double(format(pv[o,2],digits = 6)))
  			out <- cbind.data.frame(loc,pvalue)
  			message(paste('The p-values for each locus in population pair ',toString(pop1),'-',toString(pop2),' are:',sep = ''))
  			message('-------------------')
  			print(out)
  			message('-------------------')
			write.table(out,file = paste('P-values_',toString(pop1),'_',toString(pop2),'.dat',sep = ''),row.names = FALSE,sep = '\t\t')
			message(paste('The above results are saved in file: P-values_',toString(pop1),'_',toString(pop2),'.dat',sep = ''))
		}  else {
			message('multilocus estimates of differentiation in populations ',toString(pop1),' and ',toString(pop2),' are negative; cannot compute the p-values...',sep = '')	
		}
	}
}

make.2D.histogram <- function(x,a,b,n.bins) { # This function returns a 2D (relative) histogram of the data, 
	d <- (b - a) / (n.bins - 1)
	h <- array(0,n.bins)
	k <- trunc((x - a) / d) + 1  
	c <- table(apply(k,1,paste,collapse = ','))
	n <- names(c)
	coordX <- strsplit(n,',')
	coordN <- lapply(coordX,as.numeric)
	pos <- t(as.data.frame(coordN,optional = TRUE))
	r <- replace(h,pos,c)
	n <- dim(x)[1]
	ab <- matrix(c(a,b),2,2)
	list(nc = r,ab = ab)
}

cumulative.distribution.of.probabilities <- function(x) { # This function returns a list with (1) [...], and (2) the [...] (i.e., the frequency times the count)
	list <- as.data.frame(table(x))
	c1 <- as.numeric(as.vector(list[,1]))
	c2 <- as.numeric(as.vector(list[,2]))
	o <- order(sort(c1,decreasing = TRUE))
	cbind(c1[o],cumsum(c1[o] * c2[o]))
}

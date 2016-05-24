PS.Main <- function(dat, para=list())
{
	### check arguments
	cat("Checking parameters of dat...", fill=T)
	dat <- Check.Argu.Dat(dat)
	cat("Checking parameters of para...", fill=T)
	para <- Check.Argu.Para(para, dat)
	
	### filter out genes that has too small counts
	dat <- PS.Filter(dat, ct.sum=para$ct.sum, ct.mean=para$ct.mean)
	
	### begin calculation
	time_st <- Sys.time()
	
	### estimate sequencing depth
	seq.depth <- Est.Depth(dat$n)$cmeans	# This seq depth is for the original data
	#Display.Warning(dat, seq.depth)

	# doing order transformation
	if (para$trans)
	{
		dat$n.ori <- dat$n
		cat("\nDoing order transformation...\n")
		dat$n <- Trans.Data(n=dat$n, div=para$div, 
			pow.file=para$pow.file)
	} else
	{
		cat("\nWithout doing power transformation...\n")
	}
	
	# calculate the score statistics
	cat("\nCalculating score statistics and permutated score statistics...\n")
	set.seed(para$seed)
	ps.obj <- Score.Stat(dat=dat, para=para)
	
	# estimate FDR
	cat("\nCalculating FDR...\n")
	fdr.res <- PS.FDR(ps.obj=ps.obj)
	
	res <- PS.Sum(dat=dat, fdr.res=fdr.res, seq.depth=seq.depth, trans=para$trans)
	
	time_ed <- Sys.time()
	print(time_ed - time_st)

	return(res)
}

############################################################
# check arguments of dat
############################################################
Check.Argu.Dat <- function(dat)
{
	if (is.null(dat$n))
	{
		stop("the data matrix dat$n must be specified!")
	}
	
	if (is.null(dat$y))
	{
		stop("the outcome vector dat$y must be specified!")
	}
	
	if (is.null(dat$type))
	{
		stop("the outcome type dat$type must be specified!")
	}
	
	dat$n <- as.matrix(dat$n)
	if (!is.matrix(dat$n) | !is.numeric(dat$n))
	{
		stop('n must be a numeric matrix!')
	}
	
	if (sum(dat$n >= 0) != length(dat$n))
	{
		stop('All elements of n must be non-negative!')
	}
	
	if (!is.numeric(dat$y) | length(dat$y) <= 1)
	{
		stop('y must be a numeric vector with length no less than 2!')
	}
	
	if (ncol(dat$n) != length(dat$y))
	{
		stop('Number of column of n should be the same as the length of y!')
	}
	
	dat$type <- match.arg(dat$type, c('twoclass', 'multiclass', 'quant'), several.ok = FALSE)
	
	if (is.null(dat$pair))
	{
		dat$pair <- F
	}
	
	if (!is.logical(dat$pair))
	{
		stop('pair must be logical (either TRUE of FALSE)!')
	}
	
	if (is.null(dat$gname))
	{
		dat$gname <- 1 : nrow(dat$n)
	}
	
	if (length(dat$gname) != nrow(dat$n))
	{
		stop("Length of dat$gname must equal the number of genes nrow(dat$n)!")
	}
	
	return(dat)
}

############################################################
# check arguments of para
############################################################
Check.Argu.Para <- function(para, dat)
{
	if (is.null(para$trans))
	{
		para$trans <- T
	} else
	{
		if (!is.logical(para$trans))
		{
			stop('trans is not logical!')
		}
	}
	
	if (is.null(para$npermu))
	{
		para$npermu <- 100
	} else
	{
		if (!is.numeric(para$npermu) | para$npermu <= 0)
		{
			stop('npermu must be a positive integer!')
		}
	}
	
	if (is.null(para$div))
	{
		para$div <- 10
	} else
	{
		if (!is.numeric(para$div) | para$div < 5)
		{
			stop('div must be a integer no smaller than 5!')
		}
	}
	
	if (is.null(para$pow.file))
	{
		para$pow.file <- "pow.txt"
	}
	
	if (is.null(para$ct.sum))
	{
		para$ct.sum <- 5
	}
	
	if (!is.numeric(para$ct.sum))
	{
		stop("para$ct.sum must be numeric!")
	}
	
	if (is.null(para$ct.mean))
	{
		para$ct.mean <- 0.5
	}
	
	if (!is.numeric(para$ct.mean))
	{
		stop("para$ct.mean must be numeric!")
	}
	
	if (is.null(para$seed))
	{
		para$seed <- 10
	}
	
	return(para)
}

############################################################
# display the warning message
############################################################
Display.Warning <- function(dat, seq.depth)
{
	if (dat$type == 'twoclass')
	{
		if (abs(log(mean(seq.depth[dat$y == 1]) / mean(seq.depth[dat$y == 2]))) > log(2))
		{
			cat("WARNING: \nThe mean sequencing depths of the two classes are quite different, so our program may underestimate the false discovery rates!\n")
		}
	}
}

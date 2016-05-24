##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: http://www.stats.ox.ac.uk/~snijders/siena
## *
## * File: sienaprint.r
## *
## * Description: This file contains the print and summary modules for the
## * classes siena, sienaFit, sienaDependent, sienaAlgorithm, and sienaBayesFit
## *
## ****************************************************************************/
##@print.siena Methods
print.siena <- function(x, ...)
{
    ##@onlys internal print.siena; prints matrix attributes of x$depvar
	onlys <- function(x, attrib){
	# function to print the attribute attrib (character string)
	# of object x$depvars
	# This attribute must be a logical matrix
	# with named columns, and rows indicating periods.
		textattr <- deparse(substitute(attrib),  backtick=FALSE)
		uponlys <- t(as.matrix(sapply(x$depvars, function(y){attr(y,attrib)})))
		for (j in 1:dim(uponlys)[2])
		{
			if (any(uponlys[,j]))
			{
				cat(textattr)
				cat(":  ", colnames(uponlys)[j], ":    ", sep="")
				if (all(uponlys[,j]))
				{
					cat("all periods\n")
				}
				else
				{
					word <-
						ifelse((sum(uponlys[,j]) <= 1), "period ", "periods ")
					cat(word)
					cat(which(uponlys[,j]), sep=", ")
					cat("\n")
				}
			}
		}
	}
# begin main method siena.print proper
	if (!inherits(x, "siena"))
	{
        stop("not a legitimate Siena data object")
	}
	cat('Dependent variables: ', paste(names(x$depvars), collapse=", "), "\n")
	cat('Number of observations:', x$observations, "\n\n")
	if (!is.null(x$nodeSet))
	{
		tmp <- sapply(x$nodeSet, length)
		if (length(x$nodeSet) <= 1)
		{
			tmp <- c('Nodeset        ' = 'Number of nodes', tmp)
		}
		else
		{
			tmp <- c('Nodesets       ' = 'Number of nodes', tmp)
		}
		print(tmp, quote=FALSE)
		cat("\n")
	}

	for (j in (1:length(x$depvars)))
	{
		xj <- x$depvars[[j]]
		mymat <- matrix("", 7, 2)
		mymat[1,] <- c("Dependent variable", attr(xj, "name"))
		mymat[2,] <- c("Type",               attr(xj, "type"))
		mymat[3,] <- c("Observations",       attr(xj,  "netdims")[3])
		if (attr(xj, "type") == "bipartite")
		{
			mymat[4,] <- c("First nodeset ", attr(xj, "nodeSet")[1])
			mymat[5,] <- c("Second nodeset ",attr(xj, "nodeSet")[2])
			nrows <- 5
		}
		else
		{
			mymat[4,] <- c("Nodeset ", attr(xj, "nodeSet"))
			nrows <- 4
		}
		if (attr(xj, "type") == "behavior")
		{
			mymat[nrows+1,] <- c("Range",
				paste(attr(xj, "range2"),  collapse=" - "))
		}
		else
		{
			mymat[nrows+1,] <- c("Densities",
				paste(signif(attr(xj, "density"), 2),  collapse=" "))
		}
		mymat2 <- apply(mymat[0:nrows+1,], 2, format)
		write.table(mymat2, row.names=FALSE, col.names=FALSE, quote=FALSE)
		cat("\n")
	}

	if (length(x$cCovars) > 0)
	{
		cat('Constant covariates: ', paste(names(x$cCovars), collapse=", "), "\n")
	}
	if (length(x$vCovars) > 0)
	{
		cat('Changing covariates: ',
			paste(names(x$vCovars), collapse = ", "), "\n")
	}
	if (length(x$dycCovars) > 0)
	{
		cat('Constant dyadic covariates: ',
			paste(names(x$dycCovars), collapse=", "), "\n")
	}
	if (length(x$dyvCovars) > 0)
	{
		cat('Changing dyadic covariates: ',
			paste(names(x$dyvCovars), collapse=", "), "\n")
	}

	onlys(x, 'uponly')
	onlys(x, 'downonly')

	attrs <- attributes(x)
	highers <- attrs[["higher"]]
	disjoints <- attrs[["disjoint"]]
	atleastones <- attrs[["atLeastOne"]]
	if (any(highers))
	{
		cat("Higher: ", names(highers)[highers], sep="", "\n")
	}
	if (any(disjoints))
	{
		cat("Disjoint: ", names(disjoints)[disjoints], sep="",  "\n")
	}
	if (any(atleastones))
	{
		cat("atLeastOne: ", names(atleastones)[atleastones], sep="",  "\n")
	}
	invisible(x)
}
##@print.sienaGroup Methods
print.sienaGroup <- function(x, ...)
{
	if (!inherits(x, "sienaGroup"))
	{
        stop("not a legitimate Siena group data object")
	}
	att <- attributes(x)
	cat('Dependent variables: \n')
	cat(paste(att$netnames, ":", att$types),'\n')
	cat('Total number of periods:', att$observations)
	cat("\nmore to be added!\n")
	invisible(x)
}

print.sienaDependent <- function(x, ...)
{
	if (!inherits(x, "sienaDependent"))
	{
        stop("not a legitimate Siena dependent variable object")
	}
	mymat <- matrix("", 4, 2)
	mymat[1,] <- c("Type",               attr(x, "type"))
	mymat[2,] <- c("Observations",       attr(x,  "netdims")[3])
	if (attr(x, "type") == "bipartite")
	{
		mymat[3,] <- c("First nodeset ",
						paste(attr(x, "nodeSet")[1], " (", attr(x, "netdims")[1],
							" elements)", sep=""))
		mymat[4,] <- c("Second nodeset ",
						paste(attr(x, "nodeSet")[2], " (", attr(x, "netdims")[2],
							" elements)", sep=""))
		nrows <- 4
	}
	else
	{
		mymat[3,] <- c("Nodeset ",
						paste(attr(x, "nodeSet"), " (", attr(x, "netdims")[1],
							" elements)", sep=""))
		nrows <- 3
	}

	mymat <- apply(mymat, 2, format)
	write.table(mymat[1:nrows,], row.names=FALSE, col.names=FALSE, quote=FALSE)
	cat("\n")
	invisible(x)
}

##@print.sienafit Methods
print.sienaFit <- function(x, tstat=TRUE, ...)
{
	if (!inherits(x, "sienaFit"))
	{
        stop("not a legitimate Siena model fit")
	}
	if (!x$OK)
	{
		cat("Error end of estimation algorithm\n")
	}
	else if (x$termination == "UserInterrupt")
	{
		cat("User interrupted run, object possibly incomplete\n")
	}
	else
	{
		cat("Estimates, standard errors and convergence t-ratios\n\n")
		tmp <- sienaFitThetaTable(x, tstat=tstat)
		mydf <- tmp$mydf
		mymat <- as.matrix(mydf)
		mymat[, 'value'] <- format(round(mydf$value, digits=4))
		mymat[, 'se'] <- format(round(mydf$se, digits=4))
		mymat[, 'tstat'] <- paste(" ", format(round(mydf$tstat, digits=4)))
		mymat[is.na(mydf$tstat), 'tstat'] <- ' '
		mymat[, 'type'] <- format(mymat[, 'type'])
		mymat[, 'text'] <- format(mymat[, 'text'])
		mymat[mydf$row < 1, 'row'] <-
			format(mydf[mydf$row < 1, 'row'])
		mymat[mydf[,'row'] >= 1, 'row'] <-
			paste(format(mydf[mydf$row >= 1, 'row']), '.', sep='')
		mymat <- rbind(c(rep("", 4), "Estimate", "", "Standard", "",
						 "Convergence"),
					   c(rep("", 6),  "  Error", "", "  t-ratio"), mymat)
		mymat <- apply(mymat, 2, format)
		tmp1 <- apply(mymat, 1, function(x) paste(x, collapse=" "))
		addtorow <- tmp$addtorow
		for (i in 1:length(tmp1))
		{
			if (length(addtorow$command) > 0)
			{
				for (j in 1:length(addtorow$command))
				{
					ii <- match(i-1, addtorow$pos[[j]])
					if (!is.na(ii))
						if (i == 2 | addtorow$command[j] == 'Network Dynamics')
							cat( addtorow$command[j], '\n')
						else
							cat('\n', addtorow$command[j], '\n', sep='')
				}
			}
			cat(tmp1[i], '\n')
		}

		cat("\nTotal of", x$n, "iteration steps.\n\n")
		if (x$termination == "UserInterrupt")
			cat(" \n*** Warning ***",
				"Estimation terminated early at user request.\n")
	}
	invisible(x)
}

##@summary.sienaFit Methods
summary.sienaFit <- function(object, ...)
{
    if (!inherits(object, "sienaFit"))
	{
        stop("not a legitimate Siena model fit")
	}
    class(object) <- c("summary.sienaFit", class(object))
    object
}
##@print.summary.sienaFit Methods
print.summary.sienaFit <- function(x, ...)
{
	if (!inherits(x, "summary.sienaFit"))
	{
        stop("not a legitimate summary of a Siena model fit")
	}
	print.sienaFit(x)
	cat(c("Overall maximum convergence ratio: ",
		sprintf("%8.4f", x$tconv.max), "\n\n"))
	if (sum(x$test) > 0) ## we have some score tests
	{
		testn <- sum(x$test)
		if (x$maxlike)
		{
			cat("Score test <c>\n\n")
		}
		else
		{
			cat("Generalised score test <c>\n\n")
		}
		cat("Testing the goodness-of-fit of the model restricted by\n")
		j <- 0
		for (k in 1:x$pp)
		{
			if (x$test[k])
			{
				j <- j + 1
				cat(c(" (", j, ")   ",
					  format(paste(x$requestedEffects$type[k], ":  ",
								   x$requestedEffects$effectName[k],
								   sep=""),
							 width=50), " = ",
					  sprintf("%8.4f", x$theta[k]),"\n"),
					sep = "")
			}
		}
		cat("_________________________________________________\n")
		cat("                ")
		cat("   \n")
		if (testn > 1)
		{
			cat('Joint test:\n-----------\n')
		}
		cat(c('   c = ',sprintf("%8.4f", x$testresOverall),
			  '   d.f. = ',j,'   p-value '), sep='')
		pvalue <- 1 - pchisq(x$testresOverall, j)
        if (pvalue < 0.0001)
		{
            cat('< 0.0001\n')
		}
        else
		{
            cat(c('= ', sprintf("%8.4f\n", pvalue)), sep = '')
		}
        if (testn==1)
		{
            cat(c('\n   one-sided (normal variate): ',
                     sprintf("%8.4f", x$testresulto[1])), sep = '')
		}
        if (testn> 1)
        {
            cat('\n\n')
            for (k in 1:j)
            {
                cat(c('(', k, ') tested separately:\n'), sep='')
                cat('-----------------------\n')
                cat(' - two-sided:\n')
                cat(c('  c = ', sprintf("%8.4f", x$testresult[k]),
                         '   d.f. = 1  p-value '), sep = '')
                pvalue<- 1-pchisq(x$testresult[k],1)
                if (pvalue < 0.0001)
				{
                    cat('< 0.0001\n')
				}
                else
				{
                    cat(c('= ', sprintf("%8.4f", pvalue), '\n'), sep = '')
				}
                cat(c(' - one-sided (normal variate): ',
					  sprintf("%8.4f", x$testresulto[k])), sep = '')
                if (k < j)
				{
                    cat('\n\n')
				}
            }
        }
        cat('    \n_________________________________________________\n\n')
        cat('One-step estimates: \n\n')
        for (i in 1 : x$pp)
        {
            onestepest <- x$oneStep[i] + x$theta[i]
            cat(c(format(paste(x$requestedEffects$type[i], ':  ',
							   x$requestedEffects$effectName[i], sep = ''),
						 width=50),
                     sprintf("%8.4f", onestepest), '\n'), sep = "")
        }
        cat('\n')
   }
   if ((x$OK)&(!is.null(x$covtheta)))
   {
       cat("Covariance matrix of estimates (correlations below diagonal)\n\n")
       covcor <- x$covtheta
       correl <- x$covtheta / sqrt(diag(x$covtheta))[row(x$covtheta)] /
           sqrt(diag(x$covtheta))[col(x$covtheta)]
       covcor[lower.tri(covcor)] <- correl[lower.tri(correl)]
       printMatrix(format(round(t(covcor),digits=3),width=12))
       cat("\nDerivative matrix of expected statistics X by parameters:\n\n")
       printMatrix(format(round(x$dfra,digits=3),width=12))
       cat("\nCovariance matrix of X (correlations below diagonal):\n\n")
       covcor <- x$msf
       correl <- x$msf / sqrt(diag(x$msf))[row(x$msf)] /
		   sqrt(diag(x$msf))[col(x$msf)]
       covcor[lower.tri(covcor)] <- correl[lower.tri(correl)]
       printMatrix(format(round(t(covcor),digits=3),width=12))
   }
   invisible(x)
}

##@printMatrix Miscellaneous
printMatrix <- function(mat)
{
    cat(mat, sep=c(rep.int(' ', ncol(mat) - 1), '\n'))
}
##@print.sienaAlgorithm Methods
print.sienaAlgorithm <- function(x, ...)
{
    cat(' Project name:', x$projname, '\n')
    cat(' Use standard initial values:', x$useStdInits, '\n')
    cat(' Random seed:', x$randomSeed,'\n')
    cat(' Starting value of gain parameter:', x$firstg, '\n')
	cat(' Diagonalization parameter:', x$diagonalize, '\n')
	cat(' Dolby noise reduction:', x$dolby, '\n')
    if (any(x$MaxDegree > 0))
    {
        cat(' Restrictions on degree in simulations:')
        cat(x$MaxDegree,'\n')
    }
    cat(' Method for calculation of Derivatives:',
        c('Scores', 'Finite Differences')[as.numeric(x$FinDiff.method) + 1],
        '\n')
    cat(' Number of subphases in phase 2:', x$nsub, '\n')
    cat(' Number of iterations in phase 3:', x$n3, '\n')
	if (x$maxlike)
	{
		cat(" Fitting by maximum likelihood\n")
	}
	else
	{
		if (is.na(x$cconditional))
		{
			cat(" Unconditional simulation if more than one dependent",
				"variable\n")
		}
		else
		{
			if (x$cconditional)
			{
				cat(" Conditional simulation:")
				if (x$condname != '')
				{
					cat('conditioned on', x$condname, '\n')
				}
				else
				{
					if (x$condvarno > 0)
					{
						cat('conditioned on variable number', x$condvarno, '\n')
					}
				}
			}
			else
			{
				cat(" Unconditional simulation\n")
			}
		}
	}
    cat(" Model Type:", ModelTypeStrings[x$modelType], "\n")
    invisible(x)
}
##@sienaFitThetaTable Miscellaneous
sienaFitThetaTable <- function(x, tstat=FALSE)
{
    effects <- x$requestedEffects
    pp <- x$pp
    if (x$cconditional)
    {
        nrates <- length(x$rate)
    }
    else
    {
        nrates <- 0
    }
    pp <- pp + nrates
    ## mydf stores the data before formatting
    mydf <- data.frame(dummy=rep(" ", pp),
                       row=rep(0, pp),
                       type=rep("", pp),
                       text=rep(" ", pp),
                       value = rep(0, pp),
                       lParenthesis = rep("(", pp),
                       se = rep(0, pp),
                       rParenthesis = rep(")", pp),
                       tstat=rep(NA, pp),
                       stringsAsFactors =FALSE)
    ## add to row is the extra lines to put in to table if you wish
    addtorow <- list()
    addtorow$pos <- list()
    addsub <- 1
    if (nrates > 0) ## conditional
    {
        addtorow$command[addsub] <-
            'Rate parameters: '
        addtorow$pos[[addsub]] <- 2
        addsub <- addsub + 1
        if (length(x$rate) == 1)
        {
            mydf[1, 'row'] <- 0
            if (length(attr(x$f,'netnames')) == 1)
            {
                mydf[1, 'text'] <- 'Rate parameter'
            }
            else
            {
                mydf[1, 'text'] <- 'Rate parameter of conditioning variable'
            }
            mydf[1, 'value'] <- x$rate[1]
            mydf[1, 'se'] <- x$vrate[1]
        }
        else ## observations > 2
        {
            nn <- length(x$rate)
            nnstr <- as.numeric(paste('0.', as.character(1:nn), sep=""))
            mydf[1:nn, 'row'] <- nnstr
            mydf[1:nn, 'value'] <- x$rate
            mydf[1:nn, 'se'] <- x$vrate
            if (length(x$f$types) == 1)
            {
                mydf[1:nn, 'text'] <- paste('Rate parameter period', 1:nn)
            }
            else
            {
                mydf[1:nn, 'text'] <-
                    paste('Rate parameter cond. variable period', 1:nn)
            }
            addtorow$command[addsub] <-
                'Other parameters: '
            addtorow$pos[[addsub]] <- nn + 2
            addsub <- addsub + 1
        }
    }
    nBehavs <- sum(x$f$types == "behavior")
    nNetworks <- length(x$f$types) - nBehavs
    if (nBehavs > 0 && nNetworks > 0)
    {
        addtorow$command[addsub] <-
            "Network Dynamics"
        addtorow$pos[[addsub]] <- nrates + 2
        addsub <- addsub + 1
    }
	if (!is.null(x$covtheta))
	{
		ses <- sqrt(diag(x$covtheta))
		ses[x$fixed] <- NA
	}
    theta <- x$theta
	if (!is.null(x$covtheta))
	{
		theta[diag(x$covtheta) < 0.0] <- NA
	}

    if (nBehavs > 0)
    {
        behEffects <- effects[effects$netType == 'behavior',]
        behNames <- unique(behEffects$name)
    }
    if (nBehavs > 1)
    {
        behEffects$effectName <- paste('<',
                                       (1:nBehavs)[match(behEffects$name,
                                                         behNames)],
                                       '> ', behEffects$effectName,
                                       sep='')
        effects$effectName[effects$netType=='behavior'] <-
            behEffects$effectName
    }
    mydf[nrates + (1:x$pp), 'row'] <-  1:x$pp
    mydf[nrates + (1:x$pp), 'type' ] <- ifelse(effects$type == "creation",
                                               "creat", effects$type)
    mydf[nrates + (1:x$pp), 'text' ] <- effects$effectName
    mydf[nrates + (1:x$pp), 'value' ] <- theta
	if (exists("ses"))
	{
		mydf[nrates + (1:x$pp), 'se' ] <- ses
	}
    if (!is.null(x$tstat))
    {
        mydf[1:nrates, "tstat"] <- NA
        mydf[nrates + (1:x$pp), 'tstat' ] <- x$tstat
    }

    if (nBehavs > 0 && nNetworks > 0)
    {
        nNetworkEff <- nrow(effects) - nrow(behEffects)
        addtorow$command[addsub] <-
            'Behavior Dynamics'
        addtorow$pos[[addsub]] <- nrates + 2 + nNetworkEff
        addsub <- addsub + 1
    }
    return(list(mydf=mydf, addtorow=addtorow))
}

##@sienaFitCovarianceCorrelation Miscellaneous
sienaFitCovarianceCorrelation <- function(x)
{
    covcor <- x
    correl <- x/sqrt(diag(x))[row(x)]/ sqrt(diag(x))[col(x)]
    covcor[lower.tri(covcor)] <- correl[lower.tri(correl)]
    return(covcor)

}

##@xtable fake in case package not loaded
xtable <- function(x, ...)
{
	xtable::xtable(x, ...)
}

##@xtable.sienaFit Methods
xtable.sienaFit <- function(x, caption = NULL, label = NULL, align = NULL,
                            digits = NULL, display = NULL, ...)
{
    tmp <- sienaFitThetaTable(x)
    mydf <- tmp$mydf
    addtorow <- tmp$addtorow
    ## find out whether the type is html or latex
    dots <- substitute(list(...))[-1] ##first entry is the word 'list'
    if (!is.null(dots[["type"]]))
    {
        type <- dots[["type"]]
    }
    else
    {
        type <- "latex"
    }
    if (!is.null(addtorow$command))
    {
        if (type =="latex")
        {
            use <- addtorow$command != 'Network Dynamics'
            addtorow$command <- paste('\\multicolumn{4}{l}{', addtorow$command,
                                      '} \\\\ \n')
            use[1] <- FALSE
            addtorow$command[use] <- paste('\\\\ ', addtorow$command[use])
        }
        else ##html
        {
			## use <- addtorow$command != 'Network Dynamics'
            addtorow$command <- paste("<TR> <TD colspan=9 align=left>",
                                      addtorow$command,
                                      "</TD> </TR> <TR> </TR> \n")
			##  use[1] <- FALSE
			##  addtorow$command[use] <- paste('\\\\ ', addtorow$command[use])
        }
    }
    else
    {
        addtorow <- NULL
    }
    mydf[mydf$row < 1, 'row'] <-
        format(mydf[mydf$row < 1, 'row'])
    mydf[mydf[,'row'] >= 1, 'row'] <- paste(format(mydf[mydf$row >= 1,
             'row']), '.', sep='')
    tmp <- list(xtable::xtable(mydf, caption=caption, label=label, align=align,
							   digits=digits, display=display),
				add.to.row=addtorow,
                include.colnames=FALSE, include.rownames=FALSE, ...)
    class(tmp) <- c("xtable.sienaFit", "xtable")
    tmp
}
##@print.xtable.sienaFit Methods
print.xtable.sienaFit <- function(x, ...)
{
    addtorow <- x[["add.to.row"]]
    if (!is.null(addtorow))
    {
        x$add.to.row$pos <- lapply(x$add.to.row$pos, function(x)x-2)
        do.call("print", x)
    }
    else
    {
        do.call("print", x[-2])
    }
    invisible(x)
}

##@print.chains.data.frame Methods
print.chains.data.frame <- function(x, ...)
{
    NextMethod(x)
	initState <- attr(x, "initialStateDifferences")
	endState <- attr(x, "endStateDifferences")
	if (nrow(initState) > 0)
	{
		cat("\n Initial State Differences:\n")
		print(initState)
	}
	if (nrow(endState) > 0)
	{
		cat("\n End State Differences:\n")
		print(endState)
	}
}

##@maketemp Methods
makeTemp <- function(x, ...)
{
# This reproduces a part of sienaFit but with Bayesian headings;
# for use in print.sienaBayesFit and summary.sienaBayesFit
		name1 <- x$requestedEffects$groupName[1]
		x$requestedEffects <-
			x$requestedEffects[x$requestedEffects$groupName==name1,]
		x$pp <- length(x$theta)
		x$fixed <- x$fixed[x$effects$groupName==name1]
		tmp <- sienaFitThetaTable(x)
		mydf <- tmp$mydf
		mymat <- as.matrix(mydf[,names(mydf)!="tstat"])
		mymat[, 'value'] <- format(round(mydf$value, digits=4))
		mymat[, 'se'] <- format(round(mydf$se, digits=4))
		mymat[, 'type'] <- format(mymat[, 'type'])
		mymat[, 'text'] <- format(mymat[, 'text'])
		mymat[mydf$row < 1, 'row'] <-
			format(mydf[mydf$row < 1, 'row'])
		mymat[mydf[,'row'] >= 1, 'row'] <-
			paste(format(mydf[mydf$row >= 1, 'row']), '.', sep='')
		mymat <- rbind(c(rep("", 4), "Post.   ", "", "Post.   ", ""),
					   c(rep("", 4),  "mean    ", "", "s.d.    ", ""),
						mymat)
		mymat <- apply(mymat, 2, format)
		tmp1 <- apply(mymat, 1, function(x) paste(x, collapse=" "))
		list(tmp, tmp1)
}

##@print.sienaBayesFit Methods
print.sienaBayesFit <- function(x, ...)
{
	printmat <- function(mat)
	{
		cat(sprintf("%8.4f",mat), sep=c(rep.int(' ', ncol(mat) - 1), '\n'))
	}
	if (!inherits(x, "sienaBayesFit"))
	{
        stop("not a legitimate Siena Bayes model fit")
	}
	else
	{
		cat("Note: this summary does not contain a convergence check.\n\n")
		if (length(x$f$groupNames) > 1)
		{
			cat("Groups:\n")
			s <- ""
			ml <- max(nchar(x$f$groupNames))
			fmt <- paste("%-", ml+2, "s", sep="")
			for (i in 1:dim(x$ThinParameters)[2])
			{
				s <- paste(s, gettextf(fmt,x$f$groupNames[i]))
				if ((nchar(s) + ml > 80) | (i == dim(x$ThinParameters)[2]))
				{
					cat(paste(s,"\n"))
					s <- ""
				}
			}
			cat("\n")
		}

		cat("Posterior means and standard deviations ")
		# Make temporary changes to make x look like a sienaFit object
		# so that sienaFitThetaTable can be applied.
		# This is done in function makeTemp.
		if (length(x$f$groupNames) <= 1)
		{
			cat("\n\n")
			x$theta <- colMeans(x$ThinParameters[,1,])
			x$covtheta <- cov(x$ThinParameters[,1,])
		}
		else
		{
			cat("for global mean parameters\n\n")
			nmain <- dim(x$ThinPosteriorMu)[1]
			if (x$frequentist)
			{
			first <- nmain - x$lengthPhase3 + 1
			}
			else
			{
			first <- 1
			}
			x$theta <- colMeans(x$ThinPosteriorMu[first:nmain,])
			x$covtheta <- cov(x$ThinPosteriorMu)
		}
		tmps <- makeTemp(x)
		tmp <- tmps[[1]]
		tmp1 <- tmps[[2]]
		addtorow <- tmp$addtorow
		for (i in 1:length(tmp1))
		{
			if (length(addtorow$command) > 0)
			{
				for (j in 1:length(addtorow$command))
				{
					ii <- match(i-1, addtorow$pos[[j]])
					if (!is.na(ii))
						if (i == 2 | addtorow$command[j] == 'Network Dynamics')
							cat( addtorow$command[j], '\n')
						else
							cat('\n', addtorow$command[j], '\n', sep='')
				}
			}
			cat(tmp1[i], '\n')
		}
		if (length(x$f$groupNames) > 1)
		{
			cat("\n")
			mean.Sigma <-
				apply(x$ThinPosteriorSigma[first:nmain,,], c(2,3), mean)
			sd.Sigma <- apply(x$ThinPosteriorSigma[first:nmain,,], c(2,3), sd)
			cat("Posterior mean of global covariance matrix\n")
			printmat(mean.Sigma)
			cat("\nPosterior standard deviations of ")
			cat("elements of global covariance matrix\n")
			printmat(sd.Sigma)
		}
		cat("\nTotal of", dim(x$ThinPosteriorMu)[1], "samples.\n\n")
	}
	invisible(x)
}

##@summary.sienaBayesFit Methods
summary.sienaBayesFit <- function(object, ...)
{
    if (!inherits(object, "sienaBayesFit"))
	{
        stop("not a legitimate Siena Bayes model fit")
	}
    class(object) <- c("summary.sienaBayesFit", class(object))
    object
}
##@print.summary.sienaBayesFit Methods
print.summary.sienaBayesFit <- function(x, ...)
{
	if (!inherits(x, "summary.sienaBayesFit"))
	{
        stop("not a legitimate summary of a Siena Bayes model fit")
	}
	nmain <- dim(x$ThinPosteriorMu)[1]
	if (x$frequentist)
	{
		cat("Frequentist estimation.\n")
		first <- nmain - x$lengthPhase3 + 1
	}
	else
	{
		first <- 1
		cat("Prior distribution:\n")
		cat("\nMu      ")
		for (i in seq(along=x$priorMu))
		{
			cat(sprintf("%8.4f", x$priorMu[i,1]),"\n        ")
		}
		cat("\nSigma   ")
		for (i in seq(along=x$priorMu))
		{
			cat(sprintf("%8.4f", x$priorSigma[i,]),"\n        ")
		}
		cat("\nDf       ",sprintf("%1d", x$priorDf),"\n")
		if (length(x$f$groupNames) >= 2)
		{
			cat("\nKappa  ",sprintf("%8.4f", x$priorKappa),"\n")
		}
		cat("\nFor the basic rate parameters, ")
		cat("the prior is on the square root scale.\n\n")
	}
	print.sienaBayesFit(x)

	cat("Posterior means and standard deviations per group\n")
	for (i in 1:length(x$f$groupNames))
	{
		cat("\n", x$f$groupNames[i], "\n")
		# Make temporary changes to make x look like a sienaFit object
		# so that sienaFitThetaTable can be applied.
		# This is done in function makeTemp.
		x$theta <- colMeans(x$ThinParameters[first:nmain,i,])
		x$covtheta <- cov(x$ThinParameters[first:nmain,i,])
		tmps <- makeTemp(x)
		tmp <- tmps[[1]]
		tmp1 <- tmps[[2]]
		addtorow <- tmp$addtorow
		for (i in 1:length(tmp1))
		{
			if (length(addtorow$command) > 0)
			{
				for (j in 1:length(addtorow$command))
				{
					ii <- match(i-1, addtorow$pos[[j]])
					if (!is.na(ii))
						if (i == 2 | addtorow$command[j] == 'Network Dynamics')
							cat( addtorow$command[j], '\n')
						else
							cat('\n', addtorow$command[j], '\n', sep='')
				}
			}
			cat(tmp1[i], '\n')
		}
	}
	invisible(x)
}


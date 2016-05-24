#' @nord
setMethod('print', 'mfboot', function(x, ...){
	cat(x$nboot, "bootstrap samples")
    cat("\n", paste(100 * (1 - x$alpha), "%", sep = ""), "confidence interval\n\n")
    cat('\nComparing', x$compare[2], 'to', x$compare[1], '\n')
    print(x$stat)
    cat('\n')
	})
	
	
#' @nord
setMethod('show', 'mfboot', function(object){print(object)})

#' @nord
setMethod('print', 'mfhlboot', function(x, ...){
    cat("\n", x$nboot, " bootstrap samples", sep='')
    cat("\n", paste(100 * (1-x$alpha), "%", sep = ""), " confidence intervals", sep='')
    cat('\nComparing',x$compare[2],'to',x$compare[1],'\n\n')
    cat("\nMitigated Fraction\n\n")
    print(x$MFstat)
    cat("\n\nHodges-Lehmann\n\n")
    print(x$HLstat)
    cat('\n\nQuartile Differences (quartiles of ',x$compare[2],' - quartiles of ',x$compare[1],')\n\n',sep='')
    print(x$QDIFstat)
    cat("\n\nQuartiles of", x$compare[1], "\n")
    print(x$QXstat)
    cat("\n\nQuartiles of", x$compare[2], "\n")
    print(x$QYstat)
    cat('\n')
	})

#' @nord
setMethod('show', 'mfhlboot', function(object){print(object)})

#' @nord
setMethod('print', 'mfmp', function(x, ...){
	cat(x$what, '\n')
	print(x$ci)
	})

#' @nord
setMethod('show', 'mfmp', function(object){print(object)})

#' @nord
setMethod('print', 'mfbootcluster', function(x, ...){
    cat('\n\n', x$what, sep = '')
    cat('\nComparing', x$compare[2], 'to', x$compare[1], '\n')
    cat("\n", paste(100 * (1 - x$alpha), "%", sep = ""), "confidence interval\n\n")
    print(x$stat)
	if(!is.null(x$excludedClusters)){
		exc <- paste(x$excludedClusters, collapse = ', ')
	} else {
		exc <- 'None'
	}
    cat('\nExcluded Clusters\n', exc, '\n')
    cat('\n')
	})
	
#' @nord
setMethod('show', 'mfbootcluster', function(object){print(object)})

#' @nord
setMethod('print', 'mfcluster', function(x, ...){
    cat('\nComparing', x$compare[2], 'to', x$compare[1],'\n')
    cat('\nMF =', x$All$mf, '\n')
    cat('\nBy Cluster\n')  
    byclus <- na.omit(x$byCluster)
    attributes(byclus)$na.action <- NULL
    print(byclus)
    cat('\nAll\n')
    print(x$All)
	if(!is.null(x$excludedClusters)){
		exc <- paste(x$excludedClusters, collapse = ', ')
	} else {
		exc <- 'None'
	}
    cat('\nExcluded Clusters\n', exc, '\n')
    cat('\n')
})

#' @nord
setMethod('show', 'mfcluster', function(object){print(object)})

#' @nord
setMethod('print', 'mfcomponents', function(x, ...){
    cat('\nMF =', x$mf, 'comparing', x$compare[2], 'to', x$compare[1], '\n')
    subj <- data.frame(x$subj)
    mfju <- unique(subj$mf.j)
    nu <- length(mfju)
    subtab <- data.frame(mf.j = mfju, freq = rep(NA, nu), min.y = rep(NA, nu),
		max.y = rep(NA, nu))
    subtab <- subtab[rev(order(subtab$mf.j)), ]
    for(i in 1:nrow(subtab)){
        ys <- subj$y[subj$mf.j == subtab$mf.j[i]]
        subtab$freq[i] <- length(ys)
        subtab$min.y[i] <- min(ys)
        subtab$max.y[i] <- max(ys)
    }
    cat('\nMF Subject Components\n\n')
    print(subtab, row.names = F)
    cat('\n')
	})

#' @nord
setMethod('show', 'mfcomponents', function(object){print(object)})
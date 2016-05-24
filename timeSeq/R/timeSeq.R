# Main function.

# Load needed library
# library(gss)
# library(doParallel)
# library(foreach)
# library(reshape)
# library(edgeR)

timeSeq = function(data, group.label, gene.names, reads = NULL, exon.length = NULL, 
	gene.length = NULL, exon.level = TRUE, gene.level = FALSE, p.values = FALSE, 
	n_cores = NULL, iterations = 10, offset = TRUE) {

	mkran = function(formula, data) {
		with(data, 
            {
			form.wk <- terms.formula(formula)[[2]]
			if (!("|" %in% strsplit(deparse(form.wk), "")[[1]])) 
				stop("gss error in mkran: missing | in grouping formula")
			term.wk <- strsplit(deparse(form.wk), " \\| ")[[1]]
			z2.wk <- eval(parse(text = term.wk[2]))
			if (!is.factor(z2.wk)) 
				stop(paste("gss error in mkran: ", term.wk[2], " should be a factor"))
			z <- NULL
			lvl.z2 <- levels(z2.wk)
			for (i in lvl.z2) z <- cbind(z, leng[1:length(lvl.z2)][as.numeric(i)] * as.numeric(z2.wk == i))
			if (term.wk[1] == "1") {
				init <- 0
				env <- length(levels(z2.wk))
				fun <- function(zeta, env) diag(10^(-zeta), env)
				sigma <- list(fun = fun, env = env)
			} else {
				z1.wk <- eval(parse(text = term.wk[1]))
				if (!is.factor(z1.wk)) 
					stop(paste("gss error in mkran: ", term.wk[1], " should be a factor"))
				ind <- lvl.wk <- NULL
				nz <- length(lvl.z2)
				nsig <- length(levels(z1.wk))
				for (i in levels(z1.wk)) {
					zz.wk <- z2.wk[z1.wk == i, drop = TRUE]
					ind <- c(ind, list((1:nz)[lvl.z2 %in% levels(zz.wk)]))
					lvl.wk <- c(lvl.wk, levels(zz.wk))
				}
				if (max(table(lvl.wk) > 1)) 
					stop("gss error in mkran: ", term.wk[2], " should be nested under ", term.wk[1])
				init <- rep(0, length(levels(z1.wk)))
				env <- list(size = nz, nsig = nsig, ind = ind)
				fun <- function(zeta, env) {
					wk <- rep(0, env$size)
					for (i in 1:env$nsig) wk[env$ind[[i]]] <- 10^(-zeta[i])
					diag(wk)
				}
				sigma <- list(fun = fun, env = env)
			}
			list(z = z, sigma = sigma, init = init)
		    }
        )
	}

    group.length = length(group.label)
    group.1 = length(group.label[group.label == 1])
    group.2 = length(group.label[group.label == 2])
    List = unique(gene.names)

    norm.factor = calcNormFactors(data)
    if (is.null(reads)) effective.lib.size = colSums(data) * norm.factor else effective.lib.size = reads * norm.factor

	desired.cores = n_cores
	if (length(desired.cores) == 0) {
	    max.cores = detectCores()
	    desired.cores = max.cores - 1		
	} 
    cl = makeCluster(desired.cores)
    registerDoParallel(cl)

	fits = foreach(genenames = List, .combine = rbind, .packages = c("reshape", "gss")) %dopar% {

		if (exon.level) {
			gene = matrix(data[gene.names == genenames, ], ncol = group.length)
			n = dim(gene)[1]
			if (n < 1) return(sprintf("Gene: %s has no data!", genenames))
			colnames(gene) = group.label
			rownames(gene) = NULL

			m.dt = melt(gene)
			offset.lib = rep(effective.lib.size, each = n)
			if (!is.null(exon.length)) {
				ndt = data.frame(id = as.factor(m.dt[, 1]), g = as.factor(m.dt[, 2]), t = c(rep(1:group.1, each = n), rep(1:group.2, each = n)), leng = rep(exon.length[gene.names == genenames], group.length), count = m.dt[, 3])
				random_effect = mkran(~1 | id, ndt)
		        if(offset){
        			fit1 = try(gssanova1(count ~ t + g + t:g, offset = log(offset.lib), data = ndt, family = "nbinomial", random = random_effect), silent = TRUE)
		            fit2 = try(gssanova1(count ~ t + g, data = ndt, offset = log(offset.lib), family = "nbinomial", random = random_effect), silent = TRUE)
		        } else {
          			fit1 = try(gssanova1(count ~ t + g + t:g, data = ndt, family = "nbinomial", random = random_effect), silent = TRUE)
			        fit2 = try(gssanova1(count ~ t + g, data = ndt, family = "nbinomial", random = random_effect), silent = TRUE)
        		}
			
			} else {
				ndt <- data.frame(id = as.factor(m.dt[, 1]), g = as.factor(m.dt[, 2]), t = c(rep(1:group.1, each = n), rep(1:group.2, each = n)), 
				count = m.dt[, 3])
		        if(offset){
        			fit1 = try(gssanova1(count ~ t + g + t:g, offset = log(offset.lib), data = ndt, family = "nbinomial", random = ~1 | id), silent = TRUE)
		            fit2 = try(gssanova1(count ~ t + g, offset = log(offset.lib), data = ndt, family = "nbinomial", random = ~1 | id), silent = TRUE)
        		} else {
		            fit1 = try(gssanova1(count ~ t + g + t:g, data = ndt, family = "nbinomial", random = ~1 | id), silent = TRUE)
		            fit2 = try(gssanova1(count ~ t + g, data = ndt, family = "nbinomial", random = ~1 | id), silent = TRUE)
        		}				
			}

			if (("try-error" %in% class(fit1))) NPDE.ratio = NA
			else {
				pj.tc.temp = try(project(fit1, c("t", "g")), silent = TRUE)
				if (("try-error" %in% class(pj.tc.temp))) NPDE.ratio = NA else NPDE.ratio = pj.tc.temp$ratio
			}
			if (("try-error" %in% class(fit2))) PDE.ratio = NA
			else {
				pj.tg.temp = try(project(fit2, c("t")), silent = TRUE)
				if (("try-error" %in% class(pj.tg.temp))) PDE.ratio = NA else PDE.ratio = pj.tg.temp$ratio
			}
		} else if (gene.level) {
			gene = matrix(data[gene.names == genenames, ], ncol=group.length)
			n = dim(gene)[1]
			if (n < 1) return(sprintf("Gene: %s has no data!", genenames))
			colnames(gene) = group.label
			rownames(gene) = NULL

			m.dt = melt(gene)
			offset.lib = rep(effective.lib.size, each = n)
			if (!is.null(gene.length)) {
				ndt = data.frame(id = as.factor(m.dt[, 1]), g = as.factor(m.dt[, 2]), t = c(rep(1:group.1, each = n), rep(1:group.2, each = n)), leng = rep(gene.length[gene.names == genenames], group.length), count = m.dt[, 3])
				random_effect = mkran(~1 | id, ndt)
	        	if(offset){
    	 			fit1 = try(gssanova1(count ~ t + g + t:g, offset = log(offset.lib), data = ndt, family = "nbinomial", random = random_effect), silent = TRUE)
            		fit2 = try(gssanova1(count ~ t + g, data = ndt, offset = log(offset.lib), family = "nbinomial", random = random_effect), silent = TRUE)
				} else {
            		fit1 = try(gssanova1(count ~ t + g + t:g, data = ndt, family = "nbinomial", random = random_effect), silent = TRUE)
	            	fit2 = try(gssanova1(count ~ t + g, data = ndt, family = "nbinomial", random = random_effect), silent = TRUE)
        		}
			} else {
				ndt = data.frame(id = as.factor(m.dt[, 1]), g = as.factor(m.dt[, 2]), t = c(rep(1:group.1, each = n), rep(1:group.2, each = n)), count = m.dt[, 3])
	        	if(offset){
	  	        	fit1 = try(gssanova1(count ~ t + g + t:g, offset = log(offset.lib), data = ndt, family = "nbinomial", random = ~1 | id), silent = TRUE)
		            fit2 = try(gssanova1(count ~ t + g, offset = log(offset.lib), data = ndt, family = "nbinomial", random = ~1 | id), silent = TRUE)
        		} else {
		            fit1 = try(gssanova1(count ~ t + g + t:g, data = ndt, family = "nbinomial", random = ~1 | id), silent = TRUE)
          			fit2 = try(gssanova1(count ~ t + g, data = ndt, family = "nbinomial", random = ~1 | id), silent = TRUE)
        		}
			}

			if (("try-error" %in% class(fit1))) NPDE.ratio = NA
			else {
				pj.tc.temp = try(project(fit1, c("t", "g")), silent = TRUE)
				if (("try-error" %in% class(pj.tc.temp))) NPDE.ratio = NA else NPDE.ratio = pj.tc.temp$ratio
			}
			if (("try-error" %in% class(fit2))) PDE.ratio = NA
			else {
				pj.tg.temp = try(project(fit2, c("t")), silent = TRUE)
				if (("try-error" %in% class(pj.tg.temp))) PDE.ratio = NA else PDE.ratio = pj.tg.temp$ratio
			}
		}
		return(c(genenames, NPDE.ratio, PDE.ratio))
	}

	stopCluster(cl)
	fits <- matrix(fits, ncol = 3)
	
	if (p.values) {
		if (iterations <= 1) {
			cat("Iterations must be a positive integer greater than 1.\n")
			return(cat("ERROR!"))
		}	
		
		pvalue = array(0, c(length(List), 2))
		for (i in 1 : length(List)){
			genenames = List[i]
			gene = matrix(data[gene.names == genenames, ], ncol=group.length)
			out = timeSeq.permutate(gene, genenames, gene.names, group.length, group.1, group.2, group.label, fits[fits[, 1] == genenames], 
									n_cores, iterations, effective.lib.size, exon.length, exon.level, gene.length, gene.level, offset)
			#return(out)						
			pvalue[i, ] = c(out[[1]], out[[2]])
		}
		pvalue = data.frame(List, pvalue)
		colnames(pvalue) = c("genename", "NPDE", "PDE")
	} else pvalue = NULL


	count = numeric(length(List))
	max_length = 0
	for (i in c(1:length(List))) {
		count[i] = dim(matrix(data[gene.names == List[i], ], ncol=group.length))[1]
		if (count[i] > max_length) max_length = count[i]
	}

	table = array(0, c(length(List), max_length, group.length))
	for (i in c(1 : length(List))) {
		m = unlist(data[gene.names == List[i], 1:group.length])
		m = t(matrix(m, c(count[i], group.length)))
		#print(m)
		#cat("\n")
		for (j in c(1:count[i])) table[i, j, ] = m[, j]
	}

	NPDE.ratio = numeric(length(List))
	PDE.ratio = numeric(length(List))
	m = fits
	for (i in 1 : length(List)) {
		NPDE.ratio[i] = as.numeric(m[as.character(m[, 1]) == List[i], 2])
		PDE.ratio[i] = as.numeric(m[as.character(m[, 1]) == List[i], 3])
	}

	###Sort by ratios
	sorted = timeSeq.sort(List, NPDE.ratio, PDE.ratio, table, count)

	out = list(sorted = sorted, 
	           count = count,
	           NPDE.ratio = NPDE.ratio,
               PDE.ratio = PDE.ratio, 
               genenames = List, 
               table = table, 
               data = data, 
               gene.names = gene.names, 
		       effective.lib.size = effective.lib.size, 
		       exon.length = exon.length,
		       group.label = group.label,
               group.length = group.length, 
               group1.length = group.1, 
               group2.length = group.2, 
               exon.level = exon.level, 
		       p.values = p.values, 
               pvalue = pvalue)
               
	class(out) = "timeSeq.obj"
	return(out)
}

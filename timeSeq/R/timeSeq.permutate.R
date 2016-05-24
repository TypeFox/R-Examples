timeSeq.permutate = function(gene_initial, gene_name, gene.names, group.length, group1.length, group2.length, group.label, 
	ratio_initial, n_cores = NULL, iterations = 10, effective.lib.size, exon.length, exon.level, gene.length, gene.level, offset) {
  
	if (iterations <= 1) {
		cat("Iterations must be a positive integer greater than 1.\n")
		return(cat("ERROR!"))
	}	
  
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
	
	desired.cores = n_cores
	if (length(desired.cores) == 0) {
	    max.cores = detectCores()
	    desired.cores = max.cores - 1		
	} 	
	
	cl.perm = makeCluster(desired.cores)
	registerDoParallel(cl.perm)	
	perm.fits = foreach(k = 1 : iterations, .combine = rbind, .packages = c("reshape", "gss")) %dopar% {    
    	if(exon.level){
			##Don't use order as the name of variable. 
			order.perm = 1 : ncol(gene_initial)
        	dim.gene = dim(gene_initial)[1]
 			if (dim.gene < 1) return(cat(paste("Gene has no data!\n")))
            order.perm[1 : group.length] = sample(1 : group.length, group.length)    	
            gene.perm = matrix(gene_initial[, order.perm], ncol = group.length)
        	colnames(gene.perm) = group.label
        	rownames(gene.perm) = NULL
		      
        	m.dt = melt(gene.perm)
        	offset.lib = rep(effective.lib.size[order.perm], each = dim.gene)
        	if(!is.null(exon.length)){
  				ndt = data.frame(id = as.factor(m.dt[,1]), g = as.factor(m.dt[,2]), t = c(rep(1:group1.length,each=dim.gene), rep(1 : group2.length,each = dim.gene)), leng = rep(exon.length[gene.names == gene_name], group.length), count = m.dt[,3])
				random_effect = mkran(~1|id, ndt)
		        if(offset){
          			fit1 = try(gssanova1(count ~ t + g + t: g, offset = log(offset.lib), data = ndt, family = "nbinomial",random = random_effect), silent = TRUE)
          			fit2 = try(gssanova1(count ~ t + g, data = ndt, offset = log(offset.lib), family = "nbinomial",random = random_effect), silent = TRUE) 
		        } else {
          			fit1 = try(gssanova1(count ~ t + g + t: g, data = ndt, family = "nbinomial",random = random_effect), silent = TRUE)
          			fit2 = try(gssanova1(count ~ t + g, data = ndt, family = "nbinomial",random = random_effect), silent = TRUE) 
        		}                               
      		} else {
        		ndt = data.frame(id = as.factor(m.dt[,1]), g = as.factor(m.dt[,2]), t = c(rep(1 : group1.length, each = dim.gene), rep(1 : group2.length, each = dim.gene)), count = m.dt[,3])
	            if(offset){
              		fit1 = try(gssanova1(count ~ t + g + t: g, offset = log( offset.lib), data = ndt, family = "nbinomial", random= ~1|id), silent = TRUE)
              		fit2 = try(gssanova1(count ~ t + g, offset = log( offset.lib), data = ndt, family = "nbinomial",random=~1|id), silent = TRUE)
            	} else {
              		fit1 = try(gssanova1(count ~ t + g + t: g, data = ndt, family = "nbinomial", random= ~1|id), silent = TRUE)
              		fit2 = try(gssanova1(count ~ t + g, data = ndt, family = "nbinomial",random=~1|id), silent = TRUE)
            	}			
            }      
      
			if (( 'try-error' %in% class(fit1) )) NPDE.ratio = NA
        	else{
        		pj.tc.temp = try(project(fit1, c("t","g")), silent = TRUE)
        		if (('try-error' %in% class(pj.tc.temp))) NPDE.ratio = NA else NPDE.ratio = pj.tc.temp$ratio
        	}
      
        	if (( 'try-error' %in% class(fit2) )) PDE.ratio = NA
        	else {
        		pj.tg.temp = try(project(fit2, c("t")), silent = TRUE)
        		if (('try-error' %in% class(pj.tg.temp))) PDE.ratio = NA else PDE.ratio = pj.tg.temp$ratio        
        	}
        } else if(gene.level){
      
        	order.perm = 1 : ncol(gene_initial)
        	dim.gene = dim(gene_initial)[1]
 			if (dim.gene < 1) return(cat(paste("Gene has no data!\n")))
	        order.perm[1 : group.length] = sample(1 : group.length, group.length)      
        	gene.perm = matrix(gene_initial[, order.perm], ncol = group.length)

      		colnames(gene.perm) = group.label
      	  	rownames(gene.perm) = NULL
			m.dt = melt(gene.perm)
			offset.lib = rep(effective.lib.size[order.perm], each = dim.gene)
			if(!is.null(gene.length)){
				ndt = data.frame(id = as.factor(m.dt[,1]), g = as.factor(m.dt[,2]), t = c(rep(1 : group1.length, each = dim.gene), rep(1 : group2.length, each = dim.gene)), leng = rep(gene.length[gene.names == gene_name], group.length), count = m.dt[,3])
                random_effect = mkran(~1|id, ndt)
            	if(offset){
            		fit1 = try(gssanova1(count ~ t + g + t: g, offset = log(offset.lib), data = ndt, family = "nbinomial",random = random_effect), silent = TRUE)
            		fit2 = try(gssanova1(count ~ t + g, data = ndt, offset = log(offset.lib), family = "nbinomial",random = random_effect), silent = TRUE)
            	} else {
            		fit1 = try(gssanova1(count ~ t + g + t: g, data = ndt, family = "nbinomial",random = random_effect), silent = TRUE)
            		fit2 = try(gssanova1(count ~ t + g, data = ndt, family = "nbinomial",random = random_effect), silent = TRUE)
            	}        		                                
			} else {
        		ndt = data.frame(id = as.factor(m.dt[,1]), g = as.factor(m.dt[,2]), t = c(rep(1 : group1.length, each = dim.gene), rep(1 : group2.length, each = dim.gene)), count = m.dt[,3])
	            if(offset){
      	        	fit1 = try(gssanova1(count ~ t + g + t: g, offset = log( offset.lib), data = ndt, family = "nbinomial",random= ~1|id), silent = TRUE)
              		fit2 = try(gssanova1(count ~ t + g, offset = log( offset.lib), data = ndt, family = "nbinomial",random= ~1|id), silent = TRUE)
	            } else {
                	fit1 = try(gssanova1(count ~ t + g + t: g, data = ndt, family = "nbinomial",random= ~1|id), silent = TRUE)
              		fit2 = try(gssanova1(count ~ t + g, data = ndt, family = "nbinomial",random= ~1|id), silent = TRUE)
            	}		
     	    }      
      
			if (( 'try-error' %in% class(fit1) )) NPDE.ratio = NA
			else {
				pj.tc.temp = try(project(fit1, c("t","g")), silent = TRUE)
				if (('try-error' %in% class(pj.tc.temp))) NPDE.ratio = NA else NPDE.ratio = pj.tc.temp$ratio
			}      
			if (( 'try-error' %in% class(fit2) )) PDE.ratio = NA
      	  	else {
        		pj.tg.temp = try(project(fit2, c("t")), silent = TRUE)
				if (('try-error' %in% class(pj.tg.temp))) PDE.ratio = NA else PDE.ratio = pj.tg.temp$ratio
      	  	}      
        }		
		return(c(NPDE.ratio, PDE.ratio))    
    }
  
	stopCluster(cl.perm)

    ratio.initial = as.matrix(ratio_initial)
    NPDE.state = as.numeric(ratio.initial[2,])
    PDE.state = as.numeric(ratio.initial[3,])
    perm.fits = as.matrix(perm.fits, byrow = TRUE)
    
    if ((!is.na(NPDE.state)) && (sum(!is.na(as.numeric(perm.fits[,1]))) >= 1)) 
		pvalue.NPDE = sum(as.numeric(perm.fits[,1][!is.na(as.numeric(perm.fits[,1]))]) > NPDE.state) / sum(!is.na(as.numeric(perm.fits[,1])))
    else pvalue.NPDE = NA
   
    if ((!is.na(PDE.state)) && (sum(!is.na(as.numeric(perm.fits[,2]))) >= 1))
		pvalue.PDE = sum(as.numeric(perm.fits[,2][!is.na(as.numeric(perm.fits[,2]))]) > PDE.state) / sum(!is.na(as.numeric(perm.fits[,2])))
    else pvalue.PDE = NA
		
    out = list(pvalue.NPDE, pvalue.PDE)
    return(out)
}


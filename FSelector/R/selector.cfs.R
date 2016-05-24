### CFS
# classification and regression
# continous and discrete data
cfs <- function(formula, data) {
	cont_correlation <- function(a, b) {
		result = 0
		if(!is.factor(a) && !is.factor(b)) { # both continous
			complete = complete.cases(a) & complete.cases(b)
			if(!any(complete))
				return(0)
			vec1 = a[complete]
			vec2 = b[complete]
			if(sd(vec1) == 0 || sd(vec2) == 0)
				return(0)
			result = cor(vec1, vec2)
		} else if(is.factor(a) && is.factor(b)) { # both discrete
			tab = table(a, b)
			alevels = rownames(tab)
			blevels = colnames(tab)
			
			result = sum(sapply(alevels, function(avalue) {
					avec = as.numeric(a == avalue)
					complete_a = complete.cases(a)
					return(sum(sapply(blevels, function(bvalue) {
							bvec = as.numeric(b == bvalue)
							complete_b = complete.cases(b)
							complete = complete_a & complete_b
							avec_complete_data = avec[complete]
							bvec_complete_data = bvec[complete]
							if(sd(avec_complete_data, na.rm=TRUE) == 0 || sd(bvec_complete_data, na.rm=TRUE) == 0)
								return(0)
							return(tab[avalue, bvalue] / length(a) * cor(avec_complete_data, bvec_complete_data))
						})))
				}))
		} else { # continous and discrete
			cont = NULL;
			disc = NULL;
			if(is.factor(a)) {
				cont = b
				disc = a
			} else {
				cont = a
				disc = b
			}
			
			cont_complete = complete.cases(cont)
			disc_table = table(disc)
			disc_levels = names(disc_table)
			
			if(length(disc_levels) == 0) {
				result = 0
			} else {
				result = sum(sapply(disc_levels, function(lev) {
						disc_vec = as.numeric(disc == lev)
						disc_vec_complete = complete.cases(disc_vec)
						complete = cont_complete & disc_vec_complete
						disc_vec_complete_data = disc_vec[complete]
						cont_complete_data = cont[complete]
						if(sd(cont_complete_data) == 0)
							return(0)
						return(disc_table[lev] / length(disc) * cor(disc_vec_complete_data, cont_complete_data))
					}))
			}
		}
		return(result)
	}

	# uses parent.env (correlations)
	get_correlation <- function(attr1, attr2, classification, new_data, entropies) {
		#lazy evaluation
		if(!is.na(correlations[attr1, attr2])) {
			return(correlations[attr1, attr2])
		}
		
		tmp_res = NA
		if(classification) { #discrete class
			tmp_res = 2.0 * (entropies[attr1] + entropies[attr2] - entropyHelper(data.frame(cbind(new_data[[attr1]], new_data[[attr2]])))) / (entropies[attr1] + entropies[attr2])
		} else { #continous class
			tmp_res = cont_correlation(new_data[[attr1]], new_data[[attr2]])
		}
        if(is.nan(tmp_res)) {
            # all entropies (individual + joint) are 0
            tmp_res = 0
        }

		correlations[attr1, attr2] <<- tmp_res
		correlations[attr2, attr1] <<- tmp_res
		
		return(tmp_res)
	}

	# uses parent.env
	evaluator <- function(attrs) {
		ff_sum = 0
		ff_count = 0
		fc_sum = 0
		attr_count = length(attrs)
		
		if(attr_count <= 0)
			stop("Attributes not specified")
		
		for(i in 1:attr_count) {
			attr1 = attrs[i]
			
			# feature-class correlation
            cor = get_correlation(attr1, 1, classification, new_data, entropies)
            fc_sum = fc_sum + cor
			
			# feature-feature correlation
			if(i == attr_count) {
				next()
			}
			for(j in (i+1):attr_count) {
				attr2 = attrs[j]
                cor = get_correlation(attr1, attr2, classification, new_data, entropies)
                ff_count = ff_count + 1
                ff_sum = ff_sum + cor
			}
		}

		ff_cor = ff_sum / ff_count
		fc_cor = fc_sum / attr_count
		
		if(attr_count == 1)
			return(fc_cor)
		else
			return(attr_count * fc_cor / sqrt(attr_count + attr_count * (attr_count - 1) * ff_cor))
	}
	
	new_data = get.data.frame.from.formula(formula, data)
	
	# prepare correlation matrix
	classification = is.factor(new_data[[1]])
	attr_count = dim(new_data)[2]
	attr_names = colnames(new_data)
	correlations = matrix(rep(NA, attr_count ^ 2), nrow = attr_count, ncol = attr_count,
		dimnames = list(attr_names, attr_names))
	
	entropies = NULL
	if(classification) {
		new_data = supervised.discretization(formula, data = new_data)
		new_data = get.data.frame.from.formula(formula, new_data)
		entropies = sapply(new_data, entropyHelper)		
	}

	result = best.first.search(names(new_data)[-1], evaluator)
		
	return(result)
}

### RELIEF
# classification and regression
# continous and discrete data
relief <- function(formula, data, neighbours.count = 5, sample.size = 10) {
	# uses parent.env
	find_neighbours <- function(instance_idx) {
		instance = new_data[instance_idx,, drop = FALSE]

		# for every other instance
		for(current_idx in 1:instances_count) {
			if(instance_idx == current_idx)
				next()
			current_instance = new_data[current_idx,, drop = FALSE]
			if(is.na(current_instance[1, 1]))
				next()
			
			dist = instance_distance(instance, current_instance)
			
			if(classification)
				class_no = which(classes == current_instance[[1]])
			else
				class_no = 1
			if(nn_stored_count[class_no] < neighbours.count) {
				nn_stored_count[class_no] <<- nn_stored_count[class_no] + 1
				n_array[class_no, nn_stored_count[class_no], ] <<- c(dist, current_idx)
			} else {
				max_idx = which.max(n_array[class_no, , 1])
				max_value = n_array[class_no, max_idx, 1]
				if(dist < max_value) {
					n_array[class_no, max_idx, ] <<- c(dist, current_idx)
				}
			}
		}
	}
	
	# uses parent.env
	update_weights <- function(instance_idx) {
		instance = new_data[instance_idx,, drop = FALSE]
		instance_class = instance[1, 1]
		instance_class_no = which(classes == instance_class)
		
		if(classification) {
			# for each attribute
			for(attr_idx in 1:attributes_count) {
				col_idx = attr_idx + 1
				
				# nearest hits
				hits_sum = 0
				if(nn_stored_count[instance_class_no] > 0) {
					hits_sum = sum(sapply(1:nn_stored_count[instance_class_no], function(n_idx) {
							n_instance_idx = n_array[instance_class_no, n_idx, 2]
							n_instance = new_data[n_instance_idx,, drop = FALSE]
							return(field_distance(col_idx, instance, n_instance))
						}))
					hits_sum = hits_sum / nn_stored_count[instance_class_no]
				}
				
				# nearest misses
				misses_sum = 0
				if(class_count > 1) {
					misses_sum = sum(sapply((1:class_count)[-instance_class_no], function(class_no) {
							class_misses_sum = 0
							if(nn_stored_count[class_no] > 0) {
								class_misses_sum = sum(sapply(1:nn_stored_count[class_no], function(n_idx) {
										n_instance_idx = n_array[class_no, n_idx, 2]
										n_instance = new_data[n_instance_idx,, drop = FALSE]
										return(field_distance(col_idx, instance, n_instance))
									}))
								class_misses_sum = class_misses_sum * class_prob[class_no] / nn_stored_count[class_no]
							}
							return(class_misses_sum)
						}))
					
					
					misses_sum = misses_sum / (1 - class_prob[instance_class_no])
				}
				results[attr_idx] <<- results[attr_idx] - hits_sum + misses_sum
			}
		} else {
			if(nn_stored_count[1] > 0) {
				for(n_idx in 1:nn_stored_count[1]) {
					n_instance_idx = n_array[1, n_idx, 2]
					n_instance = new_data[n_instance_idx,, drop = FALSE]
					class_diff = field_distance(1, instance, n_instance)
					ndc <<- ndc + class_diff / nn_stored_count[1]
					for(attr_idx in 1:attributes_count) {
						col_idx = attr_idx + 1
						attr_diff_norm = field_distance(col_idx, instance, n_instance) / nn_stored_count[1]
						nda[attr_idx] <<- nda[attr_idx] + attr_diff_norm
						ndcda[attr_idx] <<- ndcda[attr_idx] + class_diff * attr_diff_norm
					}
				}
			}
		}
	}
	
	# parameters: data.frame, data.frame
	instance_distance <- function(instance1, instance2) {
		len = dim(instance1)[2]
		if(len != dim(instance2)[2])
			stop("Instances of different lengths")
		if(len <= 1)
			stop("Too few attributes")
		
		result = sapply(2:len, function(i) {
				return(field_distance(i, instance1, instance2))
			})
		#return(sqrt(sum(result ^ 2))) #sqrt not needed
		res = sum(result ^ 2)
		if(is.na(res)) {
			stop("Internal error. Distance NA.")
		}
		return(res)
	}
	
	# uses parent.env
	# parameters: index, data.frame, data.frame
	field_distance <- function(col_idx, instance1, instance2) {
		value1 = instance1[1, col_idx]
		value2 = instance2[1, col_idx]
		attr_idx = col_idx - 1 # skip class
		
		if(is.factor(value1) && is.factor(value2)) {
			if(is.na(value1) && is.na(value2)) {
				if(classification)
					return(1 - sum(p_val_in_class[[attr_idx]][, instance1[1, 1]] * p_val_in_class[[attr_idx]][, instance2[1, 1]]))
				else
					return(1 - p_same_val[[attr_idx]])
			} else if(is.na(value1) || is.na(value2)) {
				if(is.na(value1)) {
					known_value = value2
					unknown_class = instance1[1, 1]
				} else {
					known_value = value1
					unknown_class = instance2[1, 1]
				}
				if(classification)
					return(1 - p_val_in_class[[attr_idx]][known_value, unknown_class])
				else
					return(1 - p_val[[attr_idx]][known_value])
			} else if(value1 == value2) {
				return(0)
			} else { #if(value1 != value2)
				return(1)
			}
		} else if(is.numeric(value1) && is.numeric(value2)) {
			if(is.na(value1) && is.na(value2)) {
				return(1)
			} else if(is.na(value1)) {
				return(max(value2, 1 - value2))
			} else if(is.na(value2)) {
				return(max(value1, 1 - value1))
			} else {
				return(abs(value1 - value2))
			}		
		} else {
			stop("Unsupported value type")
		}
	}

	new_data = get.data.frame.from.formula(formula, data)
	new_data = normalize.min.max(new_data)
	
	# for discrete classes
	class_vector = NULL
	class_count = NULL
	class_prob = NULL
	classes = NULL
	p_val_in_class = NULL
	p_val = NULL
	p_same_val = NULL
	
	# for continous class
	ndc = NULL
	nda = NULL
	ndcda = NULL
	
	results = NULL
	n_array = NULL
	nn_stored_count = NULL
	classification = NULL
	sample_instances_idx = NULL

	instances_count = dim(new_data)[1]
	attributes_count = dim(new_data)[2] - 1
	attr_names = dimnames(new_data)[[2]][-1]
	
	if(neighbours.count < 1) {
		neighbours.count = 1
		warning(paste("Assumed: neighbours.count = ", neighbours.count))
	}

	if(sample.size < 1) {
		warning(paste("Assumed: sample.size = ", sample.size))
		sample.size = 1
		sample_instances_idx = sample(1:instances_count, 1)
	} else if(sample.size > instances_count) {
		warning(paste("Assumed: sample.size = ", sample.size))
		sample.size = instances_count
		sample_instances_idx = 1:instances_count
	} else {
		sample_instances_idx = sort(sample(1:instances_count, sample.size, replace=TRUE))
	}
	
	classification = is.factor(new_data[[1]])
	if(classification) {
		class_vector = new_data[[1]]
		class_prob = table(class_vector)
		class_prob = class_prob / sum(class_prob)
		classes = names(class_prob)
		class_count = length(classes)
		
		p_val_in_class = lapply(new_data[-1], function(vec) {
				if(!is.factor(vec) || !any(is.na(vec)))
					return(NULL)
				tab = table(vec, class_vector)
				return(apply(tab, 2, function(x) {
						s = sum(x) 
						if(s == 0)
							return(x)
						else
							return(x / s)
					}))
			})
	} else {
		class_count = 1
		ndc = 0
		nda = array(0, attributes_count)
		ndcda = array(0, attributes_count)
	
		p_val = lapply(new_data[-1], function(vec) {
				if(!is.factor(vec) || !any(is.na(vec)))
					return(NULL)
				tab = table(vec)
				if(sum(tab) != 0) {
					tab = tab / sum(tab)
				}
				return(tab)
			})
		p_same_val = lapply(p_val, function(attr) {
				if(is.null(attr))
					return(NULL)
				return(sum(attr ^ 2))
			})
	}

	n_array = array(0, c(class_count, neighbours.count, 2))
	nn_stored_count = array(0, class_count)
	results = rep(0, attributes_count)

	sapply(sample_instances_idx, function(current_instance_idx) {
		current_instance = new_data[current_instance_idx,, drop = FALSE]
		if(is.na(current_instance[[1]]))
			return(NULL)
		
		nn_stored_count[] <<- 0
		n_array[] <<- Inf
		find_neighbours(current_instance_idx)
		update_weights(current_instance_idx)
	})
	

	if(classification) {
		results = results / sample.size
		return(data.frame(attr_importance = results, row.names = attr_names))
	} else {
		results = ndcda / ndc - ((nda - ndcda) / (sample.size - ndc))
		results = data.frame(attr_importance = results, row.names = attr_names)
		#results = normalize.min.max(results)
		return(results)
	}
	
}

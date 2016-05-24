###########################################
# odebug
###########################################
# Analyzes ovariable outputs with respect to indices and dependencies. 
###########################################

odebug <- function(x, variance = FALSE) {
	if (nrow(x@output)==0) x <- EvalOutput(x)
	
	out <- list() 
	
	# Item 1 - output lenghts
	
	out[["output_rows"]] <- list()
	out[["output_rows"]][[x@name]] <- nrow(x@output)
	
	if (nrow(x@dependencies)>0) {
		for (i in x@dependencies$Name){
			if (class(get(i)) == "ovariable") {
				out[["output_rows"]][[i]] <- nrow(get(i)@output)
			}
		}
	}
	
	# Item 2 - NAs
	
	out[["NAs"]] <- list()
	out[["NAs"]][[x@name]] <- list()
	out[["NAs"]][[x@name]][["total"]] <- sum(is.na(result(x)))
	
	if (nrow(x@dependencies)>0) {
		for (i in x@dependencies$Name){
			out[["NAs"]][[i]] <- list()
			if (class(get(i)) == "ovariable") {
				out[["NAs"]][[i]][["total"]] <- sum(is.na(result(get(i))))
			}
		}
	}
	
	# Item 3 - marginals
	
	out[["marginals"]] <- list()
	out[["marginals"]][[x@name]] <- colnames(x@output)[x@marginal]
	
	if (nrow(x@dependencies)>0) {
		for (i in x@dependencies$Name){
			if (class(get(i)) == "ovariable") {
				out[["marginals"]][[i]] <- colnames(get(i)@output)[get(i)@marginal]
			}
		}
		
		common_marginals <- NULL
		all_marginals <- NULL
		for (i in out[["marginals"]]) {
			#if (length(common_marginals) == 0) {
			#	common_marginals <- i 
			#} else { 
			common_marginals <- intersect(common_marginals, i)
			#}
			all_marginals <- union(all_marginals, i)
		}
		
		matching_marginals <- NULL
		for (i in out[["marginals"]]) {
			for (j in out[["marginals"]]) {
				matching_marginals <- union(matching_marginals, intersect(all_marginals, i))
			}
		}
		
		out[["marginals"]][["all"]] <- all_marginals
		out[["marginals"]][["common"]] <- common_marginals
		out[["marginals"]][["matching"]] <- matching_marginals
		
		# Item 4 - missing locations in common marginals
		
		locs <- list()
		missing <- list()
		for (j in x@dependencies$Name){
			missing[[j]] <- list()
		}
		for (i in common_marginals) {
			locs[[i]] <- NULL
			for (j in x@dependencies$Name){
				if (class(get(j)) == "ovariable") {
					locs[[i]] <- union(locs[[i]], get(j)@output[[i]])
				}
			}
			for (j in x@dependencies$Name){
				if (class(get(j)) == "ovariable") {
					missing[[j]][[i]] <- setdiff(locs[[i]], get(j)@output[[i]])
				}
			}
		}
		
		out[["marginal_index_locations"]] <- locs
		out[["missing_locations"]] <- missing
		
		# Item 4.1 - mispelt index names
	}
	
	# Item 5 - marginal variance analysis
	
	if (variance) {
		# Create second order combinations of marginals in x
		margs <- out[["marginals"]][[1]]
		margs <- margs[margs!="Iter"] # Iteration number needs to be taken out
		combinations <- combn(margs, 2)
		combinations <- apply(combinations, 2, paste, collapse = "*")
		combinations <- paste(combinations, collapse = "+")
		
		out[["variance_analysis"]] <- 
		summary(
			aov(
				as.formula(
					paste(
						"result(x) ~ ", 
						combinations, 
						sep = ""
					)
				), 
				data = x@output
			)
		)
	}
	
	return(out)
}


#test <- data.frame(A = c("x","y","z"), B = rep(c("a","b","c"), each = 3), C = rep(c("1","2","3"), each = 3*3), Result = runif(1*3*3*3))
#test.aov <- aov(as.formula("Result ~ A*B + A*C + B*C"), data=test)
#summary(test.aov)

#test <- data.frame(A = c("x","y","z"), B = rep(c("a","b","c"), each = 3), C = rep(c("1","2","3"), each = 3*3), testResult = runif(1*3*3*3))
#test <- Ovariable(name = "test", output = test, marginal = c(T, T, T, F))
#odebug(test, variance = TRUE)

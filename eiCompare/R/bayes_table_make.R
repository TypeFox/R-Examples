bayes_table_make <-
function (ei_bayes_object, cand_vector, table_names) {

	# Used for Later Sorting/Colnames
	seq_split <- 2:length(cand_vector)
	rn <- c(insert(cand_vector, ats= seq_split,values=rep("se",length(cand_vector)-1)), "se")
	# Summarize Bayes Object to get posterior means/devs
	ei_bayes_object <- summary(ei_bayes_object)
	means <- ei_bayes_object$coef # get the estimates
	means <- data.frame(means[, "Mean"], means[, "Std. Dev."])
	means <- t(means) # Transpose it

	# Have to break apart the data to put in correct order
	list_holder <- list()

	for (i in 1:length(cand_vector)) {
		subs <- grep(cand_vector[i],colnames(means), value=T) # use grep() to collect appropriate subsetted column names
		subs_data <- means[,subs] # Then extract that data and put into list
		colnames(subs_data) <- table_names # Need to put on same column names for rbind() later
		list_holder[[i]] <- subs_data
	}
	# LDPLY puts lists together into table
	out <- ldply(list_holder, rbind)*100
	out <- data.frame(rn, out) # Add on column of names

	# Adding on Total Row	
	tot <- colSums(out[seq(1,nrow(out),2),2:ncol(out)])
	just_data <- out[,2:ncol(out)]
	add <- rbind(just_data, tot)
	add <- data.frame(1:nrow(add), add)
	colnames(add) <- c("Candidate", table_names)
	add[,1] <- c(as.character(out[,1]), "Total")
	out <- add 
	
	return(out)

}

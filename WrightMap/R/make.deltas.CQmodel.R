make.deltas.CQmodel <-
function(item.params, item.table = NULL, interactions = NULL, step.table = NULL, item.sign = NULL, inter.sign = NULL, 
	step.sign = NULL, ...) {
	#print("deltas")
	RMP <- item.params$RMP
	if (is.null(item.table)) {
		parts <- unique(unlist(item.params$parts))
		#print(parts)
		if (length(parts) > 2) 
			stop("Please specify an item.table as well as the interactions and/or step.table")
		tables <- names(RMP)
		interactions.at <- grep("\\*", tables)
		if(length(interactions.at) > 0) {
		interactions <- tables[interactions.at]
		tables <- tables[-interactions.at]
		}
		item.table <- tables[1]
		if (length(tables) == 2) 
			step.table <- tables[2]
		else if (length(tables) > 2) 
			stop("Please specify tables")
	}

	eqn <- item.params$equation

	if (is.null(item.sign)) {
		item.sign <- ifelse(grepl(paste("-", item.table, sep = ""), eqn), -1, 1)
		#print(item.table)
		#print(eqn)
	}
	item.name <- ifelse(item.params$imported,"Parameters",item.table)
	#print(item.params$imported)
	item.table <- RMP[[item.table]]
	throlds = item.table$est
	throlds <- throlds[!is.na(throlds)]
	
		item.names <- unlist(item.table[item.name])
	
	if (item.name == "step") 
		item.names <- item.names[item.names != 0]

	if (!is.null(step.table) || !is.null(interactions)) {
	
		if (!is.null(step.table)) {
			if (is.null(step.sign)) {
				step.sign <- ifelse(grepl(paste("-", step.table, sep = ""), eqn), -1, 1)
			}
			step.name <- step.table
			step.table <- RMP[[step.table]]
			steps <- step.table$est
			steps <- steps[!is.na(steps)]
		} else {
			cross.parts <- unlist(strsplit(interactions, "\\*"))
			step.name <- cross.parts[cross.parts != item.name]
			steps <- 0
			step.sign <- 1
		}


		if (!is.null(interactions)) {
			if (is.null(inter.sign)) {
				inter.sign <- ifelse(grepl(paste("-", interactions, sep = ""), eqn), -1, 1)
			}
			inter.name <- interactions
			interactions <- RMP[[interactions]]
			if (step.name == "step") 
				step.col = "step"
			else step.col = paste("n", step.name, sep = "_")
			if (item.name == "step") 
				item.col = "step"
			else item.col = paste("n", item.name, sep = "_")

			crosses <- reshape(interactions[c(item.col, step.col, "est")], direction = "wide", timevar = step.col, idvar = item.col)
			crosses <- crosses[colSums(!is.na(crosses)) != 0]
			crosses <- crosses[rowSums(!is.na(crosses)) > 1, ]
		} else {
			crosses <- 0
			inter.sign <- 1
		}
		throlds <- make.deltas(throlds, crosses, steps, item.sign, step.sign, inter.sign)
	
		if (!is.null(step.table)) 
			step.names <- unlist(step.table[step.name])
		else step.names <- unique(unlist(interactions[step.name]))
		if (step.name == "step") 
			step.names <- step.names[step.names != 0]
		colnames(throlds) <- step.names
		rownames(throlds) <- item.names
	}
	else
		names(throlds) <- item.names


	
	
	

	# message("Using ", item.name, ifelse(!is.null(interactions), paste(" and", inter.name), ""), ifelse(!is.null(step.table), paste(" and", 
		# step.name), ""), " tables to create delta parameters")

	return(throlds)

}

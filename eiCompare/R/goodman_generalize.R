goodman_generalize <-
function(cand_vector, race_group, total, data, table_names, ...) {
	
	# Functions
	list_extract <- function(x) x[,1:2] # sends to lapply to extract indiv column estimates
	
	# Table/Output Row Labeling
	seq_split <- 2:length(cand_vector)
	rn <- c(insert(cand_vector, ats= seq_split,values=rep("se",length(cand_vector)-1)), "se")
	
	# Remove any missing datas
	data <- na.omit(data) 

	#Loop Placeholder
	race_group_table <- list()

	# Loop over Race Vector
	for (k in 1:length(race_group)) {
	
		cand_table <- list()
		for (i in 1:length(cand_vector)) {
			form <- formula(paste(cand_vector[i], race_group[k], "+", total)) 			
			summary(res <-  lm(form, data = data, ...) )
			vote_pct <- coef(res)[1] + coef(res)[2]
			ste <- coef(summary(res))[, "Std. Error"]
			vote_ste <- ste[1] + ste[2]
			cand_table[[i]] <- c(vote_pct, vote_ste)*100
		}

		cand_table <- unlist(cand_table) # cand_table is for one racial group and all candidates
		cand_table <- data.frame(rn, cand_table) # Add in vector for labeling
		race_group_table[[k]] <- cand_table # Put candidate results into list

	}
	
	if(length(race_group) == 1) { # For when there is just % Minority vs. % White, for example
	
		race_group_table <- data.frame(race_group_table)
	
	} else{ # For when there are multiple groups (e.g., pct_hisp, pct_asian, pct_white
		race_group_table <- data.frame( race_group_table ) # Turn list into data.frame
		race_group_table <- race_group_table[,c(1,seq(2,ncol(race_group_table),2))] # clean up table

	}
  # Adding on Total Row	
	tot <- colSums(race_group_table[seq(1,nrow(race_group_table),2),2:ncol(race_group_table)])
	just_data <- race_group_table[,2:ncol(race_group_table)]
	add <- rbind(just_data, tot)
	add <- data.frame(1:nrow(add), add)
	colnames(add) <- c("Candidate", table_names)
	add[,1] <- c(as.character(race_group_table[,1]), "Total")
	race_group_table <- add 
	
	return(race_group_table)
	
}

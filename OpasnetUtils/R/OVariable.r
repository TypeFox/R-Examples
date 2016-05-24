# SETCLASS OVARIABLE ################### Defines the S4 class "ovariable" which is the basic building block in open assessments.
setClass(
	"ovariable", 
	representation(
		name			= "character",
		output			= "data.frame", 
		data			= "data.frame", 
		marginal		= "logical", 
		formula			= "function", 
		dependencies	= "data.frame",
		ddata			= "character"
	),
	prototype = prototype(
		name			= character(),
		output			= data.frame(),
		data			= data.frame(),
		marginal		= logical(),
		formula			= function(...){0},
		dependencies	= data.frame(),
		ddata			= character()
	)
)

####################
# Ovariable (constructor)
#################
# Constructs an ovariable and optionally downloads ddata and saves the variable for use in other codes on the server. 
# The point of this constructor is to simplify variable creation: where before many different functions would have to
# be used to get data into a usable format now a simple call Ovariable(<variable_name>, ddata = <page_ident>, save = TRUE)
# will do the trick.
#####################
Ovariable <- function(
		name = character(), 
		data = data.frame(), 
		formula = function(...){0}, 
		dependencies = data.frame(), 
		ddata = character(),
		output = data.frame(),
		marginal = logical(),
		subset = character(), # add subset postfix to ddata e.g. "Op_enXXXX" -> "Op_enXXXX.subset"
		getddata = TRUE, # download dynamic data immediately (as opposed to waiting until evaluation)
		save = FALSE, # save on server using objects.put
		public = TRUE, # if save = TRUE, use objects.store instead to make it publicly available
		...
) {
	if (length(subset) > 0) ddata <- paste(ddata, opbase.sanitize_subset_name(subset), sep='.')
	
	out <- new(
			"ovariable",
			name = name,
			data = data,
			formula = formula,
			dependencies = dependencies,
			ddata = ddata,
			output = output,
			marginal = marginal
	)
	if (getddata) out <- ddata_apply(out)
	if (save){
		assign(name, out)
		if (public) objects.store(list = name, ...) else objects.put(list = name, ...)
	}
	return(out)
}

# SETMETHOD MATH ################### Math defines basic mathematical operations (log, exp, abs, ...) for ovariables
setMethod(
		f = "Math", 
		signature = signature(x = "ovariable"), 
		definition = function(x) {
			test <- paste(x@name, "Result", sep = "") %in% colnames(x@output)
			rescol <- ifelse(test, paste(x@name, "Result", sep = ""), "Result")
			x@output[[rescol]] <- callGeneric(x@output[[rescol]])
			return(x)
		}
)

# SETMETHOD OPS ##################
#########################################################################################
# Arithmetic operations of ovariables: first they are merged by index columns,
# then the operation is performed for the respective Result columns.
# If one of the expressions is numeric, it is first converted to an ovariable.
#########################################################################################

setMethod(
		f = "Ops", 
		signature = signature(e1 = "ovariable", e2 = "ovariable"), 
		definition = function(e1, e2) {
			
			# EvalOutput if not done yet
			
			if(nrow(e1@output) == 0) e1 <- EvalOutput(e1)
			if(nrow(e2@output) == 0) e2 <- EvalOutput(e2)
			
			# First check presence of name specific Result-columns
			
			test1 <- "Result" %in% colnames(e1@output)
			test2 <- "Result" %in% colnames(e2@output)
			
			test3 <- paste(e1@name, "Result", sep = "") %in% colnames(e1@output)
			test4 <- paste(e2@name, "Result", sep = "") %in% colnames(e2@output)
			
			# If found take note
			
			rescol1 <- ifelse(!test1, paste(e1@name, "Result", sep = ""), "Result")
			rescol2 <- ifelse(!test2, paste(e2@name, "Result", sep = ""), "Result")
			
			if(!(test1 | test3) | !(test2 | test4)) stop("No result column found while operating mathematically with ovariables!\n")
			
			#if (!(test1 & test2)) {
			#	rescol1 <- "Result.x"
			#	rescol2 <- "Result.y"
			#}
			
			# If not change prefixless Result to Result.x/y
			exl_list <- ""
			if (test1) {
				colnames(e1@output)[colnames(e1@output)=="Result"] <- "Result.x"
				rescol1 <- "Result.x"
				exl_list <- c(exl_list, "Result.x")
			}
			
			if (test2) {
				colnames(e2@output)[colnames(e2@output)=="Result"] <- "Result.y"
				rescol2 <- "Result.y"
				exl_list <- c(exl_list, "Result.y")
			}
			
			test5 <- rescol1 == rescol2
			if (test5) {
				colnames(e2@output)[colnames(e2@output)==rescol2] <- "Result.y"
				rescol2 <- "Result.y"
				exl_list <- c(exl_list, "Result.y")
			}
			
			# Now merging should be possible without any confusion
			
			out <- merge(e1, e2)@output
			
			# Call generic function on the two Result-columns
			
			out$Result <- callGeneric(out[[rescol1]], out[[rescol2]])
			
			out <- new(
					"ovariable",
					#	dependencies = data.frame(Name = c(e1@name, e2@name)),
					output = out[
							!colnames(out) %in% exl_list | colnames(out) == "Result"
					]
			)
			out <- CheckMarginals(out, deps = list(e1, e2), verbose = FALSE)
			return(out)
		}
)

setMethod(
		f = "Ops", 
		signature = signature(e1 = "ovariable", e2 = "numeric"), 
		definition = function(e1, e2) {
			#e2 <- new("ovariable", output = data.frame(Result = e2))
			#out <- callGeneric(e1, e2) # Call above definition
			#out@name <- e1@name
			result(e1) <- callGeneric(result(e1), e2)
			return(e1)
		}
)

setMethod(
		f = "Ops", 
		signature = signature(e1 = "numeric", e2 = "ovariable"), 
		definition = function(e1, e2) {
			#e1 <- new("ovariable", output = data.frame(Result = e1))
			#out <- callGeneric(e1, e2) # Call above definition
			#out@name <- e2@name
			result(e2) <- callGeneric(e1, result(e2))
			return(e2)
		}
)

# SETMETHOD MERGE ########### merge ovariable outputs

setMethod(f = "merge", 
		signature = signature(x = "ovariable", y = "ovariable"),
		definition = function(x, y, all = FALSE, sort = FALSE, ...) {
			if (nrow(x@output) == 0) stop("X output missing!")
			if (nrow(y@output) == 0) stop("Y output missing!")
			
			x@output <- dropall(x@output)
			y@output <- dropall(y@output)
			
			temp <- fill.na.merge(x, y)
			x <- temp[[1]]
			y <- temp[[2]]
			
			#temp <- merge(x@output, y@output, all = all, sort = sort, ...)#, by = test)
			by_auto <- intersect(colnames(x@output), colnames(y@output))
			
			# Unkeep matching Result columns from y to avoid bugs (false if by_auto is NULL)
			if (any(grepl("Result$|Source$", by_auto))) {
				y@output <- y@output[!colnames(y@output) %in% by_auto[grepl("Result$|Source$", by_auto)]]
				by_auto <- by_auto[!grepl("Result$|Source$", by_auto)]
			}
			
			if (length(by_auto) == 0) {
				if (ncol(x@output) == 1) {
					a <- data.frame(rep(x@output[[1]], each = nrow(y@output)))
					colnames(a) <- colnames(x@output)
				} else {
					a <- x@output[rep(1:nrow(x@output), each = nrow(y@output)), ]
				}
				temp <- data.frame(a, y@output[rep(1:nrow(y@output), times = nrow(x@output)), ])
				if (ncol(y@output) == 1) {
					colnames(temp)[length(colnames(temp))] <- colnames(y@output)
				}
			} else {
				type <- "inner"
				if (all == TRUE) type <- "full" else  {
					args <- list() #list(...)
					if (!is.null(args$all.x)) {
						if (args$all.x) type <- "left"
					}
					if (!is.null(args$all.y)) {
						if (args$all.y) {
							if (type == "left") type <- "full" else type <- "right"
						}
					}
				}
				test <- list()
				for (i in by_auto) {
					if (!is.factor(x@output[[i]])) x@output[[i]] <- factor(x@output[[i]])
					if (!is.factor(y@output[[i]])) y@output[[i]] <- factor(y@output[[i]])
				}
				temp <- join(x@output, y@output, by_auto, type, "all")
			}
			
			temp <- new("ovariable", output = temp)
			#temp <- CheckMarginals(temp, deps = list(x,y))
			return(temp)
		}
)

setMethod(f = "merge", 
		signature = signature(x = "ovariable", y = "numeric"),
		definition = function(x, y, ...) {
			y <- new("ovariable", output = data.frame(Result = y))
			return(callGeneric(x, y, ...))
		}
)

setMethod(f = "merge", 
		signature = signature(x = "numeric", y = "ovariable"),
		definition = function(x, y, ...) {
			x <- new("ovariable", output = data.frame(Result = x))
			return(callGeneric(x, y, ...))
		}
)

setMethod(f = "merge", 
		signature = signature(x = "ovariable", y = "data.frame"),
		definition = function(x, y, ...) {
			y <- new("ovariable", output = y)
			return(callGeneric(x, y, ...))
		}
)

setMethod(f = "merge", 
		signature = signature(x = "data.frame", y = "ovariable"),
		definition = function(x, y, ...) {
			y <- new("ovariable", output = x)
			return(callGeneric(x, y, ...))
		}
)

# SETMETHOD PLOT ################ plot diagrams about ovariable data

setMethod(
		f = "plot",
		signature = signature(x = "ovariable"),
		definition = function(x) {
			plot(
					x    = x@output[, paste("Source", x@name, sep = "")], 
					y    = x@output$Result, 
					xlab = paste("Source", x@name, sep = ""), 
					ylab = x@output[x@output[, paste("Source", x@name, sep = "")] == "Data", "Unit"][1], 
					main = x@name
			)
		}
)

# SETMETHOD summary ################### Summary defines how summaries of ovariables are shown.
setMethod(
		f = "summary", 
		signature = signature(object = "ovariable"), 
		definition = function(object, function_names = character(), marginals = character(), hide_source = TRUE, ...) {
			#test <- paste(object@name, "Result", sep = "") %in% colnames(object@output)
			#rescol <- ifelse(test, paste(object@name, "Result", sep = ""), "Result")
			#object@output <- object@output[ , -grep("Description|Source", colnames(object@output))] # not a necessary line
			
			# EvalOutput if not done yet
			
			if(nrow(object@output) == 0) object <- EvalOutput(object)
			
			# If no function names are defined then use defaults which depend on whether the data is probabilistic or not
			if("Iter" %in% colnames(object@output) && !"Iter" %in% marginals) {
				if (length(function_names)==0) function_names <- c("mean", "sd", "min", "Q0.025", "median", "Q0.975", "max")
				#object@output <- object@output[object@output$Iter == 1, ]
			}
			else {
				if (length(function_names)==0) function_names <- c("mean")
			}
			function_names <- unique(function_names)
			
			functions <- list()
			for(fname in function_names) {
				functions <- c(functions, get(fname))
			}
			
			# If marginals are not defined the data is assumed probabilistic and the summary to be about the distribution
			if(length(marginals)==0) {
				marginals <- colnames(object@output)[object@marginal & colnames(object@output) != "Iter"]
				
				# Hide single source source-columns if hide_source is TRUE.
				if (hide_source == TRUE) {
					source_cols <- marginals[grep("Source$", marginals)]
					for (i in source_cols) {
						locs <- unique(object@output[[i]])
						locs <- locs[!is.na(locs)]
						if (length(locs) == 1) {
							marginals <- marginals[marginals != i]
						}
					}
				}
			}
			
			# Remove NA results to reduce problems
			object@output <- object@output[!is.na(result(object)),]
			
			# Apply the selected functions
			temp <- list()
			for(fun in functions){
				temp[[length(temp)+1]] <- tapply(result(object), object@output[marginals], fun)
			}
			#out <- data.frame()
			
			# Convert results to data.frames and remove useless rows
			for(i in 1:length(temp)){
				temp[[i]] <- as.data.frame(as.table(temp[[i]]))
				temp[[i]] <- temp[[i]][!is.na(temp[[i]][["Freq"]]),]
				colnames(temp[[i]])[colnames(temp[[i]])=="Freq"] <- function_names[i]
			}
			
			# Merging
			if(length(temp)>1) {
				out <- merge(temp[[1]], temp[[2]], all = TRUE)
				if(length(temp)>2) {
					for(i in 3:length(temp)) {
						out <- merge(out, temp[[i]], all = TRUE)
					}
				}
			}
			else {
				out <- temp[[1]]
			}
			#if(nrow(object@output) > 200) {
			#	object@output <- object@output[1:200, ]
			#}
			#return(object@output)
			return(out)
		}
)

Q0.025 <- function(x){
	return(quantile(x, probs = 0.025))
}

Q0.975 <- function(x){
	return(quantile(x, probs = 0.975))
}



####################
# result
######################
### result returns a vector that contains the result column of the
### output of a given ovariable. The vector contains the original column
### name as the attribute comment.
### e1 is the ovariable to operate with.

result <- function(e1) { # e1 must be an ovariable or a data.frame.
	
# Should we allow people to use this for data.frames as well?
#	if(class(e1) == "data.frame") e1 <- new("ovariable", name = character(), output = e1)
	
	# First check presence of name specific Result-columns
	
	test1 <- "Result" %in% colnames(e1@output)
	
	test3 <- paste(e1@name, "Result", sep = "") %in% colnames(e1@output)
	
	# If found take note
	
	rescol1 <- ifelse(test1, "Result", paste(e1@name, "Result", sep = ""))
	
	if(!(test1 | test3)) stop("No result column found while operating mathematically with ovariables!\n")
	
	out <- e1@output[[rescol1]]
	comment(out) <- rescol1 # Save the column name for later use
	
	return(out)
}

## "result<-" is a function that tells what is done if content is assigned into Getrescol(ovariable).
## e1 is the ovariable into which something is assigned.
## value is the thing to assign into the ovariable.

assign("result<-", function(e1, value) {
			e1@output[[comment(result(e1))]] <- value
			return(e1)}
)

####################
# ddata_apply
############################
# This function will download newest available data in the base according to the defined ddata link (page identifier).
# Normally if data already exists it is left alone. Replacement can be forced with a parameter.
# Remember to use ddata_tidy = FALSE for old data with "obs" as Iteration column.
#########################
ddata_apply <- function(
		ovariable, 
		ddata_tidy = TRUE, 
		force_ddata = FALSE, 
		...
) { 
	if (length(attributes(ovariable)) < 8) return(ovariable) # line for compatibility with old ovariable definitions
	if ((identical(ovariable@data, data.frame()) | force_ddata) & !identical(ovariable@ddata, character())) {
		ovariable@data <- opbase.data(ovariable@ddata)
		if (ddata_tidy) ovariable@data <- tidy(ovariable@data, ovariable@name, direction='long') # data from base should
		# always be taken as is (ddata_tidy = FALSE) or to direction 'wide'
	}
	return(ovariable)
}


# continuousOps merges two ovariables by continuous indices and performs an operation.
# O1, O2 are ovariables. O1 is of main interest, while O2 has information that links to O2 via continuous index or indices.
# All locations in these continuous indices of O1 are created for O2 assuming that the value in the previous location of cols applies.
# Note that this is asymmetric. Locations in O2 that are missing from O1 are omitted.
# continuousOps assumes that all continuous indices are in the same dimension, the first one being the main index.
# Additional indices affect the outcome only if there are (approximate) ties. Therefore, avoid using this with several continuous indices.
# However, if continuous indices are NOT shared by both O1 and O2, they cause no trouble.
# fun is the name of a function that is performed after merge. Typically it is '*', '+' or some other Ops.
# cols is the vector of continuous indices that are merged. It is only needed if there are > 1 indices and the order is critical.
# Otherwise, shared continuous indices are identified automatically.

continuousOps <- function(O1, O2, fun, cols = NULL)
{
	rescol <- paste(O2@name, "Result", sep = "")
#	O1 <- unkeep(O1, paste(O1@name, "Source", sep = "")) # Remove these because they will create unpredictable results 
#	O2 <- unkeep(O2, paste(O2@name, "Source", sep = "")) # with orbind and merge. You should consider removing also other redundant columns.

	# From O1, we only want the shared continuous indices. The rest will return in the end with Ops.
	commons <- intersect(colnames(O1@output)[O1@marginal], colnames(O2@output)[O2@marginal])
	contcommon <- commons[sapply(O1@output[commons], FUN = is.numeric)]
	contcommon <- unique(c(cols, contcommon)) # The user-defined order is implemented if given.

	out <- 	merge( # Take the common continuous indices from the main ovariable and combine with additional data. Fill gaps.
		O1@output[contcommon], 
		O2@output[O2@marginal | colnames(O2@output) == rescol], 
		all = TRUE
	)

	marginals <- colnames(out)[colnames(out) %in% union(colnames(O2@output)[O2@marginal], colnames(O1@output)[O1@marginal])]
	discretes <- marginals[!sapply(out[marginals], FUN = is.numeric)] # Find continuous indices among marginals.

	if(!all(marginals %in% contcommon)) out <- fillna(out, setdiff(marginals, contcommon))
	
	if(length(discretes) > 0) out <- orbind(out, unique(out[discretes])) # Add a NA between each cell defined by the non-continuous indices.

	# Sort along the continuous indices, each group separately. Groups are unique combinations of discrete index locations.
	out <- out[do.call(order, (out[unique(c(discretes, contcommon, marginals))])) , ] 
	
	out[[rescol]] <- ifelse(
		is.na(out[[rescol]] & !is.na(out[[contcommon[1]]])), # If value is missing and this is not a group border.
		c(NA, out[1:(nrow(out) - 1) , rescol]), # value from the previous row.
		out[[rescol]]
	)
	
	out <- unique(out[!is.na(out[[rescol]]) , ]) # removes all rows that are before the first location in O2.
	out <- Ovariable(name = sub("Result$", "", rescol), output = out, marginal = colnames(out) %in% marginals)
	out <- do.call(fun, list(O1, out)) # Perform an Ops or other function with the original main ovariable and the data-enhanced ovariable.
	return(out)
}

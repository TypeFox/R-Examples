## Definition

setClass("mRMRe.Network", representation(topologies = "list", mi_matrix = "matrix", causality_list = "list",
                sample_names = "character", feature_names = "character", target_indices = "integer"))

## Wrappers
`mRMR.network` <- function(...)
{
    return(new("mRMRe.Network", ...)) 
}


## FIXME: Add wrappers for network

## initialize

setMethod("initialize", signature("mRMRe.Network"), function(.Object, data, prior_weight, target_indices, levels,
                layers, ..., mi_threshold = -Inf, causality_threshold = Inf)
{
    if (missing(layers))
        layers <- 1L
    
    
    #.Object@causality_vector <- as.numeric(sapply(seq(featureCount(data)), function(i){ return(NA) }))

    .Object@mi_matrix <- matrix(nrow = featureCount(data), ncol = featureCount(data), dimnames = list(featureNames(data), featureNames(data)))
    .Object@sample_names <- sampleNames(data)
    .Object@feature_names <- featureNames(data)
    .Object@target_indices <- as.integer(target_indices)
    .Object@topologies <- list()
	.Object@causality_list <- list()
    
    for (i in 1:layers)
    {
        filter <- new("mRMRe.Filter", data = data, prior_weight = prior_weight, target_indices = target_indices,
                levels = levels, ...)

        solutions <- solutions(filter, mi_threshold = mi_threshold, causality_threshold = causality_threshold)
       	causality <- causality(filter)
		
        lapply(names(solutions), function(i) { 
					.Object@topologies[[i]] <<- solutions[[i]]
					.Object@causality_list[[i]] <<- causality[[i]]
				})
		screen <- which(!is.na(mim(filter, method="cor")))
        .Object@mi_matrix[screen] <- mim(filter, method="cor")[screen]

        new_target_indices <- unique(unlist(solutions))
        new_target_indices <- new_target_indices[!is.na(new_target_indices)]
        target_indices <- new_target_indices[!as.character(new_target_indices) %in% names(.Object@topologies)]
        
        if (length(target_indices) == 0)
            break()
    }

    if (length(target_indices) == 0)
        return(.Object)

    ## Perform last-layer linking  

    filter <- new("mRMRe.Filter", data = data, prior_weight = prior_weight, target_indices = target_indices,
            levels = levels, ...)
    solutions <- solutions(filter, mi_threshold = mi_threshold, causality_threshold = causality_threshold)
    
    lapply(target_indices, function(target_index)
    {
        solution <- solutions[[as.character(target_index)]]
        new_solutions <- apply(solution, c(1, 2), function(feature_index)
                    ifelse(as.character(feature_index) %in% names(.Object@topologies), feature_index, NA))
        
        if (sum(is.na(new_solutions)) > 0)
            .Object@topologies[[as.character(target_index)]] <<- new_solutions
    })
    
    return(.Object)
})

## show

setMethod("show", signature("mRMRe.Network"), function(object)
{
    str(object)
})

## sampleNames

setMethod("sampleCount", signature("mRMRe.Network"), function(object)
		{
			return(length(object@sample_names))
		})

## sampleNames

setMethod("sampleNames", signature("mRMRe.Network"), function(object)
{
    return(object@sample_names)
})

## featureCount

setMethod("featureCount", signature("mRMRe.Network"), function(object)
		{
			return(length(object@feature_names))
		})

## featureNames

setMethod("featureNames", signature("mRMRe.Network"), function(object)
{
    return(object@feature_names)
})

## solutions

setMethod("solutions", signature("mRMRe.Network"), function(object)
{
    # filters[[target]][solution, ] is a vector of selected features
    # in a solution for a target; missing values denote removed features
            
    return(object@topologies)
})

## mim

setMethod("mim", signature("mRMRe.Network"), function(object)
{
    # mi_matrix[i, j] contains the biased correlation between
    # features i and j (i -> j directionality)
            
    return(object@mi_matrix)
})

## causality

setMethod("causality", signature("mRMRe.Network"), function(object)
{
    # causality_matrix[[target]][feature] contains the causality coefficient
    # between feature and target (feature -> target directionality)
            
    return(object@causality_list)
})

## adjacencyMatrix

setMethod("adjacencyMatrix", signature("mRMRe.Network"), function(object)
{
    adjacency_matrix <- matrix(0, nrow = length(object@feature_names), ncol = length(object@feature_names), dimnames=list(object@feature_names, object@feature_names))
    
    lapply(names(object@topologies), function(target_index)
    {
        connected_indices <- as.vector(object@topologies[[target_index]])
        connected_indices <- unique(connected_indices[!is.na(connected_indices)])
        if(length(connected_indices) > 0) {
            adjacency_matrix[connected_indices, as.integer(target_index)] <<- 1
            if(length(causality(object)) == 0)
                adjacency_matrix[as.integer(target_index), connected_indices] <<- 1
        }
    })

    # adjacency matrix: parents (seletected features) in rows, children (target features) in columns
    return(adjacency_matrix)
})

## adjacencyMatrixSum

setMethod("adjacencyMatrixSum", signature("mRMRe.Network"), function(object)
{
    adjacency_matrix <- matrix(0, nrow = length(object@feature_names), ncol = length(object@feature_names), dimnames=list(object@feature_names, object@feature_names))
    
    lapply(names(object@topologies), function(target_index)
    {
        connected_indices <- as.vector(object@topologies[[target_index]])
        connected_indices <- sort(connected_indices[!is.na(connected_indices)])
        connected_indices_count <- table(connected_indices)
        connected_indices <- unique(connected_indices)
        if(length(connected_indices) > 0) {
            adjacency_matrix[connected_indices, as.integer(target_index)] <<- connected_indices_count
            if(length(causality(object)) == 0)
                adjacency_matrix[as.integer(target_index), connected_indices] <<- connected_indices_count
        }
    })

    # adjacency matrix: parents (seletected features) in rows, children (target features) in columns
    return(adjacency_matrix)
})

## visualize

setMethod("visualize", signature("mRMRe.Network"), function(object)
{
    ## FIXME : Cannot find a way to display vertex names...
    
    adjacency <- adjacencyMatrix(object)
	used_rows <- apply(adjacency, 1, sum) > 0
	used_cols <- apply(adjacency, 2, sum) > 0
	adjacency <- adjacency[used_rows | used_cols, used_cols | used_rows]
    graph <- graph.adjacency(adjacency, mode = "undirected", add.rownames = 'label')
	
    return(plot.igraph(graph))
})

## target

setMethod("target", signature("mRMRe.Network"), function(object)
{
	return(object@target_indices)
})

setMethod("scores", signature("mRMRe.Network"), function(object)
{
	mi_matrix <- mim(object)
	targets <- names(solutions(object))

	scores <- lapply(targets, function(target) {
				apply(solutions(object)[[target]], 2, function(solution) {
							sapply(1:length(solution), function(i) {

										feature_i <- solution[i]
										if(is.na(feature_i))
											return(NA)
										if(i == 1)
											return(mi_matrix[as.numeric(target), feature_i])
								
										ancestry_score <- mean(sapply((i-1):1, function(j) mi_matrix[feature_i, solution[j]]))
										return(mi_matrix[as.numeric(target), feature_i] - ancestry_score)
									})							
						})
			})
	names(scores) <- targets
	return(scores)
})

make.matrix <- function(vertices, interactions) {

	# remove extra comma at the beginning
	vertices <- sub("^,","",vertices);
	interactions <- sub("^,","",interactions);

	vertices.list <- unlist(strsplit(vertices,","));
	interactions.list <- unlist(strsplit(interactions,","));

	# let retrieve the gene ids i-e vertices names
    captions <- rep("NA", length(vertices.list));
	for (i in seq(1, length(vertices.list), 1)) {
		captions[i] <- gsub("\"", "", vertices.list[i]);
		}

	adjacency.matrix <- matrix(
		data = 0,
		nrow = length(vertices.list),
		ncol = length(vertices.list),
		dimnames = list(
			captions,
			captions
			)
		);

	# lets populate the adjacency matrix
	for (i in 1:length(interactions.list)) {
		interactors <- unlist(strsplit(gsub("\"","", interactions.list[i]), ":"));
		adjacency.matrix[interactors[1], interactors[2]] <- 1;
		adjacency.matrix[interactors[2], interactors[1]] <- 1;
		}

	return(adjacency.matrix);

	}

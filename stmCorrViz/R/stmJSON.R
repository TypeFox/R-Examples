stmJSON <-
function(mod, documents_raw=NULL, documents_matrix=NULL,
                    title="STM Model", clustering_threshold=1.5,
                    labels_number=7, verbose){

#   The stmJSON function is a top-layer function that outputs a nested JSON
#   structure for use in D3 visualizations of an STM model. The overall procedure
#   is as follows (for greater detail, refer to the comments in the individual
#   subroutines):
#     1. The routine retrieves the theta matrix from the STM object, then computes
#        the correlation of topics based on it, hence uses correlation to compute
#        distances among topics, and finally performs a hierarchical clustering of
#        the topics (this is achieved by calling hclust).
#     2. The algorithm then finds all binary splits in the middle of the hierarchical
#        clustering tree whose clustering height measure is below a user-specified
#        clustering threshold. All these splits are marked as "aggregation points".
#     3. The routine retrieves the merge matrix from the output of hclust, and
#        produces a new merge list (with the same conventions) by getting rid of
#        all the splits performed at aggregation points along the tree. While the
#        hclust merge matrix only contains binary splits, the new merge list can
#        (and presumably will) contain non-binary cluster splits. The subroutine
#        collapseMergedList() is used to complete this step in the procedure.
#     4. The merge list is turned into a structure of nested lists with a recursive
#        function call. Each level of this nested structure corresponds to a
#        node in the hierarchical representation of the STM model. The terminal
#        nodes are topics. Topics are children to clusters. Each terminal node is
#        represented as a list of two elements: name and size. Each non-terminal
#        node is a list of two elements: name, and a list of children.
#     5. Lastly, the nested-list stucture is turned into a JSON object by using the
#        rjson library.
#     6. The assignClusterNames() subroutine is used for assigning labels to
#        clusters. This is currently non-functional.

  # Require libraries
  #require(stm)
  #require(slam)
  #require(jsonlite)

  # Generate baseline topic list
  out <- list()

  if(verbose==TRUE)
    cat("Performing hierarchical topic clustering ... \n")

  # Run hclust subroutine
  clust <- clusterAnalysis(mod, labels_number=labels_number)

  # Obtain document-term matrix
  # documents_raw <- iconv(documents_raw, to='utf-8', sub="")
  # temp <- textProcessor(documents_raw, metadata, verbose=FALSE)
  # meta<-temp$meta
  # vocab<-temp$vocab
  # docs<-temp$documents_raw
  # prep_output <- prepDocuments(docs, vocab, meta, verbose=FALSE)
  # docs<-prep_output$documents_raw

  if(verbose==TRUE)
    cat("Generating JSON representation of the model ... \n")

  # Extra data objects
  thoughts <- stm::findThoughts(mod, documents_raw)
  full_labels <- stm::labelTopics(mod)
  topic_proportions <- colMeans(mod$theta)
  #exclusivity_scores <- exclusivity(mod)
  #semcoh_scores <- semanticCoherence(mod, docs, 10)

  # Find aggregation points
  K <- mod$settings$dim$K
  topic_to_topic_splits <- c()
  for(i in seq(K-1))
    if(clust$merge[i,1] <0 && clust$merge[i,2] < 0)
      topic_to_topic_splits <- c(topic_to_topic_splits, i)
  aggregate <- setdiff(which(clust$height <= clustering_threshold),
                         topic_to_topic_splits)

  # Produce merge list
  merge_list <- list()
  for(i in seq(K-1))
    merge_list[[i]] <- clust$merge[i,]
  names(merge_list) <- 1:(K-1)

  # Collapse merge list
  merge_list <- collapseMergeList(merge_list, clust$merge, aggregate, K)

  # Build collapsed-clusters data structure
  top_layer <- merge_list[paste(K-1)][[1]]
  out$children <- buildClusters(list(), current = top_layer, merge_list,
                                labels=clust$labels, full_labels,
                                thoughts, topic_proportions)

  # Implementation with diagnostics
  #   out$children <- buildClusters(list(), current = top_layer, merge_list,
  #                                 labels=clust$labels, full_labels,
  #                                 thoughts, exclusivity_scores,
  #                                 semcoh_scores)

  # Get beta weights for model
  beta_weights <- getBetaWeights(mod, documents_matrix)

  # Assign cluster names
  out <- assignClusterNames(out, labels_number, beta_weights, mod$vocab)

  # Root Information
  out$name <- title
  out$this_root <- TRUE
  out$summary <- utils::capture.output(mod)
  out$proportions <- topic_proportions

  # Convert structure to JSON
  out_JSON <- jsonlite::toJSON(out, force=TRUE)
  return(out_JSON)
}

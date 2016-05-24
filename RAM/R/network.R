network_data <- function (data, is.OTU=TRUE, metadata) {

  if (!requireNamespace("reshape2")) {
    stop("Package reshape2 is required for this function")
  }


  if ( is.OTU ) {
    if ( ! ("taxonomy" %in% names(data)) ) {
      stop("Are you sure this is an OTU table?")
    }
    d <- LCA.OTU(data)
    rownames(d) <- paste(rownames(d), d$LCA, sep="_")
    d <- d[, -ncol(d)]
    d <- as.data.frame(t(d))
  } else {
    if ( "taxonomy" %in% names(data) ) {
      stop("Are you sure this is NOT an OTU table?")
    }
    d <- data
  }

  m <- metadata[match(rownames(d), rownames(metadata)),]
  if ( !identical(rownames(d), rownames(m))) {
    stop("Samples in otu table and metadatadata do not match")
  }

  ## nodes  
  col <- ncol(m) + 4
  df1 = data.frame(matrix(vector(), nrow(d), col, dimnames=list(c(), c("node_name", "ntype", "degree", "weighted_degree", names(m)))), stringsAsFactors=F)
  df1$node_name <- rownames(d)
  df1$ntype <- "user_node"
  # node degree (number of edges)
  df1$degree <- rowSums(decostand(d, "pa"))
  # weighted degree (total counts of edges)
  df1$weighted_degree <- rowSums(d)
  for ( i in names(m) ) {
    df1[[i]] <- m[[i]]
    class(df1[[i]]) <- class(m[[i]])
  }
  
  d.t <- as.data.frame(t(d))
  df2 = data.frame(matrix(vector(), nrow(d.t), col, dimnames=list(c(), c("node_name", "ntype", "degree", "weighted_degree", names(m)))), stringsAsFactors=F)
  df2$node_name <- rownames(d.t)
  df2$ntype <- "tax_node"
  df2$degree <- rowSums(decostand(d.t, "pa"))
  df2$weighted_degree <- rowSums(d.t)
  for ( i in names(m) ) {
    df2[[i]] <- ""
    if ( class(m[[i]]) == "Date" ) {
      class(df2[[i]]) <- class(m[[i]])
    }
  }
  
 
  if ( identical(names(df1), names(df2)) ) {
    net_node <- rbind(df1, df2)
  }

  # edges
  d$id <- rownames(d)
  d.melt <- reshape2::melt(d)
  names(d.melt)[1] <- "from"
  names(d.melt)[2] <- "to"
  names(d.melt)[3] <- "eweight"
  d.melt <- d.melt[d.melt$eweight != 0,]
  net_edge <- merge(d.melt, m, by.x="from", by.y="row.names", all=TRUE)
   
  return(list(net_node=net_node, net_edge=net_edge))

}

  

  


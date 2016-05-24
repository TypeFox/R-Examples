
# 0.9.6 add; assigns vertices randomly to contacts, assuming that all vertices 
# are equally likely participate in a contact
vertex_randomization <- function(edges)
{
  vertex_columns <- c("VertexFrom", "VertexTo")
  vertices <- unique(unlist(edges[, vertex_columns]))
  edges$VertexFrom <- edges$VertexTo    
  repeat
  {
    invalid <- edges$VertexFrom == edges$VertexTo
    invalid_count <- sum(invalid)
    if (invalid_count == 0) break
    edges[invalid, vertex_columns] <- sample(vertices, invalid_count * 2, replace=T)
  }
  return(edges)
}

# 0.9.6 add; assigns vertices randomly to contacts, but limits how often a 
# vertex participates in a contact to how often it originally participated
contact_randomization <- function(edges)
{
  vertex_columns <- c("VertexFrom", "VertexTo")
  vertices <- unlist(edges[, vertex_columns])
  edges[, vertex_columns] <- sample(vertices, length(vertices))
  edge_count <- nrow(edges)
  repeat
  {
    invalid <- edges$VertexFrom == edges$VertexTo
    invalid_count <- sum(invalid)
    if (invalid_count == 0) break
    for (i in which(invalid)) edges <- swap(edges, i, sample(vertex_columns, 1), sample(edge_count, 1), sample(vertex_columns, 1)) 
  }
  return(edges)
}

# 0.9.6 add, see Holme & Saramaki, Physics Reports 519 (2012), p. 118
time_reversal <- function(edges)
{
  max_time <- max(edges$TimeStop)
  min_time <- min(edges$TimeStart)
  tmp <- edges$TimeStart
  edges$TimeStart <- max_time - edges$TimeStop + min_time
  edges$TimeStop <- max_time - tmp + min_time
  return(edges)
}

# 0.9.6 add, see Holme & Saramaki, Physics Reports 519 (2012), p. 117
randomly_permuted_times <- function(edges)
{
  edge_count <- nrow(edges)
  new_order <- sample(edge_count, edge_count)
  edges$TimeStart <- edges$TimeStart[new_order]
  edges$TimeStop <- edges$TimeStop[new_order]
  return(edges)
}

# 0.9.6 add, see Holme & Saramaki, Physics Reports 519 (2012), p. 117
random_times <- function(edges)
{
  delta <- edges$TimeStop - edges$TimeStart
  max_time <- max(edges$TimeStop)
  range <- min(edges$TimeStart):max_time
  edges$TimeStop <- max_time + 1
  repeat
  {
    invalid <- edges$TimeStop > max_time
    invalid_count <- sum(invalid)
    if (invalid_count == 0) break
    edges$TimeStart[invalid] <- sample(range, invalid_count, replace=T)
    edges$TimeStop[invalid] <- edges$TimeStart[invalid] + delta[invalid]    
  }
  return(edges)
}

# 0.9.6 add, see Holme & Saramaki, Physics Reports 519 (2012), p. 117
randomized_contacts <- function(edges)
{
  unique_edges <- unique(edges[, c("VertexFrom", "VertexTo")])
  unique_edge_count <- nrow(unique_edges)
  map <- sample(unique_edge_count, nrow(edges), replace=T)
  map[sample(length(map), unique_edge_count)] <- 1:unique_edge_count
  edges$VertexFrom <- unique_edges$VertexFrom[map]
  edges$VertexTo <- unique_edges$VertexTo[map] 
  return(edges)
}

# 0.9.6 add, see Holme & Saramaki, Physics Reports 519 (2012), p. 117
edge_randomization <- function(edges)
{
  return(randomize_edges_helper(edges, F))
}

# 0.9.6 add, see Holme & Saramaki, Physics Reports 519 (2012), p. 116
randomized_edges <- function(edges)
{  
  return(randomize_edges_helper(edges, T))
}

# 0.9.6 add
randomize_edges_helper <- function(edges, randomize_vertices)
{
  
  # list all edges 
  vertex_columns <- c("VertexFrom", "VertexTo")
  unique_edges <- unique(edges[, vertex_columns])
  unique_edge_count <- nrow(unique_edges)
  
  # map each edge to a randomly chosen edge 
  edge_map <- cbind(unique_edges, unique_edges[sample(unique_edge_count, unique_edge_count), ])  
  new_vertex_columns <- c("NewVF", "NewVT")
  colnames(edge_map) <- c(vertex_columns, new_vertex_columns)
  
  # randomize which vertices that are connected by the randomly chosen edges
  if (randomize_vertices)
  {
    edge_map[, new_vertex_columns] <- sample(unlist(edge_map[, new_vertex_columns]), unique_edge_count * 2)    
    repeat
    {
      invalid <- (edge_map$NewVF == edge_map$NewVT) | (duplicated(edge_map[, new_vertex_columns]))
      if (sum(invalid) == 0) break
      for (i in which(invalid)) edge_map <- swap(edge_map, i, sample(new_vertex_columns, 1), sample(unique_edge_count, 1), sample(new_vertex_columns, 1)) 
    }
  }
  
  # replace vertices of original edges with the new ones 
  original_colnames <- colnames(edges)
  attribute_columns <- original_colnames[!(original_colnames %in% vertex_columns)]
  edges <- merge(edges, edge_map)
  edges <- edges[, c(new_vertex_columns, attribute_columns)]
  colnames(edges)[1:length(new_vertex_columns)] <- vertex_columns
  
  # done
  return(edges)
  
}

# 0.9.6 add
swap <- function(df, r1, c1, r2, c2)
{
  tmp <- df[r1, c1]
  df[r1, c1] <- df[r2, c2]
  df[r2, c2] <- tmp
  return(df)
}
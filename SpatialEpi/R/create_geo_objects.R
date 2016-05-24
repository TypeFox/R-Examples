create_geo_objects <-
function(max.prop, population, centroids, sp.obj){

# Number of areas
n <- nrow(centroids)


#-----------------------------------------------------------------------------
# Set up single zones
#-----------------------------------------------------------------------------
if(max.prop < max(normalize(population))){
  print(paste("max.prop needs to be at least", max(normalize(population))))
}

# Get basic single zone info
geoInfo <- zones(centroids, population, max.prop)
nearest.neighbors <- geoInfo$nearest.neighbors
cluster.coords <- geoInfo$cluster.coords
n.zones <- nrow(cluster.coords)

# Create list of length n.zones indicating the component areas for each zone
cluster.list <- vector(mode="list", length=n.zones)
counter <- 1
for(i in 1:n) {
  nn <- nearest.neighbors[[i]]
  for(j in 1:length(nn)) {
    cluster.list[[counter]] <- nn[1:j] 
    counter <- counter + 1  
  } 
}

#-----------------------------------------------------------------------------
# Generate overlap object which tracks the overlap between single zones
#-----------------------------------------------------------------------------
# For each area, get "presense":  list all unit zones that it is included in
presence <- vector(mode="list",length=n)
for(i in 1:n){
  presence[[i]] <- which(sapply(cluster.list, function(x){is.element(i, x)}))
}


#-------------------------------------------------------------------------------
# Compute "buffer" of areas between single zones.
#-------------------------------------------------------------------------------
# Compute adjacency matrix
adj <- poly2nb(sp.obj, queen=TRUE)
adj <- nb2mat(adj, zero.policy=TRUE)[1:n, 1:n]

# Convert Adjacency Matrix To A List
adj.new <- vector("list",length=n)
for(i in 1:n){
  adj.new[[i]] <- which(adj[i,]!=0)
}
adj <- adj.new
rm(adj.new)

# Update presence list to incorporate buffer.  We need to preserve original 
# presence list all the way thru, so temporarily store results here
presence_temp <- presence

# Loop thru all areas
for(i in 1:n){
  current.zones <- presence[[i]]
  adjacent.areas <- adj[[i]]
  adjacent.zones <- unique(unlist(presence[adjacent.areas]))
  presence_temp[[i]] <- sort(unique(c(adjacent.zones, current.zones)))
}

# restore presence list
presence <- presence_temp


#----------------------------------------------------------------
# Return list
#----------------------------------------------------------------
overlap <- list(presence=presence, cluster.list=cluster.list)
return(list(
  overlap=overlap, 
  cluster.coords=cluster.coords)
)
}

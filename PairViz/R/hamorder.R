

order_tsp <- function(d, method = "nearest", cycle=TRUE,improve=FALSE,path_dir = path_cor,...)
{   #returns SHORTEST CYCLE or PATH via TSP
	#method must be one of"nearest_insertion", "farthest_insertion", "cheapest_insertion",         
	#"arbitrary_insertion" "nn", "repetitive_nn"
	# "2-opt" and if concorde package is loaded, "concorde"
	# Sometimes "2-opt" gets stuck and does not return!!!
#	require(TSP)
	tsp <- TSP(d)
	if (cycle != TRUE) tsp <- insert_dummy(tsp, label= "cut")
	tour <- solve_TSP(tsp,method=method,...)
	if (improve && method != "2-opt")
	  tour <- solve_TSP(tsp,method="2-opt",control=list(tour=tour))
	if (cycle != TRUE) tour <- cut_tour(tour,"cut")
    tour <- as.integer(tour) 
    if (is.function(path_dir))
       tour <-  best_orientation(tour,d,cycle, path_dir)
    return(tour)
	}
	
order_best <- function(d, maxexact=9,nsamples=50000,path_weight=sum,cycle=FALSE,path_dir = path_cor,...) {
     #returns SHORTEST PATH or approximate 
  #   require(gtools)
	if (class(d)=="dist") {
	   d <- as.matrix(d)
	   dnames <- labels(d)}
	else 
         dnames <- colnames(d)
      n <- nrow(d)
      if (n <= maxexact)
     	 perms <<- permutations(n,n)
      else perms <- t(sapply(1:nsamples, function(i) sample(n,n)))
      pathlens <- apply(perms,1, function(h) path_weight(path_values(h,d,cycle)))
      o <- perms[which.min(pathlens),]
      if (is.function(path_dir))
         o <-  best_orientation(o,d,cycle, path_dir) 
      return(o)
	}

#--------------------------------------------

	
# path_values <- function(path,d,cycle=FALSE)  {
	# n <- length(path)
	# o <- cbind(path[-n],path[-1])
	# if (cycle) path <- rbind(path,c(path[n],path[1]))
	# return(d[o])	
	# }

path_values <- function(path,d,cycle=FALSE)  {
	n <- length(path)
	o <- cbind(path[-n],path[-1])
	if (cycle) o <- rbind(o,c(path[n],path[1]))
	return(d[o])	
	}


best_orientation <- function(path,d, cycle=FALSE, path_dir= path_cor,from=NULL){
	# Finds the best cycle/path preserving adjacencies in path.
	# For cycles, If from is NULL, the  best start is found first

  vecshift <- function(vec,s) {
 	v <- vec[s:length(vec)]
 	if (s > 1) v <- c(v,vec[1:(s-1)])
 	return(v)
  	}
  d <- as.matrix(d)
  n <- length(path)
  v <- path_values(path,d, cycle)
  if (cycle) {
     if (is.null(from)) {
    	dirs <- sapply(1:length(v),function(s) { 
    		          vs <- vecshift(v,s) 
    		          return(c(path_dir(vs),path_dir(rev(vs))))})
    	if (max(dirs[1,]) > max(dirs[2,]))
           path <- vecshift(path,which.max(dirs[1,]))
         else
           path <- rev(vecshift(path,which.max(dirs[2,])))
         } else {
        path <- vecshift(path,match(from,path))
         v <- path_values(path,d, cycle)
        if (path_dir(v) < path_dir(rev(v))) path <-  c(path[1], rev(path[-1]))   
       }} else 
        if (path_dir(v) < path_dir(rev(v))) path <- rev(path)       
   return(path)
}
	
	
path_cor <- function(edgew,method="kendall")
	cor(1:length(edgew),edgew,method=method)
	
	
	
weighted_hpaths <- function(d, path1 = NULL,paths=NULL,
           matrix=TRUE,cycle=NULL,path_weight=sum,path_dir = path_cor,...)	{
           
	# The first path is given by path- if not provided, path/cycle with smallest path_weight.
	# Using path_dir find best start and orientation for path1, and use to relabel other rows of paths.
	# Using path_dir re-orientation rows 2..k of paths
	# Permute successive paths in order of path length, as given by pathfn.
	# and permute rows of paths using total path_weight.	 #If matrix is TRUE, returns a matrix where each row is a hamiltonian 
     # path, otherwise concatenates the rows into a vector.
 
	d <- as.matrix(d)
	n <- nrow(d)
	if (is.null(cycle)) cycle <- n %% 2 != 0    
	if (is.null(path1)) path1 <- order_tsp(d, cycle=cycle,...)
	if (is.null(paths)) paths <- hpaths(n,cycle=cycle)
	path1 <- best_orientation(path1,d,cycle, path_dir) 
	paths <- permute_hpaths(path1,paths,matrix=TRUE)
	for (i in 2: nrow(paths)) 
	   paths[i,] <- best_orientation(paths[i,],d,cycle,path_dir,path1[1])
	hlen <- apply(paths,1, function(h) path_weight(path_values(h,d,cycle)))
	roword <- c(1,order(hlen[-1])+1)
	paths <- paths[roword,]
	if (!matrix)  {
		paths <- as.vector(t(paths))
 		if (n%%2 != 0) paths <- close_path(paths)
 		}
 	return(paths)
	}
	

# set class value for cluster scatter list

clv.paste <- function(...) paste(..., sep="")
cls.class <- function() { return("cls.list") }

# important functions to check validity of inputs in many functions in this package 

data.validity <- function(data, name)
{
	if( !is.matrix(data) ) 
		if( !is.data.frame(data) )
			stop(paste("Bad usage: input '", name, "' should be matrix or data.frame type.", sep=""))
	
	data = as.matrix(data)
	if( !is.double(data) )
		stop(paste("Bad usage: input '", name, "' is not double type.", sep=""))
	
	return( data )
}

data.validity.int <- function(data, name)
{
	if( !is.matrix(data) ) 
		if( !is.data.frame(data) )
			stop(paste("Bad usage: input '", name, "' should be matrix or data.frame type.", sep=""))
	
	if( !is.integer(data) )
		stop(paste("Bad usage: input '", name, "' is not integer type.", sep=""))
	
	return( as.matrix(data) )
}


# this function checks integer vectors that contains numbers of clusters 
# to which data will be clustered

cls.num.vect.validity <- function(clust, obj.num, name)
{
	# informtaion about object id's should be a vector type 
	if( is.vector(clust) == FALSE ) 
		stop(paste("Bad usage: input '", name ,"' should be vector type.", sep=""))

	# cluster id should be a integer type, if not, it have to be coerced 
	if( !is.integer(clust) )
	{
		clust = as.integer(clust)
		warning(paste("Vector '", name,"' should be an integer so it was coerced to integer type.", sep=""))
		if( NA %in% clust )
			stop(paste("Coercion of vector '", name, "' produced NA values ", sep=""))
		if( TRUE %in% ( clust > 100 ) )
			warning("Cluster number is very big (more than 100).")
	}
	
	# cluster id can't be less than 2
	if( length( clust[ clust < 2 ] ) != 0 || length( clust[ clust > obj.num ] ) != 0 )
		stop(paste("Bad input data: vector '", name, "' contains numbers less than 2 
			  but every cluster number should be a number from 2 to object number.", sep=""))
	
	return(clust)
}

# this function checks integer vectors that contains cluster id's for each clustered object 

cls.id.vect.validity <- function(clust, name)
{
	# informtaion about object id's should be a vector type 
	if( is.vector(clust) == FALSE ) 
		stop(paste("Bad usage: input '", name ,"' should be vector type.", sep=""))

	# cluster id should be a integer type, if not, it have to be coerced 
	if( !is.integer(clust) )
	{
		clust = as.integer(clust)
		warning(paste("Vector '", name,"' should be an integer so it was coerced to integer type.", sep=""))
		if( NA %in% clust )
			stop(paste("Coercion of vector '", name, "' produced NA values ", sep=""))
		if( TRUE %in% ( clust > 100 ) )
			warning("Cluster number is very big (more than 100).")
	}
	# cluster id can't be less than 1
	if( length( clust[ clust < 1 ] ) != 0 )
		stop(paste("Bad input data: vector '", name, "' contains numbers less than 1 
			  but cluster id should be a number from 1 to max(clust).", sep=""))

	# check empty clusters
	if( FALSE %in% (1:max(clust) %in% clust) )
		warning(paste("Vector '", name, "' contains empty clusters, one or more values 
				from 1 to max(clust) do not appear in this vector.", sep=""))

	return(clust)
}

dist.validity <- function(dist)
{
	dist <- pmatch(dist, c("euclidean","manhattan","correlation"))
	if(is.na(dist) || dist == -1)
		stop(paste("Bad usage: input '", dist, "' is not supported, 
			 available are: euclidean | manhattan | correlation. ", sep=""))

	return(as.integer(dist))
}

# functions which check validity of cluster measure methods choosen by the user

check.intracls.method <- function(intracls)
{
	methods = c("complete", "average", "centroid")
	choosen.methods = methods %in% intracls
	if( !(TRUE %in% choosen.methods) )
		stop("Bad input data: 'intracls' vector does not contain any supported intracluster distance, 
			  supported are: complete | average | centroid." )
	return(choosen.methods)
	
}

check.intercls.method <- function(intercls)
{
	methods = c("single", "complete", "average", "centroid", "aveToCent", "hausdorff")
	choosen.methods = methods %in% intercls
	if( !(TRUE %in% choosen.methods) )
		stop("Bad input data: 'intercls' vector does not contain any supported intracluster distance, 
			  supported are: single | complete | average | centroid | aveToCent | hausdorff" )
	return(choosen.methods)	
}

check.intracls.diss.mx.method <- function(intracls)
{
	methods = c("complete", "average")
	choosen.methods = methods %in% intracls
	if( !(TRUE %in% choosen.methods) )
		stop("Bad input data: 'intracls' vector does not contain any supported intracluster distance, 
			  supported are: complete | average." )
	
	return(choosen.methods)
}

check.intercls.diss.mx.method <- function(intercls)
{
	methods = c("single", "complete", "average", "hausdorff")
	choosen.methods = methods %in% intercls
	if( !(TRUE %in% choosen.methods) )
		stop("Bad input data: 'intercls' vector does not contain any supported intracluster distance, 
			  supported are: single | complete | average | hausdorff" )
	return(choosen.methods)
	
}

check.avail.methods <- function(user.methods, vec.name, supp.methods)
{	
	choosen.methods = supp.methods %in% user.methods
	if( !(TRUE %in% choosen.methods) )
		stop(paste( "Bad input data: \"", vec.name, "\" vector does not contain any supported method name, supported are: ", paste( supp.methods , collapse=" | " ), sep="" ) )
	return(choosen.methods)
	
}

# functions very usefull in visualization of "clusterScatterMeasures results
cut.matrix <- function(matrix, not.empty.cls, names)
{
	if( length(names) != length(not.empty.cls) )
		matrix = as.matrix(matrix[not.empty.cls,not.empty.cls])
	colnames(matrix) = names
	rownames(matrix) = names
	return(matrix)
}

cut.vector <- function(vect, not.empty.cls, names)
{
	vect = vect[not.empty.cls]
	vect = matrix(vect,1,length(vect))
	colnames(vect) = names
	return(vect)
}

cluster.size <- function(clust, cl.num=0)
{
	clust_num = cl.num
	if( is.numeric(clust_num) == FALSE || cl.num <= 0) clust_num = max(clust)

	idx <- .Call( "clv_clustersSizeExt",
			  as.integer(clust),
			  as.integer(clust_num),
			  PACKAGE="clv"
			)

	return(idx)
}

cls.attrib <- function(data, clust)
{
	data = data.validity(data, "data")
	clust = cls.id.vect.validity(clust, "clust")
	
	if(dim(data)[1] != length(clust))
		stop("Bad input data: number of 'data' objects do not agree with length of vector 'clust'.")

	clust_num = as.integer(max(clust))

    idx <- .Call( "clusterAttrib",
				  data,
				  clust,
				  clust_num,
				  PACKAGE="clv"
				)

	# check which clusters are empty
	not.empty.cls = (idx$cluster.size > 0)
	names = paste( "c", 1:clust_num, sep="")[not.empty.cls]

	if( length(names) > 1 )
	{
		idx$cluster.center = as.matrix(idx$cluster.center[not.empty.cls,])
		rownames(idx$cluster.center) = names
	}
	else
	{
		idx$cluster.center = t(as.matrix(idx$cluster.center[not.empty.cls,]))
		rownames(idx$cluster.center) = names
	}

	return(idx)
}

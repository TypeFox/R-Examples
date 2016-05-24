
# starnard internal cluster scatter measures: 
# - distances between clusters
# - density of each cluster


cls.scatt.data <- function( data, clust, dist="euclidean" )
{
	data = data.validity(data, "data")
	clust = cls.id.vect.validity(clust, "clust")
	dist = dist.validity(dist)
	
	if(dim(data)[1] != length(clust))
		stop("Bad input data: number of 'data' objects do not agree with length of vector 'clust'.")

	clust_num = as.integer(max(clust))

    idx <- .Call( "clusterScatterMeasures",
				  data,
				  clust,
				  clust_num,
				  dist,
				  PACKAGE="clv"
				)
	
	# check which clusters are empty
	not.empty.cls = (idx$cluster.size > 0)
	names = paste( "c", 1:clust_num, sep="")[not.empty.cls]
	idx$intracls.complete = cut.vector(idx$intracls.complete,not.empty.cls,names)
	idx$intracls.average = cut.vector(idx$intracls.average,not.empty.cls,names)
	idx$intracls.centroid = cut.vector(idx$intracls.centroid,not.empty.cls,names)

	if( length(names) > 1 )
	{
		# cut all rows and columns which represents empty clusters
		idx$intercls.single = cut.matrix(idx$intercls.single, not.empty.cls, names)
		idx$intercls.complete = cut.matrix(idx$intercls.complete, not.empty.cls, names)
		idx$intercls.average = cut.matrix(idx$intercls.average, not.empty.cls, names)
		idx$intercls.centroid = cut.matrix(idx$intercls.centroid, not.empty.cls, names)
		idx$intercls.ave_to_cent = cut.matrix(idx$intercls.ave_to_cent, not.empty.cls, names)
		idx$intercls.hausdorff = cut.matrix(idx$intercls.hausdorff, not.empty.cls, names)
	
		idx$cluster.center = as.matrix(idx$cluster.center[not.empty.cls,])
		rownames(idx$cluster.center) = names
	}
	else
	{
		m = NA
		idx$intercls.single = m
		idx$intercls.complete = m
		idx$intercls.average = m
		idx$intercls.centroid = m
		idx$intercls.ave_to_cent = m
		idx$intercls.hausdorff = m
		idx$cluster.center = t(as.matrix(idx$cluster.center[not.empty.cls,]))
		rownames(idx$cluster.center) = names
	}
	class(idx) <- cls.class()

    return(idx)
}

# wersion with dissimilarity matrix 

cls.scatt.diss.mx <- function( diss.mx, clust )
{
	diss.mx = data.validity(diss.mx, "diss.mx")
	clust = cls.id.vect.validity(clust, "clust")

	if(dim(diss.mx)[1] != dim(diss.mx)[2])
		stop("Bad input data: 'diss.mx' should be a square (symetric) matrix.")
	if(dim(diss.mx)[1] != length(clust))
		stop("Bad input data: number of 'diss.mx' objects do not agree with length of vector 'clust'.")

	clust_num = as.integer(max(clust))

    idx <- .Call( "clusterScatterMeasuresDissMx",
				  diss.mx,
				  clust,
				  clust_num,
				  PACKAGE="clv"
				)
	
	# check which clusters are empty
	not.empty.cls = (idx$cluster.size > 0)
	names = paste( "c", 1:clust_num, sep="")[not.empty.cls]
	idx$intracls.complete = cut.vector(idx$intracls.complete,not.empty.cls,names)
	idx$intracls.average = cut.vector(idx$intracls.average,not.empty.cls,names)

	if( length(names) > 1 )
	{
		# cut all rows and columns which represents empty clusters
		idx$intercls.single = cut.matrix(idx$intercls.single, not.empty.cls, names)
		idx$intercls.complete = cut.matrix(idx$intercls.complete, not.empty.cls, names)
		idx$intercls.average = cut.matrix(idx$intercls.average, not.empty.cls, names)
		idx$intercls.hausdorff = cut.matrix(idx$intercls.hausdorff, not.empty.cls, names)
	}
	else
	{
		m = NA
		idx$intercls.single = m
		idx$intercls.complete = m
		idx$intercls.average = m
		idx$intercls.hausdorff = m
	}
	class(idx) <- cls.class()

    return(idx)
}

# Dunn index:
# 
# 		min( dist(Ck, Cl) )/max(diam(Cm))
#
# k,l,m numbers of clusters which come from the same partitioning 

clv.Dunn <- function( index.list, intracls, intercls )
{

	# sprawdzenie czy ind.list jest obiektem pochodz±cym z funckji cls.scatt.measures
	if( class(index.list) != cls.class() ) 
		stop("Bad input data: 'index.list' is not an object created by function 'cls.scatt.measures(..)'
				or 'cls.scatt.measures.diss.mx(..)' .")

	idx = index.list

	if( length(index.list) == 7 )
	{
		intra.bool = check.intracls.diss.mx.method(intracls)
		inter.bool = check.intercls.diss.mx.method(intercls)
		intra.name = c("comp", "ave")
		inter.name = c("sin", "comp", "ave", "haus")
		intra.list = list( idx$intracls.complete, idx$intracls.average )[intra.bool]
		inter.list = list( idx$intercls.single, idx$intercls.complete, idx$intercls.average,
					   idx$intercls.hausdorff )[inter.bool]
	}
	else if( length(index.list) == 11 )
		{
			intra.bool = check.intracls.method(intracls)
			inter.bool = check.intercls.method(intercls)
			intra.name = c("comp", "ave", "cent")
			inter.name = c("sin", "comp", "ave", "cent", "aveto", "haus")
			intra.list = list( idx$intracls.complete, idx$intracls.average, idx$intracls.centroid )[intra.bool]
			inter.list = list( idx$intercls.single, idx$intercls.complete, idx$intercls.average,
								idx$intercls.centroid, idx$intercls.ave_to_cent, idx$intercls.hausdorff )[inter.bool]
		}
		else stop("Bad input data: 'index.list' is not an object created by function 'cls.scatt.measures(..)'
					or 'cls.scatt.measures.diss.mx(..)' .")

	intra.num = length(intra.bool[intra.bool])
	inter.num = length(inter.bool[inter.bool])
	
	result = matrix(0, inter.num, intra.num)
	
	not.empty.cls = (index.list$cluster.size > 0)

	for( i in 1:intra.num )
		intra.list[[i]] = max(intra.list[[i]])

	# find minimum in each chosen matrix but only in columns and rows which represents not empty cluster
	# and only in lower tringle part of matrix
	
	boolmx = lower.tri( inter.list[[1]][ not.empty.cls, not.empty.cls ] )
	for( i in 1:inter.num )
	{
		inter.list[[i]] = inter.list[[i]][ not.empty.cls, not.empty.cls ]
		inter.list[[i]] = min( inter.list[[i]][ boolmx ] )	
	}

	for( i in 1:inter.num ) 
		for( j in 1:intra.num ) 
			result[i,j] = inter.list[[i]]/intra.list[[j]]

	rownames(result) = inter.name[inter.bool]
	colnames(result) = intra.name[intra.bool]

	return(result)
}

# Davies-Bouldin index:
# 	sum(forall Ck) [ max(Ck != Cl)(diam(Ck)+diam(Cl))/dist(Ck,Cl)) ]

clv.DB.ind <- function( vect, mx, clust.num )
{
	result = 0
	for( i in 1:clust.num )
	{
		max = 0
		for( j in 1:clust.num)
		{
			if( i != j)
			{
				if( i > j ) tmp = (vect[i]+vect[j])/mx[i,j]
				else tmp = (vect[i]+vect[j])/mx[j,i]
				if( max < tmp ) max = tmp
			}
		}
		result = result + max
	}
	return(result/clust.num)
}


# Davies-Bouldin index:
# 
# 		sum(forall Ck) [ max(Ck != Cl)(diam(Ck)+diam(Cl))/dist(Ck,Cl)) ]
#
# k,l numbers of clusters which come from the same partitioning
# as output we have matrix ( dim = c(lenght(intercls),lenght(intracls)) ) of DB indicies
# matrix dimension depends how many "diam" and "dist" measures will be chosen by the user
# ( dim(matrix) = c(lenght(intercls),lenght(intracls)) ) 


clv.Davies.Bouldin <- function( index.list, intracls, intercls )
{
	# sprawdzenie czy ind.list jest obiektem pochodz±cym z funckji cls.scatt.measures
	if( class(index.list) != cls.class() ) 
		stop("Bad input data: 'index.list' is not an object created by function 'cls.scatt.measures(..)'
				or 'cls.scatt.measures.diss.mx(..)' .")

	idx = index.list

	if( length(index.list) == 7 )
	{
		intra.bool = check.intracls.diss.mx.method(intracls)
		inter.bool = check.intercls.diss.mx.method(intercls)
		intra.name = c("comp", "ave")
		inter.name = c("sin", "comp", "ave", "haus")
		intra.list = list( idx$intracls.complete, idx$intracls.average )[intra.bool]
		inter.list = list( idx$intercls.single, idx$intercls.complete, idx$intercls.average,
					   idx$intercls.hausdorff )[inter.bool]
	}
	else if( length(index.list) == 11 )
		{
			intra.bool = check.intracls.method(intracls)
			inter.bool = check.intercls.method(intercls)
			intra.name = c("comp", "ave", "cent")
			inter.name = c("sin", "comp", "ave", "cent", "aveto", "haus")
			intra.list = list( idx$intracls.complete, idx$intracls.average, idx$intracls.centroid)[intra.bool]
			inter.list = list( idx$intercls.single, idx$intercls.complete, idx$intercls.average,
								idx$intercls.centroid, idx$intercls.ave_to_cent, idx$intercls.hausdorff )[inter.bool]
		}
		else stop("Bad input data: 'index.list' is not an object created by function 'cls.scatt.measures(..)'
					or 'cls.scatt.measures.diss.mx(..)' .")
	
	intra.num = length(intra.bool[intra.bool])
	inter.num = length(inter.bool[inter.bool])
	
	not.empty.cls = (idx$cluster.size > 0)
	clust_num = length(not.empty.cls[not.empty.cls])

	# cut all those rows and columns which represents empty cluetrs 
	for( i in 1:intra.num ) 
		intra.list[[i]] = intra.list[[i]][ not.empty.cls ]

	for( i in 1:inter.num )
		inter.list[[i]] = inter.list[[i]][ not.empty.cls, not.empty.cls ]

	# create result matrix
	result = matrix(0, inter.num, intra.num)
	# nad compute all indicies
	for( i in 1:inter.num )
		for( j in 1:intra.num )
			result[i,j] = clv.DB.ind(intra.list[[j]],inter.list[[i]], clust_num)

	rownames(result) = inter.name[inter.bool]
	colnames(result) = intra.name[intra.bool]

	return(result)
}

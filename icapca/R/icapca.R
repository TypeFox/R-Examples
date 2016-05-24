ica_pca <- function(x, inf_crit='unc', components=0, center=TRUE, subgaussian_range=NULL, supergaussian_range=NULL, gaussian_range=NULL, hinted_subgaussian_sources=NULL, hinted_supergaussian_sources=NULL, hinted_unspecified_sources=NULL, seed=0, offset_random=0, fold=0, xval_epsilon=0.0, desired_initialization=0, sample_order=NULL, sample_offset=0, samples=0)
{
	if(!is.matrix(x)) stop("x must be a matrix")
	mrows<-dim(x)[1]
	ncols<-dim(x)[2]

	if(!is.character(inf_crit)) stop("inf_crit must be one of 'unc', 'aic', 'bic', 'xval', 'caicj'")
	if(inf_crit=='unc') inf_crit_numeric=0
	else if(inf_crit=='aic') inf_crit_numeric=1
	else if(inf_crit=='bic') inf_crit_numeric=2
	else if(inf_crit=='xval') inf_crit_numeric=3
	else if(inf_crit=='caicj') inf_crit_numeric=4
	else stop(paste(inf_crit, "is not a valid value for inf_crit"))
	
	if(components==0) components<-mrows
	else{
		components=as.integer(components)
		if(is.na(components) || !is.finite(components) || components<=0){
			stop("unable to  coerce components to a finite positive integer")
		}
	}

	subgaussian_range_min=0
	subgaussian_range_max=components
	if(!is.null(subgaussian_range)){
		subgaussian_range=as.integer(subgaussian_range)
		if(length(subgaussian_range)<1 || length(subgaussian_range)>2) stop("subgaussian_range must be one or two non-negative integers")
		subgaussian_range_min=as.integer(subgaussian_range[1]);
		if(is.na(subgaussian_range_min) || !is.finite(subgaussian_range_min) || subgaussian_range_min<0) stop("subgaussian_range must be one or two non-negative integers")
		if(length(subgaussian_range)==1) subgaussian_range_max=subgaussian_range_min
		else{
			subgaussian_range_max=as.integer(subgaussian_range[2]);
			if(is.na(subgaussian_range_max) || !is.finite(subgaussian_range_max) || subgaussian_range_max<0) stop("subgaussian_range must be one or two non-negative integers")
			if(subgaussian_range_max<subgaussian_range_min){
				warning("exchanging requested minimum and maximum for subgaussian range")
				temp=subgaussian_range_min
				subgaussian_range_min=subgaussian_range_max
				subgaussian_range_max=temp
			}
		}
	}
	supergaussian_range_min=0
	supergaussian_range_max=components
	if(!is.null(supergaussian_range)){
		supergaussian_range=as.integer(supergaussian_range)
		if(length(supergaussian_range)<1 || length(supergaussian_range)>2) stop("supergaussian_range must be one or two non-negative integers")
		supergaussian_range_min=as.integer(supergaussian_range[1]);
		if(is.na(supergaussian_range_min) || !is.finite(supergaussian_range_min) || supergaussian_range_min<0) stop("supergaussian_range must be one or two non-negative integers")
		if(length(supergaussian_range)==1) supergaussian_range_max=supergaussian_range_min
		else{
			supergaussian_range_max=as.integer(supergaussian_range[2]);
			if(is.na(supergaussian_range_max) || !is.finite(supergaussian_range_max) || supergaussian_range_max<0) stop("supergaussian_range must be one or two non-negative integers")
			if(supergaussian_range_max<supergaussian_range_min){
				warning("exchanging requested minimum and maximum for supergaussian range")
				temp=supergaussian_range_min
				supergaussian_range_min=supergaussian_range_max
				supergaussian_range_max=temp
			}
		}
	}
	gaussian_range_min=0
	gaussian_range_max=components
	if(!is.null(gaussian_range)){
		gaussian_range=as.integer(gaussian_range)
		if(length(gaussian_range)<1 || length(gaussian_range)>2) stop("gaussian_range must be one or two non-negative integers")
		gaussian_range_min=as.integer(gaussian_range[1]);
		if(is.na(gaussian_range_min) || !is.finite(gaussian_range_min) || gaussian_range_min<0) stop("gaussian_range must be one or two non-negative integers")
		if(length(gaussian_range)==1) gaussian_range_max=gaussian_range_min
		else{
			gaussian_range_max=as.integer(gaussian_range[2]);
			if(is.na(gaussian_range_max) || !is.finite(gaussian_range_max) || gaussian_range_max<0) stop("gaussian_range must be one or two non-negative integers")
			if(gaussian_range_max<gaussian_range_min){
				warning("exchanging requested minimum and maximum for gaussian range")
				temp=gaussian_range_min
				gaussian_range_min=gaussian_range_max
				gaussian_range_max=temp
			}
		}
	}
	if(subgaussian_range_max + supergaussian_range_max + gaussian_range_max==0)
	 stop("ranges for subgaussian, supergaussian and gaussian sources cannot all be zero")
	ranges<-c(subgaussian_range_min, subgaussian_range_max, supergaussian_range_min , supergaussian_range_max, gaussian_range_min, gaussian_range_max)

	source_rows=0
	hinted_subgaussian_source_count=0
	if(!is.null(hinted_subgaussian_sources)){
		if(!is.matrix(hinted_subgaussian_sources)) stop("hinted_subgaussian_sources must be a matrix")
		if(dim(hinted_subgaussian_sources)[2]!=ncols) stop("number of columns in hinted_subgaussian sources must match number of columns in x")
		hinted_subgaussian_source_count=hinted_subgaussian_source_count+dim(hinted_subgaussian_sources)[1]
		source_rows=source_rows+hinted_subgaussian_source_count
	}
	hinted_supergaussian_source_count=0
	if(!is.null(hinted_supergaussian_sources)){
		if(!is.matrix(hinted_supergaussian_sources)) stop("hinted_supergaussian_sources must be a matrix")
		if(dim(hinted_supergaussian_sources)[2]!=ncols) stop("number of columns in hinted_supergaussian sources must match number of columns in x")
		hinted_supergaussian_source_count=hinted_supergaussian_source_count+dim(hinted_supergaussian_sources)[1]
		source_rows=source_rows+hinted_supergaussian_source_count
	}
	hinted_unspecified_source_count=0
	if(!is.null(hinted_unspecified_sources)){
		if(!is.matrix(hinted_unspecified_sources)) stop("hinted_unspecified_sources must be a matrix")
		if(dim(hinted_unspecified_sources)[2]!=ncols) stop("number of columns in hinted_unspecified sources must match number of columns in x")
		hinted_unspecified_source_count=hinted_unspecified_source_count+dim(hinted_unspecified_sources)[1]
		source_rows=source_rows+hinted_unspecified_source_count
	}
	if(source_rows>mrows) stop("number of combined rows in hinted sources cannot exceed number of rows in x")
	
	source_matrix<-0
	if(source_rows>0){
		source_matrix=matrix(nrow=source_rows, ncol=ncols)
		start=1
		if(hinted_subgaussian_source_count>0){
			source_matrix[start:(start+hinted_subgaussian_source_count-1),]=hinted_subgaussian_sources
			start=start+hinted_subgaussian_source_count
		}
		if(hinted_supergaussian_source_count>0){
			source_matrix[start:(start+hinted_supergaussian_source_count-1),]=hinted_supergaussian_sources
			start=start+hinted_supergaussian_source_count
		}
		if(hinted_unspecified_source_count>0){
			source_matrix[start:(start+hinted_unspecified_source_count-1),]=hinted_unspecified_sources
			start=start+hinted_unspecified_source_count
		}
	}
	
	ica_s<-vector(mode="numeric", length=components*ncols)
	loglikelihood<-0
	distribution<-vector(mode="numeric", length=components)
	variance<-vector(mode="numeric", length=components)
	probability<-vector(mode="numeric", length=components+1)
	leave_rows_uncentered<-as.integer(!center)
	error=0
	
	seed=as.integer(seed)
	offset_random=as.integer(offset_random)
	xval_epsilon=as.double(xval_epsilon)
	desired_initialization=as.integer(desired_initialization)
	

	if(!is.null(sample_order)){
		if(!is.vector(sample_order) || !is.integer(sample_order)) stop("sample_order must be a vector of integers")
		if(length(sample_order)!=ncols) stop("sample_order must be a vector with length matching the number of columns in the input data")
		sample_order<-as.integer(sample_order);
	}
	samples=as.integer(samples)
	if(samples>0){
		sample_offset=as.integer(sample_offset)
		if(is.null(sample_order)) stop("sample_order cannot be null when samples>0")
	}	
	fold=as.integer(fold)
	if(fold>0){
		if(inf_crit_numeric!=3) stop("fold can only be non-zero when inf_crit=='xval'")
		if(fold==1) stop("fold cannot be equal to one");
		if(ncols%%fold!=0) stop("fold must be a divisor of the number columns in the input matrix");
	}
	
	if(seed>0){
		set.seed(seed, kind="default", normal.kind="default")
	}
	
	result<-.C("ica_pca_R", as.integer(error), as.integer(seed), as.integer(offset_random),
	 as.integer(sample_order), length(sample_order), as.integer(sample_offset), as.integer(samples),
	 as.integer(inf_crit_numeric), as.integer(fold), as.double(xval_epsilon), as.integer(ranges), 
	 as.double(x), as.integer(mrows), as.integer(ncols), 
	 as.integer(components), as.double(ica_s), as.double(loglikelihood), as.integer(desired_initialization),
	 as.integer(distribution), as.double(variance), 
	 as.double(probability),
	 as.integer(hinted_subgaussian_source_count),
	 as.integer(hinted_supergaussian_source_count),
	 as.integer(hinted_unspecified_source_count),
	 as.double(t(source_matrix)), as.integer(leave_rows_uncentered),
	 PACKAGE="icapca")
	 	
	if(result[[1]]>0) return(NULL)
	ica_s<-matrix(result[[16]], nrow=result[[15]], ncol=ncols)
	loglikelihood<-result[[17]]
	distribution<-as.vector(result[[19]][1:result[[15]]])
	variance<-as.vector(result[[20]][1:result[[15]]])
	probability<-as.vector(result[[21]][1:(result[[15]]+1)])
	
	subgaussian_range<-as.vector(result[[11]][1:2])
	supergaussian_range<-as.vector(result[[11]][3:4])
	gaussian_range<-as.vector(result[[11]][5:6])
	result<-list(s=ica_s, loglikelihood=loglikelihood, distribution=distribution, variance=variance, probability=probability, subgaussian_range=subgaussian_range, supergaussian_range=supergaussian_range, gaussian_range=gaussian_range)
	return(result)
}
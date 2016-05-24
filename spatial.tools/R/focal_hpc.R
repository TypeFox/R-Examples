

# focal_hpc helper functions

focal_hpc_precheck <- function(x,window_dims,window_center,processing_unit,verbose)
{
	if(verbose) message("Performing pre-checks...")
	
	if(length(window_dims)==1) window_dims=c(window_dims,window_dims)	
	if(length(window_center)==1) window_center <- c(window_center,window_center)
	if(is.na(window_center[2])) window_center[2] <- ceiling(window_dims[2]/2)
	
	if(verbose) { message(paste("window_dims:",window_dims,sep="")) }
	if(verbose) { message(paste("window_center:",window_center,sep="")) }
	
	window_rows=window_dims[2]
	window_cols=window_dims[1]
	
	if(is.list(x))
	{
		layer_names <- sapply(x,names,simplify=FALSE)
	} else
	{
		layer_names=names(x)
	}
	
	if(any(window_dims>1))
	{
		if(verbose) message("Focal processing mode...")
		processing_mode="focal"
		if(is.null(processing_unit)) processing_unit <- "single"
	} else
	{
		if(verbose) message("Pixel processing mode...")
		processing_mode="pixel"
		if(is.null(processing_unit)) processing_unit <- "chunk"
	}
	
	startrow_offset=-(window_center[2]-1)
	endrow_offset=window_rows-window_center[2]
	
	return(list(window_dims=window_dims,window_center=window_center,
					window_rows=window_rows,window_cols=window_cols,
					layer_names=layer_names,
					processing_mode=processing_mode,processing_unit=processing_unit,
					startrow_offset=startrow_offset,endrow_offset=endrow_offset))
}

focal_hpc_test <- function(x,fun,window_center,window_dims,args,
		layer_names,
		startrow_offset,endrow_offset,
		processing_mode,processing_unit,chunk_format,
		verbose)
{

	window_index <- NULL # Why is this here?
	
	if(verbose) { message("Checking the function on a small chunk of data.") }
	
	# Add additional info to the args.
	if(!is.null(args)) {
		args$window_center=window_center
		args$window_dims=window_dims
		args$layer_names=layer_names
	} else
	{
		args=list(window_center=window_center)
		args$window_dims=window_dims
		args$layer_names=layer_names
	}
	
	# We are going to pull out the first row and first two pixels to check the function...
	
	
	if(processing_unit=="single")
	{
		if(verbose) { message("processing_unit=window...")}
		if(class(x)=="list")
		{
			r_check <- sapply(X=x,
					FUN=function(X,window_dims,chunk_format)
					{
						getValuesBlock_enhanced(X, r1=1, r2=window_dims[2], c1=1,c2=window_dims[1],
								format=chunk_format)
						
					},window_dims=window_dims,chunk_format=chunk_format,
					simplify=FALSE)
		} else
		{
			r_check <- getValuesBlock_enhanced(x, r1=1, r2=window_dims[2], c1=1,c2=window_dims[1],
					format=chunk_format)
		}
	} else
	{
		# The function works on the entire chunk.
		if(verbose) { message("processing_unit=chunk...")}
		if(class(x)=="list")
		{
			r_check <- sapply(X=x,
					FUN=function(X,window_dims,chunk_format)
					{
#						getValuesBlock_enhanced(X, r1=1, r2=1, c1=1,c2=2,
#								format=chunk_format)
						
						test_chunk <- getValuesBlock_enhanced(X, 
								r1=1, r2=window_dims[2], c1=1,c2=(window_dims[1]+1),
								format=chunk_format)
						
						if(processing_mode=="focal")
						{
							x_off <- (1:window_dims[1])-ceiling(window_dims[1]/2)
							y_off <- (1:window_dims[2]) 
							#						z_off <- 1:dim(test_chunk[[1]])[3]
							z_off <- 1:dim(test_chunk)[3]
							xyz_off <- as.data.frame(t(expand.grid(x_off,y_off,z_off)))
							#						central_index <- (1:(dim(test_chunk[[1]])[1]
							central_index <- (1:(dim(test_chunk)[1]
											-(window_dims[1]-1)))+(ceiling(window_dims[1]/2))-1
							
#							focal_shifted_array <- mapply(
#									FUN=function(chunk,layer_names,window_index,window_dims,xyz_off,central_index)
#									{
#										focal_shifted_array <- sapply(X=xyz_off,
#												FUN=function(X,chunk,central_index)
#												{						
#													central_index_off <- central_index+X[1]
#													outchunk <- chunk[central_index_off,X[2],X[3],drop=FALSE]
#													outchunk <- aperm(outchunk,c(2,1,3))
#													return(outchunk)
#												},chunk=chunk,central_index=central_index)
#										dim(focal_shifted_array) <- c(length(central_index),prod(window_dims),dim(chunk)[3])
#										dimnames(focal_shifted_array) <- vector(mode="list",length=3)
#										if(!is.null(layer_names)) dimnames(focal_shifted_array)[[3]]=layer_names # [as.numeric(xyz_off[3,])]
#										return(focal_shifted_array)
#									},chunk=test_chunk,layer_name=layer_names,
#									MoreArgs=list(window_index=window_index,window_dims=window_dims,
#											xyz_off=xyz_off,central_index=central_index),
#									SIMPLIFY=FALSE)
							
							focal_shifted_array <- sapply(X=xyz_off,FUN=
											function(X,chunk,central_index)
									{
										central_index_off <- central_index+X[1]
										outchunk <- chunk[central_index_off,X[2],X[3],drop=FALSE]
										outchunk <- aperm(outchunk,c(2,1,3))
										return(outchunk)
									},chunk=test_chunk,central_index=central_index)	
							dim(focal_shifted_array) <- c(length(central_index),
									prod(window_dims),dim(test_chunk)[3])
							dimnames(focal_shifted_array) <- vector(mode="list",length=3)
							if(!is.null(layer_names)) { dimnames(focal_shifted_array)[[3]] <- names(X) }
							#layer_names }	
							return(focal_shifted_array)
						} else
						{
							return(test_chunk)
						}
					},window_dims=window_dims,chunk_format=chunk_format,simplify=FALSE)
		} else
		{
			test_chunk <- getValuesBlock_enhanced(x, 
					r1=1, r2=window_dims[2], c1=1,c2=(window_dims[1]+1),
					format=chunk_format)
			
			if(processing_mode=="focal")
			{
				x_off <- (1:window_dims[1])-ceiling(window_dims[1]/2)
				y_off <- (1:window_dims[2]) 
				z_off <- 1:dim(test_chunk)[3]
				xyz_off <- as.data.frame(t(expand.grid(x_off,y_off,z_off)))
				central_index <- (1:(dim(test_chunk)[1]
								-(window_dims[1]-1)))+(ceiling(window_dims[1]/2))-1
				
				focal_shifted_array <- sapply(X=xyz_off,FUN=
								function(X,chunk,central_index)
						{
							central_index_off <- central_index+X[1]
							outchunk <- chunk[central_index_off,X[2],X[3],drop=FALSE]
							outchunk <- aperm(outchunk,c(2,1,3))
							return(outchunk)
						},chunk=test_chunk,central_index=central_index)	
				dim(focal_shifted_array) <- c(length(central_index),
						prod(window_dims),dim(test_chunk)[3])
				dimnames(focal_shifted_array) <- vector(mode="list",length=3)
				if(!is.null(layer_names)) { dimnames(focal_shifted_array)[[3]] <- layer_names }
				
				r_check <- focal_shifted_array
				
			} else
			{
				r_check <- test_chunk
			}
		}
	}
	
	# Add additional info to the args.
	r_check_args=args
	r_check_args$x=r_check
	r_check_function <- do.call(fun, r_check_args)
	
	# Coerce to list
	if(class(r_check_function) != "list") r_check_function <- list(r_check_function)
	
	r_check_output_classes <- sapply(r_check_function,class)
	
	if(processing_unit=="single")
	{
		if(!all(r_check_output_classes=="numeric"))
		{
			stop("window processing units require numeric vector outputs.  Please check your function.")
		} else 
		{
			outbands=sapply(r_check_function,length)
			outfiles=length(r_check_function)
		}
	}
	
	if(processing_unit=="chunk")
	{
		r_check_output_dims <- sapply(r_check_function,dim)
		r_check_output_dims_col_check <- all(r_check_output_dims[1,]==2)
		r_check_output_dims_row_check <- all(r_check_output_dims[2,]==1)
		r_check_output_class_check <- all(r_check_output_classes=="array")
		
#		if(class(r_check_function)!="array" || 
#				dim(r_check_function)[1] != 2 ||
#				dim(r_check_function)[2] != 1)
		
		if(!r_check_output_class_check || !r_check_output_dims_col_check || !r_check_output_dims_row_check)
		{
			message("chunk processing units require array vector outputs.  Please check your function.")
			stop(dim(r_check_function))
		} else 
		{
			outbands <- r_check_output_dims[3,]
			outfiles=length(r_check_function)
		}
	}
#	}
	if(verbose) { message(paste("Number of output bands determined to be:",outbands,sep=" ")) }
	if(verbose) { message(paste("Number of output files determined to be:",outfiles,sep=" ")) }
	
	return(list(outbands=as.list(outbands),outfiles=outfiles))
}

focal_hpc_chunk_setup <- function(x,window_dims,window_center,
		chunk_nrows,startrow_offset,endrow_offset,minblocks,blocksize,verbose)
{
	nodes <- getDoParWorkers() 
	if(minblocks=="max" ) minblocks <- nodes
	
#	tr=blockSize(x,chunksize=(chunk_nrows*nodes+(window_dims[2]-1))*ncol(x))
	
	if(is.list(x))
	{
		nlayers_x <- sum(sapply(x,nlayers))
		x <- x[[1]]
	} else
	{
		nlayers_x <- nlayers(x)
	}
	
	blocksize_reduction=4
	
	if(is.null(blocksize))
		tr=blockSize(x,n=nlayers_x*blocksize_reduction,minrows=window_dims[2],minblocks=minblocks)
	else
		tr=blockSize(x,chunksize=ncol(x)*blocksize*nlayers_x,
				n=nlayers_x,minrows=window_dims[2],minblocks=minblocks)
	
	if (tr$n < nodes) nodes <- tr$n
	
	if(verbose) message(paste("Total number of blocks to process:",tr$n))
	
	tr$row2 <- tr$row + tr$nrows - 1
	
	tr$focal_row=tr$row+startrow_offset
	tr$focal_row2=tr$row2+endrow_offset
	
	tr$focal_row[tr$focal_row<1]=1
	tr$focal_row2[tr$focal_row2>nrow(x)]=nrow(x)
	
	tr$startrow_offset = startrow_offset
	tr$endrow_offset = endrow_offset
	
#	bottom_right_buffer <- window_dims-window_center
#	top_left_buffer <- window_center-c(1,1)
#	buffers <- c(top_left_buffer,bottom_right_buffer)
#	names(buffers) <- c("left","top","right","bottom")
#	tr$buffers <- buffers
	
	texture_tr=list(rowcenters=((tr$row[1]:tr$row2[1])-startrow_offset))
	texture_tr$row=texture_tr$rowcenters+startrow_offset
	texture_tr$row2=texture_tr$rowcenters+endrow_offset
	
	return(list(tr=tr,texture_tr=texture_tr))
}

focal_hpc_focal_getChunk <- function(x,tr,format,r,i,r_old,chunkArgs)
{
	# Create some blank variables:
	window_center <- NULL
	window_dims <- NULL
	datatype <- NULL
	
	list2env(chunkArgs,envir=environment())
	
	startrow_offset <- tr$startrow_offset
	endrow_offset <- tr$endrow_offset
	
	chunk_format=format
	
	if(i==1)
	{
		r <- getValuesBlock_enhanced(x, r1=tr$focal_row[i], r2=tr$focal_row2[i], c1=1, c2=ncol(x),
				format=chunk_format)
	} else
	{
		r <- getValuesBlock_enhanced(x, r1=(tr$focal_row2[(i-1)]+1), r2=tr$focal_row2[i], 
				c1=1, c2=ncol(x),format=chunk_format)
	}
	
	if(i==1)
	{
		# Add top cap
		if((1-(tr$row[1]+startrow_offset))>0)
			r=abind(
					array(data=NA,dim=c(ncol(x),(1-(tr$row[1]+startrow_offset)),nlayers(x))),
					r,
					along=2)
	}
	
	if(i==tr$n)
	{
		# Add bottom cap
		if(nrow(x)-tr$row2[tr$n]+endrow_offset>0)
			r=abind(r,
					array(data=NA,dim=c(ncol(x),(nrow(x)-tr$row2[tr$n]+endrow_offset),nlayers(x))),
					along=2)
	}
	
	# TODO: WE NEED TO BE ABLE TO SUBTRACT STUFF HERE ALSO (if center is outside)
	left_cap=window_center[2]-1
	right_cap=window_dims[2]-window_center[2]
	
	if(left_cap>0)
	{
		# Add left cap.
		r=abind(
				array(data=NA,dim=c(left_cap,dim(r)[2],dim(r)[3])),
				r,
				along=1)
	}
	
	if(right_cap>0)
	{
		# Add right cap.
		r=abind(
				r,
				array(data=NA,dim=c(right_cap,dim(r)[2],dim(r)[3])),
				along=1)
	}
	
	if(i>1 && window_dims[2]>1)
	{
		r <- abind(r_old,r,along=2)
	}
	return(r)	
}

focal_hpc_focalChunkFunction <- function(chunk,chunkArgs)
{	
	# TODO: as.numeric() the outputs of the functions.
	
	
	# Create some blank variables (to avoid R CMD CHECK errors):
	x <- NULL
	layer_names <- NULL
	fun <- NULL
	window_dims <- NULL
	# window_center <- NULL
	outbands <- NULL
	window_center <- NULL
	processing_unit <- NULL
	datatype <- NULL
	verbose <- NULL
	
	#
	e <- list2env(chunkArgs,envir=environment())
	
	# Add additional info to the args.
	if(!is.null(args)) {
		args$window_center=window_center
		args$window_dims=window_dims
		args$layer_names=layer_names
	} else
	{
		args=list(window_center=window_center)
		args$window_dims=window_dims
		args$layer_names=layer_names
	}
	
	if(is.list(x))
	{
		ncol_x=ncol(x[[1]])
		image_dims=dim(x[[1]])
	} else
	{
		ncol_x=ncol(x)
		image_dims=dim(x)
	}
	window_index=1:ncol_x
	
	# Needs to be sped up
#	system.time(
	
	if(processing_unit=="single")	
	{
		r_out <-
				mapply(
						function(window_index,chunk,args,window_dims)
						{
							if(is.list(chunk))
							{
								x_array <- mapply(FUN=
												function(chunk,layer_names,window_index,window_dims)
										{
											x_array=chunk[(window_index:(window_index+window_dims[2]-1)),,,drop=FALSE]									
											dimnames(x_array) <- vector(mode="list",length=3)
											if(!is.null(layer_names)) dimnames(x_array)[[3]]=layer_names
											return(x_array)
										},chunk=chunk,layer_name=layer_names,
										MoreArgs=list(window_index=window_index,window_dims=window_dims),
										SIMPLIFY=FALSE)	
								
							} else
							{
								x_array=chunk[(window_index:(window_index+window_dims[2]-1)),,,drop=FALSE]									
								dimnames(x_array) <- vector(mode="list",length=3)
								if(!is.null(layer_names)) dimnames(x_array)[[3]]=layer_names
							}
							fun_args=args
							fun_args$x=x_array
							r_out <- do.call(fun, fun_args)
							return(unlist(r_out))
						}
						,
						window_index,
						MoreArgs=list(chunk=chunk$processing_chunk,args=args,window_dims=window_dims),SIMPLIFY=TRUE
				)
		
		outbands_numeric <- unlist(outbands)
		outbands_end <- cumsum(outbands_numeric)
		outbands_start <- c(1,(outbands_end+1)[-length(filename)])
		
		if(class(r_out)=="numeric") dim(r_out) <- c(1,length(r_out))
		
		r_out <- mapply(function(outbands_start,outbands_end,r_out,ncol_x) 
				{
					subarray <- array(t(r_out[(outbands_start:outbands_end),,drop=FALSE]),dim=c(ncol_x,1,(outbands_end-outbands_start+1)))
					return(subarray)
				}
				,
				outbands_start=as.list(outbands_start),
				outbands_end=as.list(outbands_end),
				MoreArgs=list(r_out=r_out,ncol_x=ncol_x),SIMPLIFY=FALSE)
		
	} else
	{
		#	processing_chunk=chunk$processing_chunk
		if(verbose) message("Chunk mode...")
		if(is.list(chunk$processing_chunk))
		{
			x_off <- (1:window_dims[1])-ceiling(window_dims[1]/2)
			y_off <- (1:window_dims[2]) 
			z_off <- 1:dim(chunk$processing_chunk[[1]])[3]
			xyz_off <- as.data.frame(t(expand.grid(x_off,y_off,z_off)))
			
			central_index <- (1:(dim(chunk$processing_chunk[[1]])[1]-(window_dims[1]-1)))+(ceiling(window_dims[1]/2))-1
			
			focal_shifted_array <- mapply(
					FUN=function(chunk,layer_names,window_index,window_dims,xyz_off,central_index)
					{
						focal_shifted_array <- sapply(X=xyz_off,
								FUN=function(X,chunk,central_index)
								{						
									central_index_off <- central_index+X[1]
									outchunk <- chunk[central_index_off,X[2],X[3],drop=FALSE]
									outchunk <- aperm(outchunk,c(2,1,3))
									return(outchunk)
								},chunk=chunk,central_index=central_index)
						dim(focal_shifted_array) <- c(length(central_index),prod(window_dims),dim(chunk)[3])
						dimnames(focal_shifted_array) <- vector(mode="list",length=3)
						if(!is.null(layer_names)) dimnames(focal_shifted_array)[[3]]=layer_names # [as.numeric(xyz_off[3,])]
						return(focal_shifted_array)
					},chunk=chunk$processing_chunk,layer_name=layer_names,
					MoreArgs=list(window_index=window_index,window_dims=window_dims,
							xyz_off=xyz_off,central_index=central_index),
					SIMPLIFY=FALSE)
		} else
		{
			x_off <- (1:window_dims[1])-ceiling(window_dims[1]/2)
			y_off <- (1:window_dims[2]) 
			z_off <- 1:dim(chunk$processing_chunk)[3]
			xyz_off <- as.data.frame(t(expand.grid(x_off,y_off,z_off)))
			
			central_index <- (1:(dim(chunk$processing_chunk)[1]-(window_dims[1]-1)))+(ceiling(window_dims[1]/2))-1
			
			focal_shifted_array <- sapply(X=xyz_off,FUN=
							function(X,chunk,central_index)
					{
						central_index_off <- central_index+X[1]
						outchunk <- chunk[central_index_off,X[2],X[3],drop=FALSE]
						outchunk <- aperm(outchunk,c(2,1,3))
						return(outchunk)
					},chunk=chunk$processing_chunk,central_index=central_index)	
			dim(focal_shifted_array) <- c(length(central_index),prod(window_dims),dim(chunk$processing_chunk)[3])
			dimnames(focal_shifted_array) <- vector(mode="list",length=3)
			if(!is.null(layer_names)) { dimnames(focal_shifted_array)[[3]] <- layer_names }
		}
		fun_args=args
		fun_args$x=focal_shifted_array
		r_out <- do.call(fun, fun_args)
		if(!is.list(r_out)) r_out <- list(r_out)
		
#		r_out <- mapply(
#				function(r_out,outbands,ncol_x)
#				{
#					dim(r_out) <- c(ncol_x,1,outbands)
#				},
#				r_out=r_out,outbands=outbands,MoreArgs=list(ncol_x=ncol_x))
		
		#	dim(r_out) <- c(ncol_x,1,outbands)
	}
	
	
	
	out_image_dims=mapply(
			function(filenum,image_dims,outbands)
			{
				c(image_dims[2],image_dims[1],outbands)
			},
			filenum=as.list(seq(filename)),
			outbands=outbands,
			MoreArgs=list(image_dims=image_dims),SIMPLIFY=FALSE)

	if(!is.list(r_out)) r_out <- list(r_out)
	
#	chunk_position=list(
#			1:ncol_x,
#			chunk$row_center,
#			1:outbands
#	)
#	writeSuccess <- FALSE
#	while(!writeSuccess)
#	{
#		writeSuccess=TRUE
#		tryCatch(
#				binary_image_write(filename=filename,mode=real64(),image_dims=image_dims,
#						interleave="BSQ",data=r_out,data_position=chunk_position)
#				,
#				error=function(err) writeSuccess <<- FALSE)	
#	}	
	
	
	chunk_position <- lapply(X=outbands,
			FUN=function(X,chunk,ncol_x)
			{
				chunk_position=list(
						1:ncol_x,
						chunk$row_center,
						1:X
				)
				return(chunk_position)
			},
			chunk=chunk,ncol_x=ncol_x)
	
	writeMapply <- mapply(function(filename,image_dims,r_out,chunk_position)
			{
#				print(filename)
#				print(image_dims)
#				print(r_out)
				writeSuccess=FALSE
				
				while(!writeSuccess)
				{
					writeSuccess=TRUE
					tryCatch(
							binary_image_write(filename=filename,mode=datatype,image_dims=image_dims,
									interleave="BSQ",data=r_out,data_position=chunk_position)
							,
							error=function(err) writeSuccess <<- FALSE)	
				}
			},filename=as.list(filename),r_out=r_out,chunk_position=chunk_position,image_dims=out_image_dims)
}


focal_hpc_focal_processing <- function(tr,texture_tr,chunkArgs)
{
	# Create some blank variables:
	verbose <- NULL
	x <- NULL
	chunk_format <- NULL
	chunk <- NULL
	window_dims <- NULL
	datatype <- NULL
	.packages <- NULL
	
	list2env(chunkArgs,envir=environment())
	
	if(is.list(x))
	{
		r_old <- vector(mode="list",length=length(x))
	} else
	{
		r_old <- NULL
	}
	
	for(i in 1:tr$n)
	{
		if(verbose) cat("Iteration: ",i," of ",tr$n,"\n")
		
		if(is.list(x))
		{
			r <- mapply(FUN=function(X,tr,format,i,r_old,chunkArgs)
					{
						chunkArgs$x <- X
						# spatial.tools:::
						focal_hpc_focal_getChunk(x=X,tr=tr,format=format,i=i,r_old=r_old,
								chunkArgs=chunkArgs)	
					},X=x,r_old=r_old,
					MoreArgs=list(tr=tr,format=chunk_format,i=i,chunkArgs=chunkArgs),
					SIMPLIFY=FALSE
			)
		} else
		{
			r <- 
			#		spatial.tools:::
			focal_hpc_focal_getChunk(x=x,tr=tr,format=chunk_format,i=i,r_old=r_old,
							chunkArgs=chunkArgs)
			
		}
		# We need to divide up the chunks here.
		# This is going to cause memory issues if we aren't careful...
		j=1:tr$nrows[i]
		row_centers=tr$row[i]:tr$row2[i]
		
		#(tr$row[i]:tr$row2[i])-(tr$startrow_offset)
		chunkList=mapply(function(j,r,texture_tr,row_centers,chunk_format)
				{
					if(is.list(r))
					{
						processing_chunk <- sapply(X=r,FUN=function(X,texture_tr,j,chunk_format)
								{
									if(chunk_format=="array")
									{
										# Something is wrong with j
									#	print(j)
										processing_chunk=X[,texture_tr$row[j]:texture_tr$row2[j],,drop=FALSE]
									}
									if(chunk_format=="raster")
									{
										processing_chunk=getValuesBlock_enhanced(X,
												r1=texture_tr$row[j],r2=texture_tr$row2[j],
												format="raster")
									}
									return(processing_chunk)
								},texture_tr=texture_tr,j=j,chunk_format=chunk_format,
								simplify=FALSE)
					} else
					{
						if(chunk_format=="array")
						{
							processing_chunk=r[,texture_tr$row[j]:texture_tr$row2[j],,drop=FALSE]
						}
						if(chunk_format=="raster")
						{
							processing_chunk=getValuesBlock_enhanced(r,
									r1=texture_tr$row[j],r2=texture_tr$row2[j],
									format="raster")
						}
					}
					out_chunk=list(row_center=row_centers[j],
							processing_chunk=processing_chunk)
					return(out_chunk)
				}
				,j,MoreArgs=list(r=r,texture_tr=texture_tr,row_centers=row_centers,chunk_format=chunk_format),
				SIMPLIFY=FALSE)
		
		foreach(chunk=chunkList, .packages=c("rgdal","raster","spatial.tools","mmap",.packages),
						.verbose=verbose) %dopar% 
		#		spatial.tools:::
		focal_hpc_focalChunkFunction(chunk,chunkArgs)

		if(i<tr$n && window_dims[2] > 1)
			if(is.list(r))
			{
				r_old <- sapply(X=r,FUN=function(X,tr,i)
						{
							r <- X
							r_old <- array(data=r[,(dim(r)[2]-(tr$focal_row2[i]-tr$focal_row[i+1])):dim(r)[2],],
									dim=c(
											dim(r)[1],
											length((dim(r)[2]-(tr$focal_row2[i]-tr$focal_row[i+1])):dim(r)[2]),
											dim(r)[3])
							)
							return(r_old)
						},
						tr=tr,i=i,
						simplify=FALSE)
			} else
			{
				r_old <- array(data=r[,(dim(r)[2]-(tr$focal_row2[i]-tr$focal_row[i+1])):dim(r)[2],],
						dim=c(
								dim(r)[1],
								length((dim(r)[2]-(tr$focal_row2[i]-tr$focal_row[i+1])):dim(r)[2]),
								dim(r)[3])
				)
			}
	}
}

focal_hpc_pixelChunkFunction <- function(chunkID,tr,x,
		chunk_format,fun,fun_args,layer_names,outbands,filename,datatype)
{
	
	# Seeing some memory creep, hopefully this will help:
	# gc()
	# Read the chunk
	if(is.list(x))
	{
		r <- sapply(X=x,FUN=function(X,chunkID,chunk_format)
				{
					
					getValuesBlock_enhanced(X,r1=tr$row[chunkID],r2=tr$row2[chunkID],
							c1=1,c2=ncol(X),format=chunk_format)
				},chunkID=chunkID,chunk_format=chunk_format,simplify = FALSE)
		
		fun_args$x=r
		
#		if(chunk_format=="array")
#		{
#			# FIX THIS
#			for(i in seq(length(x)))	
#			{
#				dimnames(fun_args$x[[i]]) <- vector(mode="list",length=3)
#				if(!is.null(layer_names[[i]])) dimnames(fun_args$x[i])[[3]]=layer_names[[i]]
#			}
#		}
	} else
	{
		r <- getValuesBlock_enhanced(x,r1=tr$row[chunkID],r2=tr$row2[chunkID],
				c1=1,c2=ncol(x),format=chunk_format)
		
		fun_args$x=r
		
		if(chunk_format=="array")
		{
			dimnames(fun_args$x)=vector(mode="list",length=3)
			if(!is.null(layer_names)) dimnames(fun_args$x)[[3]]=layer_names
		}
	}
	
	# Execute the function.
	r_out <- do.call(fun, fun_args)
	
	if(!is.list(r_out)) r_out <- list(r_out)
	
	# Write the output
	
	
	if(is.list(x))
	{
		image_dims=dim(x[[1]])
		image_dims=c(image_dims[2],image_dims[1],image_dims[3])
	} else
	{
		image_dims=dim(x)
		image_dims=c(image_dims[2],image_dims[1],image_dims[3])
	}
	
	chunk_position <- lapply(X=outbands,
			FUN=function(X,x,tr)
			{
				if(is.list(x))
				{
#					image_dims=dim(x[[1]])
#					image_dims=c(image_dims[2],image_dims[1],image_dims[3])
					chunk_position=list(
							1:ncol(x[[1]]),
							tr$row[chunkID]:tr$row2[chunkID],
							1:X
					)
				} else
				{
#					image_dims=dim(x)
#					image_dims=c(image_dims[2],image_dims[1],image_dims[3])
					chunk_position=list(
							1:ncol(x),
							tr$row[chunkID]:tr$row2[chunkID],
							1:X)
				}
				return(chunk_position)
			},
			x=x,tr=tr)
	
	writeMapply <- mapply(function(filename,image_dims,r_out,chunk_position)
			{
#				print(filename)
#				print(image_dims)
#				print(r_out)
				writeSuccess=FALSE
				
#				browser()
				
				while(!writeSuccess)
				{
					writeSuccess=TRUE
					tryCatch(
							binary_image_write(filename=filename,mode=datatype,image_dims=image_dims,
									interleave="BSQ",data=r_out,data_position=chunk_position)
							,
							error=function(err) writeSuccess <<- FALSE)	
				}
			},filename=as.list(filename),r_out=r_out,chunk_position=chunk_position,MoreArgs=list(image_dims=image_dims)
	)
	
	
	return(NULL)
}

focal_hpc_pixel_processing <- function(tr,chunkArgs)
{
	# Create some blank variables:
	x <- NULL
	chunk_format <- NULL
	fun <- NULL
	layer_names <- NULL
	outbands <- NULL
	verbose <- NULL
	datatype <- NULL
	.packages <- NULL
	
	list2env(chunkArgs,envir=environment())
	chunkID <- seq(tr$n)
	foreach(chunkID=chunkID, .packages=c("rgdal","raster","spatial.tools","mmap",.packages),.verbose=verbose) %dopar% 
	#		spatial.tools:::
	focal_hpc_pixelChunkFunction(chunkID,tr,x,chunk_format,fun,args,layer_names,outbands,
					filename,datatype)
}

#' Engine for performing fast, easy-to-develop pixel and focal raster calculations with parallel processing capability.
#' @param x Raster*. A Raster* used as the input into the function.  Multiple inputs should be stack()'ed or list()'ed together.
#' @param fun function. A focal function to be applied to the image. See Details.
#' @param args list. Arguments to pass to the function (see ?mapply).  Note that the 'fun' should explicitly name the variables.
#' @param window_dims Vector. The size of a processing window in col x row order.  Be default, a single pixel (c(1,1).
#' @param window_center Vector. The local coordinate of the center of a processing window.  By default the middle of the processing window.  UNSUPPORTED.
#' @param chunk_format Character. The format to send the chunk to the function.  Can be "array" (default) or "raster".
#' @param minblocks Numeric. The minimum number of chunks to divide the raster into for processing.  Defaults to 1.
#' @param blocksize Numeric. The size (in rows) for a block of data.  If unset, focal_hpc will attempt to figure out an optimal blocksize.
#' @param filename Character. Filename(s) of the output raster.
#' @param outformat Character. Outformat of the raster. Must be a format usable by hdr(). Default is 'raster'. CURRENTLY UNSUPPORTED.
#' @param overwrite Logical. Allow files to be overwritten? Default is FALSE.
#' @param processing_unit Character. ("single"|"chunk") Will be auto-set if not specified ("chunk" for pixel-processing, "single" for focal processing).  See Description.
#' @param outbands Numeric. If known, how many bands in each output file?  Assigning this and outfiles will allow focal_hpc to skip the pre-check.
#' @param outfiles Numeric. If known, how many output files?  Assigning this and outbands will allow focal_hpc to skip the pre-check.
#' @param setMinMax Logical. Run a setMinMax() on each output file after processing (this will slow the processing down). Default is FALSE.
#' @param additional_header Character. Create additional output headers for use with other GIS systems (see \code{\link{hdr}}). Set to NULL to suppress.  Default is "ENVI".
#' @param datatype Character.  Output number type.  See ?dataType.  Default is "FLT8S".  
#' @param debugmode Logical or Numeric.  If TRUE or 1, the function will enter debug mode during the test phase.  If debugmode equals 2, the function will stop after the test phase, but won't explicitly enter debug mode.  This is useful if the user function has a browser() statement within it.  Note the inputs will be an array of size 2 columns, 1 row, and how ever many input bands.
#' @param .packages Character. A character vector of package names needed by the function (parallel mode only).
#' @param clearworkers Logical. Force the workers to clear all objects upon completing (releasing memory)?  Default=TRUE.
#' @param verbose logical. Enable verbose execution? Default is FALSE.  
#' @param ... Additional parameters (none at present).
#' @author Jonathan A. Greenberg (\email{spatial.tools@@estarcion.net})
#' @seealso \code{\link{rasterEngine}}, \code{\link{foreach}}, \code{\link{mmap}}, \code{\link{dataType}}, \code{\link{hdr}} 
#' @details focal_hpc is designed to execute a function on a Raster* object using foreach, to
#' achieve parallel reads, executions and writes. Parallel random writes are achieved through the use of
#' mmap, so individual image chunks can finish and write their outputs without having to wait for
#' all nodes in the cluster to finish and then perform sequential writing.  On Windows systems,
#' random writes are possible but apparently not parallel writes.  focal_hpc solves this by trying to
#' write to a portion of the image file, and if it finds an error (a race condition occurs), it will
#' simply retry the writes until it successfully finishes.  On Unix-alikes, truly parallel writes
#' should be possible.
#' 
#' Note that \code{\link{rasterEngine}} is a convienence wrapper for focal_hpc and, in general, should be used instead
#' of focal_hpc directly.  
#'
#' focal_hpc operates in two modes, which have different input and outputs to the function:
#' 
#' Pixel based processing:
#' 
#' 1) If chunk_format=="array" (default), the input to the function should assume an array of dimensions 
#' x,y,z where x = the number of columns in a chunk, y = the number of rows in the chunk, and 
#' z = the number of bands in the chunk.  If chunk_format=="raster", the input to the function
#' will be a raster subset.
#' Note that we are ordering the array using standards for geographic data, (columns, rows, bands), 
#' not how R usually thinks of arrays (rows, columns, bands).
#' 
#' 2) The output of the function should always be an array with the x and y dimensions matching
#' the input, and an arbitrary number of band outputs.  Remember to order the dimensions as
#' columns, rows, bands (x,y,z).
#' 
#' Local window processing:
#' 
#' 1) The function should be written to process a SINGLE window at a time, given the dimensions
#' of window_dims, so the input to the function should assume a window of dimensions window_dims 
#' with a local center defined by window_center.  As with before, the input can be passed to 
#' the function as an array (suggested) or a small raster.
#' 
#' 2) The output should be a single pixel value, so can either be a single value, or a vector
#' (which is assumed to be multiple bands of a single pixel).
#' 
#' The speed of the execution when running in parallel will vary based on the specific setup, 
#' and may, indeed, be slower than a sequential execution (e.g. with calc() ), 
#' particularly on smaller files.  Note that by simply running sfQuickStop(), focal_hpc
#' will run in sequential mode.
#' 
#' @examples

#'  tahoe_highrez <- brick(system.file("external/tahoe_highrez.tif", package="spatial.tools"))
#' # Pixel-based processing:
#'ndvi_function <- function(x)
#'{
#'	# Note that x is received by the function as a 3-d array:
#'	red_band <- x[,,2]
#'	nir_band <- x[,,3]
#'	ndvi <- (nir_band - red_band)/(nir_band + red_band)
#'	# The output of the function should also be a 3-d array,
#'	# even if it is a single band:
#'	ndvi <- array(ndvi,dim=c(dim(x)[1],dim(x)[2],1))
#'	return(ndvi)
#'}
#' 
#'sfQuickInit(cpus=2)
#'tahoe_ndvi <- focal_hpc(x=tahoe_highrez,fun=ndvi_function)
#'sfQuickStop()
#' 
#' \dontrun{ 
#'# Focal-based processing:
#'local_smoother <- function(x)
#'{
#'	# Assumes a 3-d array representing
#'	# a single local window, and return
#'	# a single value or a vector of values.
#'	smoothed <- apply(x,3,mean)
#'	return(smoothed)
#'}
#' 
#' # Apply the function to a 3x3 window:
#' sfQuickInit(cpus=2)
#' tahoe_3x3_smoothed <- focal_hpc(x=tahoe_highrez,fun=local_smoother,window_dims=c(3,3))
#' sfQuickStop()
#' 
#' # Example with 7 x 7 window in full parallel mode:
#' sfQuickInit()
#' tahoe_7x7_smoothed <- focal_hpc(x=tahoe_highrez,fun=local_smoother,window_dims=c(7,7))
#' sfQuickStop()
#' }
#' 
#' @importFrom abind abind
#' @import iterators
#' @import foreach
#' @import mmap
#' @import raster
#' @import compiler
#' 
#' @export

focal_hpc <- function(x,
		fun,args=NULL, 
		window_dims=c(1,1), 
		window_center=c(ceiling(window_dims[1]/2),ceiling(window_dims[2]/2)),
		filename=NULL, overwrite=FALSE,
		outformat="raster",additional_header="ENVI",
		datatype="FLT8S",
		processing_unit=NULL,chunk_format="array",
		minblocks="max",blocksize=NULL,
		outbands=NULL,outfiles=NULL,
		setMinMax=FALSE,
		debugmode=FALSE,
		.packages=NULL,
		clearworkers=TRUE,
		verbose=FALSE,
		...) 
{
	# Required libraries:
#	require("raster")
#	require("foreach")
#	require("rgdal")
#	require("mmap")
#	require("abind")
	
	# Create some blank variables to avoid warnings:
	
	layer_names <- NULL
	startrow_offset <- NULL
	endrow_offset <- NULL
#	processing_unit <- NULL
	chunk_nrows <- NULL
	tr <- NULL
	processing_mode <- NULL
	texture_tr <- NULL
	
# 	Make sure all the packages are loaded.
	loaded_packages <- lapply(.packages, require, character.only=T)
	
#	if(!is.list(x)) { x <- list(x) }
	
	# Register a sequential backend if one is not already registered:
	if(!getDoParRegistered()) 
	{
		if(verbose) { warning("No parallel backend registered.  Operating in sequential mode.") }
		registerDoSEQ()
	}
	# Prechecks.
	list2env(
	#		spatial.tools:::
	focal_hpc_precheck(x,window_dims,window_center,processing_unit,verbose),envir=environment())
	
	# Add in new formals to the function if args doesn't match the function:
	
	
	# Fix missing ellipses in function.  Thanks to Ista Zahn for the solution.
	# http://r.789695.n4.nabble.com/Checking-for-and-adding-arguments-to-a-function-tp4685450p4685452.html
	
	base_formals <- formals(fun)
	base_formals_names <- names(base_formals)
	# Add in args if missing
	missing_formals_names <- setdiff(names(args),base_formals_names)
	missing_formals <- args[names(args) %in% missing_formals_names]
	
	
	new_formals <- c(
			unlist(base_formals[base_formals_names != "..."],recursive=FALSE),
	#		unlist(missing_formals,recursive=FALSE),
			missing_formals,
			alist(...=)
	)
	
	formals(fun) <- new_formals
	
	
	# Debug mode:
	if(debugmode==TRUE) debug(fun)
#	else(undebug(fun))
	
	# Test focal_hpc and determine the number of outbands and files.
	if(is.null(outbands) || is.null(outfiles))
	{
		outbands_and_files <- 
		#		spatial.tools:::
		focal_hpc_test(x,fun,window_center,window_dims,args,layer_names,
						startrow_offset,endrow_offset,
						processing_mode,processing_unit,chunk_format,
						verbose)
		outbands=outbands_and_files$outbands
		outfiles=outbands_and_files$outfiles
	}
	
	if(debugmode==TRUE || debugmode==2) 
	{
		undebug(fun)
		stop("Debug mode finished.  Stopping run.")
	}
	
	# Set up chunking parameters.
	list2env(
	#		spatial.tools:::
	focal_hpc_chunk_setup(
					x=x,window_dims=window_dims,window_center=window_center,
					chunk_nrows=chunk_nrows,startrow_offset=startrow_offset,
					endrow_offset=endrow_offset,
					minblocks=minblocks,
					blocksize=blocksize,
					verbose=verbose),
			envir=environment())
	
	# print(class(x))
	
	# Create blank image file.
	if(!is.list(x)) reference_raster <- x else reference_raster <- x[[1]]
	
	# This should be moved to pre-checks:
	if(is.null(filename)) filename <- lapply(seq(outfiles),function(x) NULL)
	
	out <- mapply(function(filename,outbands,reference_raster,overwrite,verbose)
			{
				create_blank_raster(filename=filename,
						format="raster",datatype=datatype,bandorder="BSQ",
						nlayers=outbands,
						create_header=TRUE,reference_raster=reference_raster,
						additional_header=additional_header,
						overwrite=overwrite,verbose=verbose)
			},filename=filename,outbands=outbands,
			MoreArgs=list(reference_raster=reference_raster,overwrite=overwrite,verbose=verbose)
	)
	
	# Create chunk arguments.
	if(verbose) { message("Loading chunk arguments.") }
	chunkArgs = list(fun=fun,x=x,x_ncol=ncol(x),tr=tr,
			window_dims=window_dims,window_center=window_center,
			layer_names=layer_names,
			args=args,filename=out,
			outbands=outbands,processing_unit=processing_unit,
			verbose=verbose,layer_names=layer_names,
			chunk_format=chunk_format,
			.packages=.packages,
			datatype=datatype)
	
	# Processing:
	if(processing_mode=="focal")
	{
		#spatial.tools:::
		focal_hpc_focal_processing(tr,texture_tr,chunkArgs)
	} else
	{
		#spatial.tools:::
		focal_hpc_pixel_processing(tr,chunkArgs)
	}
	
# browser()
	
# Flush the workers.

	if(clearworkers && getDoParName() != "doSEQ")
	{
		cleared <- foreach(nworkers=seq(getDoParWorkers())) %dopar%
				{
					rm(list=ls(all.names=TRUE))
					return(TRUE)
				}
	}

	focal_out <- sapply(X=as.list(out),FUN=function(X,additional_header,setMinMax)
			{
				focal_out <- brick(X)
				if(setMinMax) focal_out <- setMinMax(focal_out) else focal_out@data@haveminmax <- FALSE
				if(!is.null(additional_header))
					hdr(focal_out,additional_header)
				return(focal_out)
			},additional_header=additional_header,setMinMax=setMinMax)
	
	if(length(focal_out)==1) return(focal_out[[1]]) else return(focal_out)
	
}


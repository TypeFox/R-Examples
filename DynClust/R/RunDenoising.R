
#' @name RunDenoising
#' @title Denoising step of a dynamical image sequence
#' @description Performs the denoising step of a dynamic sequence
#' of images. It is also the first step of the clustering.
#'
#' The denoising procedure is iteratively applied on each voxel.
#'
#' The denoised version of Fx is obtained with a three stages procedure:
#' 1) Selection of time-homogeneous voxels in the sub-mask around the voxel of interest;
#' 2) Growth of spatial neighborhoods made of time-homogeneous voxels obtained at stage 1 with sizes growing geometrically---each neiborhood is associated to a denoised version by averaging over its members;
#' 3) Selection of the largest spatial neighborhood such that its associated denoised version is time-homogeneous with all the previous ones.
#' 
#' Time homogeneity is tested with function \code{MultiTestH0}.
#'
#' Further details about the denoising method and the statistical test of homogeneity can be found in the references.
#'
#' @seealso \code{\link{GetDenoisingResults}}, \code{\link{MultiTestH0}}.
#' 
#' @param data.array a (2D or 3D)+T array containing the dynamic
#' sequence of images (the dataset). The last dimension is the time.
#' @param data.var a numeric indicating the variance of the dataset (default 1).
#' If set to NULL, the variance is computed using a baseline image. See \code{enhStart} parameter.
#' @param depth a numeric indicating the depth of a voxel.
#' @param alpha a numeric value indicating the global level of the
#' multitest.
#' @param mask.size a vector indicating the size of the spatial
#' hypercube defined around voxels used to search for neighbors.
#'
#' If NA (default): sqrt(dim(data.array)[1:length(dim(data.array))-1]).
#'
#' If NULL (complete image): dim(data.array)[1:length(dim(data.array))-1]
#' 
#' @param nproc a numeric value indicating the number of processors
#' used for parallel computation.
#' @param enhStart an integer, if larger than 1, a baseline is
#' computed as a median image obtain from time indexes between 1 and enhStart-1.
#' Default value \code{ifelse(is.null(var),2,1)}.
#4 This baseline is removed from each image of the sequence.
#' @return a list containing:
#' \itemize{
#' \item \code{info.den}, a list of list whose length is the number of voxels,
#' each sub-list contains the result of buildEstimate for one voxel.
#' \item \code{data.proj}, the projections of the dynamics.
#' a list containing a
#' denoised version of the dataset as an array, as well as a list for
#' which each element contains a list with the voxel index, the
#' indexes of its neighboors, the resulting denoised signal, and the variance of the
#' denoised signal
#' \item \code{var}, a numeric providing the known variance
#' }
#' @references Rozenholc, Y. and Reiss, M. (2012) \emph{Preserving time structures while denoising a dynamical image}, Mathematical Methods for Signal and Image Analysis and Representation (Chapter 12),  Florack, L. and Duits, R. and Jongbloed, G. and van~Lieshout, M.-C. and Davies, L. Ed., Springer-Verlag, Berlin
#' @references Lieury, T. and Pouzat, C. and Rozenholc, Y. (submitted) \emph{Spatial denoising and clustering of dynamical image sequence: application to DCE imaging in medicine and calcium imaging in neurons}  
#' @author Tiffany Lieury, Christophe Pouzat, Yves Rozenholc
#' @example ../R/Denoising-example.R
#' @export
RunDenoising <- function(data.array,                         
                         data.var=1,                                       
                         depth=1,
                         alpha=.05,
                         mask.size=NA,
                         nproc=1,
                         enhStart=ifelse(is.null(var),2,1)) {
    
     ## sets the sizes of the crowns
     ## sizes increase geometrically as a serie with a0=1 and q=1.5
     incr.size <- c(1,4,8,16,36,92,212,477) #increasing crown sizes 
     ball.size <- cumsum(incr.size) #ball sizes
     
     ## data dimensions
     dim <- dim(data.array)
     ## number of dimensions
     ndim <- length(dim)
     ## special status for time
     ntime <- dim[ndim]
     ## spatial dimensions
     coord <- dim[-ndim]
     ## number of voxels
     nvox <- prod(coord)
     ## transforms the kD+T data set into a matrix ntime x nvox
     data.array <- matrix(aperm(data.array,c(ndim,1:(ndim-1))),nrow=ntime,ncol=nvox)
     ## number of partitions = number of tests = iter + 1
     iter <- floor(log2(ntime))-1      ### partition sizes 2^(0:iter)
     ## size of the finest partition
     Dmax <- 2^iter
     ## number of intervals + 1 in all partitions     
     from = 2^(iter+1)
     ## test thresholds with Bonferroni correction adapted to the partition number
     thrs = qchisq(1-alpha/(iter+1),2^(0:iter))
     
     
     ## mask initialization
     if ((!is.null(mask.size))&&(is.na(mask.size))) {
         p=(1000/prod(coord))^(1/length(coord))
         mask.size = ceiling(p*coord)
         print(paste('mask size :'))
         print(mask.size)
     }
     if (is.null(mask.size)) mask.size=coord
     
     ## ################### TESTING FEASIBILITY #############################
     ## tests can the data be analyzed
     if(nproc<1) stop("nproc must be at least equal to 1")
     if(!is.null(dim(data.var)) | length(data.var)>1) stop("data.var should be a numeric")
     if(length(dim)<=2) stop("number of dimensions should not be smaller than 3")
     if(!is.null(mask.size)){
         if (any(mask.size<=1)) stop('coef of mask.size should be larger than 1')
         mask.size <- round(mask.size)}
     
    ## ############# BASELINE and VARIANCE COMPUTATION  ####################
     ## baseline removal and computation of the variance if needed
     baseline = NULL
     if (enhStart>1) {
         baseline = apply(data.array[1:(enhStart-1),],2,median)
         data.array = sweep(data.array,2,baseline)
         if (is.null(data.var)||(data.var<0)) {
             data.var <- ifelse(is.null(data.var),1,-data.var)*var(as.vector(data.array[1:(enhStart-1),]))
             print(paste('estimated variance:',data.var))
         } 
     }
     
     ## ##################### DATA NORMALIZATION ##############################
     ## standardized the data with respect to the common variance
     keep.sd = sqrt(data.var)
     if (keep.sd!=1) data.array = data.array/keep.sd
     
     data.var = 1
     
        
    ## ########################### CREATING LOCAL FUNCTIONS #############################
        
    RetProj <- function(){      
        ## Inialize with the finest partition
        ## locations of time indexes in the finest partition
        loc <- ceiling((1:ntime)*Dmax/ntime)
        ## lengths of the finest partition intervals // returned as a vector
        num <- c(tapply(rep(1,ntime),loc,sum))
        ## storage locations 
        to <- from-1; 
        from <- to-length(num)+1
        ## projections 
        data.proj[from:to,] <<- apply(data.array,2,function(col) tapply(col,loc,sum))
        ## loop from finer to thicker partitions
        for (idx in (iter-1):0) {
            ## new lengths of the partition intervals
            num.new <- num[seq(1,length(num),by=2)]+num[seq(2,length(num),by=2)]
            ## new storage locations 
            to.new <- from-1; 
            from.new <- to.new-length(num.new)+1
            ## new projections
            data.proj[from.new:to.new,] <<- 
                data.proj[seq(from,to,by=2),]+data.proj[seq(from+1,to,by=2),]
            ## normalize old projections
            data.proj[from:to,] <<- data.proj[from:to,]/sqrt(num)
            ## update
            num <- num.new
            to <- to.new
            from <- from.new
        }
        ## last normalization is not done in the loop 
        data.proj[1,] <<- data.proj[1,]/sqrt(num)
      }    
            
    ## local function returning a the denoised estimation of the time dynamics and the current voxel neighboors
    ## used to compute it. 
    ## Parameters :
    ## a numeric value indicating the index of the current voxel to denoise
    buildEstimate <- function(pix.idx){
        ## local function used to find neighbors of a center voxel indexed by pix.idx
        ## neighbors are voxels that are time-homogeneous with the center       
        ## neighbors indexes are returned ordered with respect to the euclidian distance to the center
        findNeighbors <- function(pix.idx){
            ## Find the neighboors of the current voxel pix.idx
            ## build the mask hypercube around pix.idx
            inf <- pmax(1,data.coord[pix.idx,]-mask.size)
            sup <- pmin(coord,data.coord[pix.idx,]+mask.size)
            tmp <- paste(paste(inf,':',sup,sep=''),collapse=',')
            mask <- eval(parse(text=paste("c(arr.coord[",tmp,"])")))
                                    
            ## ##### THIS IS NOT EFFICIENT AS ALL MASK VOXELS ARE TESTED ##########
            ## test in mask voxels which are homogenous with fd.idx.pix 
            bool = MultiTestH0(data.proj[,mask,drop=FALSE]-data.proj[,pix.idx],2,thrs)
            ## ##### SELECTION SHOULD BE APPLIED IN THE RING CONSTRUCTION #########
                            
            ## retrieve voxels from mask which are time homogenous with pix.idx
            neighbors.idx <- mask[bool]
            
            ## computes the euclidean distance (in space) between pix.idx and its neighboors
            ## the euclidian distance from pix.idx to voxels of the neiborhoods
            dist <- colSums((t(data.coord[neighbors.idx,,drop=F])-data.coord[pix.idx,])^2)
            
            ## voxel indexes ordered with respect to the increasing euclidian distance
            dist.tri <- order(dist)
            
            ## returns the sorted neighborhood indexes 
            neighbors.idx <- neighbors.idx[dist.tri]
            neighbors.idx
        }
        ## indexes of the neighbors sorted by euclidian distance
        ## pix.idx is its first neighbors as distance=0
        neighbors.idx <- findNeighbors(pix.idx) 
                      
        ## projection of the dynamics of the neighbors
        neighbors.proj <- data.proj[,neighbors.idx,drop=F]
      
        ## Denoising procedure
        ## V1 à V8 are the neighboorhoods of the voxel x time x of increasing size. Here we want to find the biggest Vx 
        ## which is still coherent with the smaller ones.
        ## number of possible balls 
        nV <- sum(ball.size<=length(neighbors.idx))
        ## Iv.neighb is a list of the neighbors in each ball
        Iv.neighb <- list()
        ## Iv is the matrix of the projection estimates build over the successive balls
        Iv <- matrix(NA,nrow=nproj,ncol=nV)
        ## data.varIv vector of the associated variances 
        data.varIv <- rep(NA,l=nV)
        ## Jv is a matrix of the projection estimates build over the successive ring 
        Jv <- matrix(NA,nrow=nproj,ncol=nV)
        ## data.varJv vector of the associated variances 
        data.varJv <- rep(NA,l=nV)
      
        ## Initialize
        limits = kV = 1
        Iv[,kV] = neighbors.proj[,1] # pix.idx is at the first place
    
        data.varIv[kV] = data.var
        while (kV<nV){
            kV = kV+1
            limits.new <- ball.size[kV]
            ring.idx = (limits+1):limits.new
            
            Jv[,kV] <- rowMeans(neighbors.proj[,ring.idx,drop=F])
            data.varJv[kV] <- data.var/length(ring.idx)
            
            ## test thresholds with Bonferroni correction adapted to both partition number and interior balls
            thrs = qchisq(1-alpha/(iter+1)/(kV-1),2^(0:iter))
      
            ## test time coherence
            testcoh <- MultiTestH0(Iv[,1:(kV-1),drop=F]-Jv[,kV],
                                      data.varIv[1:(kV-1)]+data.varJv[kV],
                                      thrs)
            ## if no time coherence with previous estimates // BREAK
            if(sum(testcoh)< kV-1) { ## one test at least is false
                kV = kV-1
                break
            }
            ## otherwise update projection estimates
            limits = limits.new
            Iv[,kV] <- rowMeans(neighbors.proj[,1:limits,drop=F])
            data.varIv[kV] <- data.var/limits
     
        }
        
        ## the denoised dynamics rescaled
        Ix = keep.sd*rowMeans(data.array[,neighbors.idx[1:limits],drop=F])
        if (!is.null(baseline)) Ix = Ix+baseline[pix.idx]
        
        list(Lx=data.coord[neighbors.idx[1:limits],],Cx=data.coord[pix.idx,],Px=Iv[,kV],
             Ix=Ix,Vx= neighbors.idx[1:limits]) 
        
    #### returns a list containing:
    #### 'Vx' a vector containing all the neighbors indexes used to build the denoised dynamic
    #### 'Ix' a vector containing the denoised dynamic
    #### 'Px' a vector containing the denoised projection
    #### 'Lx' a matrix containing the original coordinates of the neighbors in data.array 
    #### 'Cx' a vector containing the original coordinates of the center
    }            
                
    ## Define the storage matrix
    data.proj <- matrix(nrow=from-1,ncol=nvox)
    
    RetProj()       
    
    ## ########################### Build geometry #############################    
    ## matrix of all the 3D coordinates at column indexes of data.matrix corresponds to row indexes in data.coord
    tmp <- paste(paste("1:dim[",1:(ndim-1),"]"),collapse=',')
    data.coord <- eval(parse(text=paste("expand.grid(",tmp,")")))
    data.coord <- as.matrix(data.coord)
    arr.coord <- array(1:nvox,dim=coord)
    
    ## transforms the data into a matrix of projections
    nproj <- dim(data.proj)[1]
    
    ## ########################## LAUNCH ANALYSIS #############################     
    ## for each voxel x time to visit do
    res.visited <- mclapply(1:nvox,buildEstimate,mc.preschedule=FALSE,mc.cores=nproc)
    
    list(info.den=res.visited,data.proj=data.proj,var=keep.sd^2,baseline=baseline)    
}

#' NIfTI-to-R conversion
#' 
#' Reads in a NIfTI (.nii) file and puts the data in a 4-dimensional array.
#' 
#' 
#' @param niifilename the path for the .nii file.
#' @param which.vols which volumes (images) to include.  In terms of the 4D
#' array, this refers to subsetting in the fourth dimension. If \code{NULL}
#' (the default), all volumes are included.
#' @param savename if non-\code{NULL}, the name of the .RData file to which the
#' 4D array will be saved.
#' @param remove.zero optional when maskname is not provided. If \code{TRUE}, a
#' binary array indicating the voxels with nonzero measures based on the first
#' three dimension of the nii file will be provided. If \code{FALSE}, a 3D
#' array with \code{TRUE} everywhere will be provided.
#' @param maskname name of a .nii file providing a "mask", a 3D binary array
#' indicating which voxels to include.
#' @param ind,ind.auto \code{ind} is an optional list saying which indices
#' (which slices of the image) to include in each of the three dimensions.  If
#' \code{NULL}, this will be all slices with nonzero data if \code{ind.auto =
#' TRUE}, and all slices otherwise.
#' @param coord coordinates of the first three dimensions of the 4D array
#' created.
#' @return a 4-dimensional array.
#' @author Lei Huang \email{huangracer@@gmail.com}, Philip Reiss
#' \email{phil.reiss@@nyumc.org} and Rong Jiao \email{jiaorong007@@gmail.com}
#' @seealso \code{\link{R2nii}}
#' @export
nii2R <-
function(niifilename, which.vols=NULL, savename=NULL, remove.zero=TRUE, maskname=NULL, ind=NULL, 
         ind.auto=TRUE, coord=NULL)   {
    
    bigobj = oro.nifti::readNIfTI(niifilename)[]
	  if (length(dim(bigobj))==3) {
		  tmp = array(dim=c(dim(bigobj), 1))
		  tmp[ , , , 1] = bigobj[]
		  bigobj = tmp
		  rm(tmp)
	  }
    
    if (is.null(ind) & ind.auto)    ind = get.ind(bigobj[,,,1])
    else if (is.null(ind) & !ind.auto)    for (i in 1:3)  ind[[i]]=1:dim(bigobj[])[i]
    x.ind = ind[[1]]
    y.ind = ind[[2]]
    z.ind = ind[[3]]
    if (is.null(which.vols)) which.vols = 1:(dim(bigobj)[4])
    
    if (is.null(coord)) {
  		coord = list()
  		for (l in 1:3) coord[[l]] = 1:(dim(bigobj[])[l])
    }
    
    coord[[1]] = coord[[1]][x.ind]
    coord[[2]] = coord[[2]][y.ind]
    coord[[3]] = coord[[3]][z.ind]
        
    dir.create("./temp")
    cat("Saving each image to temp directory...\n")
    for (l in which.vols)   {
        cat("Image #",l,"\n")
        temp = bigobj[x.ind, y.ind, z.ind, l]
        save(temp, file=paste("./temp/Sep-", l, ".RData", sep=""))
    }
    dim.nii = dim(bigobj[,,,1])
    rm(temp, bigobj)
    d4 = array(NA, c(length(x.ind), length(y.ind), length(z.ind), length(which.vols)))
    cat("Adding each image 4D array...\n")
    i=0
    for (l in which.vols)   {
        i = i+1
        cat("Image #",i,"\n")
        load(paste("./temp/Sep-", l, ".RData", sep=""))
        d4[,,,i]=temp
    }
    
    if (!is.null(maskname))   has.data = oro.nifti::readNIfTI(maskname)[x.ind,y.ind,z.ind]!=0
    else if (remove.zero)   has.data =  apply(d4, 1:3, function(mat) !all(mat==0|is.infinite(mat)|is.na(mat)))
         else   has.data = array(TRUE, dim(d4)[1:3])
    
    attr(d4, "x.ind") = x.ind
    attr(d4, "y.ind") = y.ind
    attr(d4, "z.ind") = z.ind
    attr(d4, "which.vols") = which.vols
    attr(d4, "dim.nii") = dim.nii
    attr(d4, "coord") = coord
    attr(d4, "has.data") = has.data
    if (!is.null(savename)) save(d4, file=paste(savename,".RData",sep=""))
    return(d4)
}


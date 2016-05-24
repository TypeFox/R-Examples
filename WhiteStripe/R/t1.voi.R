#' @title t1.voi.hist
#'
#' @description Histogram of VOI of T1 template image
#' @name t1.voi.hist
#' @docType data
#' @aliases t1.voi.hist
#' @usage t1.voi.hist
#' @format A volume of interest histogram from a T1 image for smoothing
#' @keywords datasets
#' @examples
#' \dontrun{
#' if (download_img_data()){
#' t1 = readNIfTI(system.file("T1Strip.nii.gz", package="WhiteStripe"))
#' t1.voi = make_img_voi(t1)
#' any(is.na(t1.voi))
#' # FALSE
#' t1.voi.hist = hist(t1.voi, 
#' breaks=2000, 
#' plot=FALSE) 
#' #save(t1.voi.hist, file="data/t1.voi.hist.rda", compress = TRUE,
#' # compression_level=9)
#' }
#' } 
NULL

#' @title t2.voi.hist
#'
#' @description Histogram of VOI of T2 template image
#' @name t2.voi.hist
#' @docType data
#' @aliases t2.voi.hist
#' @usage t2.voi.hist
#' @format A histogram volume of interest from a T2 image for smoothing
#' @keywords datasets
#' @examples
#' \dontrun{
#' if (download_img_data()){
#' t2 = readNIfTI(system.file("T2Strip.nii.gz", package="WhiteStripe"))
#' t2.voi = make_img_voi(t2)
#' any(is.na(t2.voi))
#' # FALSE 
#' t2.voi.hist = hist(t2.voi, 
#' breaks=2000, 
#' plot=FALSE)  
#' #save(t2.voi.hist, file="data/t2.voi.hist.rda", compress = TRUE,
#' # compression_level=9) 
#' }
#' } 
NULL

#' @title smoothed_histogram
#'
#' @description Smoothed histogram of image
#' @name s.hist
#' @docType data
#' @aliases s.hist
#' @usage s.hist
#' @format A GAM from mgcv for x and y from histograms
#' @keywords datasets
#' @examples
#' \dontrun{ 
#' data(t2.voi.hist)
#' y = t2.voi.hist$counts
#' x = t2.voi.hist$mids
#' x = x[!is.na(y)];
#' y = y[!is.na(y)]
#' # 70 used for speed of example
#' s.hist = smooth_hist(x, y, k=70)
#' }
NULL

#' @title midpoints of VOI histogram
#'
#' @description Points from VOI histogram
#' @name xvals
#' @docType data
#' @aliases xvals
#' @usage xvals
#' @format x values from histogram for VOI
#' @keywords datasets
NULL

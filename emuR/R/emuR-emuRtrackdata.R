requireNamespace("data.table", quietly = T)
requireNamespace("dplyr", quietly = T)

##' create emuRtrackdata object
##' 
##' Joins \code{\link{emuRsegs}} and \code{\link{trackdata}} objects
##' to create an \code{\link{emuRtrackdata}} object that is a sub-class of
##' a \code{\link{data.table}} (and \code{\link{data.frame}}) object. This object 
##' can be viewed as a flat version of a \code{\link{trackdata}} object that also 
##' contains all the information of a \code{\link{emuRsegs}} object. It is meant to
##' ease integration with other packages as it is based on the well known 
##' \code{\link{data.table}} and \code{\link{data.frame}} objects.
##' @param sl seglist of class \code{\link{emuRsegs}}
##' @param td \code{\link{trackdata}} object generated from sl
##' @return emuRtrackdata object
##' @import data.table
##' @export
##' @examples
##' \dontrun{
##' 
##' ##################################
##' # prerequisite: loaded ae emuDB 
##' # (see ?load_emuDB for more information)
##' 
##' # query emuDB (to get object of class emuRsegs)
##' sl = query(emuDBhandle = ae, 
##'            query = "Phonetic == i:")
##'            
##' # get formats for SEGMENTs in sl (to get object of class trackdata)
##' td = get_trackdata(emuDBhandle = ae, 
##'                    seglist = sl,
##'                    onTheFlyFunctionName = "forest")
##' 
##' # create emuRtrackdata object
##' create_emuRtrackdata(sl = sl, td = td)
##' 
##' }
create_emuRtrackdata <- function(sl, td){
  
  ########################
  # check parameters
  # check correct classes
  if(!inherits(sl, "emuRsegs") || !inherits(td, "trackdata")){
    stop("emuRtrackdata could not be created: sl is not of class 'emuRsegs' or td arguments is not of class 'trackdata'")
  }
  
  # check same number of items
  if(dim(td$index)[1] != nrow(sl)){
    stop("emuRtrackdata could not be created: td and sl objects don't have the same number of elements (dim(td$index)[1] != nrow(sl))")
  }
  
  
  nframes = 1 + apply(td$index, 1, diff)
  inds = rep(1:nrow(td), nframes)
  # expand seglist
  expSl = sl[inds,]
  
  times = tracktimes(td)
  start.time = rep(start(td), nframes)
  n.time = times - start.time
  rownames(td$data) = NULL
  res = data.table(sl_rowIdx = inds, expSl, times_rel = n.time, times_orig = times, td$data)
  class(res) <- c("emuRtrackdata", class(res))
  return(res)
}


# ##' Function to extract an emuRtrackdata at a single time point of to create another EMU-trackdata object between two times
# ##' 
# ##' this function is intended to mimic some of the dcut behaviour
# ##' 
# ##' @param emuRtrackdata emuRtrackdata object
# ##' @param leftTime bla
# ##' @param rightTime bli
# ##' @param single blup
# ##' @param average blap
# ##' @param prop woooord
# cut.emuRtrackdata = function(emuRtrackdata, leftTime, rightTime, 
#                              single = TRUE, average = TRUE, prop = FALSE){
#   
#   if (missing(rightTime)) {
#     if(length(leftTime == 1)){
#       res = emuRtrackdata %>% 
#         dplyr::group_by("sl_rowIdx") %>% 
#         dplyr::mutate(times_propDiff = 1 + ("times_orig" - min("times_orig")) / (max("times_orig") - min("times_orig")) - leftTime) %>%
#         dplyr::filter(times_propDiff == min("times_propDiff")) %>%
#         dplyr::select(-times_propDiff)
#     }else{
#       
#     }
#   }
#   # return(as.data.table(res))
# }


#######################
# FOR DEVELOPMENT
# library('testthat')
# test_file('tests/testthat/test_emuRtrackdata.R')

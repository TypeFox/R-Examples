#'Generate a gapdata class object from a dendrogram object
#'
#'This function takes a dendrogram class object as an input, and generate a gapdata class object as an output.
#'By parsing the dendrogram object based on parameters for gaps, gaps between leaves in a dendrogram are introduced,
#'and the coordinates of the leaves are adjusted. The gaps can be based on the a height (or distance) threshold to
#'to introduce the gaps of the same width, or quantitative mapping of distance values mapped linearly or exponentially.
#'
#'
#' @param d dendrogram class object
#' @param mode gap mode, either "threshold" or "quantitative"
#' @param mapping in case of quantitative mode, either "linear" or "exponential" mapping
#' @param ratio the percentage of width allocated for the sum of gaps.
#' @param scale the sclae log base for the exponential mapping
#' @param threshold the height at which the dendrogram is cult to infer clusters
#' @param verbose logical for whether in verbose mode or not
#' @param ... ignored
#' @export gap_data
#' @aliases gap_data
#' @return a list of data frames that contain coordinates for drawing a gapped dendrogram
#'
#'

gap_data <- function(d, mode=c("quantitative", "threshold"), mapping=c("exponential", "linear"), ratio= 0.2, scale = 0.5, threshold = 0, verbose=FALSE,  ...){
  #arguments
  mode <- match.arg(mode)
  #number of nodes
  N = attr(d, "members")
  if(verbose) print(paste0("total number of nodes = ", N))
  #allocate gap space (default 20%)  
  gap_total = ratio*N
  if(verbose)print(paste0("total length of gap = ", gap_total))
  
  #annotate gap to dendrogram
  if(mode == "quantitative"){
    mapping <- match.arg(mapping)
    max_height = attr(d, "height")
    #calculate the sum of distance
    sum = sum_distance(d, sum =0, mapping = mapping, scale=scale, max_height=max_height)
    #recursively calculate gap
    if(verbose)print("calculate_gap() -----------")
    d = calculate_gap(d=d, sum=sum, gap_total=gap_total, mode= mode, mapping=mapping, scale=scale,max_height= max_height, verbose=verbose)
  }else if(mode == "threshold"){
    #count total number of branches above the threshold
    count = count_gap(d=d, count=0, threshold=threshold)
    #calculate gap_size
    gap_size = gap_total/count
    #recursively calculate gap
    if(verbose)print("calculate_gap() -----------")
    d = calculate_gap(d=d, mode= mode, gap_size=gap_size, threshold=threshold, verbose=verbose)
  }
  
  #re-evaluate the positions for each leaves
  if(verbose)print("assign_positions() -----------")
  d = assign_positions(d, runningX = 1, verbose=verbose)
  
  #re-evaluate the branch positions
  d = assign_branch_positions(d)
  
  #extract a list
  l = extract_list(d, type="rectangle")
  #add column names
  names(l$segments) <- c("x", "y", "xend", "yend")
  names(l$labels) <- c("x", "y", "label")

  #compose a gapdata class object
  output <- as.gapdata(d = d, segments = l$segments, labels = l$labels)
  
  output #return
}





#'Calculate the gaps based on distance
#'
#'This function takes a dendrogram class object and other attributes to calcuate the 
#'size of gaps between leaves. The gap is stored in a  leaf to its left. The function
#'is called recursively.
#'
#'
#' @param d dendrogram class object
#' @param sum the sum of distance
#' @param gap_total the total width allocated for gaps
#' @param mode gap mode, either "threshold" or "quantitative"
#' @param mapping in case of quantitative mode, either "linear" or "exponential" mapping
#' @param scale the sclae log base for the exponential mapping
#' @param max_height the highest distance value, which is the value of the first dendrogram branch
#' @param threshold the threshold value for threshold mode
#' @param gap_size the size of gap for threshold mode
#' @param verbose logical for whether in verbose mode or not
#' @param ... ignored
#' @export calculate_gap
#' @aliases calculate_gap
#' @return the annotated dendrogram class object
#' @keywords internal
#' 

calculate_gap <- function(d, sum, gap_total, mode=c("quantitative", "threshold"), mapping=c("exponential", "linear"), scale = 0.2, max_height=0,  threshold=2, gap_size=0, verbose=FALSE, ...){
  a = attributes(d) #attributes
  left = d[[1]]
  right = d[[2]]
  #distance of this branch
  h = a$height
  # map this distance to allocated space
  mode <- match.arg(mode)
  gap = 0;
  if(mode == "quantitative"){
    mapping <- match.arg(mapping)
    if(mapping == "linear"){
      # linear mapping
      gap = map(h, 0, sum, 0, gap_total)
    }else if(mapping == "exponential"){
      #exponential mapping      
      h.exp = map.exp(h, 0, max_height, 0, 1, scale = scale)
      gap = map(h.exp, 0, sum, 0, gap_total)
    }
  }else if(mode =="threshold"){
    if(h > threshold){
      gap = gap_size
    }
  }
  
  #print(paste0("verbose: h=", format_number(h), " gap_total=", format_number(gap_total), " sum=", format_number(sum), " gap=", gap))
  
  if(is.leaf(left) && is.leaf(right)){
    #save the gap on the left
    attr(left, "right_gap")=gap  
    if(verbose) print(paste0("gap of ", attr(left, "label"), " is ", format_number(gap)))
  }else if(!is.leaf(left) && is.leaf(right)){
    #find the most right leaf
    most_right = get_most_right_leaf(left)
    #assign position for right
    attr(most_right, "right_gap") = gap
    #assign the value
    if(verbose) print(paste0("gap of ", attr(most_right, "label"), " is ", format_number(gap), " sub_leaf"))
    left = set_most_right_leaf(left, most_right)   
    #go deeper on left next
    left = calculate_gap(d=left, sum=sum, gap_total=gap_total, mode= mode, mapping=mapping, scale=scale,max_height= max_height, threshold=threshold, gap_size=gap_size, verbose=verbose)
  }else if(is.leaf(left) && !is.leaf(right)){
    #save the gap on the left
    attr(left, "right_gap")=gap  
    if(verbose) print(paste0("gap of ",attr(left, "label"), " is ", format_number(gap)))
    #go into the next layer
    right = calculate_gap(d=right, sum=sum, gap_total=gap_total, mode= mode, mapping=mapping, scale=scale,max_height= max_height, threshold=threshold, gap_size=gap_size, verbose=verbose)
  }else{
    #find the most right leaf
    most_right = get_most_right_leaf(left)
    attr(most_right, "right_gap") = gap
    if(verbose) print(paste0("gap of ", attr(most_right, "label"), " is ", format_number(gap), " sub_sub"))
    #assign the value
    left = set_most_right_leaf(left, most_right)    
    #go into the next layer
    left = calculate_gap(d=left, sum=sum, gap_total=gap_total, mode= mode, mapping=mapping, scale=scale,max_height= max_height, threshold=threshold, gap_size=gap_size, verbose=verbose)
    right = calculate_gap(d=right, sum=sum, gap_total=gap_total, mode= mode, mapping=mapping, scale=scale,max_height= max_height, threshold=threshold, gap_size=gap_size, verbose=verbose)
  } 
  d[[1]] = left
  d[[2]] = right
  d
}
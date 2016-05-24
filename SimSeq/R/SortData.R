SortData <-
function(counts, treatment, replic = NULL, sort.method, norm.factors = NULL){
  #Trims and sorts matrix of counts into appropriate format to be used
  #in data.sim
  
  #Check inputs
  counts <- as.matrix(counts)
  n.col <- dim(counts)[2]
  
  if (!is.null(replic)){
    if(length(replic) != n.col) 
      stop("Error: Length of replic vector must equal number of 
           columns in counts matrix")
    if (any(tabulate(replic) > 2)) 
      stop("Error: Number of observations per replicate in counts 
           matrix must not be greater than 2.")
  }
  
  if(length(treatment) != n.col)
    stop("Error: Length of treatment vector must equal number of 
           columns in counts matrix")
  if (length(unique(treatment)) != 2) 
    stop("Error: Number of treatment groups in counts matrix 
         must be equal to two.")
  if (!sort.method %in% c("paired", "unpaired")) 
    stop("Error: sort.method must be set to either 'paired' or 'unpaired'.")
  if (sort.method == "paired" && is.null(replic)) 
    stop("Error: Must specify replic vector when sort.method equals 'paired'.")
  
  if(!is.null(norm.factors)){
    if (!is.numeric(norm.factors)) 
      stop("Error: norm.factors must be a positive numeric vector with 
           length equal to the number of columns in the counts matrix.")
    if (is.numeric(norm.factors)){
      if (any(norm.factors <= 0) || length(norm.factors) != n.col) 
        stop("Error: norm.factors must be a positive numeric vector with 
             length equal to the number of columns in the counts matrix.")
    }
  }

  
  if(sort.method == "paired"){ 
    #sort inputs by replic then treatment
    sorting <- order(replic, treatment, decreasing = FALSE)
    replic <- replic[sorting]
    treatment <- treatment[sorting]
    counts <- counts[, sorting, drop = FALSE]
    norm.factors <- norm.factors[sorting]
    
    #remove replicates in inputs that are non-paired
    counts <- counts[, replic %in% unique(replic)[tabulate(replic) == 2], drop = FALSE]
    norm.factors <- norm.factors[replic %in% unique(replic)[tabulate(replic) == 2]]
    treatment <- factor(treatment[replic %in% unique(replic)[tabulate(replic) == 2]])
    sorting <- sorting[replic %in% unique(replic)[tabulate(replic) == 2]]
    replic <- factor(replic[replic %in% unique(replic)[tabulate(replic) == 2]])
  }
  
  if (sort.method == "unpaired"){
    #sort inputs by treatment
    sorting <- order(treatment)
    replic <- replic[sorting]
    treatment <- treatment[sorting]
    counts <- counts[, sorting, drop = FALSE]
    norm.factors <- norm.factors[sorting]
    
    #no need to remove any columns (replicates)
  }
  
  return(list(counts = counts, replic = replic, treatment = treatment, 
              norm.factors = norm.factors, sorting = sorting))
  }
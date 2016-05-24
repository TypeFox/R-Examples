#DAGstructure object contains all relevant information of a certain DAG structure
setClass("DAGstructure",
         representation(
           parents = "list",
           children = "list",
           sets = "list",
           twoway = "logical"
         )
)

#function that construct a DAG structure from a given list of sets. Output is a DAGstructure object. 
construct <- function(sets) {
  if(length(sets) == 0)
  {
    stop("sets is empty!")
  }
  
  sets <- removeMultiples(sets)
  p_and_c <- make_p_and_c(sets)
  twoway <- check(p_and_c$parents, p_and_c$children, sets)
  new("DAGstructure", parents = p_and_c$parents,
      children = p_and_c$children,
      sets = sets,
      twoway = twoway)
}

#function that removes duplicate sets
removeMultiples <- function(sets)
{
  unique_sets <- list()
  for(i in 1:length(sets))
  {
    ispresent <- FALSE
    
    if(length(unique_sets) > 0)
      for(j in 1:length(unique_sets))
        #if(all(sets[[i]] %in% unique_sets[[j]]) && all(unique_sets[[j]] %in% sets[[i]]))
        if(length(sets[[i]]) == length(unique_sets[[j]]) && all(sets[[i]] %in% unique_sets[[j]]) && all(unique_sets[[j]] %in% sets[[i]]))
        {
          ispresent <- TRUE
          break
        }
    
    if(!ispresent)
    {
      unique_sets[[length(unique_sets) + 1]] <- sets[[i]]
      names(unique_sets)[length(unique_sets)] <- names(sets)[i]
    }
  }
  if(length(sets) != length(unique_sets))
    warning(paste("Given sets contained", length(sets) - length(unique_sets), "duplicate(s)."))
  
  return(unique_sets)
}

#given sets: makes list of parents and children, for each set, the indices of the sets that are its parents or children are stored
make_p_and_c <- function(sets)
{
  
  # order all sets from biggest to smallest 
  n <- length(sets)
  lengths <- rep(0,n)
  for(i in 1:n)
    lengths[i] <- length(sets[[i]])
  perm <- order(lengths, decreasing = TRUE)
  
  # makes empty lists to store parents and children
  parents <- vector("list", n)
  children <- vector("list", n)
  
  for(i in perm)
    for(j in perm)
      if (i != j && length(sets[[i]]) <= length(sets[[j]]) && all(sets[[i]] %in% sets[[j]]))
        #if(i != j && all(sets[[i]] %in% sets[[j]]))
      {
        # test whether i is not a child of a child of j
        ischild <- TRUE
        
        for(k in children[[j]])
        {
          if(all(sets[[i]] %in% sets[[k]]))
          {
            ischild <- FALSE
            break
          }
        }
        
        if(ischild)
        {
          children[[j]] <- c(children[[j]], i)
          parents[[i]] <- c(parents[[i]], j)
        }
      }
  return(list(children=children,parents=parents))
  
}

#check whether each node is the intersection of its child noded (if so: twoway property holds)
check <- function(parents, children, sets)
{
  for(i in 1:length(parents))
  {
    unionchild <- NULL
    for(k in children[[i]])
    {
      unionchild <- union(unionchild,sets[[k]])
    }
    
    #check only for non leaf nodes
    if(length(children[[i]])>0)
    {
      if(!(all(sets[[i]] %in% unionchild)))
        return(FALSE)    
    }
  }
  return(TRUE)
}


setGeneric("istwoway", function(object) standardGeneric("istwoway"))
setMethod("istwoway", "DAGstructure", function(object) object@twoway)

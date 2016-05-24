#======================================================================
# List utilities
#======================================================================
# combine the two lists, replacing duplicates in list1 with values from list2
merge.lists <- function (list1, list2)
{
  if(length(list2)==0){return(list1)}
  if(length(list1)==0){return(list2)}
  allNames <- unique(c(names(list1), names(list2)))
  merged <- list1 # we will copy over/replace values from list2 as necessary
  for (x in allNames)
  {
    a <- NULL ; b <- NULL
    if(x %in% names(list1)){ a <- list1[[x]]}
    if(x %in% names(list2)){ b <- list2[[x]]}
    if (is.null(a)) {
      # only exists in list2, copy over
      merged[[x]] <- b
    } else if (is.list(a) && is.list(b)) {
      # recurse
      merged[[x]] <- merge.lists(a, b)
    } else if (!is.null(b)) {
      # replace the list1 value with the list2 value (if it exists)
      merged[[x]] <- b
    }
  }
  return(merged)
}

# transform an S4 class instance to a list
S4.as.list <- function(x)
{
  mapply(function(y) slot(x,y),
         slotNames(class(x)[1]),
         SIMPLIFY=FALSE)
}

# as for lapply but apply names and remove any NULLs before return
lapply_ex <- function(x, fun, names=NULL, merge=F,...)
{
  nms <- names(x)
  res <- lapply(x, fun, ...)
  if(!is.null(names))
  {
    names(res)<- names
  }
  else
  {
    names(res)<- nms
  }
  res <- delete.NULLs(res)
  if(merge)
  {
    res <- do.call(rbind,res)
  }
  return(res)
}


# reduce dimensions of array by one by coercing top level index elements to a list
multi.array.to.list <- function(x, nms=NULL, index=3)
{
  # get a list from elements in (third) array index
  ls <- apply(x, MARGIN=index, function(x){list(x)})
  ls <- lapply(ls,
               function(x)
               {
                 x <- x[[1]]
                 x[!is.finite(x)]<-0
                 return(x)
               }
  )

  names(ls) <- nms
  return(ls)
}

# replicate an object n times and return as a list
lrep <- function(x, n)
{
  lapply(1:n, function(i){return(x)})
}

# remove null values from a list
delete.NULLs  <-  function (x.list)
{
  x.list[unlist(lapply(x.list, length) != 0)]
}



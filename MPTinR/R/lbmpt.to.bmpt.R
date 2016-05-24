### lbmpt.to.mpt ###

# .isNode
.isNode <- function(op) grepl("[[:alpha:]]",op)

# .isLeaf
.isLeaf <- function(op){
  v <- grepl("^[[:digit:]]+$",op)
  if(sum(!v)>0){
    out <- FALSE
  }else{
    out <- TRUE
  }
  if(is.list(op)){
    out <- FALSE
  }
  return(out)
}

# .parse.rec
.parse.rec <-  function(cfl){
 if(.isNode(cfl[1]) & .isLeaf(cfl[2]) & .isLeaf(cfl[3])){
    if(is.na(cfl[4])){
      l <- list(list(c(cfl[1], cfl[2]),(c(paste("(1-",cfl[1],")", sep=""), cfl[3]))),"")# here is something wrong!!
      return(l)
    }else{
    l <- list(list(c(cfl[1], cfl[2]),(c(paste("(1-",cfl[1],")", sep=""), cfl[3]))),cfl[4:length(cfl)])# here is something wrong!!
    return(l)
  }}
  if(.isNode(cfl[1]) & .isLeaf(cfl[2])){
    secondcall <- .parse.rec(cfl[3:length(cfl)])
    l <- list(list(c(cfl[1], cfl[2]), c(paste("(1-",cfl[1],")", sep=""), secondcall[1])),secondcall[[2]])
    return(l)
  }
  if(.isNode(cfl[1])){
    firstcall  <-  .parse.rec(cfl[2:length(cfl)])
    if(.isLeaf(firstcall[[2]])&length(firstcall[[2]])==1){
      l <- list(list(list(cfl[1], firstcall[[1]]), c(paste("(1-",cfl[1],")", sep=""), firstcall[[2]])),cfl[2:length(cfl)]) 
      return(l)
     } else{
       if(.isLeaf(firstcall[[2]][1])){
        secondcall <- .parse.rec(firstcall[[2]][-1])
        l <- list(list(list(cfl[1], firstcall[[1]]), c(paste("(1-",cfl[1],")", sep=""), firstcall[[2]][1])), firstcall[[2]][-1])
        return(l)
       } else{
      secondcall <- .parse.rec(firstcall[[2]])
      l <- list(list(list(cfl[1], firstcall[[1]]), c(paste("(1-",cfl[1],")", sep=""), secondcall[1])),secondcall[[2]])
      return(l)
    }
    }
  }
  }

# .LinearizeNestedList
.LinearizeNestedList <- function(NList, LinearizeDataFrames=FALSE,
                                NameSep="/", ForceNames=FALSE) {
  # LinearizeNestedList:
  #
  # https://sites.google.com/site/akhilsbehl/geekspace/
  # articles/r/linearize_nested_lists_in_r
  #
  # Akhil S Bhel
  #
  # Implements a recursive algorithm to linearize nested lists upto any
  # arbitrary level of nesting (limited by R's allowance for recursion-depth).
  # By linearization, it is meant to bring all list branches emanating from
  # any nth-nested trunk upto the top-level trunk s.t. the return value is a
  # simple non-nested list having all branches emanating from this top-level
  # branch.
  #
  # Since dataframes are essentially lists a boolean option is provided to
  # switch on/off the linearization of dataframes. This has been found
  # desirable in the author's experience.
  #
  # Also, one'd typically want to preserve names in the lists in a way as to
  # clearly denote the association of any list element to it's nth-level
  # history. As such we provide a clean and simple method of preserving names
  # information of list elements. The names at any level of nesting are
  # appended to the names of all preceding trunks using the `NameSep` option
  # string as the seperator. The default `/` has been chosen to mimic the unix
  # tradition of filesystem hierarchies. The default behavior works with
  # existing names at any n-th level trunk, if found; otherwise, coerces simple
  # numeric names corresponding to the position of a list element on the
  # nth-trunk. Note, however, that this naming pattern does not ensure unique
  # names for all elements in the resulting list. If the nested lists had
  # non-unique names in a trunk the same would be reflected in the final list.
  # Also, note that the function does not at all handle cases where `some`
  # names are missing and some are not.
  #
  # Clearly, preserving the n-level hierarchy of branches in the element names
  # may lead to names that are too long. Often, only the depth of a list
  # element may only be important. To deal with this possibility a boolean
  # option called `ForceNames` has been provided. ForceNames shall drop all
  # original names in the lists and coerce simple numeric names which simply
  # indicate the position of an element at the nth-level trunk as well as all
  # preceding trunk numbers.
  #
  # Returns:
  # LinearList: Named list.
  #
  # Sanity checks:
  #
  stopifnot(is.character(NameSep), length(NameSep) == 1)
  stopifnot(is.logical(LinearizeDataFrames), length(LinearizeDataFrames) == 1)
  stopifnot(is.logical(ForceNames), length(ForceNames) == 1)
  if (! is.list(NList)) return(NList)
  #
  # If no names on the top-level list coerce names. Recursion shall handle
  # naming at all levels.
  #
  if (is.null(names(NList)) | ForceNames == TRUE)
    names(NList) <- as.character(1:length(NList))
  #
  # If simply a dataframe deal promptly.
  #
  if (is.data.frame(NList) & LinearizeDataFrames == FALSE)
    return(NList)
  if (is.data.frame(NList) & LinearizeDataFrames == TRUE)
    return(as.list(NList))
  #
  # Book-keeping code to employ a while loop.
  #
  A <- 1
  B <- length(NList)
  #
  # We use a while loop to deal with the fact that the length of the nested
  # list grows dynamically in the process of linearization.
  #
  while (A <= B) {
    Element <- NList[[A]]
    EName <- names(NList)[A]
    if (is.list(Element)) {
      #
      # Before and After to keep track of the status of the top-level trunk
      # below and above the current element.
      #
      if (A == 1) {
        Before <- NULL
      } else {
        Before <- NList[1:(A - 1)]
      }
      if (A == B) {
        After <- NULL
      } else {
        After <- NList[(A + 1):B]
      }
      #
      # Treat dataframes specially.
      #
      if (is.data.frame(Element)) {
        if (LinearizeDataFrames == TRUE) {
          #
          # `Jump` takes care of how much the list shall grow in this step.
          #
          Jump <- length(Element)
          NList[[A]] <- NULL
          #
          # Generate or coerce names as need be.
          #
          if (is.null(names(Element)) | ForceNames == TRUE)
            names(Element) <- as.character(1:length(Element))
          #
          # Just throw back as list since dataframes have no nesting.
          #
          Element <- as.list(Element)
          #
          # Update names
          #
          names(Element) <- paste(EName, names(Element), sep=NameSep)
          #
          # Plug the branch back into the top-level trunk.
          #
          NList <- c(Before, Element, After)
        }
        Jump <- 1
      } else {
        NList[[A]] <- NULL
        #
        # Go recursive! :)
        #
        if (is.null(names(Element)) | ForceNames == TRUE)
          names(Element) <- as.character(1:length(Element))
        Element <- .LinearizeNestedList(Element, LinearizeDataFrames,
                                       NameSep, ForceNames)
        names(Element) <- paste(EName, names(Element), sep=NameSep)
        Jump <- length(Element)
        NList <- c(Before, Element, After)
      }
    } else {
      Jump <- 1
    }
    #
    # Update book-keeping variables.
    #
    A <- A + Jump
    B <- length(NList)
  }
  return(NList)
}

# .append.helper
.append.helper <- function(ap, dct){
 
  
  for(key in unique(names(dct))){
    if(!is.list(dct[[key]])){
      for(q in which(names(dct)==key)){
      dct[[q]] <- c(dct[[q]], ap)
    }} else{
      dct <- .LinearizeNestedList(dct)
      names(dct)[1:length(names(dct))] <- substr(names(dct)[1], 1, 1)
    dct <- lapply(dct, function(x) c(x, ap))
    }
  }
  return(dct)
}


# .rendEq.rec
.rendEq.rec <- function(parsed){
  left <- parsed[[1]]
  right <- parsed[[2]]
  lp <- list()
  rp <- list()
  if(.isLeaf(left[2])){
    lp[[left[2]]] <- left[[1]]
  } else{
    lp <- .append.helper(left[[1]], .rendEq.rec(left[[2]]))
  }
  if(.isLeaf(right[2])){
    rp[[right[2]]] <- right[[1]]
  } else{
    rp <- .append.helper(right[[1]], .rendEq.rec(right[[2]]))
  }

  return(c(lp,rp))
}
  


# .parse
.parse <- function(cfl) .parse.rec(cfl)[[1]]



.renderEquation <- function(parsed, category.names= TRUE){
  trees <- .rendEq.rec(parsed)
m.vec <- vector("character", length(unique(names(trees))))
for(key in unique(names(trees))){
f.category <- vector("character", length(which(names(trees)== key)))
z <- 1
for(i in which(names(trees)== key)){
text <- vector("character", length(seq_along(trees[[i]])))
if(length(trees[[i]])==1){
  f.category[z] <- trees[[i]]
  
}else{
  trees[[i]] <- trees[[i]][length(trees[[i]]):1]
for(y in seq_along(trees[[i]])[-length(trees[[i]])]){
  text[y] <- paste(trees[[i]][y],"*", sep="")
   }
text[length(trees[[i]])] <- trees[[i]][length(trees[[i]])]
  

f.category[z] <- paste(text, collapse="")
}
z <- z+1
}
if(length(f.category)==1){
  m.vec[as.numeric(key)] <- f.category
  if(category.names){
  m.vec[as.numeric(key)] <- paste(m.vec[as.numeric(key)], "# category", key)
  } 
}else{
b.category <- vector("character", length(f.category))
for(z in seq_along(f.category)[-length(f.category)]){
  b.category[z] <- paste(f.category[[z]],"+", sep= "")
  b.category[length(f.category)] <- f.category[[length(f.category)]]
}
m.vec[as.numeric(key)] <- paste(b.category, collapse="")
if(category.names){
m.vec[as.numeric(key)] <- paste(m.vec[as.numeric(key)], "# category", key)
}
}
}
return(m.vec)
}



# lbmpt.to.mpt
lbmpt.to.mpt <- function(model.list, outfile = NULL, category.names = TRUE){
    vec <- c()
    y <- 1
    for(l in seq_along(model.list)){
      parsed <- .parse(model.list[[l]])
      raw.model <- .renderEquation(parsed, category.names)
        
      for(q in seq_along(raw.model)){
        vec[y] <- raw.model[q]
        y <- y+1
      }
      
      vec[y] <- ""
      y <- y+1
    }
    if (is.null(outfile)) return(writeLines(vec))
    else writeLines(vec, con = outfile)
  }


# test
# lbmpt.to.mpt(list(c("a", "1", "b", "2", "3"),c("x", "b", "1", "2", "c", "x", "1", "3", "e", "4", "3"),c("a", "b","1", "2", "c","3", "4") ))
# 
# lbmpt.to.mpt(list(c("a", "1", "b", "2", "3"),c("x", "b", "1", "2", "c", "x", "1", "3", "e", "4", "3"),c("a", "b","1", "2", "c","3", "4") ), outfile="test1.txt")





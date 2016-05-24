# This code plots Venn Diagrams for up to 5 sets. The
# function getVennCounts is passed a list of vectors.
# This is transformed into a table indicating the
# number of intersections for each intersection. This table
# is generated for any number of sets.

# The function drawVennDiagram plots circles (up to three
# sets) or ellipses (4 and 5 sets) to depict the sets.
# The sum of values placed is the number of entries of
# each set.

# Function to determine values of a venn diagram
# It works for an arbitrary large set of input sets.
#

getVennCounts <- function(l, universe, verbose=F, ...)
  UseMethod("getVennCounts")

getVennCounts.data.frame <- function(l, universe=NA, verbose=F, ...)
  {
    if (verbose) cat("Interpreting data as data.frame.\n")
    if( !all(unique(unlist(l)) %in% c(0,1))  )
      stop("Only indicator columns permitted")

    l <- lapply( l, function(x) which(as.logical(x)))
    getVennCounts.list(l, universe=universe, verbose=verbose)
  }

getVennCounts.matrix <- function(l, universe=NA, verbose=F, ...)
{
  getVennCounts.data.frame(as.data.frame(l), universe=NA, verbose=F, ...)
}

# l offers a list of arrays, their values are to
# be tested for the size of their intersects.
getVennCounts.list<-function(l, universe=NA, verbose=F, intersections=TRUE) {
  if (verbose) cat("Interpreting data as list.\n")
  numSets<-length(l)
  result.table<-NULL
  result.table.names<-NULL

  memberList <- list()

  # Iteration over all possible intersections involving all sets
  # or the complement (negation) of those sets.
  for (i in 0:(-1 + 2^numSets)) {
    # i2 is a binary representation of that number
    i2<-baseOf(i,2,numSets)

    # some debug output
    #print(paste(i,":",paste(i2,collapse="",sep="")))

    # p.pos determines the position in number
    #       which is also the set that is inspected

    sel<-universe

    # positive selection first
    for (p.pos in which(1 == i2) ) {
      current.set<-l[[p.pos]]
      if (!is.null(dim(current.set))) {
        # circumventing strange experiences with data.frames
        warning(paste("List element [[",p.pos,"]] has dimensions, but all elements are considered.\n",sep=""))
        current.set<-as.character(as.matrix(current.set))
        dim(current.set)<-NULL
      }
      #print(paste("set ",p.pos,", val=1: ",paste(current.set,collapse=",")))
      if (is.null(sel)) {
        #print("Sel is null")
      } else if (1 == length(sel) && is.na(sel)) {
        sel<-current.set
      }
      else {
        w<-which(sel %in% current.set)
        if (length(w)>0) {
          sel<-sel[w]
        }
        else {
          sel<-NULL
        }
      }
    }

    # something should be in sel now, otherwise
    # the number will be 0

    # negative selection
    for (p.pos in which(0 == i2) ) {
      if (is.null(sel) || ( 1 == length(sel) && is.na(sel))) {
        # The complement is not known, hence no checks done
      }
      else {
        current.set<-l[[p.pos]]
        if (!is.null(dim(current.set))) {
          warning(paste("List element [[",p.pos,"]] has dimensions, but all elements are considered.\n",sep=""))
          current.set<-as.character(as.matrix(current.set))
          dim(current.set)<-NULL
        }
        w<-which( ! sel %in% current.set)
        #print(paste("set ",p.pos,", val=1: ",paste(current.set,collapse=",")))
        if (length(w)>0) {
          sel<-sel[w]
        }
        else {
          sel<-NULL
        }
      }
    }
    #print(paste("sel:",paste(sel,collapse=",")))

    if(is.null(sel) || (1 == length(sel) && is.na(sel))) {
      sel<-NULL
    }

    r.name<-paste(i2,collapse="")
    if (intersections) {
      memberList[[r.name]] <- sel
    }

    r<-length(sel)
    result.row<-c(r,i2)
    dim(result.row)<-c(1,length(result.row))
    rownames(result.row)<-c(r.name)
    #print(paste("Adding ",r.name))
    if (is.null(result.table)) {
      result.table<-result.row
    }
    else {
      result.table<-rbind(result.table,result.row)
    }
    #if (is.null(result.table)) {
    #	result.table<-r
    #	result.table.names<-r.name
    #}
    #else {
    #	result.table<-c(result.table,r)
    #	result.table.names<-c(result.table.names,r.name)
    #}
  }
  #names(result.table)<-result.table.names
  if (is.null(names(l))) {
    colnames(result.table)<-c("num",LETTERS[1:numSets])
  }
  else{
    colnames(result.table)<-c("num",names(l))
  }
  if (intersections) {
    attr(result.table,"intersections") <- memberList
  }
  class(result.table) <- "venn"
  return(result.table)
}

#print(getVennCounts(list(A,B,C,D)))
#print(getVennCounts(list(a=A,b=B,c=C,d=D)))

venn <- function(data,
                 universe=NA,
                 small=0.7,
                 showSetLogicLabel=FALSE,
                 simplify=FALSE,
                 show.plot=TRUE,
                 intersections=TRUE,
                 names,
                 ...
                 )
{
  counts <- getVennCounts(data,
                          universe=universe,
                          intersections=intersections
                          )

  if(show.plot)
    drawVennDiagram(data=counts,
                    small=small,
                    showSetLogicLabel=showSetLogicLabel,
                    simplify=simplify,
                    ...
                    )

  # use VennMemberNames to properly label and order the 'intersection' table
  if(intersections)
    attr(counts, "intersections") <- vennMembers(l=data,
                                                 universe=universe,
                                                 names=names
                                                 )

  invisible(counts)
}


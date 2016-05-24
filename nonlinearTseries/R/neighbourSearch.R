# neighbour Search in phase space
# Reference: Efficient neighbor searching in nonlinear time series analysis, schreiber

#wrapper for the boxAssistant C function.
#This function computes the accumulated histogram of the "takens" takens'
#vectors using a wrapped grid of "number.boxes" per dimension. Each box
#has a size of "radius".
#Further details may be seen in the C function comments
boxAssistant=function(takens, radius,number.boxes=NULL){
  #estimates number.boxes if it has not been specified
  if (is.null(number.boxes)) number.boxes = estimateNumberBoxes(takens,radius)
  #parameters for the C function
  numberTakens=nrow(takens)
  embedding.dim=ncol(takens)
  boxes=rep(0,number.boxes*number.boxes+1)
  possibleNeighbours=rep(0,numberTakens)
  #call C
  ret=.C("boxAssistant",takens=as.double(takens),numberTakens=as.integer(numberTakens),
                    embeddingD=as.integer(embedding.dim),eps= as.double(radius),
                    numberBoxes=as.integer(number.boxes),boxes=as.integer(boxes),
                    possibleNeighbours=as.integer(possibleNeighbours),
                    PACKAGE="nonlinearTseries" )
  #format solution and return it
  sol=list(boxes=ret$boxes,possibleNeighbours=ret$possibleNeighbours)
  sol = propagateTakensAttr(sol, takens)
  
  return(sol)
}

################################################################################
#' neighbour search
#' @description
#' This function finds all the neighbours of a given Takens' vector. The neighbours 
#' are found using a box assisted algorithm that creates a wrapped grid with a given
#' number of  boxes per dimension. 
#' @param takens The matrix containing all the Takens' vectors (see \link{buildTakens}).
#' @param positionTakens Integer denoting the Takens' vector whose neighbours will be searched.
#' @param radius Distance in which the algorithm will search for neighbours.
#' @param number.boxes Integer denoting the number of boxes per dimension that will be used
#' to construct a wrapped grid (see Schreiber). If the user does not specify a number of boxes,
#' this function estimates a proper number.
#' @return A containing all the neighbours of the 
#' \emph{positionTakens-th} Takens' vector. If the list is empty, that means that there is
#' no neighbour of the \emph{positionTakens-th} Takens' vector in the given radius.
#' @seealso \code{\link{findAllNeighbours}}.
#' @references  Schreiber, T. Efficient neighbor searching in nonlinear time series analysis. Int. J. Bifurcation and Chaos, 
#' 5, p. 349, (1995).
#' @author Constantino A. Garcia
#' @export neighbourSearch
#' @useDynLib nonlinearTseries
neighbourSearch=function(takens,positionTakens,radius,number.boxes=NULL){
  #estimates number.boxes if it has not been specified
  if (is.null(number.boxes)) number.boxes = estimateNumberBoxes(takens, radius)
  #arguments for the C function
  #position in the vector was specified in "R-form" (starts at 1)
  #Translate to "C-form" (starts at 0)
  positionTakens=positionTakens-1
  numberTakens=nrow(takens)
  embedding.dim=ncol(takens)
  lengthBoxes=number.boxes*number.boxes+1
  boxes=rep(0,lengthBoxes)
  possibleNeighbours=rep(0,numberTakens)
  neighList=rep(-1,numberTakens)
  nfound=0
 
  #get histogram using boxes
  sol=boxAssistant(takens, radius,number.boxes)
  boxes=sol$boxes;
  possibleNeighbours=sol$possibleNeighbours
 # call c code
  cneighs=.C("neighbourSearchFromBoxes",takens=as.double(takens),positionTakens=as.integer(positionTakens),
                          numberTakens=as.integer(numberTakens), embeddingD=as.integer(embedding.dim),
                          eps=as.double(radius), numberBoxes=as.integer(number.boxes),
                          boxes=as.integer(boxes),possibleNeighbours=as.integer(possibleNeighbours),
                          neighList=as.integer(neighList),
                          nfound=as.integer(nfound),
                          PACKAGE="nonlinearTseries")
  #arrange solutions
  if (cneighs$nfound==0){
    finalNeighs=list(nfound=0,neighList=list())
  }else{
    #eliminate -1 and repeated numbers. convert c-vector positions to R positions
    neighList=cneighs$neighList[1:(cneighs$nfound)]+ 1
    finalNeighs=list(nfound=cneighs$nfound,neighList=neighList)
  }
  finalNeighs = propagateTakensAttr(finalNeighs, takens)
  # remember to translate the C-index into R-index
  attr(finalNeighs,"takens.index") = positionTakens + 1
  attr(finalNeighs,"radius") = radius
  return(finalNeighs)
}

################################################################################
#' neighbour search
#' @description
#' This function finds all the neighbours of all the vectors from Takens'
#' vector array. The neighbours are found using a box assisted algorithm that
#' creates a wrapped grid of a given number of boxes per dimension. 
#' @param takens The matrix containing all the Takens' vectors (see \link{buildTakens}).
#' @param radius Distance in which the algorithm will search for neighbours.
#' @param number.boxes Integer denoting the number of boxes per dimension that will be used
#' to construct a wrapped grid (see Schreiber). If the user does not specify a number of boxes,
#' this function estimates a proper number.
#' @return A list in which the n-th position contains another list with all the neighbours of the 
#' n-th Takens' vector. If the list is empty, that means that there is no neighbour of the n-th Takens'
#' vector in the given radius.
#' @references Schreiber, T. Efficient neighbor searching in nonlinear time series analysis. Int. J. Bifurcation and Chaos, 
#' 5, p. 349, (1995).
#' @author Constantino A. Garcia
#' @examples 
#' \dontrun{
#' # Find all the neighbours Takens' vectors build from the Henon time
#' # series. The size of the neighbourhood is set to 0.1.
#' h=henon(start = c(0.63954883, 0.04772637), do.plot = FALSE)
#' takens = buildTakens(h$x,embedding.dim=2,time.lag=1)
#' neighbours=findAllNeighbours(takens,0.1)
#' }
#' @seealso \code{\link{neighbourSearch}}.
#' @export findAllNeighbours
findAllNeighbours=function(takens,radius,number.boxes=NULL){
  #estimates number.boxes if it has not been specified
  if (is.null(number.boxes)) number.boxes = estimateNumberBoxes(takens, radius)
  #arguments for the C function
  numberTakens=nrow(takens)
  embedding.dim=ncol(takens)
  lengthBoxes=number.boxes*number.boxes+1
  boxes=rep(0,lengthBoxes)
  possibleNeighbours=rep(0,numberTakens)
  neighList=rep(-1,numberTakens)
  nfound=0
  
  #get histogram using boxes
  sol=boxAssistant(takens, radius,number.boxes)
  boxes=sol$boxes;
  possibleNeighbours=sol$possibleNeighbours
  
  allneighs=list()
  for (i in 1:numberTakens){
    #position in the vector was specified in "R-form" (starts at 1)
    #Translate to "C-form" (starts at 0)
    positionTakens=i-1
    #call the c function
    cneighs=.C("neighbourSearchFromBoxes",takens=as.double(takens),positionTakens=as.integer(positionTakens),
               numberTakens=as.integer(numberTakens), embeddingD=as.integer(embedding.dim),
               eps=as.double(radius), numberBoxes=as.integer(number.boxes),
               boxes=as.integer(boxes),possibleNeighbours=as.integer(possibleNeighbours),
               neighList=as.integer(neighList),
               nfound=as.integer(nfound),
               PACKAGE="nonlinearTseries")
    #arrange solutions
    if (cneighs$nfound==0){
      allneighs[[i]]=list()
    }else{
      #eliminate -1 and repeated numbers. convert c-vector positions to R positions
      auxiliarNeighList=(cneighs$neighList[1:(cneighs$nfound)]+ 1)
      allneighs[[i]]=as.vector(auxiliarNeighList)
    }
  }
  
  allneighs = propagateTakensAttr(allneighs, takens)
  attr(allneighs,"radius") = radius
  return (allneighs)
}




#brute force algorithm that computes all distances  to find
#each neighbourhood. This function will be used to test the
#box assisted algorithm
findAllNeighboursDist=function(data,radius){
  dst=as.matrix(dist(data,method="maximum"))
  dst[col(dst)==row(dst)]=Inf
  sol=list()
  l=nrow(data)
  for (i in 1:l){
    aux=as.vector(which(dst[i,]<radius))
    if (length(aux)==0){
      sol[[i]]=list()
    }else{
      sol[[i]]=aux
    }
  }
  return(sol)
}

# finds at least N neighbours for each takens' vector stored in the takens array 
# using firstly neigbourhoods of size radius. If there are not found enough neighbours, the
# neighbourhood size is increased by +radius.increment.
find.at.least.N.neighboursBIS=function(takens,radius,radius.increment,N,number.boxes=NULL){
  n.takens=nrow(takens)
  if (N > n.takens ) stop(paste("Can not find ",N," neighbours if there exist only ",n.takens," takens' vectors\n"))
  if (is.null(number.boxes)) number.boxes = estimateNumberBoxes(takens, radius)
  # data structure to store the solution
  all.neighs=list()
  # let's find at least N neighbours for each vector
  original.radius=radius
  current.radius=radius
  for (i in 1:n.takens){
     neighs.aux= find.at.least.N.neighbours.for.vectorBIS(takens,i,current.radius,radius.increment,N,number.boxes)
     vectors.and.distances=getDistancesAndSort(i,neighs.aux,takens)
     which.vec = vectors.and.distances$dist<original.radius
     if (sum(which.vec)>=N){
       all.neighs[[i]] = vectors.and.distances$neigh[which.vec]
       current.radius=original.radius
     }else{
       all.neighs[[i]] = vectors.and.distances$neigh#[1:5]
       current.radius = vectors.and.distances$dist[[5]] +0.001#add a small increment to ensure that the fith element would be in the neighbourhood
     } 
  }
  return(all.neighs)
}

getDistancesAndSort=function(reference,neighs.positions,takens){
  l = length(neighs.positions)
  dst=rep(0,l)
  for (i in 1:l){
    dst[[i]]=dist(rbind(takens[reference,], takens[neighs.positions[[i]], ]),method="maximum")
  }
  new.pos=order(dst,decreasing=FALSE)
  return(list(neighs=neighs.positions[new.pos],dist=dst[new.pos]))
}

find.at.least.N.neighbours.for.vectorBIS=function(takens,positionTakens,radius,radius.increment,N,number.boxes){
  n.found = 0
  while (n.found < N){
    neighbours = neighbourSearch(takens,positionTakens,radius,number.boxes)
    n.found = neighbours$nfound
    if (n.found < N)  radius = radius + radius.increment
  }
  return (neighbours$neighList)
}

##########3
find.at.least.N.neighbours=function(takens,radius,radius.increment,N,number.boxes=NULL){
  n.takens=nrow(takens)
  if (N > n.takens ) stop(paste("Can not find ",N," neighbours if there exist only ",n.takens," takens' vectors\n"))
  if (is.null(number.boxes)) number.boxes = estimateNumberBoxes(takens, radius)
  # data structure to store the solution
  all.neighs=list()
  # let's find at least N neighbours for each vector
  for (i in 1:n.takens){
    all.neighs[[i]] = find.at.least.N.neighbours.for.vector(takens,i,radius,radius.increment,N,number.boxes)
  }
  return(all.neighs)
}

find.at.least.N.neighbours.for.vector=function(takens,positionTakens,radius,radius.increment,N,number.boxes){
  n.found = 0
  while (n.found < N){
   neighbours = neighbourSearch(takens,positionTakens,radius,number.boxes)
   n.found = neighbours$nfound
   radius = radius + radius.increment
  }
  return (neighbours$neighList)
}

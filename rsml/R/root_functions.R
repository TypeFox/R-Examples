

######################################################################

#' Root constructor
#' @param nodes = the nodes composing the root. Can be null
#' @param parent = the identifier of the root parent. Can be null
#' @param children = vector of children roots (root objects)
#' @param id = root unique identifier
#' @param insertion = insertion position of the root on its parent. Can be null is no parent
#' @param insertion_angle = insertion angle of the root on its parent. Can be null is no parent
#' @return the root
#' @author Guillaume Lobet - guillaume.lobet(at)ulg.ac.be
#' @export
#' @examples
#' r <- root()
#' 
#' n <- node(1, 1)
#' r <- root(n)
root = 
  function(nodes = NULL, 
           parent = "", 
           children = NULL, 
           id="", 
           insertion=NULL, 
           insertion_angle=NULL)
  {
    rts = list(nodes = nodes, 
               parent = parent, 
               children = children, 
               id=id,
               insertion = insertion,
               insertion_angle = insertion_angle)
    class(rts) = "root"
    rts
  }

######################################################################

#' Compute the insertion of the root on its parent
#' @param parent = the parent root
#' @param current = the current root
#' @return the insertion position
#' @export
#' @examples
#' data(lupin)
#' r <- lupin$roots[[1]]
#' r1 <- r$children[[1]]
#' getInsertionPosition(r, r1)

getInsertionPosition = 
  function(parent, current)
  {
    if(class(parent) != "root" || class(current) != "root" )
      stop("need objects of class root")
    insertion <- 0
    minDist <- 1000
    for(n in parent$nodes){
      dx <- current$nodes[[1]]$x - n$x
      dy <- current$nodes[[1]]$y - n$y
      dist = sqrt(dx^2 + dy^2)
      if(dist < minDist){
        minDist <- dist
        insertion <- n$bLength
      }
    }
    insertion
  }

######################################################################
#' Compute the insertion angle of the root on its parent
#' @param parent = the parent root
#' @param current = the current root
#' @return the insertion angle, in degree
#' @export
#' @examples
#' data(lupin)
#' r <- lupin$roots[[1]]
#' r1 <- r$children[[1]]
#' getInsertionAngle(r, r1)
getInsertionAngle = 
  function(parent, current)
  {
    if(class(parent) != "root" || class(current) != "root" )
      stop("need objects of class root")
    parentAngle <- 0
    minDist <- 1000
    for(n in parent$nodes){
      dx <- current$nodes[[1]]$x - n$x
      dy <- current$nodes[[1]]$y - n$y
      dist = sqrt(dx^2 + dy^2)
      if(dist < minDist){
        minDist <- dist
        parentAngle <- as.numeric(n$orientation)
      }
    }
    count <- 4
    ang <- 0
    if(length(current$nodes) < count) count <- length(current$nodes)
    for(i in 1:count){
      ang <- ang + as.numeric(current$nodes[[i]]$orientation)
    }
    ang = ang / count
    
    insertAngl <- 0
    if (ang > parentAngle){
      insertAngl = ang - parentAngle
    }else{
      insertAngl = parentAngle - ang
    }
    if (insertAngl > pi ) insertAngl = (2 * pi) - insertAngl
    
    insertAngl * (180 / pi)
  }

######################################################################

#' Add a node to an existing root
#' @param ro = the current root
#' @param no = the current node
#' @return the root, with the added node
#' @author Guillaume Lobet - guillaume.lobet(at)ulg.ac.be
#' @keywords rsml
#' @export
#' @examples
#' n <- node(1, 1)
#' r <- root()
#' r <- addNodeToRoot(r, n)
addNodeToRoot = 
  function(ro, no)
  {
      if(class(no) != "node" || class(ro) != "root" )
        stop("need objects of class node and root")
     ro$nodes[[length(ro$nodes) + 1]] <- no
     ro
  }

######################################################################
#' Get the number of nodes in a root
#' @param obj of class root
#' @return the number of nodes in the root
#' @author Guillaume Lobet - guillaume.lobet(at)ulg.ac.be
#' @keywords rsml
#' @export
#' @examples
#' data(lupin)
#' r <- lupin$roots[[1]]
#' nNode(r)
nNode = function(obj) length(obj$nodes)

######################################################################

#' Compute the xrange of the root
#' @param obj of class root
#' @keywords rsml
#' @return c(x1,x2) where x1 and x2 are the x limits of the root
#' @export
#' @examples
#' data(lupin)
#' r <- lupin$roots[[1]]
#' xrange(r)
xrange = 
  function(obj)
  {
    xmin = 1e9;
    xmax = -1e9;
    for(i in 1:nNode(obj)){
      n <-  obj$nodes[[i]]
      if(n$x < xmin) xmin <- n$x
      if(n$x > xmax) xmax <- n$x
    }   
    c(xmin, xmax)
  }

######################################################################

#' Compute the xrange of the root
#' @param obj of class root
#' @keywords rsml
#' @return c(y1,y2) where y1 and y2 are the y limits of the root
#' @export
#' @examples
#' data(lupin)
#' r <- lupin$roots[[1]]
#' yrange(r)
yrange = 
  function(obj)
  {
    ymin = 1e9;
    ymax = -1e9;
    for(i in 1:nNode(obj)){
      n <-  obj$nodes[[i]]
      if(n$y < ymin) ymin <- n$y
      if(n$y > ymax) ymax <- n$y
    }   
    c(ymin, ymax)
  }


######################################################################

#' Compute the xrange of the root
#' @param obj of class root
#' @keywords rsml
#' @return c(y1,y2) where y1 and y2 are the y limits of the root
#' @export
#' @examples
#' data(lupin)
#' r <- lupin$roots[[1]]
#' zrange(r)
zrange = 
  function(obj)
  {
    zmin = 1e9;
    zmax = -1e9;
    for(i in 1:nNode(obj)){
      n <-  obj$nodes[[i]]
      if(n$y < zmin) zmin <- n$z
      if(n$y > zmax) zmax <- n$z
    }   
    c(zmin, zmax)
  }



######################################################################

#' Compute the length of the root based on the coordinatres of its nodes
#' @param x object of class root
#' @return the length of the root
#' @author Guillaume Lobet - guillaume.lobet(at)ulg.ac.be
#' @keywords rsml
#' @export
#' @examples
#' data(lupin)
#' r <- lupin$roots[[1]]
#' length(r)
length.root = 
  function(x)
  {
    obj <- x
    obj$nodes[length(obj$nodes)][[1]]$bLength
  }

######################################################################

#' Add a child root (lateral) to an existing root
#' @param current = the current root
#' @param child = the child root to attach
#' @return the current root, with the additional child attached
#' @keywords rsml
#' @author Guillaume Lobet - guillaume.lobet(at)ulg.ac.be
#' @export
#' @examples
#' data(lupin)
#' current <- lupin$roots[[1]]
#' child <- current$children[[1]]
#' current <- addChildToRoot(current, child)
addChildToRoot = 
  function(current, child)
  {
    if(class(child) != "root" || class(current) != "root" )
      stop("need objects of class root")
    current$children[[length(current$children) + 1]] <- child
    current
  }

######################################################################

#' Get the number of children in a root
#' @param obj of class root
#' @return the number of child root in the current root
#' @author Guillaume Lobet - guillaume.lobet(at)ulg.ac.be
#' @keywords rsml
#' @export
#' @examples
#' data(lupin)
#' r <- lupin$roots[[1]]
#' nChild(r)
nChild = function(obj) length(obj$children)

######################################################################

#' Compute the mean insertion angle of the children (lateral) roots
#' @param obj of class root
#' @return the mean lateral angle
#' @author Guillaume Lobet - guillaume.lobet(at)ulg.ac.be
#' @keywords rsml
#' @export
#' @examples
#' data(lupin)
#' r <- lupin$roots[[1]]
#' meanInsertionAngle(r)
meanInsertionAngle = 
  function(obj)
  {
    ang <- 0
    for(i in 1:nChild(obj)){
      ang <- ang + obj$children[[i]]$insertion_angle
    }
    mean <- ang / nChild(obj)
    if(length(mean) == 0) mean=0
    mean
  }

######################################################################

#' Compute the mean interbranch distance of the children (lateral) roots
#' @param obj of class root
#' @param allroot if true, compute the interbanch distance on the whole root
#' @return the mean interbranch distance
#' @author Guillaume Lobet - guillaume.lobet(at)ulg.ac.be
#' @export
#' @keywords rsml
#' @examples
#' data(lupin)
#' r <- lupin$roots[[1]]
#' meanInterbranch(r)
meanInterbranch = 
  function(obj, allroot=F)
  {
    min <- 10000
    max <- 0
    if(!allroot){
      for(i in 1:nChild(obj)){
        if(obj$children[[i]]$insertion < min) min <- obj$children[[i]]$insertion
        if(obj$children[[i]]$insertion > max) max <- obj$children[[i]]$insertion
      }
      nChild(obj) / (max-min)
    }else{ 
      nChild(obj) / length(obj)
    }
  }

######################################################################

#' Compute the length of the root and its children based on the coordinates of its nodes
#' @param obj of class root
#' @return the total lenght of the root and children
#' @author Guillaume Lobet - guillaume.lobet(at)ulg.ac.be
#' @export
#' @keywords rsml
#' @examples
#' data(lupin)
#' r <- lupin$roots[[1]]
#' totalLength(r)
totalLength = 
  function(obj)
  {
    l = length(obj)
    for(i in 1:nChild(obj)){
      l <- l + length(obj$children[[i]])
      for(j in 1:nChild(obj$children[[i]])){
        l <- l + length(obj$children[[i]]$children[j])
      }
    }
    l
  }

######################################################################

#' Get the coordinates of the root nodes
#' @param obj of class root
#' @return a dataframe containing the node coordinates (x, y, z)
#' @author Guillaume Lobet - guillaume.lobet(at)ulg.ac.be
#' @keywords rsml
#' data(lupin)
#' r <- lupin$roots[[1]]
#' coords(r)
coords = 
  function(obj)
  {
    coords <- data.frame(x=numeric(nNode(obj)), y=numeric(nNode(obj)), z=numeric(nNode(obj)))
    for(i in 1:nNode(obj)){
      coords$x[i] <- obj$nodes[[i]]$x
      coords$y[i] <- obj$nodes[[i]]$y
      coords$z[i] <- obj$nodes[[i]]$z   
    }
    coords
  }

###################################

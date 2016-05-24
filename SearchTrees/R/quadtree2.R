setClass("SearchTree", representation = list(ref = "externalptr", numNodes = "integer", dataNodes = "integer", maxDepth = "integer", maxBucket = "integer", totalData = "integer", dataType="character", "VIRTUAL"))
setClass("QuadTree", contains="SearchTree")

setGeneric("rectLookup",
           function(tree, ptOne, ptTwo, xlims, ylims)
           standardGeneric("rectLookup")
           )

setMethod("rectLookup", "QuadTree",
          function(tree, ptOne, ptTwo, xlims, ylims)
          {
            if(missing(xlims))
              xlims = sort(c(ptOne[1], ptTwo[1]))
            if(missing(ylims))
              ylims = sort(c(ptOne[2], ptTwo[2]))
            .Call("R_Rectangle_Lookup", tree, as.numeric(xlims), as.numeric(ylims))
          }
          )


findMaxDepth = function(maxDepth, minNodeArea, xlim, ylim)
  {

    if (!missing(minNodeArea))
      {
        totArea = (xlim[2] - xlim[1]) * (ylim[2] - ylim[1])
        areas = which(totArea / (4^(1:10)) <= minNodeArea)
        
        if(!length(areas))
          {
            warning("The minNodeArea selected lead to a maximum depth > 10, which is very memory intensive for negligable benefit. Using maximum depth of 10.")
           
            maxDepth = 10
          } else {
            maxDepth = areas[1]
          }

      }
    as.integer(maxDepth)
  }
    

createTree = function(data, treeType = "quad", dataType = "point", columns = if (dataType=="point") 1:2 else 1:4, ...)
  {
    if (tolower(treeType) == "quad")
      {
        if(dataType == "point")
          {
            if(length(columns) != 2)
              stop("wrong number of columns for this index type.")
            x = data[,columns[1]] 
            y = data[,columns[2]]
            ret = quadTree2(x,y, ... )
          } else if (dataType == "rect") {
            ret = rectTree(data[, columns[1]], data[, columns[2]], data[, columns[3]], data[, columns[4]])
            dataType = "rectangle"
          }
        ret@dataType = dataType
        ret@totalData = length(data[,columns[1]])
        ret
      }
  }
        
quadTree2 = function(x, y, maxDepth = 7, minNodeArea, ...)
  { 
    xrange = range(x)
    yrange = range(y)
    maxDepth = findMaxDepth(maxDepth, minNodeArea, xrange, yrange)

    x = as.numeric(x)
    y = as.numeric(y)
    .Call("R_Build_Quadtree_Pt", x, y, max(x), min(x), max(y), min(y), maxDepth)

  }
  
rectTree = function(x1, x2, y1, y2, maxDepth = 7, minNodeArea, ...)
  {
    x1 = as.numeric(x1)
    x2 = as.numeric(x2)
    y1 = as.numeric(y1)
    y2 = as.numeric(y2)
    xlim = c(min(x1), max(x2))
    ylim = c(min(y1), max(y2))
    maxDepth = findMaxDepth(maxDepth, minNodeArea, xlim, ylim)
    .Call("R_Build_Quadtree_Rect", x1, x2, y1, y2, max(x2), min(x1), max(y2), min(y1), maxDepth)
  }
setGeneric("knnLookup",
           function(tree, newx, newy, newdat, columns=1:2,   k = 5)
           standardGeneric("knnLookup")
           )
setMethod("knnLookup", "QuadTree",
          function(tree, newx, newy, newdat, columns,  k)
          {

            k = as.integer(k)
            if(missing(newx))
              newx = as.numeric(newdat[,columns[1]])
            if(missing(newy))
              newy = as.numeric(newdat[,columns[2]])
            inds = .Call("R_Find_Neighbors_Pts", tree, newx, newy, k )
            matrix(inds, byrow = TRUE, ncol = k)
          }
          )

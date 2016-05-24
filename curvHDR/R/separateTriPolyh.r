########## R function: separateTriPolyh ##########

# For separating a triangular mesh into contiguous
# components.

# Last changed: 07 MAR 2016

separateTriPolyh <- function(contour3dObj)
{
   # Set function for compute for each triangle the indices of
   # triangles that share an edge with it (this function was
   # developed by Luke Tierney in June 2008, and could be
   # made more efficient).

   triangleNeighbors <- function(tris) 
   {
      ve <- t2ve(tris)
      vt <- vertexTriangles(ve)
      ib <- ve$ib
      n.tri <- ncol(ib)
      tn <- vector("list",n.tri)
      for (i in 1 : n.tri)
      {
         v1 <- unique(vt[[ib[1,i]]])
         v2 <- unique(vt[[ib[2,i]]])
         v3 <- unique(vt[[ib[3,i]]])
         i12 <- intersect(v1,v2)
         i23 <- intersect(v2,v3)
         i31 <- intersect(v3,v1)
         u <- union(union(i12,i23),i31)
         tn[[i]] <- u[u!=i]
      }
      return(tn)
   }
   
   # Set function for implementing Dijkstra's version of Rem's algorithm
   # for computing equivalence classes based on a number of vertices 1:nvert
   # and a set of N edges provided as an N x 2 matrix
   # (also developed by Luke Tierney in June, 2008).

   getPatches <- function(nvert,edges) 
   {
      f <- 1:nvert

      if (!(is.vector(edges)) && dim(edges)[1] != 0)
      {
         nedge <- nrow(edges)

         for (e in 1:nedge) 
         {
            p0 <- edges[e, 1]
            q0 <- edges[e, 2]
            p1 <- f[p0]
            q1 <- f[q0]
            while (p1 != q1)
            {
               if (q1 < p1)
               {
                  f[p0] <- q1
                  p0 <- p1
                  p1 <- f[p1]
              }
              else
              {
                 f[q0] <- p1
                 q0 <- q1
                 q1 <- f[q1]
              }
            }
         }
      }
      if(is.vector(edges))
      {
        if(edges[1] < edges[2])
            f[edges[2]] <- edges[1]
        else  f[edges[1]] <- edges[2]
      }

      for (v in 1:nvert)
        f[v] <- f[f[v]]

      return(split(1:nvert,f))
   }

   # Set up function for computing the edges to indicate which triangles
   # share an edge (also developed by Luke Tierney in June, 2008).
   
   triangleNeighborEdges <- function(tn) 
   {
       edges <- function(i)
       {
          v <- tn[[i]]
          if (length(v) > 0) cbind(i,v)
          else numeric(0)
       }
       do.call(rbind, lapply(1:length(tn), edges))
   }

   # Create Tierney function objects:
   
   T1Obj <- triangleNeighbors(contour3dObj)
   T2Obj <- getPatches(length(T1Obj),triangleNeighborEdges(T1Obj))
   numPolyhedra <- length(T2Obj)

   # Obtain the list of triangular mesh objects (class=Triangle3D):

   outputPolyhedra <- list(contour3dObj,numPolyhedra)

   # Add fix for apparent problem that arises when the
   # input consists of a single polyhedron:

   if ((length(outputPolyhedra)==2)&(outputPolyhedra[[2]]==1))
      outputPolyhedra <- list(outputPolyhedra[[1]])
   
   
   for (j in 1:numPolyhedra) 
   {
      outputPolyhedra[[j]] <- contour3dObj
      outputPolyhedra[[j]]$v1 <- contour3dObj$v1[T2Obj[[j]],]
      outputPolyhedra[[j]]$v2 <- contour3dObj$v2[T2Obj[[j]],]
      outputPolyhedra[[j]]$v3 <- contour3dObj$v3[T2Obj[[j]],]
   }

   return(outputPolyhedra)
}

########## End of separateTriPolyh ##########



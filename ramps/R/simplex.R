makeFirstSimplex <- function(kappa, control)
## kappa:   current kappa vector
{
   vertices <- diag(length(kappa))

   ## the multiplier to get the size of the random simplex
   mult <- control$sigma2.e$tuning[1]
   if (mult < 1) {
      ## get the coordinates of the simplex shrinked towards vertex 1
      vertices[,-1] <- vertices[,-1] + (1 - mult) * (vertices[,1] - vertices[,-1])

      ## first, generate a point uniformly in appropriate size triangle
      ## then translate vertices to make x correspond to random point
      rbaryc <- runif.simplex(length(kappa))
      vertices <- vertices + kappa - as.vector(vertices %*% t(rbaryc))
   }

   vertices
}


shrinkSimplex <- function(bx, bc, cx, cc, vertices)
## bx: coefficients of vertices to produce x (barycentric coordinates of x)
## bc: coefficients of vertices to produce cand (barycentric coords of cand)
## cx, cc:  Cartesian coords of x and cand
## vertices: vertices of the current simplex
{
   for (i in seq(bc)[bc < bx]) {
      vertices[,-i] <- vertices[,-i] + bc[i] * (vertices[,i] - vertices[,-i])
      bc <- solve(vertices, cc)
   }

   vertices
}


sliceSimplex <- function(theta, mpdfun, log=FALSE, f.theta, control, ...)
## theta:    current point
## mpdfun:   marginalized posterior density up to a scale or additive constant
## log:      work on log fun or not
## f.theta:  mpdfun(theta)
{
   y <- ifelse(log, f.theta - rexp(1), runif(1, 0, f.theta))

   ## interval for generating phi
   phi.old <- params2phi(theta, control)
   width <- control$phi$tuning
   l <- phi.old - runif(length(phi.old), 0, width)
   r <- l + width

   ## simplex for generating kappa
   kappa.old <- params2kappa(theta, control)
   vertices <- makeFirstSimplex(kappa.old, control)
   kappa.old.b <- solve(vertices, kappa.old)

   newevals <- 0
   repeat {
      ## sample phi.cand and kappa.cand
      phi.cand <- runif(length(phi.old), l, r)
      kappa.cand.b <- runif.simplex(length(kappa.old))
      kappa.cand <- as.vector(vertices %*% t(kappa.cand.b))

      # out of boundary indicator
      ob.phi <- control$phi$f(phi.cand) <= 0
      ob.kappa <- min(kappa.cand) < 0 || max(kappa.cand) > 1
      ob <- ob.phi || ob.kappa

      if (!ob) { ## cand is in the boundary
         cand <- c(phi.cand, kappa.cand)
         evalresult <- mpdfun(cand, control = control, ...)

         if (evalresult$value > -Inf) {  ## valid x0 to mpdfun
            newevals <- newevals + 1
         }

         ## here is the only exit of the loop
         if (y < evalresult$value) {
            evalresult$theta <- cand
            evalresult$newevals <- newevals
            return(evalresult)
         }
      }

      if (!ob || ob.phi) {     ## shrink for phi
         idx <- phi.cand < phi.old
         l <- ifelse(idx, phi.cand, l)
         r <- ifelse(idx, r, phi.cand)
      }
      if (!ob || ob.kappa) {   ## shrink for kappa
         vertices <- shrinkSimplex(kappa.old.b, kappa.cand.b, kappa.old,
                                   kappa.cand, vertices)
         kappa.old.b <- solve(vertices, kappa.old)
      }
   }
}

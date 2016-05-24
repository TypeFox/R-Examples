"plotParallel" <-
function(parallel,
         eig    = NA,
         x      = eig,
         model  = "components",
         legend = TRUE,
         ylab   = "Eigenvalues",
         xlab   = "Components",
         main   = "Parallel Analysis",
         ...
                           ) {                       
  if (any(!is.na(x))) eig <- eigenComputes(x, ...)
  if (!inherits(parallel, "parallel")) stop("Method is only for parallel objects")
  if (model == "factors") xlab <- "Factors"
  var        <- length(parallel$eigen$qevpea)
  if (length(eig) == 1) {  
   Component <- var:1
   Location  <- seq(from = 0, to = max(parallel$eigen$qevpea)*3, length.out = var)
   plot.default(as.numeric(Component),
                as.numeric(Location), 
                type = "n", 
                main = main, 
                xlab = xlab, 
                ylab = ylab) 
    }
    
  if (length(eig) > 1) {plotuScree(eig, main = main, xlab = xlab, ylab = ylab) }
  lines(1:var, parallel$eigen$qevpea , col = "green", type = "p", pch = 2)
  lines(1:var, parallel$eigen$mevpea,  col = "red")
  if (legend == TRUE) {
   if (length(eig) == 1) { 
     leg <-  c("Mean Eigenvalues", "Centiles of the Eigenvalues")
     tco <-  c("red", "green")
     co  <-  c("red", "green")
     pc  <-  c(NA, 2)
   }
   if (length(eig) > 1) { 
     leg <-  c("Eigenvalues", "Mean Eigenvalues", "Centiles of the Eigenvalues")
     tco <-  c("black", "red", "green")
     co  <-  c("black", "red", "green")
     pc  <-  c(1, NA, 2)
   }
   legend("topright", 
          legend   = leg,
          text.col = tco, 
          col      = co,
          pch      = pc 
          )
    }      
  }


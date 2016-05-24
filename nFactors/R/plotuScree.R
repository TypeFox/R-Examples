"plotuScree" <-
function(Eigenvalue, x=Eigenvalue, model  = "components",
         ylab   = "Eigenvalues",
         xlab   = "Components",
         main   = "Scree Plot" ,
         ...) {
 Eigenvalue  <- eigenComputes(x, ...)
 if (!inherits(Eigenvalue, "numeric")) stop("use only with \"numeric\" objects")
 if (model == "factors") xlab <- "Factors"
 par(mfrow = c(1,1))
 nk          <- length(Eigenvalue)
 Component   <- 1:nk
 plot.default(as.numeric(Component),
              as.numeric(Eigenvalue),
              type = 'b',col = "black", pch = 1,
              ylab = ylab,
              xlab = xlab,
              main = main
              ) 
 }


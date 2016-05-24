######################################
######    Irr calculation      #######
######################################

setGeneric("Irr", 
           def = function(assemblages,
                          W,
                          abundance = F,
                          Wmin = min(W),
                          Wmax = max(W)
           )
           {
             standardGeneric( "Irr" )
           }
)

setMethod("Irr",
          signature(assemblages = "vector", W = "vector"),
          function(assemblages, W, abundance = F, Wmin = min(W), Wmax = max(W))
          {
            if(any(!(names(assemblages) %in% names(W))))
            {
              stop(paste("Certain species names in assemblages have no equivalent in W:",
                         paste("\n", names(assemblages)[which(!(names(assemblages) %in% names(W)))], collapse = " ")))
            }
            # Calculating richness first because species with NAs will be removed later
            richness <- length(assemblages[assemblages > 0])
            
            if(any(is.na(W)))
            {
              cat("W contains NA values, they will be removed\n")
            }
            if(any(is.na(Wmin)))
            {
              Wmin = min(W, na.rm = TRUE)
            }
            if(any(is.na(Wmax)))
            {
              Wmax = max(W, na.rm = TRUE)
            }
            
            # Matching species names in assemblages and W
            W <- W[match(names(assemblages), names(W))]
            
            if(sum(assemblages) == 0)
            {
              IrrValue <- 0
            } else
            {
              IrrValue <- .irr(A = assemblages, W = W, abundance = abundance, Wmin = Wmin, Wmax = Wmax)
            }
            return(c(Irr = IrrValue, Richness = richness))
          }
)

setMethod("Irr",
          signature(assemblages = "vector", W = "matrix"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), Wmax = apply(W, 2, max))
          {
            if(any(!(names(assemblages) %in% rownames(W))))
            {
              stop(paste("Certain species names in assemblages have no equivalent in W:",
                         paste("\n", names(assemblages)[which(!(names(assemblages) %in% rownames(W)))], collapse = " ")))
            }
            # Calculating richness first because species with NAs will be removed later
            richness <- length(assemblages[assemblages > 0])
            
            if(any(is.na(W)))
            {
              cat("W contains NA values, they will be removed\n")
            }
            if(any(is.na(Wmin)))
            {
              Wmin = apply(W, 2, min, na.rm = TRUE)
            }
            if(any(is.na(Wmax)))
            {
              Wmax = apply(W, 2, max, na.rm = TRUE)
            }
            
            if(any(colnames(W) %in% c("Q", "R", "cut.off", paste("Q", 1:1000, sep = ""), paste("R", 1:1000, sep = ""),  paste("cut.off", 1:1000, sep = ""))))
            {
              W <- W[, -which(colnames(W) %in% c("Q", "R", "cut.off", 
                                                 paste("Q", 1:1000, sep = ""), 
                                                 paste("R", 1:1000, sep = ""),  
                                                 paste("cut.off", 1:1000, sep = ""))), drop = T]
            }
            if(any(names(Wmin) %in% c("Q", "R", "cut.off", paste("Q", 1:1000, sep = ""), paste("R", 1:1000, sep = ""),  paste("cut.off", 1:1000, sep = ""))))
            {
              Wmin <- Wmin[-which(names(Wmin) %in% c("Q", "R", "cut.off", 
                                                     paste("Q", 1:1000, sep = ""), 
                                                     paste("R", 1:1000, sep = ""),  
                                                     paste("cut.off", 1:1000, sep = "")))]
            }
            if(any(names(Wmax) %in% c("Q", "R", "cut.off", paste("Q", 1:1000, sep = ""), paste("R", 1:1000, sep = ""),  paste("cut.off", 1:1000, sep = ""))))
            {
              Wmax <- Wmax[-which(names(Wmax) %in% c("Q", "R", "cut.off", 
                                                     paste("Q", 1:1000, sep = ""), 
                                                     paste("R", 1:1000, sep = ""),  
                                                     paste("cut.off", 1:1000, sep = "")))]
            }

            if(is.vector(W))
            {
              Irr(assemblages = assemblages, W = W, abundance = abundance, Wmin = Wmin, Wmax = Wmax)
            } else
            {
              
              # Matching species names in assemblages and W
              W <- W[match(names(assemblages), rownames(W)), ]
              
              IrrValue <- sapply(colnames(W), function(x, A. = assemblages, W. = W, abundance. = abundance, Wmin. = Wmin, Wmax. = Wmax)
              {
                .irr(A = A., W = W.[, x], abundance = abundance., Wmin = Wmin.[x], Wmax = Wmax.[x])
              })
              names(IrrValue) <- paste("Irr_", colnames(W), sep = "")
              
              IrrValue <- c(IrrValue,
                            Richness = richness)
              return(IrrValue)
            }
          }
)

setMethod("Irr",
          signature(assemblages = "vector", W = "data.frame"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), Wmax = apply(W, 2, max))
          {
            W <- as.matrix(W)
            Irr(assemblages, W, abundance = abundance, Wmin = Wmin, Wmax = Wmax)
          }
)

setMethod("Irr",
          signature(assemblages = "matrix", W = "vector"),
          function(assemblages, W, abundance = F, Wmin = min(W), Wmax = max(W))
          {
            if(any(!(rownames(assemblages) %in% names(W))))
            {
              stop(paste("Certain species names in assemblages have no equivalent in W:",
                         paste("\n", rownames(assemblages)[which(!(rownames(assemblages) %in% names(W)))], collapse = " ")))
            }
            # Calculating richness first because species with NAs will be removed later
            richness <- apply(assemblages, 2, function(x) length(x[x > 0]))
            
            if(any(is.na(W)))
            {
              cat("W contains NA values, they will be removed\n")
            }
            if(any(is.na(Wmin)))
            {
              Wmin = min(W, na.rm = TRUE)
            }
            if(any(is.na(Wmax)))
            {
              Wmax = max(W, na.rm = TRUE)
            }
            
            # Matching species names in assemblages and W
            W <- W[match(rownames(assemblages), names(W))]
            
            
            IrrValue <- sapply(colnames(assemblages), 
                               function(x, 
                                        A. = assemblages, 
                                        W. = W, 
                                        abundance. = abundance, 
                                        Wmin. = Wmin, 
                                        Wmax. = Wmax)
            {
              .irr(A = A.[, x], W = W., abundance = abundance., Wmin = Wmin., Wmax = Wmax.)
            })
            
            IrrValue <- cbind(Irr = IrrValue,
                              Richness = richness)
            rownames(IrrValue) <- colnames(assemblages)
            return(IrrValue)
          }
)

setMethod("Irr",
          signature(assemblages = "matrix", W = "matrix"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), Wmax = apply(W, 2, max))
          {
            if(any(!(rownames(assemblages) %in% rownames(W))))
            {
              stop(paste("Certain species names in assemblages have no equivalent in W:",
                         paste("\n", rownames(assemblages)[which(!(rownames(assemblages) %in% rownames(W)))], collapse = " ")))
            }
            # Calculating richness first because species with NAs will be removed later
            richness <- apply(assemblages, 2, function(x) length(x[x > 0]))
            
            if(any(is.na(W)))
            {
              cat("W contains NA values, they will be removed\n")
            }
            if(any(is.na(Wmin)))
            {
              Wmin = apply(W, 2, min, na.rm = TRUE)
            }
            if(any(is.na(Wmax)))
            {
              Wmax = apply(W, 2, max, na.rm = TRUE)
            }
            
            if(any(colnames(W) %in% c("Q", "R", "cut.off", paste("Q", 1:1000, sep = ""), paste("R", 1:1000, sep = ""),  paste("cut.off", 1:1000, sep = ""))))
            {
              W <- W[, -which(colnames(W) %in% c("Q", "R", "cut.off", 
                                                 paste("Q", 1:1000, sep = ""), 
                                                 paste("R", 1:1000, sep = ""),  
                                                 paste("cut.off", 1:1000, sep = ""))), drop = T]
            }
            if(any(names(Wmin) %in% c("Q", "R", "cut.off", paste("Q", 1:1000, sep = ""), paste("R", 1:1000, sep = ""),  paste("cut.off", 1:1000, sep = ""))))
            {
              Wmin <- Wmin[-which(names(Wmin) %in% c("Q", "R", "cut.off", 
                                                     paste("Q", 1:1000, sep = ""), 
                                                     paste("R", 1:1000, sep = ""),  
                                                     paste("cut.off", 1:1000, sep = "")))]
            }
            if(any(names(Wmax) %in% c("Q", "R", "cut.off", paste("Q", 1:1000, sep = ""), paste("R", 1:1000, sep = ""),  paste("cut.off", 1:1000, sep = ""))))
            {
              Wmax <- Wmax[-which(names(Wmax) %in% c("Q", "R", "cut.off", 
                                                     paste("Q", 1:1000, sep = ""), 
                                                     paste("R", 1:1000, sep = ""),  
                                                     paste("cut.off", 1:1000, sep = "")))]
            }
            
            if(is.vector(W))
            {
              Irr(assemblages, W, abundance = abundance, Wmin = Wmin, Wmax = Wmax)
            } else
            {
              
              # Matching species names in assemblages and W
              W <- W[match(rownames(assemblages), rownames(W)), ]
              
              IrrValue <- sapply(colnames(assemblages), 
                                 function(x, 
                                          A. = assemblages, 
                                          W. = W, 
                                          abundance. = abundance, 
                                          Wmin. = Wmin, 
                                          Wmax. = Wmax)
                                 {
                                   sapply(colnames(W.), function(y, 
                                                                 A.. = A.[, x], 
                                                                 W.. = W., 
                                                                 abundance.. = abundance., 
                                                                 Wmin.. = Wmin., 
                                                                 Wmax.. = Wmax.)
                                   {
                                     .irr(A = A.., W = W..[, y], abundance = abundance.., Wmin = Wmin..[y], Wmax = Wmax..[y])
                                   })
                                 })
              rownames(IrrValue) <- paste("Irr_", colnames(W), sep = "")
              
              IrrValue <- cbind(t(IrrValue),
                                Richness = richness)
              return(IrrValue)
            }
          }
)

setMethod("Irr",
          signature(assemblages = "matrix", W = "data.frame"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), Wmax = apply(W, 2, max))
          {
            W <- as.matrix(W)
            Irr(assemblages, W, abundance = abundance, Wmin = Wmin, Wmax = Wmax)
          }
)

setMethod("Irr",
          signature(assemblages = "data.frame", W = "vector"),
          function(assemblages, W, abundance = F, Wmin = min(W), Wmax = max(W))
          {
            assemblages <- as.matrix(assemblages)
            Irr(assemblages, W, abundance = abundance, Wmin = Wmin, Wmax = Wmax)
          }
)

setMethod("Irr",
          signature(assemblages = "data.frame", W = "matrix"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), Wmax = apply(W, 2, max))
          {
            assemblages <- as.matrix(assemblages)
            Irr(assemblages, W, abundance = abundance, Wmin = Wmin, Wmax = Wmax)
          }
)

setMethod("Irr",
          signature(assemblages = "data.frame", W = "data.frame"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), Wmax = apply(W, 2, max))
          {
            assemblages <- as.matrix(assemblages)
            W <- as.matrix(W)
            Irr(assemblages, W, abundance = abundance, Wmin = Wmin, Wmax = Wmax)
          }
)


# Function to calculate Irr at the vector level
# This functions assumes that A and W are already of the same length
# and that the names of A and W are already matched
.irr <- function(A, W, Wmin, Wmax, abundance)
{
  if(any(is.na(W)))
  {
    NA.pos <- which(is.na(W))
    A <- A[-NA.pos]
    W <- W[-NA.pos]
  }
  if(!abundance)
  {
    A[A > 0] <- 1
  } 
  return((sum(A * W) / length(A[A > 0]) - Wmin) / (Wmax - Wmin))
}
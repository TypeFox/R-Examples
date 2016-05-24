######################################
######    Isr calculation      #######
######################################

setGeneric("Isr", 
           def = function(assemblages,
                          W,
                          abundance = FALSE,
                          Wmin = min(W),
                          normalise = FALSE
           )
           {
             standardGeneric( "Isr" )
           }
)

setMethod("Isr",
          signature(assemblages = "vector", W = "vector"),
          function(assemblages, W, abundance = F, Wmin = min(W), normalise = F)
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

            full.W <- W
            
            # Calculating richness first because species with NAs will be removed later
            richness <- length(assemblages[assemblages > 0])
            
            # Matching species names in assemblages and W
            W <- W[match(names(assemblages), names(W))]
            
            if(sum(assemblages) == 0)
            {
              IsrValue <- 0
            } else
            {
              IsrValue <- .isr(A = assemblages, W = W, abundance = abundance, Wmin = Wmin, normalise = normalise, full.W = full.W)
            }
            return(c(Isr = IsrValue, Richness = richness))
          }
)

setMethod("Isr",
          signature(assemblages = "vector", W = "matrix"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), normalise = F)
          {
            if(any(!(names(assemblages) %in% rownames(W))))
            {
              stop(paste("Certain species names in assemblages have no equivalent in W:",
                         paste("\n", names(assemblages)[which(!(names(assemblages) %in% rownames(W)))], collapse = " ")))
            }
            # Calculating richness first because species with NAs will be removed later
            richness <- length(assemblages[assemblages > 0])
            
            full.W <- W
            
            if(any(is.na(W)))
            {
              cat("W contains NA values, they will be removed\n")
            }
            if(any(is.na(Wmin)))
            {
              Wmin = apply(W, 2, min, na.rm = TRUE)
            }
            
            # Removing non-weight columns (outputs from rWeights())
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
            if(any(names(full.W) %in% c("Q", "R", "cut.off", paste("Q", 1:1000, sep = ""), paste("R", 1:1000, sep = ""),  paste("cut.off", 1:1000, sep = ""))))
            {
              full.W <- full.W[-which(colnames(full.W) %in% c("Q", "R", "cut.off", 
                                                              paste("Q", 1:1000, sep = ""), 
                                                              paste("R", 1:1000, sep = ""),  
                                                              paste("cut.off", 1:1000, sep = ""))), ]
            }
            
            if(is.vector(W))
            {
              Isr(assemblages, W, abundance = abundance, Wmin = Wmin, normalise = normalise)
            } else
            {
              # Matching species names in assemblages and W
              W <- W[match(names(assemblages), rownames(W)), ]
              
              IsrValue <- sapply(colnames(W), function(x, A. = assemblages, W. = W, abundance. = abundance, 
                                                       Wmin. = Wmin, full.W. = full.W, normalise. = normalise)
              {
                .isr(A = A., W = W.[, x], abundance = abundance., Wmin = Wmin.[x], full.W = full.W.[, x], normalise = normalise.)
              })
              names(IsrValue) <- paste("Isr_", colnames(W), sep = "")
              
              IsrValue <- c(IsrValue,
                            Richness = richness)
              return(IsrValue)
            }
          }
)

setMethod("Isr",
          signature(assemblages = "vector", W = "data.frame"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), normalise = F)
          {
            W <- as.matrix(W)
            Isr(assemblages, W, abundance = abundance, Wmin = Wmin, normalise = normalise)
          }
)

setMethod("Isr",
          signature(assemblages = "matrix", W = "vector"),
          function(assemblages, W, abundance = F, Wmin = min(W), normalise = F)
          {
            if(any(!(rownames(assemblages) %in% names(W))))
            {
              stop(paste("Certain species names in assemblages have no equivalent in W:",
                         paste("\n", rownames(assemblages)[which(!(rownames(assemblages) %in% names(W)))], collapse = " ")))
            }
            # Calculating richness first because species with NAs will be removed later
            richness <- apply(assemblages, 2, function(x) length(x[x > 0]))
            
            full.W <- W
            
            if(any(is.na(W)))
            {
              cat("W contains NA values, they will be removed\n")
            }
            if(any(is.na(Wmin)))
            {
              Wmin = min(W, na.rm = TRUE)
            }
            
            # Matching species names in assemblages and W
            W <- W[match(rownames(assemblages), names(W))]
            
            IsrValue <- sapply(colnames(assemblages), 
                               function(x, 
                                        A. = assemblages, 
                                        W. = W, 
                                        abundance. = abundance, 
                                        Wmin. = Wmin, 
                                        full.W. = full.W,
                                        normalise. = normalise)
                               {
                                 .isr(A = A.[, x], W = W., abundance = abundance., 
                                      Wmin = Wmin., full.W = full.W., normalise = normalise.)
                               })
            IsrValue <- cbind(IsrValue,
                              Richness = richness)
            return(IsrValue)
          }
)

setMethod("Isr",
          signature(assemblages = "matrix", W = "matrix"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), normalise = F)
          {
            if(any(!(rownames(assemblages) %in% rownames(W))))
            {
              stop(paste("Certain species names in assemblages have no equivalent in W:",
                         paste("\n", rownames(assemblages)[which(!(rownames(assemblages) %in% rownames(W)))], collapse = " ")))
            }
            # Calculating richness first because species with NAs will be removed later
            richness <- apply(assemblages, 2, function(x) length(x[x > 0]))
                        
            full.W <- W
            
            if(any(is.na(W)))
            {
              cat("W contains NA values, they will be removed\n")
            }
            if(any(is.na(Wmin)))
            {
              Wmin = apply(W, 2, min, na.rm = TRUE)
            }
            
            # Removing non-weight columns (outputs from rWeights())
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
            if(any(names(full.W) %in% c("Q", "R", "cut.off", paste("Q", 1:1000, sep = ""), paste("R", 1:1000, sep = ""),  paste("cut.off", 1:1000, sep = ""))))
            {
              full.W <- full.W[-which(colnames(full.W) %in% c("Q", "R", "cut.off", 
                                                              paste("Q", 1:1000, sep = ""), 
                                                              paste("R", 1:1000, sep = ""),  
                                                              paste("cut.off", 1:1000, sep = ""))), ]
            }
            
            if(is.vector(W))
            {
              Isr(assemblages, W, abundance = abundance, Wmin = Wmin, normalise = normalise)
            } else
            {
              # Matching species names in assemblages and W
              W <- W[match(rownames(assemblages), rownames(W)), ]
              
              IsrValue <- sapply(colnames(assemblages), 
                                 function(x, 
                                          A. = assemblages, 
                                          W. = W, 
                                          abundance. = abundance, 
                                          Wmin. = Wmin, 
                                          full.W. = full.W,
                                          normalise. = normalise)
                                 {
                                   sapply(colnames(W.), function(y, 
                                                                 A.. = A.[, x], 
                                                                 W.. = W., 
                                                                 abundance.. = abundance., 
                                                                 Wmin.. = Wmin., 
                                                                 full.W.. = full.W.,
                                                                 normalise.. = normalise.)
                                   {
                                     .isr(A = A.., W = W..[, y], abundance = abundance..,
                                          Wmin = Wmin..[y], full.W = full.W..[, y], normalise = normalise..)
                                   })
                                 })
              rownames(IsrValue) <- paste("Isr_", colnames(W), sep = "")
              
              IsrValue <- cbind(t(IsrValue),
                                Richness = richness)
              return(IsrValue)
            }
          }
)

setMethod("Isr",
          signature(assemblages = "matrix", W = "data.frame"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), normalise = F)
          {
            W <- as.matrix(W)
            Isr(assemblages, W, abundance = abundance, Wmin = Wmin, normalise = normalise)
          }
)

setMethod("Isr",
          signature(assemblages = "data.frame", W = "vector"),
          function(assemblages, W, abundance = F, Wmin = min(W), normalise = F)
          {
            assemblages <- as.matrix(assemblages)
            Isr(assemblages, W = W, abundance = abundance, Wmin = Wmin, normalise = normalise)
          }
)

setMethod("Isr",
          signature(assemblages = "data.frame", W = "matrix"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), normalise = F)
          {
            assemblages <- as.matrix(assemblages)
            Isr(assemblages, W, abundance = abundance, Wmin = Wmin, normalise = normalise)
          }
)

setMethod("Isr",
          signature(assemblages = "data.frame", W = "data.frame"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), normalise = F)
          {
            assemblages <- as.matrix(assemblages)
            W <- as.matrix(W)
            Isr(assemblages, W, abundance = abundance, Wmin = Wmin, normalise = normalise)
          }
)


# Function to calculate Isr at the vector level
# This functions assumes that A and W are already of the same length
# and that the names of A and W are already matched
.isr <- function(A, W, Wmin, abundance, normalise, full.W)
{
  if(any(is.na(W)))
  {
    NA.pos <- which(is.na(W))
    A <- A[-NA.pos]
    W <- W[-NA.pos]
  }
  if(any(is.na(full.W)))
  {
    full.W <- full.W[-which(is.na(full.W))]
  }
  if(!abundance)
  {
    A[A > 0] <- 1
  }
  if (normalise)
  {
    if (abundance)
    {
      stop("The index cannot be standardised between 0 and 1 when abundance = TRUE. Set abundance = FALSE or normalise = FALSE")
    }
    return(sum(A * W - Wmin) / (sum(full.W - Wmin))) # Formula with presence-absence data 
  } else
  {
    return(sum(A * W)) # Formula with abundance
  }
}

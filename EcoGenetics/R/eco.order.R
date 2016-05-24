#' Ordering the rows of the data frames contained in an ecogen object
#' 
#' @details This program generates an ecogen object with the rows of all
#' the data frames ordered in reference to the row names of the XY data frame.
#' This is useful when the data frames are loaded into the ecogen object,
#' but were not ordered previously. Also, this tool can be useful 
#' for reorder rows when is needed. 
#' First, the reference data frame in the slot XY will be in the desired row order. 
#' This program then aligns all the data frames by coincidence of row names
#' with those in the slot XY.  
#' 
#' @param eco Object of class "ecogen".
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' data(eco.test)
#' eco1 <- eco
#' eco1[["P"]] <- eco[["P"]][sample(1:225), ]  #object with shuffled rows
#' eco1[["E"]] <- eco[["E"]][sample(1:225), ]
#' ordered <- eco.order(eco1)
#' head(ordered[["P"]]); head(eco[["P"]])
#' 
#' }
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' 
#' @export

setGeneric("eco.order", 
           function(eco) {
             
             if(dim(eco@XY)[1] == 0) {
               stop("XY has 0 rows")
             }
             
             posit0 <- seq(along = eco@XY[, 1])
             
             #function to order rows
             int.order <- function(X0, X, p0) {
               
               p <- seq(along = rownames(X))
               temp <- merge(data.frame(rownames(X0), p0), 
                             data.frame(rownames(X), rownames(X), p), 
                             by = 1)
               temp <- temp[order(temp$p0), ]
               c.names <- colnames(X)
               r.names <- temp[, 3]
               X <- data.frame(X[temp$p, ])
               colnames(X) <- c.names
               rownames(X) <- r.names  
               X
             }
             
             eco@P <- int.order(eco@XY, eco@P, p0 = posit0)
             eco@G <- int.order(eco@XY, eco@G, p0 = posit0)
             eco@A <- int.order(eco@XY, eco@A, p0 = posit0)
             eco@E <- int.order(eco@XY, eco@E, p0 = posit0)
             eco@S <- int.order(eco@XY, eco@S, p0 = posit0)
             eco@C <- int.order(eco@XY, eco@C, p0 = posit0)
             
             #ordering missing cells
             dum <- as.matrix(eco@A - eco@A)
             dum[eco@INT@missing.cells] <- NA
             dum <- int.order(eco@XY, dum, p0 = posit0)
             eco@INT@missing.cells <- which(is.na(dum))
             
             eco
             
           })

#' Combining the columns of ecogen object
#' 
#' @param eco1 Object of class "ecogen".
#' @param eco2 Object of class "ecogen".
#' @param ... Other "ecogen" objects to combine and the specification of 
#' the data frames to combine. 
#' Can be any of the following(s): P","G", "E", "S", "C", or "ALL" (default). 
#' If a "G" data frame is provided, the program also generates 
#' the A slot coding the missing data as "0" in default option (see the
#' argument "missing").
#' The XY slot is generated automatically if present.
#' @param missing Missing data manipulation.
#' It can take three values ("0" ,"NA" or "MEAN"- i.e, the mean frequency
#' of the corresponding allele). 
#' Missing elements are coded as 0 in the default option.
#' @examples
#' 
#' \dontrun{
#' 
#' data(eco.test)
#' eco.example <- eco.cbind(eco,eco,"ALL")
#' eco.example
#' eco.example2 <- eco.cbind(eco, eco,"P", "G", missing="NA")
#' eco.example2
#' 
#' }
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @export


setGeneric("eco.cbind", 
           function(eco1, eco2, ..., 
                    missing = c("0", "MEAN", "NA")) {
             
             
             #--GENERAL CONFIGURATION----------------------#
             
             missing <- match.arg(missing)
             
             # unlist dots
             u <- unlist(list(...))
             
             # ecogen objects
             u.ecogen <- u[sapply(u, is.ecogen)]
             # all ecogen objects
             u.ecogen <- c(eco1, eco2, u.ecogen)
             
             
             # character objects-----------------
             u.char <- u[sapply(u, is.character)]
             vec <- c("P", "G", "E", "S", "C", "ALL")
             m <- vec %in% u.char
             
             if((m[6] == TRUE) | !any(m)) { # if ALL 
               m <- rep(TRUE, 5)
             }
             
             
             #--INT.CBIND FUNCTION----------------------#
             #------------------------------------------#
             
             int.cbind <- function(e1, e2) {
              
             z <- ecogen()
             # create separed lists-----------------------
             z1 <- list(e1@P, e1@G, e1@E, e1@S, e1@C)
             z2 <- list(e2@P, e2@G, e2@E, e2@S, e2@C)
             
             tem <- list()
             
             for(i in 1:5) {
               
               if(m[[i]]) {
                 # check row number-------
                 a <- nrow(z1[[i]])
                 b <- nrow(z2[[i]])
                 
                 # if any of both data frames is empty...
                 if(any(a,b) == 0) 
                 {
                   # both data frames empty, of only one of both non empty.
                   if(a == 0 && b == 0) {
                     tem[[i]] <- data.frame()
                   } else if(a == 0 && b != 0) {
                     tem[[i]] <- z2[[i]]
                   } else if(a != 0 && b == 0) {
                     tem[[i]] <- z1[[i]]
                   }
                   # both non empty. 
                 } else {
                   # check first row names consistency.
                   # if different names present, the program generates 
                   # an empty data frame.
                   if(any(rownames(z1[[i]]) != rownames(z2[[i]]))) {
                     warning(paste("Individuals in",
                                   paste("<", vec[i], ">", sep = ""),  
                                   "data frame do not have the same rownames.
                                   This will generate an empty slot."))
                     next
                   }
                   # bind both data frames
                   tem[[i]] <- cbind(z1[[i]], z2[[i]])
                 }
                 }
             }
             
             # fill XY slot--------------------------------
             # this rule retains the first non empty matrix.
             # XY consistency should be checked across paired objects?
             # one posibility is to use row names for checks in the future
             # and standard (unique) column names (X, Y, Z)
             if(nrow(e1@XY) != 0) {
               z@XY <- e1@XY
             } else {
               if(nrow(e2@XY) != 0) {
               z@XY <- e2@XY
               }
             } 
             
             
             # fill P slot-----
             if(m[1] == TRUE) {
               z@P <- tem[[1]]
             }
             
             
             # fill G and A slots-----------------------------
             G.cond <- m[2] == TRUE && all(dim(tem[[2]]) != 0) # fill G condition

             if(G.cond) { 

                 # ploidy and ncod control
                 cont1 <- e1@INT@ploidy != e2@INT@ploidy 
                 cont2 <- e1@INT@ncod != e2@INT@ncod
                 cont3 <- e1@INT@type != e2@INT@type
                 cont <- cont1 || cont2 || cont3
                 if(cont) {
                   warning("incongruence in the ploidy, number of digits
                           per allele or type of data (dominant/codominant).
                           Genetic slots will be empty. Please check your
                           data.")
                   
                   z@G <- data.frame()
                 
               } else {  # e1 and e2 are consistent
                 
                 # fill G slot--
                 z@G <- tem[[2]]
                 
                 # create an int.genind temporal object
                 tempo <- int.df2genind(tem[[2]], 
                                        missing = missing,
                                        ncod = e1@INT@ncod,
                                        ploidy = e1@INT@ploidy,
                                        type = e1@INT@type)
                 
                 
                 # fill A and the internal slot INT-------
                 z@A <- data.frame(tempo@tab)
                 z@INT@loc.fac <- tempo@loc.fac
                 z@INT@all.names <- tempo@all.names
                 z@INT@ploidy <- tempo@ploidy
                 z@INT@type <- tempo@type
                 z@INT@NA.char <- ifelse(e1@INT@NA.char == e2@INT@NA.char,
                                         e1@INT@NA.char, "NA")
                 z@INT@sep <- ifelse(e1@INT@sep == e2@INT@sep,
                                     e1@INT@sep, "")
                 z@INT@ncod <- tempo@ncod
                 z@INT@missing <- tempo@missing
                 
                 # missing data position is additive for cbind
                 z@INT@missing.cells <- as.integer(c(e1@INT@missing.cells, 
                                          e2@INT@missing.cells + 
                                            length(e1@A)))
                 z@INT@removed.image <- tempo@removed.image
                 
               }
             }
               
               # fill E, S and C slots-------
               if(m[3] == TRUE) {
                 z@E <- tem[[3]]
               }
               
               if(m[4] == TRUE) {
                 z@S <-tem[[4]]
               }
               
               if(m[5] == TRUE) {
                 z@C <- tem[[5]]
               }
             
             return(z)
             }
             #-----END INT.CBIND----------------------#
             
             #-----OUTPUT CREATION--------------------#
             # bind multiple objects using recursion
             
             len.eco <- length(u.ecogen)
             out <- u.ecogen[[1]]
             i <- 2
             while(i <= len.eco) {
             out <- int.cbind(out, u.ecogen[[i]])
             i <- i + 1
             }
               
               return(out)
             
               })

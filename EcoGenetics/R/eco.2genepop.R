#' Exporting an ecogen genetic data frame into Genepop format
#' 
#' @description This function converts the genetic 
#' data of an ecogen object into a Genepop input file. 
#' @param eco Object of class "ecogen".
#' @param grp The name of the S slot column with groups in which the sample
#' must be divided (e.g., populations). If groups are not given (grp = NULL),
#' all individuals will be assigned to a single one.
#' @param nout Number of digits in the output file
#' @param name The name of the output file.
#' @param sep Character separating alleles.
#' @return A Genepop file in the working directory.
#' @examples 
#' 
#' \dontrun{
#' 
#' data(eco.test)
#' eco.2genepop(eco, grp = "pop", name = "infile.genepop.txt")
#' # an output file "infile.genepop.txt" is generated in the working directory
#' 
#' }
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @export


setGeneric("eco.2genepop", 
           function(eco, name = "infile.genepop.txt",
                    grp = NULL, nout = 3, sep = "") {
             
             nout <- ceiling(nout)
             nout <- c(1,2,3)[c(1,2,3) %in% nout]
             if(length(nout) != 1)  {
               stop("nout must be 1, 2 or 3")
             }
             #check codes
             ncod <- int.check.ncod(eco@G, ploidy = eco@INT@ploidy, sep = sep)
             #check group consistency
             structures <- int.check.group(eco@S, grp = grp, exp.l = nrow(eco@G))
             structures <- as.factor(as.numeric(structures)) #recoding levels
             
             
             X <- eco.format(eco@G, ncod = eco@INT@ncod, 
                             nout = nout,
                             ploidy = eco@INT@ploidy,
                             fill.mode = "first")
             
             X <- cbind(rep(0, nrow(X)), X)
             X[, 1] <- paste(rownames(X), ",")
             
             lista <- list()
             grp <- rep(" ", ncol(eco@G) + 1)
             grp <- as.matrix(t(grp))
             grp[1] <- "POP"
             maxf <- max(levels(structures))
             matriz <- matrix(nrow = 0, ncol = ncol(eco@G)+1)
             for(i in 1:maxf) {
               lista[[i]] <- X[structures == i, ]
               lista[[i]] <- rbind(grp, lista[[i]])
               matriz <- rbind(matriz, lista[[i]])
             }
             matriz <- rbind(matriz, rep("", ncol(matriz)))
             
             matriz[, 1] <- as.character(matriz[, 1])
             nombres <- rep("", (ncol(X)) ^ 2)
             nombres <- as.data.frame(matrix(nombres, ncol(X), ncol(X)))
             nombres[, 1] <- c("Data exported from EcoGenetics", colnames(X[, -1]))
             colnames(nombres) <- colnames(matriz)
             matriz <- rbind(as.matrix(nombres), matriz) 
             
             write.table(matriz, name, row.names = FALSE,
                         col.names = FALSE, quote = FALSE)
             
           })

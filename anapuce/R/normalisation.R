normalisation <-
function (fileIN=NULL,Red = "F635.Median", Green = "F532.Median", n = 3, 
    flag = -100, graph = TRUE,  filter.function = filterByDefault, 
    filter.fic = NULL, filter.var = NULL,sep.write = "\t", 
    dec.write = ".", header =  TRUE, sep = "\t", skip = 0, ...) 
{
    cat("####################################### \n")
    cat("NORMALIZATION STEP \n")
    nom <- FileChoice(fileIN)
    write.table(nom, "AnalyzedArrays.txt", row.names = FALSE, 
        col.names = "Arrays", sep = sep.write, dec = dec.write)
  
    cat("Variables used in the analysis :", Red, " and ", Green,"\n")
    for (i in 1:length(nom)) {
        cat(as.character(nom[i]), "\n")
        if (graph == TRUE) 
            pdf(file = paste(as.character(nom[i]), ".pdf", sep = ""))
        if (length(grep(".gpr", as.character(nom[i]))) == 1) {
                skip <- as.integer(read.table(as.character(nom[i]), 
                nrow = 1, skip = 1,sep=sep)[1]) + 
                2
        }
       array1 <- read.table(as.character(nom[i]), header = header, sep=sep, skip=skip, ...)

        array1 <- array1[, which(colnames(array1) %in% c(colnames(array1)[1:n], 
            Red, Green, "Flags"))]
             if (i == 1) {
                    resLowess <- array1[, 1:n]
                    is.Block=TRUE   
                    if(length(grep("Block",names(resLowess)))==0)
                       is.Block=FALSE
             }       
        if (!is.null(filter.function)) {
            if (all.equal(filter.function, filterByDefault) == 
                TRUE) {             
                array1 <- filter.function(array1, flag, n, filter.fic = filter.fic, 
                  filter.var = filter.var, sep = sep, sep.write = sep.write, 
                  dec.write = dec.write, ...)#modif
            } else {
                array1 <- filter.function(array1, flag, n)
            }
        } else {
            if (!is.null(flag)) {
                cat("Flags with value(s)", flag, "will not be removed.\n Verify the filter.function argument \n")
            }
          }   
# mednull is renamed ind.signal0 
        ind.signal0 <- which((array1[, match(Red, colnames(array1))] == 
            0) | (array1[, match(Green, colnames(array1))] == 0))
        if (length(ind.signal0) > 0) {
            cat("Number of genes with at least one signal value equal to 0 : ", 
                length(ind.signal0), "\n")
            print(array1[ind.signal0,])      
            array1 <- array1[-ind.signal0, ]
        }
        A <- 0.5 * (log2(array1[, match(Red, colnames(array1))]) + 
            log2(array1[, match(Green, colnames(array1))]))
        M <- log2(array1[, match(Red, colnames(array1))]/array1[, 
            match(Green, colnames(array1))])
        ordre <- order(A)
        yfit1 <- lowess(A, M, f = 0.3)$y[order(ordre)]
        Mcorrec<-M-yfit1
        tmp <- data.frame(array1[,1:n],A, Mcorrec, yfit1)
        names(tmp)<-c(names(array1[,1:n]),"A","M.norm","corLowess")
        if (graph) {
            plot(A, M, xlab = "A", ylab = "M", main = "M-A plot : raw data", 
                pch = 20, cex = 0.5)
            lines(lowess(A, M, f = 0.3), col = 2)
            plot(A, (M - yfit1), main = "M-A plot after a global lowess normalization", 
                xlab = "A", ylab = "M", pch = 20, cex = 0.5)
            lines(lowess(A, (M - yfit1), f = 0.3), col = 3)            
            if(is.Block){
               boxplot(M ~ array1$Block, main = "Boxplot : raw data", 
                xlab = "Block", pch = 20, cex = 0.5)
               boxplot(Mcorrec ~ array1$Block, main = "Boxplot after a global lowess normalization", 
                xlab = "Block", ylab = "M", pch = 20, cex = 0.5)
                 }
           }         
         if (is.Block){ 
                 tmp <- data.frame(tmp, corBlock=rep(0, dim(tmp)[1]))
                 for (bloc in 1:max(tmp$Block)) tmp$corBlock[tmp$Block == 
                 bloc] <- median(tmp$M.norm[tmp$Block == bloc])
                 tmp$M.norm <- tmp$M.norm - tmp$corBlock
                 if (graph){
                    plot(tmp$A, tmp$M.norm, xlab = "A", ylab = "M", main = "M-A plot after the normalization step", 
                    pch = 20, cex = 0.5)
                    lines(lowess(tmp$A, tmp$M.norm, f = 0.3), col = 2)
                    boxplot(tmp$M.norm ~ tmp$Block, main = "Boxplot after the normalization step", 
                    xlab = "Block", ylab = "M", pch = 20, cex = 0.5)
                 }
          }   
        if (is.Block) {names(tmp) <- c(names(array1)[1:n], paste(c("A","M.norm", "corLowess", "corBlock"), i, sep = ""))
        } else { 
        names(tmp) <- c(names(array1)[1:n], paste(c("A","M.norm", "corLowess"), i, sep = ""))
        }
        resLowess <- merge(resLowess, tmp, all = TRUE, by = names(tmp)[1:n], sort = FALSE)  
        if (graph)  dev.off()
    }
    write.table(resLowess, file = "resNorm.txt", quote = TRUE, 
        sep = sep.write, row.names = FALSE, dec = dec.write)  
    # (c) 2007 Institut National de la Recherche Agronomique

}


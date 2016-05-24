print.MacroTxChrono <-
function (x, file = NULL, sep = ";", ...) 
{
    res<- x
    if (!inherits(res, "MacroTxChrono")) 
        stop("non convenient data")
    cat("**Results for the MacroTxChrono**\n")
    indice <- 9
    res <- array("", c(indice, 2), list(1:indice, c("name", "description")))
    res[1, ] <- c("$SentenceList", "dataset with  sentences")   
    res[2, ] <- c("$Homo.Groups", "Description homogeneous group")
    res[3, ] <- c("$Corpus", "Description of corpus")
    res[4, ] <- c("$Correlation", "Correlation between chronology and dimensions")
    res[5, ] <- c("$res.TxCA", "Results of correspondence analysis")
    res[6, ] <- c("$res.chcpc", "Results for the Constrained hierarchical clustering")
    res[7, ] <- c("$HierWord", "Characteristic words for every node of the hierarchy")
    res[8, ] <- c("$HierSegment", "Characteristic Segments for every node of the hierarchy")
    res[9, ] <- c("$VocIndex", "Vocabulary index")
    print(res[1:indice, ])
    if (!is.null(file)) {
        write.infile(MacroTxChrono, file = file, sep = sep)
        print(paste("All the results are in the file", file))
    }
}




    
      
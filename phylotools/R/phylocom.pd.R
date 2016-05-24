#### Function phylocom.pd as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010

phylocom.pd <-
function (sample = "sample", phylo = "phylo") 
{
    files <- list.files()
    if (!sample %in% files) {
        stop(paste("Sample file named", "\"", sample, "\" 
                 \n  can not be found in the working directory."))
    }
    
    if (!phylo %in% files) {
        stop(paste("Phylo file named", "\"", phylo, "\" 
                  \n can not be found in the working directory."))
    }
    
    command <- paste("phylocom pd -s", sample, " -f", phylo)
    res <- system(Sys.which(command),intern = TRUE,show.output.on.console = FALSE, ignore.stderr = TRUE)
    res2 <- unlist(lapply(strsplit(res[2:length(res)], "\t"), 
                       FUN = function(x)gsub(" ", "", x)))
                       
    res.names <- unlist(strsplit(res[1], "\t"))
    dim(res2) <- c(5, length(res2)/5)
    res2 <- t(as.data.frame(res2))
    
    if(nrow(res2) == 1){
       res22 <- t(as.data.frame(as.numeric(res2[,-1])))
       rownames(res22) <- res2[,1]
       colnames(res22) <- res.names[-1]
       return(res22)
    }
    
    res22 <- res2[,2:ncol(res2)]
    res22 <- apply(res22, 2, as.numeric)
    rownames(res22) <- res2[,1]
    colnames(res22) <- res.names[-1]
    return(res22)
}


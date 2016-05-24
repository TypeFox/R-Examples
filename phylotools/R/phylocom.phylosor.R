#### Function phylocom.phylosor as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010


phylocom.phylosor <- function(sample.file = "sample", phylo = "phylo")
{
    # sample.file = "sample"; 
    # phylo = "phylo"
    ## 
    samp <- read.table(sample.file, header = FALSE)
    cat("sample data has been successfully read into memory.\n")
    samp[,1] <- as.character(samp[,1])
    samp[,3] <- as.character(samp[,3])
    
    names.all <- unique(as.character(samp[,1]))    ## Names of the plots
    res <- rep(0, length(names.all) * length(names.all))
    dim(res) <- c(length(names.all) ,length(names.all))
    snames.all <- sort(names.all)
    colnames(res) <- snames.all
    rownames(res) <- snames.all
   ## cat("Making combinations names.\n")     
    
   for(i in 2:length(names.all)){
       for(j in 1:(i-1)){
           selected <- c(snames.all[i], snames.all[j])
           loc1 <- samp[samp[,1] %in% selected, ]
           
           ## obtain PD for each site
           write.table(loc1, "temp_pd_ind", row.names = FALSE, quote = FALSE, col.names = FALSE)
           pd.ind <- phylocom.pd(sample = "temp_pd_ind", phylo = phylo)
           unlink("temp_pd_ind")
           pd.sum <- sum(pd.ind[,2])
           
           ## obtain the sum of PD
           loc <- data.frame(rep(paste(snames.all[i], snames.all[j], sep = "_"),nrow(loc1)), loc1[,-1])
           colnames(loc) <- c("location", "abund", "species")
           loc <- unique.data.frame(loc)
           write.table(loc, "temp_pd_merged", row.names = FALSE, quote = FALSE, col.names = FALSE)
           pd.merg <- phylocom.pd(sample = "temp_pd_merged", phylo = phylo)
           unlink("temp_pd_merged")
           
           res[i,j] <-  (pd.sum -(pd.merg[1,2]))/(0.5*(pd.sum))
       }
   }
  return(as.dist(res))
}


#### Function phylocom.comstruct as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010

phylocom.comstruct <-
function(sample = "sample", phylo = "phylo", swapmethod = 2, 
          runs = 999, swaps = 1000, aw = TRUE, file = "comstruct.txt")
{
    files <- list.files()
    if(!sample%in%files){
        stop(paste("Sample file named", "\"", sample,"\" 
              can not be found in the working directory."))
    }
    if(!phylo%in%files){
        stop(paste("Phylo file named", "\"", phylo,"\" 
               can not be found in the working directory."))
    }
    a = NULL
    if(aw){
       a = "-a"
    }
    command <- paste("phylocom comstruct -r", runs, "-m", swapmethod,
                     "-w", swaps, "-s", sample, "-f", phylo, a,">" ,file)
    shell(Sys.which(command), intern = FALSE, mustWork = NA)
    cs <- read.table(file, header = TRUE, row.names = 1, skip = 1, sep="\t")
    return(cs)
}


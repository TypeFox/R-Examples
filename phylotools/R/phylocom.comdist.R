#### Function phylocom.comdist as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010

phylocom.comdist <-
function(sample = "sample", phylo = "phylo", method = "comdist",aw = TRUE, file = "comdist.txt"){
    files <- list.files()
    if(!sample%in%files){
        stop(paste("Sample file named", "\"",sample,"\" 
              can not be found in the working directory."))
    }
    if(!phylo%in%files){
        stop(paste("Phylo file named", "\"",phylo,"\" 
             can not be found in the working directory."))
    }
     
    if(aw){
        a = "-a"
    } else{
        a = NULL
    }
    
    if(method == "comdistn") {
        pc <- "comdistnn"
    } else {
        pc <- "comdist"
    }
    
    shell(Sys.which(paste("phylocom", pc, "-s",sample, "-f", phylo, a,"> ", file )), 
           intern = FALSE, mustWork = NA)
    cd <- read.table(file, header = TRUE, row.names=1, sep="\t")
    as.dist(cd)
}


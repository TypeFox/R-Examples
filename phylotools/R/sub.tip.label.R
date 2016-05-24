#### Function sub.tip.label as part of R package phylotools
#### By Jinlong Zhang  <Jinlongzhang01@gmail.com>
#### Institute of Botany, the Chinese Academy of Sciences, Beijing ,China
#### Nov- 01-2010

sub.tip.label <- 
function(tree, dat){
    if(!inherits(tree,"phylo")){
	   stop("the input tree in is not a \"phylo\" class object.")
	}
	if(!is.data.frame(dat)){
	   stop("the input dat is not a \'dataframe\'.")
	}
    tree2 <- tree
    nnn <- tree$tip.label
	if(!nrow(dat)==length(nnn)){
	   warning("The number of tip labels in phylogenetic\n
          	   tree differ from the reference table.\n")
	}
    xxx1 <- as.character(dat[,1])
    xxx2 <- as.character(dat[,2])
	if(!all(xxx1%in%nnn)){
	   unsub.dat <- xxx1[!xxx1%in%nnn]
	   cat(length(unsub.dat),"tip label(s)",unsub.dat, "can not\n
	       be found in the reference table.\n")
	}
	
	if(!all(nnn%in%xxx1)){
		unsub.tree <- nnn[nnn%in%xxx1]
	   cat(length(unsub.tree),"Names",unsub.tree, "in reference table\n
     	   can not be found in the tree.\n")
	}
	
    label <- c()
    for(i in 1:length(dat[,1])){
        for(j in 1:length(nnn)){
            if(nnn[j] ==xxx1[i]){
                label[j] <- xxx2[i]
            }
        }
    }
    tree2$tip.label <- label
    return(tree2)
}

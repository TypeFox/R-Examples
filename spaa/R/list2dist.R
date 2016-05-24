list2dist <-
function(dat){
    dat.name1 <- as.character(dat[,1])
    dat.name2 <- as.character(dat[,2])
    dat.value <- dat[,3]
    names1 <- sort(unique(as.character(dat[,1])))
    names2 <- sort(unique(as.character(dat[,2])))
    total.names <- unique(c(names1, names2))
    elements <- rep(NA, length(total.names)^2)
    dim(elements) <- c(length(total.names),length(total.names))
    rownames(elements) <- total.names
    colnames(elements) <- total.names
    
    for(i in 1:length(total.names)){
         for(j in 1:length(total.names)){
              for(k in 1:length(dat.name1)){
                  if((total.names[i] == dat.name1[k])&(total.names[j] == dat.name2[k])){
                      elements[i,j] <- dat.value[k]
                    }
              }
         }
    }
    res <- as.dist(t(elements))
    return(res)
}


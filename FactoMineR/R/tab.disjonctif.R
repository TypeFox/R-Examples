tab.disjonctif<-function (tab){
    tab<-as.data.frame(tab)
    #fonction interne permettant la réalisation d'un TDC pour un unique facteur
    modalite.disjonctif <- function(i){
        moda <- as.factor(tab[, i])
        n <- length(moda)
        x <- matrix(0, n, nlevels(moda))
        x[(1:n) + n * (unclass(moda) - 1)] <- 1
#        nom <- attributes(tab)$names[i]
#        if((ncol(tab)!=1)&(levels(moda)[1]%in%c(1:nlevels(moda),"n","N","y","Y"))) dimnames(x) <- list(attributes(tab)$row.names, paste(nom, levels(moda),sep = "."))
#        else dimnames(x) <- list(attributes(tab)$row.names, levels(moda))
        return(x)
    }
    # fin fonction interne

    if (ncol(tab)==1) {
	  res <- modalite.disjonctif(1)
	  dimnames(res) <- list(attributes(tab)$row.names, levels(tab[,1]))
	}
    else
    {
	  variable <- rep(attributes(tab)$names,sapply(tab,nlevels))
	  listModa <- unlist(sapply(tab,levels))
	  wlistModa <- which((listModa)%in%c("y","n","Y","N"))
      if (!is.null(wlistModa)) listModa[wlistModa] <- paste(variable[wlistModa],listModa[wlistModa],sep = ".")
      numlistModa <- which(unlist(lapply(listModa,is.numeric)))
      if (!is.null(numlistModa)) listModa[numlistModa] <- paste(variable[numlistModa],listModa[numlistModa],sep = ".")
      res <- lapply(1:ncol(tab), modalite.disjonctif)
      res <- as.matrix(data.frame(res, check.names = FALSE))
	  dimnames(res) <- list(attributes(tab)$row.names,listModa)
	}
    return(res)
}

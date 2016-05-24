write.linkdat = function(x, prefix="", what=c("ped", "map", "dat", "freq", "model"), merlin=FALSE) {
    generated.files = character(0)
    if(merlin) {
        pmatr = as.matrix(x)
        markerattr = attr(pmatr, 'markerattr')
    }
   
    if(any(c("map", "dat", "freq") %in% what)) 
        .map = .getMap(x, na.action=1, verbose=F)
    
    if("ped" %in% what) {
        if(!merlin) 
            pmatr = cbind(x$famid, relabel(x$pedigree, x$orig.ids), .prettyMarkers(x$markerdata, singleCol=FALSE))
        write(t(pmatr), pedname <- paste(prefix, "ped", sep="."), ncolumns=6+2*x$nMark)
        generated.files = c(generated.files, pedname)
    }
    
    if("model" %in% what && !is.null(x$model)) {
        dfreq = format(x$model$dfreq, scientific=F, decimal.mark=".")
        # if X-linked: using the female penetrances. This works in MERLIN/MINX as long as male and female penetrances are equal. 
        penetrances = if(x$model$chrom == "X") x$model$penetrances$female else x$model$penetrances
        penets = paste(format(penetrances, scientific=F, decimal.mark="."), collapse=",")
        .model = c("my_disease", dfreq, penets, "my_model")
        write(.model, modelname <- paste(prefix, "model", sep="."), sep=" \t", ncolumns=4)
        generated.files = c(generated.files, modelname)
    }
    
    if("map" %in% what) {
        write.table(.map, mapname <- paste(prefix, "map", sep="."), col.names=F, row.names=F, quote=F)
        generated.files = c(generated.files, mapname)
    }    
    
    if("dat" %in% what) {
        .dat = cbind(code=c("A", rep("M", nrow(.map))), value=c("my_disease", .map$MARKER))
        write.table(.dat, datname <- paste(prefix, "dat", sep="."), col.names=F, row.names=F, quote=F)
        generated.files = c(generated.files, datname)
    }
    
    if("freq" %in% what) {
        if(!merlin)
            markerattr = lapply(x$markerdata, attributes)
        nalls = unlist(lapply(markerattr, function(at) at$nalleles))
        L = sum(nalls) + length(nalls)
      
        cum = cumsum(c(1, nalls+1))
        length(cum) = length(nalls) #remove last
      
        col1 = rep("A", L)
        col1[cum] = "M"

        col2 = character(L)
        col2[cum] = .map$MARKER
        if(merlin)
            allalleles = unlist(lapply(nalls, seq_len))
        else
            allalleles = unlist(lapply(markerattr, function(at) at$alleles))
        col2[-cum] = allalleles

        col3 = character(L)
        allfreqs = unlist(lapply(markerattr, function(at) at$afreq))
        col3[-cum] = format(allfreqs, scientifit=F, digits=6)

        .freq = cbind(col1, col2, col3)
        write.table(.freq, freqname <- paste(prefix, "freq", sep="."), col.names=F, row.names=F, quote=F)
        generated.files = c(generated.files, freqname)
    }   
   
    invisible(generated.files)
}

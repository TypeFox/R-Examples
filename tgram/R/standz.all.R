standz.all <-
function(traq, series, wl =NULL, w.char=NULL, order = NULL, G=30){

    # "traq" is a vector with the ordered sequences of measurements for each traqueidogram
    # "series" is a vector of indicator values (i.e. a factor) with each level indicating
    #            each unique tracheidogram
    # "wl" is a vector indicating if the measurement is wall or lumen
    # "w.char" is the character indicating "wall"
    # "order" is a vector indicating the ordering of each measurement in each lumen or wall series within a tracheidogram
    
    
    which.w <- NULL
    which.l <- NULL
    
    #indication of unique series of lumen or wall series 
    #first, process vector of lumen.wall measurements (if provided)
    if(!is.null(wl)) {
if(is.null(w.char)) stop("you should indicate w.char")
if (is.na(match(w.char, unique(wl)))) stop (paste("I cant'find any '",w.char,"' in the ",deparse(substitute(wl))," vector. Please correct the w.char argument",sep="")) 
series <- paste(series,wl, sep="###_")

#trick to identify which series are lumen and which are wall
serie.unique <- unique(series)
which.w <- grep(paste("###_",w.char, sep=""),serie.unique)
                which.l <- (1:length(serie.unique))[-which.w]
series <- sub("###","",series)
     }
    # then find unique series of measuremnts
    serie.unique <- unique(series)
        
    pepe <- NULL
    #for (i in 1:length(secu)){
      for (i in 1:length(serie.unique)){
      #print(i)
tgl <-traq[series==serie.unique[i]]
        pepe <- rbind(pepe, standz(tgl1=tgl, G=G))
}
row.names(pepe)<- serie.unique
     
    result <- list(data.stdz=pepe, which.w=which.w, which.l=which.l)
    class(result)=c("standz.all", class(result))
    return(result)
    
    
}


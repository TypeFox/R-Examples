`diagnostic_window` <-
function(position, window.size,measured,modelled, use_qualV=FALSE, diff.ecdf = NA ){
    if(position %% 500==0){
           cat("At position", position, "\n")
           gc()
    }
         wind<-position:(position+window.size-1)
         w.dat<-measured[wind]
         #dim(w.dat)<-c(1,window.size)
         toReturn<-diagnostic_dawson(measured=w.dat, modelled=modelled[wind], use_qualV=use_qualV, diff.ecdf=diff.ecdf)
         return(toReturn)

}


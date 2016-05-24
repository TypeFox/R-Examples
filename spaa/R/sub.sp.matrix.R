sub.sp.matrix <-
function(spmatrix, freq = 0.5, common = NULL){  
   if (!is.null(common)){ 
        freq = NULL
        if ((!is.null(freq))&(!is.null(common))){ 
		    stop("Contradiction in argument freq and common, 
		         try to specify only one argument.\n") 
		}
        if( common > ncol(spmatrix)){ 
		    stop("There are fewer species than you want to subset.\n") 
		 }
        freq.spmatrix <- freq.calc(spmatrix)
        names.spmatrix <- colnames(spmatrix)
        sort.freq <- sort(freq.spmatrix, decreasing = TRUE)
        names.freq <- names(sort.freq)
        nn <- names.freq[1:common]
        wii <- c()
        for(i in 1:length(nn)){ 
		    wii <- append(wii, which(nn[i] == names.spmatrix)) 
	    }
        submatrix <- spmatrix[,wii]
        }else{
             if (!is.null(freq)){
                 freq.spmatrix <- freq.calc(spmatrix)
                 sel <- freq.spmatrix > freq
                 base <- 1:length(sel)
                 se <- base[sel]
                 submatrix <- spmatrix[,se]
             }
        }
return(submatrix)
}


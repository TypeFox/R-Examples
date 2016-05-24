

Lifespan <-
function (ch) 
{


Collapse.CH = function(ch = ch$ch_split) {
            ch_vec = apply(ch, 1, paste, collapse = "")
            ch = data.frame(ch = ch_vec)
            ch["ch"] = as.character(ch$ch)
            return(ch)  
}

if(is.matrix(ch)) ch = Collapse.CH(ch)

	
	if(is.vector(ch) || is.matrix(ch)) ch = data.frame(ch)
		
	if (!(is.data.frame(ch) || is.matrix(ch))) 
        stop("WARNING: ch should be a data.frame, matrix, or vector")
    if (is.data.frame(ch)) {
        if (ncol(ch) > 1) 
            stop("WARNING: if using a data.frame, it should be a single column")
        classes = apply(ch, 2, class)
        if (!(any(classes == "character"))) 
            stop("WARNING: ch in the data.frame should be a character string of ones and zeros, check to make sure it is not a factor")
    }
        

	Lifespan = c()
    for (i in 1:nrow(ch)) {
     
	 if(is.na(ch[i, 1])) {
	 Lifespan[i] = NA
	 next
	 }
		chstring = as.character(ch[i, 1])
        idx = str_locate_all(chstring, "1")
        DaysAlive = max(unlist(idx)) - min(unlist(idx))
        Lifespan[i] = DaysAlive
    }
    return(Lifespan)
}




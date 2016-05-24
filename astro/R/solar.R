solar = function(band = "r", vega = FALSE, source = FALSE){
    
    # definitions table
    tab = rbind(
        c("u", 6.38, 5.47, "Hill et al., 2010 (Table 1)"), 
        c("g", 5.15, 5.23, "Hill et al., 2010 (Table 1)"), 
        c("r", 4.71, 4.55, "Hill et al., 2010 (Table 1)"), 
        c("i", 4.56, 4.19, "Hill et al., 2010 (Table 1)"), 
        c("z", 4.54, 4.00, "Hill et al., 2010 (Table 1)"), 
        c("Y", 4.52, 3.89, "Hill et al., 2010 (Table 1)"), 
        c("J", 4.57, 3.63, "Hill et al., 2010 (Table 1)"), 
        c("H", 4.71, 3.33, "Hill et al., 2010 (Table 1)"), 
        c("K", 5.19, 3.29, "Hill et al., 2010 (Table 1)"), 
        c("U", 6.36, 5.55, "http://mips.as.arizona.edu/~cnaw/sun.html"), 
        c("B", 5.36, 5.45, "http://mips.as.arizona.edu/~cnaw/sun.html"), 
        c("V", 4.82, 4.80, "http://mips.as.arizona.edu/~cnaw/sun.html")
    )
    
    # find band
    if(band%in%tab[,1]){
        
        # find row
        row = which(tab[,1]==band)
        
        # find absolute solar magnitude
        if(vega){
            res = as.numeric(tab[row,3])
        }else{
            res = as.numeric(tab[row,2])
        }
        
        # add source
        if(source){
            res = c(res, tab[row,4])
        }
        
        # return results
        return(res)
        
    }else{
        
        # unable to find band
        cat("ERROR: requested band '", band, "' not found in 'astro' catalogue.\n", sep="")
        
    }
    
}


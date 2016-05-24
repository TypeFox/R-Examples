# Purpose        : Convert an attribute table associated with spatial data into an HTML <description> bubble
# Maintainer     : Pierre Roudier (pierre.roudier@landcare.nz);
# Contributions  : Dylan Beaudette (debeaudette@ucdavis.edu); Tomislav Hengl (tom.hengl@wur.nl);
# Status         : ready for R-forge
# Note           : required by KML writing functions;

.df2htmltable <- function(x, fix.enc = TRUE, columns = TRUE) {

    # if the user passed in TRUE, then we want all of the columns
    if(class(columns) == 'logical'){
        columns <- 1:ncol(x) # use all columns from the dataframe row
    }
    else{     # otherwise, keep only requested columns
    x <- x[, columns]
    }                  
    
    # fix encoding:
    if(fix.enc==TRUE){
    x <- data.frame(lapply(x, iconv, to = "UTF8"))
    }
    
    # get selected table data:
    att.names <- sapply(names(x), function(i) { paste('<span style="font-weight: bold; color: #000000; padding: 3px;">', as.character(i), '</span>:&nbsp;', sep = '') } )    
    att.values <- as.vector(t(sapply(x, function(i) { paste('<span style="color: #000000; padding:3px;">', as.character(i), '</span><br>', sep = '') })))
    # combine by interleaving:
    att <- matrix(paste(att.names, att.values, sep="\n"), ncol=length(names(x)), byrow=TRUE)
    html.table <- apply(att, 1, paste, collapse="\n")
    
    return(html.table) 
}

# end of script;

###########################################################################
# netClass.r
#
# Authors: 	Shannon M. Bell
#          	ORISE Fellow
#			Environmental Protection Agency 
#			National Health and Human Effects Research Laboratory
#			Research Triangle Park, NC
#
# Version 1.0	August 16, 2012
#
# Notes: SMB conceived of and wrote the algorithm.
# This algorithm is written to take an class list (ex, tree object) and convert it
# into an adjacency matrix.  
#################################################

netClass<-function(x, labels=NULL){
   
    y<-.Call("netVal", as.numeric(x), pkg="NetComp")
    if(is.null(labels)){    
        y
    }
    else{
        row.names(y)<-labels; colnames(y)<-labels
        y
    }
}
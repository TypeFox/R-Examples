# Get vice-county (VC) ID using VC name
#
# Given the name of a VC in the NBN (e.g. 'cambridgeshire'), return its ID key.
# The ID key can then be used in searchs and filters.
# This is currently an internal function

getVCid<-function(name){
    
    if(class(name) != 'character'){
        warning('VC name non-character, conversion attempted')
        name <- as.character(name)
    }
    
    VCDB<-listVCs()
    VCDB$name<-tolower(VCDB$name)
        
    name<-tolower(name)
    
    if(name %in% VCDB$name){
        VCID<-as.character(VCDB$identifier[VCDB$name==name])
        return(VCID)
    } else {
        stop(paste(name,'is not a recognised vicecounty. Use listVCs to find valid names'))
    }
}
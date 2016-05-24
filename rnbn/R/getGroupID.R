# Get group ID using group name
# 
# Given the name of a group in the NBN (e.g. 'reptile'), return its ID key.
# The ID key can then be used in searchs and filters
# This is currently an internal function

getGroupID<-function(name){
    
    if(class(name) != 'character'){
        warning('Group name non-character, conversion attempted')
        name <- as.character(name)
    }
    
    groupsDB<-listGroups()
    groupsDB$name<-tolower(groupsDB$name)
        
    name<-tolower(name)
    
    if(name %in% groupsDB$name){
        groupID<-groupsDB$key[groupsDB$name==name]
        return(groupID)
    } else {
        stop(paste(name,'is not a recognised group. Use listGroups to find valid names'))
    }
}
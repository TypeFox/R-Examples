getname <-
function(ccode,dictionary){as.character(dictionary[match(ccode,dictionary[,2]),1])}

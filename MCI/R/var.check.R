var.check <-
function(x, asdummy=FALSE) {
if (class(x) == "character")   
{
if (asdummy==FALSE) {   
return("invalid_s")
}
else {  
return("asdummy")
}
}

if (class(x) == "factor")   
{
if (asdummy==FALSE) {   
return("invalid_s")
}
else {  
return("asdummy")
}
}

if (class(x) == "integer") {  
if (length(unique(x)) > 2) { 
if (all(x>0) == TRUE) { return ("valid_n") }
else { return ("invalid_zn") }
}
else { return("valid_d") }
}

if (class(x) == "numeric") {  
if (length(unique(x)) > 2) { 
if (all(x>0) == TRUE) { return ("valid_n") }
else { return ("invalid_zn") }
}
else { return("valid_d") }
}
}

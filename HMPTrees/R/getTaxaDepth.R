getTaxaDepth <-
function(level="genus"){
if(tolower(level) == "kingdom" || tolower(level) == "k"){
return(1)
}else if(tolower(level) == "phylum" || tolower(level) == "p"){
return(2)
}else if(tolower(level) == "class" || tolower(level) == "c"){
return(3)
}else if(tolower(level) == "order" || tolower(level) == "o"){
return(4)
}else if(tolower(level) == "family" || tolower(level) == "f"){
return(5)
}else if(tolower(level) == "genus" || tolower(level) == "g"){
return(6)
}else if(tolower(level) == "species" || tolower(level) == "s"){
return(7)
}else if(tolower(level) == "subspecies" || tolower(level) == "ss"){
return(8)
}else{
if(is.na(suppressWarnings(as.numeric(level))))
stop(sprintf("%s isn't recognized.", as.character(level)))

lvl <- as.numeric(level)
if(lvl > 0)
return(lvl)

stop("'level' cannot be negative.")
}
}

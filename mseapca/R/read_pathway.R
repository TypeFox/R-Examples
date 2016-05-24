read_pathway <-
function (fullpath) {

# load XML library
library(XML)

# read XML file
a <- xmlTreeParse(fullpath)
b <- xmlRoot(a)

# metabolite set list
ID <- NaN
metabolite_set <- NaN
for (i in 1:length(b)){

# // ID
s <- b[[i]]
u <- xmlToList(s)
v <- u[names(u)=="compound"]
id <- as.character(unlist(v))

ID[i] <- list(id)
metabolite_set[i] <- as.character(xmlAttrs(b[[i]]))
}

names(ID) <- metabolite_set
return(ID)
}

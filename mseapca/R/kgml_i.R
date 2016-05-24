kgml_i <-
function (filename) {

# // load XML library
library(XML)

# // data
a <- xmlTreeParse(filename) # // load xml
b <- xmlRoot(a) # // part of xml

# // extract
d <- names(b)
f <- b[d=="entry"]

H <- NaN;N <- NaN
for (i in 1:length(f)){
g <- xmlAttrs(f[[i]])
H[i] <- g[["type"]] # check
N[i] <- g[["name"]] # name
}

# // metabolite list
K <- N[H=="compound"]

# // error checking
L <- NaN
k1 <- 1 # index
if (length(K)>0){
for (i in 1:length(K)){

a1 <- strsplit(K[i]," ") # multiple IDs
a2 <- a1[[1]]

# // ID
for (j in 1:length(a2)){
p <- unlist(strsplit(a2[j],":"))
L[k1] <- p[2] # //  KEGG ID

k1 <- k1+1
}
}
}

# metabolic pathway
m <- xmlAttrs(b)
Q <- m[["title"]]

z <- list(unique(L))
names(z) <- Q

return(z)

}

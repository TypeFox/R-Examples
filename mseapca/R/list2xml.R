list2xml <-
function (filepath, M) {

# load XML library
library(XML)

N <- xmlNode("file", attrs = c(path=filepath))

for (i in 1:length(M)){
n <- names(M[i])
c <- M[i][[1]]

# add node to XML
z = xmlNode("metabolite_set", attrs = c(name=n)) # metabolite set
for (j in 1:length(c)){
z <- append.xmlNode(z, xmlNode("compound",c[j])) # metabolite ID
}
N <- append.xmlNode(N,z)
}
return(N)
}

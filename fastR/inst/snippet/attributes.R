ddd <- data.frame(number=1:5,letter=letters[1:5])
attributes(ddd)
dim(ddd)
nrow(ddd)
ncol(ddd)
names(ddd)
row.names(ddd)
row.names(ddd) <- c("Abe","Betty","Claire","Don","Ethel")
ddd                 # row.names affects how a data.frame prints

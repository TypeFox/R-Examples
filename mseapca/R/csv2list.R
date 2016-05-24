csv2list <-
function (filepath) {

# read csv file
x <- read.csv(filepath)

metabolite_set <- unique(x[,1])

# IDs for metabolite set
id <- x[,2]
ID <- NaN
for (i in 1:length(metabolite_set)){
id_set <- id[x[,1]==metabolite_set[i]]
ID[i] <- list(unique(id_set))
}

names(ID) <- metabolite_set
return(ID)
}

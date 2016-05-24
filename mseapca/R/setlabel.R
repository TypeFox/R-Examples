setlabel <-
function (M_ID, M){

L <- NaN
for(i in 1:length(M)){

# metabolite set
m <- as.character(unlist(M[i]))
m <- unique(m)

# metabolite ID
l <- NaN
for (j in 1:length(M_ID)){

# conversion@"," ¨ ";"
a <- chartr(",",";",M_ID[j])

# split of multiple ID
b <- unlist(strsplit(a,";"))
b <- unique(b)

# matching
c <- charmatch(b,m)
l[j] <- 0
if (sum(!is.na(c))>=1){l[j] <- 1}
}
L <- cbind(L,l)
}

L <- L[,-1]

colnames(L) <- names(M)
return(L)
}

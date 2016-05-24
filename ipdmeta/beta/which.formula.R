which.formula <- 

function(f1,f2){

terms1 <- terms(f1)
terms2 <- terms(f2)

#INDEX OF WHICH TERMS OF F1 ARE SHARED BY F2

which(attr(terms1,"term.labels")%in%attr(terms2,"term.labels"))

}

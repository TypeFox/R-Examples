"form.add" <-
function(formula1, formula2){
formula( paste( deparse(formula1,width.cutoff=500), sub("~","+",deparse(formula2,width.cutoff=500))))
}


plotCat2 <-
function (X, LI, ...) 
{
    tab.cat <- data.frame(X[, sapply(X, is.factor)])
    names(tab.cat) <- names(X)[sapply(X, is.factor)]
    if (ncol(tab.cat) > 0.5) {
        dev.new()
if(ncol(tab.cat)==3){fen<-c(2,2)}
else{fen<-n2mfrow(ncol(tab.cat))
}
        par(mfrow = fen)
        for (i in 1:ncol(tab.cat)) {
            s.class(LI, tab.cat[, i], sub = names(tab.cat)[i], 
                ...)
        }

    }
}

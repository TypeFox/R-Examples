"join.kasc" <- function(pts, w)
{
    ## Verifications
    if (!inherits(w, "kasc")) stop("non convenient data")

    ## output
    sorties<-1:nrow(pts)

    ## applies the function join.asc to each map of w
    for (i in 1:length(w)) {
        carp<-getkasc(w, names(w)[i])
        fac<-join.asc(pts, carp)
        sorties<-cbind.data.frame(sorties, fac)
    }

    ## output
    sorties<-sorties[,-1]
    names(sorties)<-names(w)
    return(sorties)
}


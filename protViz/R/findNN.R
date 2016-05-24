#R
# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/R/findNN.R $
# $Id: findNN.R 6178 2014-02-27 09:33:30Z cpanse $
# $Date: 2014-02-27 10:33:30 +0100 (Thu, 27 Feb 2014) $


# TODO 
# compute score by sum of error div. by number of hits

findNN<-function(q, vec, check=FALSE) {

    if (check){ if ( is.unsorted(vec)){
        return (list(error="vec is not sorted"))
    }}


    out <- .C("findNN",
        m=as.integer(length(q)),
        n=as.integer(length(vec)),
        q=as.double(q),
        vec=as.double(vec),
        NN=as.integer(rep(-1, length(q))),
        PACKAGE = "protViz")

    return (out$NN + 1)
}

allVec <-
function( ngam){

    nlabels <- maxlabel( ngam )

    all.vecs <- t( sapply( 0:nlabels, label2vec, ngam=ngam) )
    
    rownames( all.vecs ) <- 0:nlabels 
    return( all.vecs )
}

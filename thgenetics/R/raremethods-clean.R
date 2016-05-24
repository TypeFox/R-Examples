## Cleaned up version of raremethods.R in the sims folder

zstat_perm = function(g, m=rep(1,ncol(g)),
                      aff,
                      thresh=1,
                      gsubsetmatrix=NULL,
                      use_sign=TRUE, use_weight=0, strategy=1,
                      nperm=100){
  ret_pvalue = as.double(0)
  if(is.null(gsubsetmatrix))
    gsubsetmatrix = matrix(1,nrow=1, ncol=length(m))

  .C("zstat_perm",
     as.double(g), as.integer(m), as.integer(length(m)),
     as.double(aff), as.integer(length(aff)),
     as.double(thresh),
     as.integer(gsubsetmatrix), as.integer(nrow(gsubsetmatrix)),
     as.integer(use_sign), as.integer(use_weight), as.integer(strategy),
     as.integer(nperm),
     ret_pvalue, DUP=FALSE)

  return(ret_pvalue)
}

zstat_pathway_perm = function(g, m=rep(1, ncol(g)),
                              aff,
                              thresh=1,
                              gsubsetmatrix=NULL,
                              use_sign=TRUE, use_weight=0, strategy=3,
                              nperm = 100){
  ##
  if(strategy != 3)
    stop("zstat_pathway_perm: strategy must = 3")

  ret_pvalue = as.double(0)
  if(is.null(gsubsetmatrix))
    gsubsetmatrix = matrix(1, nrow=1, ncol=length(m))

  .C("zstat_pathway_perm",
    as.double(g), as.integer(m), as.integer(length(m)),
    as.double(aff), as.integer(length(aff)),
    as.double(thresh),
    as.integer(gsubsetmatrix), as.integer(nrow(gsubsetmatrix)),
    as.integer(use_sign), as.integer(use_weight), as.integer(strategy),
    as.integer(nperm),
    ret_pvalue, DUP=FALSE)

  return(ret_pvalue)
}

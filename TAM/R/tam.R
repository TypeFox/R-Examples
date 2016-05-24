tam <-
  function(resp , irtmodel ="1PL" , formulaA=NULL, ...){
    facets <- NULL
    res <- NULL
    if (irtmodel %in% c("2PL","GPCM" , "GPCM.design", "2PL.groups") ){
      res <- tam.mml.2pl(resp, irtmodel=irtmodel,... )   
    }
    if( !is.null(formulaA) ){
      res <- tam.mml.mfr(resp, formulaA=formulaA, facets=facets, ...)
    }
    if( is.null(res) ){
      res <- tam.mml(resp, ...) 
      
    }
  }
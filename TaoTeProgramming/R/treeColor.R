"treeColor" <-
  function (branches=30, trunkColor=NULL, branchColor=NULL,
            seed=NULL) 
  {
    if(length(seed)) set.seed(seed)
    
    canvas()
    if(!length(trunkColor)) {
      trunkColor <- grep("brown", colors(), value=TRUE)
    }
    if(!length(branchColor)) {
      branchColor <- grep("green", colors(), value=TRUE)
    }
    
    # trunk
    bundle(10, orig=c(.5, .08), dest=c(.5, .95), 
           width=c(.08, .1, .02, .1), color=trunkColor)
    bundle(20, orig=c(.5, .1), dest=c(.5, .75), 
           width=c(.1, .1, .04, .1), color=trunkColor)
    bundle(20, orig=c(.5, .1), dest=c(.5, .65), 
           width=c(.1, .1, .05, .2), color=trunkColor)
    
    for(i in 1:branches) {
      height <- runif(1, .5, .95)
      extent <- (1 - height) * sign(runif(1, -1, 1))
      bundle(10, orig=c(.5, height), 
             dest=c(.5 + extent, height - runif(1, 0, .2) * 
                      abs(extent)), width=c(.01, .1, .05, .05),
             color=branchColor)
    }
    
  }

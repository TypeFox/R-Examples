build.target <-
function(x,genegroups=NULL,type) {

  targetFun <- switch(type, cor = targetCor,
                            D = targetD,
                            F = targetF,
                            G = targetG,
                            Gpos = targetGpos,
                            Gstar= targetGstar)
  res <- targetFun(x,genegroups)
}


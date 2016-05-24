tstTree <-
function(letra = character(), 
                    bandera = integer(), 
                    left = integer(), 
                    rigth = integer(), 
                    center = integer()) {
    tree <- list(ch = letra, flag = bandera, L = left, R = rigth, C = center)
    class(tree) <- append(class(tree), "tstTree")
    return(tree)
  }

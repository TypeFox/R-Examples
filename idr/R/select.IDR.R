select.IDR <-
  function(x, IDR.x, IDR.level){

    is.selected <- IDR.x < IDR.level
    x.selected <- x[is.selected,]
    n.selected <- nrow(x.selected)
    
    return(list(x=x.selected, n=n.selected, IDR.level=IDR.level))
  }

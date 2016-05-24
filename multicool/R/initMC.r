initMC = function(x){
  if(length(x) > 1){
    r = NULL

    if(is.numeric(x)){
      ## Test for whole number, with tolerance for representation
      ## From post by Tony Plate<tplate_at_acm.org>
      tolerance = .Machine$double.eps^0.5
      if(isTRUE(all(abs(x - round(x))<  tolerance))){ ## integer args
        tbl = table(x)
        elements = as.numeric(names(tbl))
        set = rep(1:length(tbl), tbl)
        
        r = new(Multicool, set)
        mcObj = list(mode = "integer", set = r$set(), elements = elements, 
                     length = r$length(), mc = r)
        
        class(mcObj) = "mc"
        return(mcObj)
      }else{ ## doubles
        tbl = table(x)
        elements = as.numeric(tbl)
        set = rep(1:length(tbl), tbl)
        
        r = new(Multicool, set)
        mcObj = list(mode = "double", set = r$set(), elements = elements,
                     length = r$length(), mc = r)
        class(mcObj) = "mc"
        return(mcObj)
      }
    }else{ ## logicals and characters and who know's whatelse
      tbl = table(x)
      elements = names(tbl)
      set = rep(1:length(tbl), tbl)
      
      r = new(Multicool, set)
      mcObj = list(mode = mode(x), set = r$set(), elements = elements,
                   length = r$length(), mc = r)
      class(mcObj) = "mc"
      return(mcObj)
    }
  }else{
    warning("The permutations of a vector of length 1 are not very interesting")
  }
}

levexp <- structure(function# Vector releveling
### Expansion or reduction of a numeric vector by matching its levels with the factor-level columns in a data frame.
(
    x, ##<<\code{numeric} vector with names of the vector representing
       ##the levels to be matched.
    levels ##<<\code{data.frame} with factor-level columns, or
           ##\code{character} vector of levels.
) {
    tx <- split(x,names(x))
    if(is.character(levels))
        dsp <- split(levels,levels)
    if(is.data.frame(levels))
        dsp <- splitFrame(levels)
    nam <- lapply(seq_len(length(tx)),
                  function(i)paste("\\b",
                                   names(tx[i]),"\\b",sep = ""))
    nms <- lapply(seq_len(length(tx)),
                  function(i)grep(nam[[i]],
                                  names(dsp),value = TRUE))
    names(nms) <- names(tx)
    nnms <- lapply(nms,length)
    nms1 <- lapply(seq_len(length(tx)),
                   function(i)rep(tx[[i]],nnms[[i]]))
    nm <- lapply(seq_len(length(tx)),
                 function(i)data.frame(nms[[i]],nms1[[i]]))
    nmd <- do.call(rbind,nm)
    nmd1 <- nmd[,2]
    names(nmd1) <- nmd[,1]
    
    nmd1 <- nmd1[names(dsp)]
    
    return(nmd1)
### numeric vector with expanded/reduced levels.
} , ex=function(){
    ##Multilevel data frame of tree-ring widths:
    data(Prings05,envir = environment())
    ## Radial increments measured on 2003:
    data(Pradii03,envir = environment())    
    
    ## Getting the factor-level names at sample level
    ntl <- names(splitFrame(Prings05,'sample'))
    ## Releveling the radii
    refs <- levexp(Pradii03,ntl)
    
})

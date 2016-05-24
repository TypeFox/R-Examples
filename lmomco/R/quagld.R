"quagld" <-
function(f,para, paracheck=TRUE) {
    if(! check.fs(f)) return()
    if(paracheck == TRUE) {
      if(! are.pargld.valid(para)) return()
    }
    La1 <- para$para[1]
    La2 <- para$para[2]
    La3 <- para$para[3]
    La4 <- para$para[4]

    x <- La1 + La2*(f**La3 - (1-f)**La4)
    names(x) <- NULL
    return(x)
}

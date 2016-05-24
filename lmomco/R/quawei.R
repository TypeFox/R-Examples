"quawei" <-
function(f,para,paracheck=TRUE) { 
    if(! check.fs(f)) return()
    if(paracheck == TRUE) {
      if(! are.parwei.valid(para)) return()
    }
    ZETA <- para$para[1]
    B    <- para$para[2]
    D    <- para$para[3]

    K  <- 1/D
    A  <- B/D
    XI <- ZETA - B 
    gev.para <- list(type = 'gev', para = c(XI,A,K))
    return(-quagev(1-f,gev.para,paracheck=paracheck))
}

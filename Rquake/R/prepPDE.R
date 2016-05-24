prepPDE <-
function(fn)
  {
    ###  extract data from the pde output
    ###  check - then return a list
    pdes = RSEIS::getPDEscreen(fn)
    pdelist = RSEIS::PDE2list(pdes)
    w1 =  which(is.na(pdelist$sec))
    zde = data.frame(pdelist)
    zde = zde[-w1, ]
    pdejsec = RSEIS::JtimL(zde)
    ipde = as.list(zde)
    ipde$jsec  = pdejsec
    any(is.na(pdejsec))
    return(ipde)

  }

"par2lmom" <-
function(para) {
    type <- para$type
       if(type == 'aep4')  {
      return(lmomaep4(para))
    }
    else if(type == 'cau') {
      return(lmomcau(para))
    }
    else if(type == 'emu') {
      return(lmomemu(para))
    }
    else if(type == 'exp') {
      return(lmomexp(para))
    }
    else if(type == 'gam') {
      return(lmomgam(para))
    }
    else if(type == 'gep') {
      return(lmomgep(para))
    }
    else if(type == 'gev') {
      return(lmomgev(para))
    }
    else if(type == 'gld') {
      return(lmomgld(para))
    }
    else if(type == 'glo') {
      return(lmomglo(para))
    }
    else if(type == 'gno') {
      return(lmomgno(para))
    }
    else if(type == 'gov') {
      return(lmomgov(para))
    }
    else if(type == 'gpa') {
      return(lmomgpa(para))
    }
    else if(type == 'gum') {
      return(lmomgum(para))
    }
    else if(type == 'kap') {
      return(lmomkap(para))
    }
    else if(type == 'kmu') {
      return(lmomkmu(para))
    }
    else if(type == 'kur') {
      return(lmomkur(para))
    }
    else if(type == 'lap') {
      return(lmomlap(para))
    }
    else if(type == 'lmrq') {
      return(lmomlmrq(para))
    }
    else if(type == 'ln3') {
      return(lmomln3(para))
    }
    else if(type == 'nor') {
      return(lmomnor(para))
    }
    else if(type == 'pe3') {
      return(lmompe3(para))
    }
    else if(type == 'ray') {
      return(lmomray(para))
    }
    else if(type == 'revgum') {
      return(lmomrevgum(para))
    }
    else if(type == 'rice') {
      return(lmomrice(para))
    }
    else if(type == 'sla') {
      return(lmomsla(para))
    }
    else if(type == 'st3') {
      return(lmomst3(para))
    }
    else if(type == 'texp') {
      return(lmomtexp(para))
    }
    else if(type == 'wak') {
      return(lmomwak(para))
    }
    else if(type == 'wei') {
      return(lmomwei(para))
    }
    else {
      stop("Did not find a valid distribution type")
    }
}


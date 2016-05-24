"par2qua" <-
function(f,para,...) {
    type <- para$type
    if(type == 'aep4') {
      return(quaaep4(f,para,...))
    }
    else if(type == 'cau') {
      return(quacau(f,para,...))
    }
    else if(type == 'emu') {
      return(quaemu(f,para,...))
    }
    else if(type == 'exp') {
      return(quaexp(f,para,...))
    }
    else if(type == 'gam') {
      return(quagam(f,para,...))
    }
    else if(type == 'gep') {
      return(quagep(f,para,...))
    }
    else if(type == 'gev') {
      return(quagev(f,para,...))
    }
    else if(type == 'gld') {
      return(quagld(f,para,...))
    }
    else if(type == 'glo') {
      return(quaglo(f,para,...))
    }
    else if(type == 'gno') {
      return(quagno(f,para,...))
    }
    else if(type == 'gov') {
      return(quagov(f,para,...))
    }
    else if(type == 'gpa') {
      return(quagpa(f,para,...))
    }
    else if(type == 'gum') {
      return(quagum(f,para,...))
    }
    else if(type == 'kap') {
      return(quakap(f,para,...))
    }
    else if(type == 'kmu') {
      return(quakmu(f,para,...))
    }
    else if(type == 'kur') {
      return(quakur(f,para,...))
    }
    else if(type == 'lap') {
      return(qualap(f,para,...))
    }
    else if(type == 'lmrq') {
      return(qualmrq(f,para,...))
    }
    else if(type == 'ln3') {
      return(qualn3(f,para,...))
    }
    else if(type == 'nor') {
      return(quanor(f,para,...))
    }
    else if(type == 'pe3') {
      return(quape3(f,para,...))
    }
    else if(type == 'ray') {
      return(quaray(f,para,...))
    }
    else if(type == 'revgum') {
      return(quarevgum(f,para,...))
    }
    else if(type == 'rice') {
      return(quarice(f,para,...))
    }
    else if(type == 'sla') {
      return(quasla(f,para,...))
    }
    else if(type == 'st3') {
      return(quast3(f,para,...))
    }
    else if(type == 'texp') {
      return(quatexp(f,para,...))
    }
    else if(type == 'tri') {
      return(quatri(f,para,...))
    }
    else if(type == 'wak') {
      return(quawak(f,para,...))
    }
    else if(type == 'wei') {
      return(quawei(f,para,...))
    }
    else {
      stop("Did not find a valid distribution type")
    }
}


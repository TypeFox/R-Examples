"are.par.valid" <-
function(para, paracheck=TRUE, ...) {
    if(paracheck) return(TRUE)
    if(is.null(para$para)) {
      warning("The parameter object is missing a para attribute")
      return()
    }
    if(is.null(para$type)) {
      warning("The parameter object is missing a type attribute")
      return()
    }
    type <- para$type
    if(type == 'aep4') {
      return(are.paraep4.valid(para,...))
    } else if(type == 'cau') {
      return(are.parcau.valid(para,...))
    }
    else if(type == 'emu') {
      return(are.paremu.valid(para,...))
    }
    else if(type == 'exp') {
      return(are.parexp.valid(para,...))
    }
    else if(type == 'texp') {
      return(are.partexp.valid(para,...))
    }
    else if(type == 'gam') {
      return(are.pargam.valid(para,...))
    }
    else if(type == 'gep') {
      return(are.pargep.valid(para,...))
    }
    else if(type == 'gev') {
      return(are.pargev.valid(para,...))
    }
    else if(type == 'gld') {
      return(are.pargld.valid(para,...))
    }
    else if(type == 'glo') {
      return(are.parglo.valid(para,...))
    }
    else if(type == 'gno') {
      return(are.pargno.valid(para,...))
    }
    else if(type == 'gov') {
      return(are.pargov.valid(para,...))
    }
    else if(type == 'gpa') {
      return(are.pargpa.valid(para,...))
    }
    else if(type == 'gum') {
      return(are.pargum.valid(para,...))
    }
    else if(type == 'kur') {
      return(are.parkur.valid(para,...))
    }
    else if(type == 'kmu') {
      return(are.parkmu.valid(para,...))
    }
    else if(type == 'lap') {
      return(are.parlap.valid(para,...))
    }
    else if(type == 'lmrq') {
      return(are.parlmrq.valid(para,...))
    }
    else if(type == 'ln3') {
      return(are.parln3.valid(para,...))
    }
    else if(type == 'kap') {
      return(are.parkap.valid(para,...))
    }
    else if(type == 'nor') {
      return(are.parnor.valid(para,...))
    }
    else if(type == 'pe3') {
      return(are.parpe3.valid(para,...))
    }
    else if(type == 'ray') {
      return(are.parray.valid(para,...))
    }
    else if(type == 'rice') {
      return(are.parrice.valid(para,...))
    }
    else if(type == 'revgum') {
      return(are.parrevgum.valid(para,...))
    }
    else if(type == 'sla') {
      return(are.parsla.valid(para,...))
    }
    else if(type == 'st3') {
      return(are.parst3.valid(para,...))
    }
    else if(type == 'tri') {
      return(are.partri.valid(para,...))
    }
    else if(type == 'wak') {
      return(are.parwak.valid(para,...))
    }
    else if(type == 'wei') {
      return(are.parwei.valid(para,...))
    }
    else {
      stop("Did not find a valid distribution type")
    }
}

"par2pdf" <-
function(f,para,...) {
    type <- para$type
    if(type == 'aep4') {
      return(pdfaep4(f,para))
    }
    else if(type == 'cau') {
      return(pdfcau(f,para))
    }
    else if(type == 'emu') {
      return(pdfemu(f,para))
    }
    else if(type == 'exp') {
      return(pdfexp(f,para))
    }
    else if(type == 'gam') {
      return(pdfgam(f,para))
    }
    else if(type == 'gep') {
      return(pdfgep(f,para))
    }
    else if(type == 'gev') {
      return(pdfgev(f,para))
    }
    else if(type == 'gld') {
      return(pdfgld(f,para,...))
    }
    else if(type == 'glo') {
      return(pdfglo(f,para))
    }
    else if(type == 'gno') {
      return(pdfgno(f,para))
    }
    else if(type == 'gov') {
      return(pdfgov(f,para))
    }
    else if(type == 'gpa') {
      return(pdfgpa(f,para))
    }
    else if(type == 'gum') {
      return(pdfgum(f,para))
    }
    else if(type == 'kap') {
      return(pdfkap(f,para))
    }
    else if(type == 'kmu') {
      return(pdfkmu(f,para))
    }
    else if(type == 'kur') {
      return(pdfkur(f,para))
    }
    else if(type == 'lap') {
      return(pdflap(f,para))
    }
    else if(type == 'lmrq') {
      return(pdflmrq(f,para))
    }
    else if(type == 'ln3') {
      return(pdfln3(f,para))
    }
    else if(type == 'nor') {
      return(pdfnor(f,para))
    }
    else if(type == 'pe3') {
      return(pdfpe3(f,para))
    }
    else if(type == 'ray') {
      return(pdfray(f,para))
    }
    else if(type == 'revgum') {
      return(pdfrevgum(f,para))
    }
    else if(type == 'rice') {
      return(pdfrice(f,para))
    }
    else if(type == 'sla') {
      return(pdfsla(f,para))
    }
    else if(type == 'st3') {
      return(pdfst3(f,para))
    }
    else if(type == 'texp') {
      return(pdftexp(f,para))
    }
    else if(type == 'tri') {
      return(pdftri(f,para))
    }
    else if(type == 'wak') {
      return(pdfwak(f,para))
    }
    else if(type == 'wei') {
      return(pdfwei(f,para))
    }
    else {
      stop("Did not find a valid distribution type")
    }
}


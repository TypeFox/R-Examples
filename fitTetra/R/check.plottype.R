check.plottype <-
function(plot.type) {
  if (plot.type=="none" || !require("grDevices", quietly=T)) "none"
  else {
    cap <- capabilities()
    if (!cap["png"] || !exists("png")) "none"
    else {
      #we can at least generate png images 
      plot.type <- tolower(plot.type)
      if (plot.type=="png") "png"
      else if (plot.type=="emf") {
        if (.Platform$OS.type=="windows" && require("devEMF", quietly=T) && exists("emf")) "emf"
        else "png"
      }
      else if (plot.type=="svg") {
        if (cap["cairo"]) "svg"
        else "png"
      }
      else if (plot.type=="pdf") {
        if (exists("pdf")) "pdf"
        else "png"                                    
      }
      else "png" #no valid type but not "none", use png
    }         
  }
}

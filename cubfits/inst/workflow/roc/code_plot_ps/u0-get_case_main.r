### Convert i.case to case.main for plotting title.

get.case.main <- function(i.case, model){
  ### cases with phi.
  if(i.case == paste(model, "_wphi_pm", sep = "")){
    i.case.main <- paste(model, ": with phi.Obs, start: posterior mean", sep = "")
  } else if(i.case == paste(model, "_wphi_scuo", sep = "")){
    i.case.main <- paste(model, ": with phi.Obs, start: SCUO", sep = "")
  } else if(i.case == paste(model, "_wphi_true", sep = "")){
    i.case.main <- paste(model, ": with phi.Obs, start: True", sep = "")
  } else if(i.case == paste(model, "_wphi_bInit", sep = "")){
    i.case.main <- paste(model, ": with phi.Obs, start: b.Init", sep = "")

  ### cases without phi.
  } else if(i.case == paste(model, "_wophi_pm", sep = "")){
    i.case.main <- paste(model, ": without phi.Obs, start: posterior mean", sep = "")
  } else if(i.case == paste(model, "_wophi_scuo", sep = "")){
    i.case.main <- paste(model, ": without phi.Obs, start: SCUO", sep = "")
  } else if(i.case == paste(model, "_wophi_true", sep = "")){
    i.case.main <- paste(model, ": without phi.Obs, start: True", sep = "")
  } else if(i.case == paste(model, "_wophi_bInit", sep = "")){
    i.case.main <- paste(model, ": without phi.Obs, start: b.Init", sep = "")

  ### cases without phi following with phi.
  } else if(i.case == paste(model, "_wphi_wophi_pm", sep = "")){
    i.case.main <- paste(model, ": without phi.Obs, start: wphi,posterior mean", sep = "")
  } else if(i.case == paste(model, "_wphi_wophi_scuo", sep = "")){
    i.case.main <- paste(model, ": without phi.Obs, start: wphi,SCUO", sep = "")

  } else{
    i.case.main <- "case.name Not Found"
  }

  ### Return.
  i.case.main <- paste(i.case.main, ", PostScale", sep = "")
  i.case.main
}

# shiny wrappers to give the shiny app a function which can be called easily



shiny_binom_single_onestage=function(){
  shiny::runApp(file.path(system.file("Shiny/binom_single_stage",package = "EurosarcBayes")))
}




shiny_binom_single_twostage=function(){
  shiny::runApp(file.path(system.file("Shiny/binom_two_stage",package = "EurosarcBayes")))
}

shiny_LINES_posterior=function(){
  shiny::runApp(file.path(system.file("Shiny/LINES_posterior_viewer",package = "EurosarcBayes")))
}

# ended

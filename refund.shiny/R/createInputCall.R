#' Create input calls for plot_shiny.fosr()
#' 
#' Internal method that constructs the input calls for plot_shiny.fosr(). The
#' variable name and values are passed as arguments, and a corresponding slider (for
#' numeric) or drop-down (for factor) input is constructed.
#' 
#' @param name variable name
#' @param variable variable values from dataset
#' 
#' @author Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu}
#' 
createInputCall = function(name, variable){
  
  if(is.numeric(variable)){
    step.width = ((max(variable) - min(variable))/50)
    call("sliderInput", inputId = name, label = name, 
         min = signif(min(variable), 2), max = signif(max(variable), 2), value = mean(variable), 
         animate = animationOptions(interval=300, loop=T))
#  } else if(is.factor(variable) & length(levels(variable)) <= 2){
#    checkboxInput(inputId = name, label = levels(variable)[2], value = FALSE)
  } else if(is.factor(variable) & length(levels(variable)) ){
    selectInput(name, label = name, choices = levels(variable), selected = 1)
  }
  
}


#' Create input calls for plot_shiny.mfpca()
#'
#' Internal method that constructs the input calls for plot_shiny.mfpca(). The
#' number of sliders to construct for each level is passed as an argument, and 
#' corresponding sliders for each FPC are constructed.
#' 
#' @param plot.npc list of 2 numeric entries giving number of sliders at each level
#' @param plotObj the mfpca object plotted in the plot_shiny.mfpca() function. 
#' @param percents the percent variance calculated for each eigen values for levels 1 and 2.
#' 
#' @return a list of numbers that indicate percent variance for selected level.
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' 
mfpcaCalls = function(plot.npc, plotObj, percents){

  calls <- PCs <- list(NA, NA)
  
  names(calls) = names(PCs) = c("level1", "level2")
  
  for(j in 1:2){
    
    #numSliders = plotObj$npc[[j]]
    
    calls[[j]] = as.list(rep(NA, plot.npc[[j]]))
    PCs[[j]] = rep(NA, plot.npc[[j]])
    
     for(i in 1:plot.npc[[j]]){
      PCnum = paste("PC ", j, ".",  i, sep="")    
      calls[[j]][[i]] =  eval(call("sliderInput", inputId= PCnum, label = paste0(PCnum, ": ", percents[[j]][[i]], "% of Level ", j,  " Variance"),
                                   min = -2, max = 2, step = .1, value = 0, post = " SD", animate = animationOptions(interval=400, loop=T)))
      PCs[[j]][i] = PCnum
    }  
  }
  ret <- list(calls, PCs); names(ret) <- c("calls", "PCs")  
  return(ret)
}

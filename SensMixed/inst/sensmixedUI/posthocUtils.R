calc.cols2 <- function(x)
{
  if(x<0.001) 
    return("p-value < 0.001")#return("red")# 
  if(x<0.01) 
    return("p-value < 0.01")#return("orange")# 
  if(x<0.05) 
    return("p-value < 0.05")#return("yellow")# 
  return("NS")#return("grey")#
}

posthocResult <- function(){
  if (is.null(Data())) {return()}
  if(input$analysis == "Consumer data"){
    if(class(Data()) == "sensmixed") {return("")} ## avoid error from xtable
    result <- Data()
    names.lsm <- "Population means for attribute"
    names.dlsm <- "Multiple comparison tests"
  }
  else{
    if(class(Data()) == "consmixed") {return("")} ## avoid error from xtable
    if(is.null(input$AttrPosthoc) || length(input$AttrPosthoc)>1)
    {return()}
    if(!("post_hoc" %in% names(Data()))) {return()}
    
    result <- Data()$step_res[[input$AttrPosthoc]]     
    
    names.lsm <- paste("Population means for attribute ", 
                       input$AttrPosthoc)
    names.dlsm <- paste("Multiple comparison tests for attribute ",
                        input$AttrPosthoc)    
  } 
  
  if(input$whichPlot == "LSMEANS"){
    ph <- result$lsmeans.table
    
    rnames <- rownames(ph)
    diffs.facs <- sapply(rnames, 
                         function(x) 
                           substring(x, 1, 
                                     substring.location(x, " ")$first[1]-1), 
                         USE.NAMES = FALSE)    
    find.fac <- diffs.facs %in% input$effsPlot
    ph <- ph[find.fac,]
    ph[, which(colnames(ph)=="p-value")] <- 
      format.pval(ph[, which(colnames(ph)=="p-value")], digits=3, eps=1e-3)
    ph_tab <- xtable(ph, align = paste(c("l", rep("c", ncol(ph))), 
                                       collapse = ""), 
                     display = c(rep("s",
                                     which(colnames(ph) == "Estimate")), 
                                 rep("f", 6), "s"))
    
    caption(ph_tab) <- names.lsm
    print(ph_tab, caption.placement="top",
          table.placement="H", 
          type = "html",
          html.table.attributes = getOption("xtable.html.table.attributes",
                                            "rules='groups' width='105%'"))
    
  }
  else{
    ph <- result$diffs.lsmeans.table
    rnames <- rownames(ph)
    diffs.facs <- sapply(rnames, 
                         function(x) 
                           substring(x, 1, 
                                     substring.location(x, " ")$first[1]-1), 
                         USE.NAMES = FALSE)    
    find.fac <- diffs.facs %in% input$effsPlot
    ph <- ph[find.fac,]
    
    ph[, 7] <- format.pval(ph[, 7], digits=3, eps=1e-3)
    
    ph_tab <- xtable(ph, align="lccccccc", 
                     display=c("s","f","f","f","f","f","f", "s"))
    caption(ph_tab) <- names.dlsm
    
    print(ph_tab, caption.placement="top",
          table.placement="H", 
          type = "html",
          html.table.attributes = 
            getOption("xtable.html.table.attributes",
                      "rules='groups' width='105%'"))
    
  }
}

posthocPlot <- function(){
  if (is.null(Data())) {return()}
  if(input$analysis == "Consumer data"){
    if(class(Data()) == "sensmixed") {return("")} ## avoid error from xtable
    plot(Data(), cex = 1.6, 
         which.plot = input$whichPlot, effs = input$effsPlot) 
  }else {
    if(!("post_hoc" %in% names(Data()))) {return()}
    
    if(is.null(input$AttrPosthoc) || length(input$AttrPosthoc)>1)
    {return()}
    
    if(input$MAM == "TRUE"){
      if(input$whichPlot == "LSMEANS")
        tab <- Data()$step_res[[input$AttrPosthoc]]$lsmeans.table
      else
        tab <- Data()$step_res[[input$AttrPosthoc]]$diffs.lsmeans.table
      plotLSMEANS(table = tab, 
                  response = Data()$step_res[[input$AttrPosthoc]]$response, 
                  which.plot = input$whichPlot, effs = input$effsPlot)
    }
    else
      plot(Data()$step_res[[input$AttrPosthoc]], cex = 1.6, 
           which.plot = input$whichPlot, effs = input$effsPlot) 
  }       
}
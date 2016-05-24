## TODO: partition the code

getNameStep <- function(){
  if(input$analysis == "Consumer data")
    return("stepoutput")
  else
    return(input$AttrStep)
}


## return table of results for the step function for the random part
## in latex or html formats using xtable
## for sensory/consumer data
stepRandResult <- function(){
  if (is.null(Data())) {return("")} 
  if(input$analysis == "Consumer data"){
    if(class(Data()) == "sensmixed") {return("")} ## avoid error from xtable
    rnd <- Data()$rand.table
    rnd[ , "p.value"] <- 
      format.pval(rnd[, "p.value"], digits=3, eps=1e-3)
    rnd_tab <- xtable(rnd, align = 
                        paste(c("l", rep("l", ncol(rnd))), collapse = ""), 
                      display= c("s","f","d","s","s"))
    caption(rnd_tab) <- paste("Likelihood ratio tests for the 
                              random-effects and their order of elimination representing Step 1 of 
                              the automated analysis") 
    
    
    
    print(rnd_tab, caption.placement="top", table.placement="H", 
          type = input$typetable,
          html.table.attributes = getOption("xtable.html.table.attributes",
                                            "rules='groups' width='100%'"))
    
  }
  else{
    if(class(Data()) == "consmixed") {return("")} ## avoid error from xtable
    if(is.null(input$AttrStep) || length(input$AttrStep)>1)
    {return()}
    st <- Data()$step_res[[input$AttrStep]] 
    
    st$rand.table[ , "p.value"] <- format.pval(st$rand.table[, "p.value"], 
                                               digits=3, eps=1e-3)
    if(ncol(st$rand.table) == 3){
      colnames(st$rand.table) <- c("Chi.sq","Chi.DF" , "p-value")
      rand.table_tv <- xtable(st$rand.table, align="llll", 
                              display=c("s","f","d","s"))
    }
    else{
      colnames(st$rand.table) <- c("Chi.sq","Chi.DF" , "elim.num", "p-value")
      rand.table_tv <- xtable(st$rand.table, align="lllll", 
                              display=c("s","f","d","d","s"))
    }     
    caption(rand.table_tv) <- paste("Likelihood ratio tests for the 
random-effects and their order of elimination representing Step 1 of 
the automated analysis for the attribute", input$AttrStep)       
    
    print(rand.table_tv, caption.placement="top", table.placement="H", 
          type = input$typetable, 
          html.table.attributes = getOption("xtable.html.table.attributes",
                                            "rules='groups' width='100%'")) 
    
  }    
}


## return table of results for the step function for the fixed part
## in latex or html formats using xtable
## for sensory/consumer data
stepFixedResult <- function(){
  if (is.null(Data())) {return()}
  if(input$analysis == "Consumer data"){
    if(class(Data()) == "sensmixed") {return("")} ## avoid error from xtable
    an <- Data()$anova.table
    an[, "Pr(>F)"] <- format.pval(an[, "Pr(>F)"], digits=3, eps=1e-3)
    if("elim.num" %in% colnames(an))
      an_tab <- xtable(an, align = paste(c("l", rep("l", ncol(an))), 
                                         collapse = ""), 
                       display = c("s","f","f","d","f","f","s", "s"))
    else
      an_tab <- xtable(an, align = paste(c("l", rep("l", ncol(an))), 
                                         collapse = ""), 
                       display = c("s","f","f","d","f","f", "s"))
    caption(an_tab) <- 
      paste("F-tests for the fixed-effects and their order of elimination representing Step 3 of the automated analysis")
    
    print(an_tab, caption.placement="top",
          table.placement="H", 
          type = input$typetable, 
          html.table.attributes = 
            getOption("xtable.html.table.attributes",
                      "rules='groups' width='100%'"))    
    
  }
  else{
    if(class(Data()) == "consmixed") {return("")} ## avoid error from xtable
    if(is.null(input$AttrStep) || length(input$AttrStep)>1)
    {return()}
    
    
    
    st <- Data()$step_res[[input$AttrStep]] 
    
    
    st$anova.table[, "Pr(>F)"] <- format.pval(st$anova.table[, "Pr(>F)"], 
                                              digits=3, eps=1e-3)
    if("dprimeav" %in% colnames(st$anova.table)){
      colnames(st$anova.table) <-
        c("Sum Sq", "Mean Sq", "NumDF", "DenDF", "F-value","d-prime", "Pr(>F)")
      anova.table_tv <- xtable(st$anova.table, align="llllllll", 
                               display=c("s","f","f","s","f","f","f", "s"))
    }else{
      colnames(st$anova.table) <-
        c("Sum Sq", "Mean Sq", "NumDF", "DenDF", "F-value", "Pr(>F)")
      anova.table_tv <- xtable(st$anova.table, align="lllllll", 
                               display=c("s","f","f","s","f","f", "s"))
    }      
    
    caption(anova.table_tv) <- 
      paste("F-tests for the fixed-effects  for the attribute",
            input$AttrStep)
    
    print(anova.table_tv, caption.placement="top",
          table.placement="H", 
          type = input$typetable,
          html.table.attributes = 
            getOption("xtable.html.table.attributes",
                      "rules='groups' width='100%'"))     
  }  
  
}
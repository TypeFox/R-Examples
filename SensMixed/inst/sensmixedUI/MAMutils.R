## return table of results for the MAManalysis function
## in  html formats using xtable
## for sensory data

## MAM ANOVA table (similar to sensmixed function)
mamanova <- function(){
  if(is.null(uploadData())) { return() }
  if(input$analysis == "Consumer data") { return() }
  if(is.null(Data()) || is.null(Data()$MAMan)){return()}
  if(class(Data()) == "consmixed") { return() }
  resmam <- Data()$MAMan[[3]][, , input$AttrMAManalysis]
  resmam[ , "Pval"] <- 
    format.pval(resmam[, "Pval"], digits=3, eps=1e-3)
  resmam <- xtable(resmam, align = 
                     paste(c("l", rep("l", ncol(resmam))), collapse = ""), 
                   display= c("s","f","f","f","f","s"))
  caption(resmam) <- paste(" MAM ANOVA table")
  
  
  
  print(resmam, caption.placement="top", table.placement="H", 
        html.table.attributes = getOption("xtable.html.table.attributes",
                                          "rules='groups' width='100%'"),
        type = "html")
}


## individual ANOVA table
mamindiv <- function(){
  if(is.null(uploadData())) { return() }
  if(input$analysis == "Consumer data") { return() }
  if(is.null(Data()) || is.null(Data()$MAMan)){return()}
  if(class(Data()) == "consmixed") { return() }
  resindiv <- Data()$MAMan[[1]][, , input$AttrMAManalysis]
  resindiv <- xtable(resindiv, align = 
                       paste(c("l", rep("l", ncol(resindiv))), collapse = ""), 
                     display= c("s",rep("f", ncol(resindiv))))
  caption(resindiv) <- paste(" Individual ANOVA table") 
  print(resindiv, caption.placement="top", table.placement="H", 
        html.table.attributes = getOption("xtable.html.table.attributes",
                                          "rules='groups' width='100%'"),
        type = "html")
}

## Performance indices (see MAM-CAP table by Caroline Peltier)
mamperf <- function(){
  if(is.null(uploadData())) { return() }
  if(input$analysis == "Consumer data") { return() }
  if(is.null(Data()) || is.null(Data()$MAMan)){return()}
  if(class(Data()) == "consmixed") { return() }
  resperf <- Data()$MAMan[[2]][, , input$AttrMAManalysis]
  resperf <- xtable(resperf, align = 
                      paste(c("l", rep("l", ncol(resperf))), collapse = ""), 
                    display= c("s",rep("f", ncol(resperf))))
  caption(resperf) <- paste(" Individual performance tests")
  print(resperf, caption.placement="top", table.placement="H", 
        html.table.attributes = getOption("xtable.html.table.attributes",
                                          "rules='groups' width='100%'"),
        type = "html")
}

## all pairwise comparisons tests + CI
mamposthoc <-function(){
  if(is.null(uploadData())) { return() }
  if(input$analysis == "Consumer data") { return() }
  if(is.null(Data()) || is.null(Data()$MAMan)){return()}
  if(class(Data()) == "consmixed") { return() }
  
  resposthoc <- Data()$MAMan[[5]][, , input$AttrMAManalysis]
  resci <- Data()$MAMan[[8]][, , input$AttrMAManalysis]
  resposthoc <- cbind(resposthoc, resci)
  resposthoc[ , "Pval"] <- 
    format.pval(resposthoc[, "Pval"], digits=3, eps=1e-3)
  resposthoc <- xtable(resposthoc, align = 
                         paste(c("l", rep("l", ncol(resposthoc))), collapse = ""), 
                       display= c("s",rep("f", ncol(resposthoc))))
  caption(resposthoc) <- paste(" Pairwise product differences")
  print(resposthoc, caption.placement="top", table.placement="H", 
        html.table.attributes = getOption("xtable.html.table.attributes",
                                          "rules='groups' width='100%'"),
        type = "html")
}

mamdiffmean <- function(){
  if(is.null(uploadData())) { return() }
  if(input$analysis == "Consumer data") { return() }
  if(is.null(Data()) || is.null(Data()$MAMan)){return()}
  if(class(Data()) == "consmixed") { return() }
  resdiffmean <- Data()$MAMan[[6]]
  ## change rownames in order not to have duplicated names
  oddrows <- seq(2, nrow(resdiffmean), by = 2)
  rnames <- rownames(resdiffmean )
  rnames[oddrows] <- unlist(lapply(oddrows, 
                                   function(x) paste(rep(" ", x), 
                                                     collapse = "")))
  rownames(resdiffmean) <- rnames
  resdiffmean <- xtable(resdiffmean, align = 
                          paste(c("l", rep("l", ncol(resdiffmean))), 
                                collapse = ""))
  caption(resdiffmean) <- paste("Post-hoc comparison of each product
with the mean of the remaining products")
  print(resdiffmean, caption.placement="top", table.placement="H", 
        html.table.attributes = getOption("xtable.html.table.attributes",
                                          "rules='groups' width='100%'"),
        type = "html")
}
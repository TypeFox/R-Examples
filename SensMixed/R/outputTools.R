#### the file contains the functions used for presenting the results in tables #
#### or plots ##################################################################


calc.cols <- function(x)
{
  if(x<0.001) 
    return("p-value < 0.001")#return("red")# 
  if(x<0.01) 
    return("p-value < 0.01")#return("orange")# 
  if(x<0.05) 
    return("p-value < 0.05")#return("yellow")# 
  return("NS")#return("grey")#
}



change.inter.symbol <- function(x, interact.symbol){
  if(grepl(":", x)){
    symb.loc <- substring.location(x, ":")
    spl.effs <- strsplit(x,":")[[1]]
    x <- paste(spl.effs, collapse=interact.symbol)
    return(x)
  }
  x
}

.changeOutput <- function(vals, pvals, isRand){
  colnames.out <- rownames(vals)
  names <- colnames(vals)
  tr <- vector("list", length(colnames.out))
  
  for(i in 1:length(colnames.out)){       
    tr[[i]] <- createTexreg(
      coef.names = names, se=vals[i,],
      coef = vals[i,],
      pvalues = pvals[i,], isRand=isRand)
  }
  
  names(tr) <- colnames.out
  return(tr)
}

# facetAdjust <- function(x, pos = c("up", "down"))
# {
#   pos <- match.arg(pos)
#   p <- ggplot_build(x)
#   gtable <- ggplot_gtable(p); #dev.off()
#   dims <- apply(p$panel$layout[2:3], 2, max)
#   nrow <- dims[1]
#   ncol <- dims[2]
#   panels <- sum(grepl("panel", names(gtable$grobs)))
#   space <- ncol * nrow
#   n <- space - panels
#   if(panels != space){
#     idx <- (space - ncol - n + 1):(space - ncol)
#     gtable$grobs[paste0("axis_b",idx)] <- list(gtable$grobs[[paste0("axis_b",panels)]])
#     if(pos == "down"){
#       rows <- grep(paste0("axis_b\\-[", idx[1], "-", idx[n], "]"), 
#                    gtable$layout$name)
#       lastAxis <- grep(paste0("axis_b\\-", panels), gtable$layout$name)
#       gtable$layout[rows, c("t","b")] <- gtable$layout[lastAxis, c("t")]
#     }
#   }
#   class(gtable) <- c("facetAdjust", "gtable", "ggplot"); gtable
# }
# 
# print.facetAdjust <- function(x, newpage = is.null(vp), vp = NULL) {
#   if(newpage)
#     grid.newpage()
#   if(is.null(vp)){
#     grid.draw(x)
#   } else {
#     if (is.character(vp)) 
#       seekViewport(vp)
#     else pushViewport(vp)
#      grid.draw(x)
#     upViewport()
#  }
#  invisible(x)
# }



.plotSensMixed <- function(val, pval, title, mult = FALSE, sep = FALSE,
                              cex = 2,                           
                              interact.symbol = ":", ylab = ""){
  ## change the interaction symbol
  if(!interact.symbol == ":")      
    rownames(pval) <- rownames(val) <-  sapply(rownames(val), change.inter.symbol, 
                                               interact.symbol) 
  
  names.effs <- LETTERS[1:nrow(val)]
  names.effs.legend <- paste(names.effs, collapse="")

  dval <- as.data.frame(val)
  dval$effs <- rownames(dval)
  dval$effs_short <- names.effs
  dval$abbreffs <- abbreviate(dval$effs)
  suppressWarnings(dval <- melt(dval))
  dpval <- as.data.frame(pval)
  dpval$effs <- rownames(dpval)
  suppressWarnings(dpval <- melt(dpval, variable_name = "pvalue"))
  dval$pvalue <- unlist(lapply(dpval$value, calc.cols))
  uniqueInitials <- unique(dval$effs_short)
  initialShapes <- unlist(lapply(uniqueInitials, utf8ToInt))
  dval$initialShapes <- unlist(lapply(dval$effs_short, utf8ToInt))
  
  
  
  
  changelabs <- function(variable, value){
    return(effsnames[value])
  }
    
  variable <- value <- pvalue <- effsnames <- effs <- NULL
  
  p1 <- ggplot(dval, aes(x=variable, y = value, fill = pvalue, group = effs)) + 
    geom_bar(position = "dodge", stat = "identity")  + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4), 
          axis.title.x=element_blank(), 
          axis.title.y = element_text(size = rel(1.4)), 
          axis.text = element_text(size = rel(1)), 
          #legend.text = element_text(size = rel(1)), 
          #legend.title = element_text(size = rel(1)))  +
          legend.position = "none") +
    scale_fill_manual(values  = c(  "NS" = "grey", "p-value < 0.01" = "orange", 
                                    "p-value < 0.05" = "yellow", 
                                    "p-value < 0.001" = "red"), 
                      name="Significance") + ylab(ylab)
  
  if(!mult)
    return(p1 + geom_point(aes(shape = effs), fill = NA, 
                           position = position_dodge(width = 0.9), 
                           size = rel(5)) + 
             scale_shape_manual(values = initialShapes, name = "Effects"))
  else
    return(p1 + facet_wrap( ~ effs, as.table = FALSE))    

}


.changeConsmixedOutputForDoc <- function(table, name.pval){  
  table[, name.pval] <- gsub("<", "&lt ", table[, name.pval])
  table  
}

## output for the sensmixed
.createDocOutputSensmixed <- function(x, file = NA, bold = FALSE, append = TRUE, 
                                      type = "html", typeEffs = 1){
  
  if(typeEffs == 1 || typeEffs ==4){
    colnames.out.rand <- rownames(x$rand$Chi)
    names <- colnames(x$rand$Chi)
    tr_rand <- vector("list", length(colnames.out.rand))
    
    for(i in 1:length(colnames.out.rand)){       
      tr_rand[[i]] <- createTexreg(
        coef.names = names, se=x$rand$Chi[i,],
        coef = x$rand$Chi[i,],
        pvalues = x$rand$pvalueChi[i,], isRand=TRUE    
      )     
    } 
    caption.rand <- "Likelihood ratio test for the random effects"
  }
  else{
    caption.rand <- ""
    colnames.out.rand <- ""
    tr_rand = NULL
  }
    
  
  
  
  ## output for the fixed effects
  if(typeEffs == 2 || typeEffs == 4){
    colnames.out.fixed <- rownames(x$fixed$Fval)
    names <- colnames(x$fixed$Fval)
    tr <- vector("list", length(colnames.out.fixed))
    
    for(i in 1:length(colnames.out.fixed)){       
      tr[[i]] <- createTexreg(
        coef.names = names, se=x$fixed$Fval[i,],
        coef = x$fixed$Fval[i,],
        pvalues = x$fixed$pvalueF[i,],
        isRand=FALSE
      )     
    }
    caption.fixed = "F-test for the fixed effects"
  }
  else{
    caption.fixed = ""
    colnames.out.fixed <- ""
    tr <- NULL
  }
        
  
  
  if(("scaling" %in% names(x)) && (typeEffs == 3 || typeEffs == 4)){
    ## output for the scaling  effects if presented
    colnames.out.scaling <- rownames(x$scaling$FScaling)
    caption.scaling <- "F-test for the scaling effects"
    names <- colnames(x$scaling$FScaling)
    tr_scal <- vector("list", length(colnames.out.scaling))
    
    for(i in 1:length(colnames.out.scaling)){       
      tr_scal[[i]] <- createTexreg(
        coef.names = names, se=x$scaling$FScaling[i,],
        coef = x$scaling$FScaling[i,],
        pvalues = x$scaling$pvalueScaling[i,],
        isRand=FALSE
      )     
    }
    if(typeEffs == 3){
      regres <- list(lscale = tr_scal)
      custom.model.names =list(
        custom.model.names.scaling = colnames.out.scaling)
      caption2 <- list(caption.scaling = caption.scaling)
    }
    else{
      regres <- list(lrand = tr_rand, lfixed = tr, lscale = tr_scal)
      custom.model.names =list(custom.model.names.rand = colnames.out.rand,
                               custom.model.names.fixed = colnames.out.fixed,
                               custom.model.names.scaling = colnames.out.scaling)
      caption2 = list(caption.rand = caption.rand,
                     caption.fixed = caption.fixed,
                     caption.scaling = caption.scaling)
    }
  }
  else{
    if(typeEffs == 3)
      stop("There is no Scaling effect in the output")
    if(typeEffs == 1){
      custom.model.names =list(
        custom.model.names.rand = colnames.out.rand)
      caption2 = list(caption.rand = caption.rand)
      regres <- list(lrand = tr_rand)
    }
    if(typeEffs == 2){
      regres <- list(lfixed = tr)
      custom.model.names =list(
        custom.model.names.fixed = colnames.out.fixed)
      caption2 = list(caption.fixed = caption.fixed)
    }
    if(typeEffs == 4){
      regres <- list(lrand = tr_rand, lfixed = tr)
    custom.model.names =list(custom.model.names.rand = colnames.out.rand,
                             custom.model.names.fixed = colnames.out.fixed)
    caption2 = list(caption.rand = caption.rand,
                   caption.fixed = caption.fixed)
    }
  }  
   
  
  if(bold)
    stars <- numeric(0)
  else
    stars <- c(0.001, 
               0.01, 0.05)
  
  
  
  if(type == "html")
    htmlreg(regres, 
          file = file, inline.css = TRUE, 
          doctype = FALSE, html.tag = TRUE, head.tag = TRUE, 
          body.tag = TRUE,
          custom.model.names = custom.model.names, 
          caption = caption2, caption.above = TRUE, bold=bold,
          stars=stars, append = append)
  if(type == "latex")
    return(texreg(regres))
    
  
}

## output for the consmixed
.createDocOutputConsmixed <- function(x, file = NA, bold = FALSE, append = TRUE){
  sink(file = file, append = append)
  
  ## tests for the random effects
  x$rand.table[, "p.value"] <- format.pval(x$rand.table[,"p.value"],
                                           digits=3, eps=1e-3)
  x$rand.table <- .changeConsmixedOutputForDoc(x$rand.table, "p.value")
  if("elim.num" %in% colnames(x$rand.table))
    xt.rand <- xtable(x$rand.table, align="lcccc", 
                      display=c("s","f","d","s","s"))
  else
    xt.rand <- xtable(x$rand.table, align="lccc", 
                      display=c("s","f","d","s"))
  caption <- NULL
  #caption(xt.rand) <- "Likelihood ratio tests for the random-effects
  #and their order of elimination"
  print(xt.rand, caption.placement="top", table.placement="H",
        sanitize.text.function=function(x){x}, size="\\small", type = "html")
  
  ## tests for the fixed effects
  x$anova.table[, "Pr(>F)"] <- format.pval(x$anova.table[,"Pr(>F)"],
                                           digits=3, eps=1e-3)
  x$anova.table <- .changeConsmixedOutputForDoc(x$anova.table, "Pr(>F)")
  if("elim.num" %in% colnames(x$anova.table)) 
    xt.anova <- xtable(x$anova.table, align="lccccccc",
                       display=c("s","f", "f", "d", "f", "f", "s", "s"))     
  else
    xt.anova <- xtable(x$anova.table, align="lcccccc",
                       display=c("s","f", "f", "d", "f", "f","s"))
  #caption(xt.anova) <- "F-tests for the fixed-effects and their order of elimination"
  
  
  print(xt.anova, caption.placement="top", table.placement="H",
        sanitize.text.function=function(x){x}, size="\\small", type = "html")
  
  ## post hoc output
  x$diffs.lsmeans.table[, "p-value"] <- 
    format.pval(x$diffs.lsmeans.table[,"p-value"], digits=3, eps=1e-3)
  x$diffs.lsmeans.table <- 
    .changeConsmixedOutputForDoc(x$diffs.lsmeans.table, "p-value")    
  xt.lsmeans <- xtable(x$diffs.lsmeans.table, align="lccccccc",
                       display=c("s","f", "f", "f", "f", "f","f", "s"))
  #caption(xt.lsmeans) <- "Differences of Least Squares Means"
  print(xt.lsmeans, caption.placement="top", table.placement="H",
        sanitize.text.function=function(x){x}, size="\\small", type = "html")
  sink()
}
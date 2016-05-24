###### Here all the general functions and the functiona that are exported 
###### are defined 

sensmixed <- function(attributes=NULL, Prod_effects, replication = NULL,
                      individual, data, product_structure = 3,
                      error_structure = "No_Rep", MAM = FALSE, 
                      mult.scaling = FALSE, oneway_rand = TRUE, MAM_PER = FALSE, 
                      adjustedMAM = FALSE, alpha_conditionalMAM = 1, 
                      calc_post_hoc = FALSE, parallel = FALSE,
                      reduce.random = TRUE, alpha.random = 0.1,
                      alpha.fixed = 0.05, interact.symbol = ":", 
                      keep.effs = NULL, ...)
{  
  result <- sensmixedFun(attributes = attributes , Prod_effects, replication, 
                         individual, data,
                         product_structure = product_structure, 
                         error_structure = error_structure,  MAM = MAM, 
                         mult.scaling = mult.scaling,  
                         oneway_rand = oneway_rand,
                         MAM_PER = MAM_PER, adjustedMAM = adjustedMAM, 
                         alpha_conditionalMAM = alpha_conditionalMAM,
                         calc_post_hoc = calc_post_hoc,
                         parallel = parallel, reduce.random = reduce.random, 
                         alpha.random = alpha.random, 
                         alpha.fixed = alpha.fixed, 
                         interact.symbol = interact.symbol,
                         keep.effs = keep.effs)
  class(result) <- "sensmixed"
  result
}


print.sensmixed <- function(x, ...)
{
 
  tr_rand <- .changeOutput(x$random$Chi, x$random$pvalueChi, isRand = TRUE)
  cat("\nTests for the random effects:\n")
  screenreg(tr_rand, custom.model.names = names(tr_rand) )
  

  
  tr_fixed <- .changeOutput(x$fixed$Fval, x$fixed$pvalueF, isRand = FALSE)
  cat("\nTests for the fixed effects:\n")
  screenreg(tr_fixed, custom.model.names = names(tr_fixed) )
  
  if("scaling" %in% names(x)){
    tr_scaling <- .changeOutput(x$scaling$FScaling, x$scaling$pvalueScaling, 
                                FALSE)
    cat("\nTests for the scaling effects:\n")
    screenreg(tr_scaling, custom.model.names = names(tr_scaling) )
  }  
}  

plot.sensmixed <- function(x, mult = FALSE, dprime = FALSE, sep = FALSE, cex = 2,  
                           interact.symbol = ":", isFixed = TRUE, 
                           isRand = TRUE, isScaling = TRUE, ...)
{
  plotSensMixed(x, mult = mult, dprime = dprime, sep = sep, cex = cex, 
                interact.symbol = interact.symbol, isFixed = isFixed, 
                isRand = isRand, isScaling = isScaling)
}

saveToDoc <- function(x, file = NA, bold = FALSE, append = TRUE, type = "html",
                      typeEffs = 1)
{
  if(!(class(x) %in% c("sensmixed", "consmixed")))
    stop("x should be of class sensmixed")
  #if(is.na(file))
  #  stop("need to specify file")
  
  if(class(x)=="sensmixed"){
   return(.createDocOutputSensmixed(x, file = file, bold = bold, append = append,
                             type = type, typeEffs = typeEffs))
  }
  if(class(x) == "consmixed"){ 
   return(.createDocOutputConsmixed(x, file = file, bold = bold, append = append))
  }
}  


consmixed <- function(response, Prod_effects, Cons_effects=NULL, Cons, data, 
                      structure = 3, alpha.random = 0.1, alpha.fixed = 0.05, ...)
{  
  result <- consmixedFun(response=response, Prod_effects, 
                         Cons_effects=Cons_effects,
                         Cons, data, structure = structure, 
                         alpha.random = alpha.random,
                         alpha.fixed = alpha.fixed)
  class(result)<-"consmixed"
  result
}

plot.consmixed <- function(x, main = NULL, cex = 1.4, 
                           which.plot = c("LSMEANS", "DIFF of LSMEANS"),
                           effs = NULL, ...)
{
  st.x <- x
  class(st.x) <- "step"  
  plot(st.x, main = main, cex = cex, 
       which.plot = which.plot,
       effs = effs)
}



print.consmixed <- function(x, ...)
{
  if(!is.null(x$rand.table))
  {
    cat("\nRandom effects:\n") 
    x$rand.table[,"p.value"] <- format.pval(x$rand.table[,"p.value"], 
                                            digits=4, eps=1e-7)
    x$rand.table[,"Chi.sq"] <- round(x$rand.table[,"Chi.sq"],2)
    print(x$rand.table)   
  }
  
  if(nrow(x$anova.table) != 0)
  {
    if(class(x$model) == "lm" | class(x$model) == "gls")
    {
      cat("\nFixed effects:\n")
      print(x$anova.table)
      cat("\nLeast squares means:\n")
      print(x$lsmeans.table) 
      cat("\nFinal model:\n")
      print(x$model)
      return()
    }
    else
    {
      cat("\nFixed effects:\n")
      x$anova.table[,"Pr(>F)"] <- format.pval(x$anova.table[,"Pr(>F)"], 
                                              digits=4, eps=1e-7)
      x$anova.table[,c("Sum Sq","Mean Sq", "F.value")] <- 
        round(x$anova.table[,c("Sum Sq","Mean Sq", "F.value")],4)
      x$anova.table[,"DenDF"] <- round(x$anova.table[,"DenDF"],2)
      print(x$anova.table) 
      
      if(!is.null(x$lsmeans.table))
      {
        cat("\nLeast squares means:\n")
        printCoefmat(x$lsmeans.table, dig.tst=3,
                     tst.ind=c(1:(which(colnames(x$lsmeans.table)=="Estimate")-1),
                               which(colnames(x$lsmeans.table)=="DF")), digits=3,
                     P.values = TRUE, has.Pvalue=TRUE)
      }
      if(!is.null(x$diffs.lsmeans.table))
      {
        cat("\n Differences of LSMEANS:\n")
        printCoefmat(x$diffs.lsmeans.table, dig.tst=1,
                     tst.ind = 
                       c(1:(which(colnames(x$diffs.lsmeans.table)=="Estimate")-1),
                               which(colnames(x$diffs.lsmeans.table)=="DF")), 
                     digits = 3 ,P.values = TRUE, has.Pvalue = TRUE)
      }
      
    }    
  }
  else
    print(x$anova.table)
  cat("\nFinal model:\n")
  print(x$model@call) 
} 


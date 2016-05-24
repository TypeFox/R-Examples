
###############################################################################
#                             SHOW METHODS
################################################################################

.printaccess <- function() {
  cat("----------------------------------------------------------------------------\n")
  cat(" Access to slots:",
      "<ecoslot.> + <name of the slot> + <(name of the object)>","\n",
      "See help(\"EcoGenetics accessors\")\n")
  cat("----------------------------------------------------------------------------\n")
  
}

#-------------------------------------------------------------------#
#' show eco.gsa
#' @keywords internal
#' @rdname eco.gsa-methods
#' @aliases show,eco.gsa-method


setMethod("show", "eco.gsa", function(object) {
  method <- (object@METHOD)[1]
  cormethod <- (object@METHOD)[2]
  
  if(length(object@MULTI) != 0) {
    cat("\n", 
        "###########################","\n",
        paste(" ", method), "\n", 
        "###########################","\n\n")
    if(method == "Mantel test" | method == "Partial Mantel test") {
      cat("  >", "Correlation coefficent used:", 
          gsub("(^[[:alnum:]{1}])","\\U\\1", cormethod, perl =TRUE), "\n")
    }
    cat(
        paste("   >", "Number of simulations:", object@NSIM), "\n",
        paste("  >", "P-values correction method:", object@ADJUST), "\n\n",
        paste("  >", "Results:"), "\n\n")
    print(object@MULTI)
    cat("\n")
    
  } else {
    cat("\n", 
        "############################","\n",
        paste(" ", method), "\n", 
        "############################","\n\n")
    if(method == "Mantel test" | method == "Partial Mantel test") {
      cat("  >", "Correlation coefficent used:", 
          gsub("(^[[:alnum:]{1}])","\\U\\1", cormethod, perl =TRUE), "\n")
    }
    cat(
        paste("  >", "Number of simulations:", object@NSIM), "\n",
        paste(" >", "Alternative:", object@ALTER), "\n", 
        paste(" >", "P-value:", object@PVAL), "\n", 
        paste(" >", "Observed value:", object@OBS), "\n",
        paste(" >", "Expected value:", object@EXP), "\n")
    cat("\n")
  }
})


#-------------------------------------------------------------------#
#' show eco.lsa
#' @keywords internal
#' @rdname eco.autol-methods
#' @aliases show,eco.autol-method

setMethod("show", "eco.lsa", function(object)  {
  cat("\n", 
      "#########################","\n",
      paste(" ", object@METHOD), "\n", 
      "#########################","\n\n",
      paste(" >", "Test:", object@TEST), "\n",
      paste(" >", "Number of simulations:", object@NSIM), "\n",
      paste(" >", "Conditional:", object@COND), "\n")
  if(object@TEST == "permutation") {
    cat(paste("  >", "P-adjust method:", object@PADJ))
  }
  cat("\n\n", 
      paste(" Results :"), "\n\n")
  print(object@OUT)
  cat("\n Results table in", aue.access("OUT",  "name of this object"), "\n")
  .printaccess()
  cat("\n")
})


#-------------------------------------------------------------------#
#' show eco.correlog
#' @keywords internal
#' @rdname eco.correlog-methods
#' @aliases show,eco.correlog-method


setMethod("show", "eco.correlog", function(object) {
  
  randtest <- object@TEST
  method <- (object@METHOD)[1]
  cormethod <- (object@METHOD)[2]
  
  if(length(randtest != 0) & randtest != "none") {
  if(randtest == "permutation") {
    cat("\n", 
        "############################","\n",
        paste(" ", method), "\n", 
        "############################","\n\n")
    if(method == "Mantel statistic" | method == "Partial Mantel statistic") {
      cat(" >", "Correlation coefficent used:", 
      gsub("(^[[:alnum:]{1}])","\\U\\1", cormethod, perl =TRUE), "\n")
    }
    cat(paste("  >", "Number of simulations:", object@NSIM), "\n",
        paste(" >", "Random test:", object@TEST), "\n",
        paste(" >", "P-adjust method:", object@PADJUST), "\n\n", 
        paste(" >", "Results:","\n\n"))
    print(object@OUT)
    cat("\n Results table(s) in", aue.access("OUT",  "name of this object"), "\n")
    .printaccess()
    cat("\n")
 
     } else if (randtest == "bootstrap") {
 
  cat("\n", 
      "############################","\n",
      paste(" ", method), "\n", 
      "############################","\n\n")
       if(method == "Mantel statistic" | method == "Partial Mantel statistic") {
         cat("Correlation coefficent used: ", 
             gsub("(^[[:alnum:]{1}])","\\U\\1", cormethod, perl =TRUE), "\n")
       }
       cat(paste("   >", "Number of simulations: ", object@NSIM), "\n",
      paste("  >", " Random test:", object@TEST), "\n",
      paste("  >", " Results: ","\n\n"))
  print(object@OUT)
  cat("\n Results table(s) in", aue.access("OUT",  "name of this object"), "\n")
  .printaccess()
  cat("\n")
  
     }
  
  } else {
    cat("\n", 
        "############################","\n",
        paste(" ", method), "\n", 
        "############################","\n\n")
    if(method == "Mantel statistic" | method == "Partial Mantel statistic") {
      cat("Correlation coefficent used: ", 
          gsub("(^[[:alnum:]{1}])","\\U\\1", cormethod, perl =TRUE),  "\n")
    }
    cat(paste(" Results: ","\n\n"))
    print(object@OUT)
    cat("\n Results table(s) in", aue.access("OUT",  "name of this object"), "\n")
    .printaccess()
    cat("\n")
    
  }
  
  })


#-------------------------------------------------------------------#
#' show eco.weight
#' @keywords internal
#' @rdname eco.weight-methods
#' @aliases show,eco.weight-method

setMethod("show", "eco.weight", function(object)  {
  cat("\n", 
      "###################","\n",
      paste(" spatial weights"), "\n", 
      "###################","\n\n",
      paste(" Method ->", object@METHOD), "\n",
      " Parameters ->", paste("(", object@PAR, " = ",  
                              object@PAR.VAL, ")", sep =""), "\n",
      paste(" Row-standardization ->", object@ROW.SD), "\n")
  if(object@METHOD == "circle" | object@METHOD == "knearest") {
    cat(paste("  Self-included ->", object@SELF), "\n",
        paste(" Number of individuals ->", nrow(object@XY)), "\n",
        paste(" Non-zero (non-self) links ->", object@NONZERO, "%"), "\n",
        paste(" Individuals with non-zero (non-self) links ->", 
              object@NONZEROIND, "%"), "\n",
        paste(" Average (non-self) links per individual ->", object@AVG), "\n")
  } else {
    cat(paste("  Non-zero (non-self) links ->", object@NONZERO, "%"), "\n",
        paste(" Number of individuals ->", nrow(object@XY)), "\n",
        paste(" Individuals with non-zero (non-self) links ->", 
              object@NONZEROIND, "%"), "\n",
        paste(" Average (non-self) links per individual ->", object@AVG), "\n\n")
  }
})


#-------------------------------------------------------------------#
#' show eco.detrend
#' @keywords internal
#' @rdname eco.detrend-methods
#' @aliases show,eco.detrend-method

setMethod("show", "eco.detrend", function(object)  {
  
  cat("\n", 
      "###################","\n",
      paste(" ", "Data detrending"), "\n", 
      "###################","\n\n",
      paste(" >", "Polynomial degree ->", object@POLY.DEG), "\n",
      paste(" >", aue.access("RES", "name of this object"), ": Residuals-detrended data:", "\n\n"))
  print(head(object@RES))
  if(nrow(object@RES)> 6L) {
  cat(paste("(more data...)", "\n"))
  }
  cat("\n", "Other data:", "\n", 
      " >", aue.access("RES", "name of this object"), ": projected coordinates", "\n",
      " >", aue.access("XY", "name of this object"), ": models", "\n",
      " >", aue.access("ANALYSIS", "name of this object"), ": eco.mlm object with details", "\n",
      "\n")
  .printaccess()
  cat("\n")
})


#-------------------------------------------------------------------#
#' show eco.lagweight
#' @keywords internal
#' @rdname eco.lagweight-methods
#' @aliases show,eco.lagweight-method

setMethod("show", "eco.lagweight", function(object)  {
  cat("\n", 
      "######################","\n",
      paste(" spatial weights list"), "\n", 
      "######################","\n\n",
      " Parameters ->", paste("(", object@PAR, " = ",  
                              round(object@PAR.VAL, 3), ")", sep =""), "\n",
      paste(" Row-standardization ->", object@ROW.SD), "\n",
      paste(" Self-included ->", object@SELF), "\n",
      paste(" Cummulative ->", object@CUMMUL), "\n",
      paste(" Number of classes ->", length(object@MEAN)), "\n",
      paste(" Method ->", object@METHOD), "\n\n")
  
})


#-------------------------------------------------------------------#
#' show eco.mlm
#' @keywords internal
#' @rdname eco.lmtree-methods
#' @aliases show,eco.mlm-method


setMethod("show", "eco.mlm", function(object) {
  cat("\n", 
      "#########################","\n",
      paste(" multiple linear model"), "\n", 
      "#########################","\n\n")
  
  cat(paste("  >", paste(aue.access("MLM",  "name of this object"), ":", sep = ""), "multiple model results"), "\n", 
      paste(" >", paste(aue.access("SUMMARY.MLM",  "name of this object"), ":", sep = ""), "summary of the results"), "\n", 
      paste(" >", paste(aue.access("ANOVA.MLM",  "name of this object"), ":", sep = ""), "analysis of variance tables"), "\n",
      paste(" >", paste(aue.access("PREDICTED",  "name of this object"), ":", sep = ""), "predicted values"), "\n",			
      paste(" >", paste(aue.access("RESIDUALS",  "name of this object"), ":", sep = ""), "residuals of the analysis"), "\n\n")
  .printaccess()
  cat("\n")
  })


#-------------------------------------------------------------------#
#' show eco.mctree
#' @keywords internal
#' @rdname eco.lmtree-methods
#' @aliases show,eco.mctree-method


setMethod("show", 	"eco.mctree", function(object) {
  cat("\n", 
      "#############################","\n",
      paste(" multiple regression trees"), "\n", 
      "#############################","\n\n")
  
  cat(paste("  >", paste(aue.access("TREES",  "name of this object"), ":", sep = ""), "trees"), "\n", 
      paste(" >", paste(aue.access("CLASSPREDICT", "name of this object"), ":", sep = ""), "predictions of the analysis"), "\n", 
      paste(" >", paste(aue.access("FREQUENCIES", "name of this object"), ":", sep = ""), "number of individuals \n\t\t\t predicted in each node"), "\n",
      paste(" >", paste(aue.access("PREDICTED", "name of this object"), ":", sep = ""), "predicted values"), "\n",			
      paste(" >", paste(aue.access("RESIDUALS", "name of this object"), ":", sep = ""), "residuals of the analysis"), "\n\n")
  .printaccess()
  cat("\n")
  })



###############################################################################
#                             SUMMARY METHODS
################################################################################

#-------------------------------------------------------------------#
#' Summary for eco.lmtree output
#' @param object Output object of \code{\link{eco.lmtree}}.
#' @return A Table with a summary of the analysis for "mlm" analysis, 
#' the plot of the trees with significant splits for "mctree" analysis.
#' @seealso \code{\link{eco.lmtree}}
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' data(eco.test)
#' #' mod <- eco.lmtree(DF1 = eco$P, DF2 = eco$E, 
#' analysis = "mlm")                                    
#' summary(mod)                                    #summary for "mlm" analysis
#' 
#' mod <- eco.lmtree(DF1 = eco$P, DF2 = eco$E,
#' analysis = "mctree", fact = eco$S$structure)               
#' summary(mod)                                    #summary for "mctree" analysis
#' 
#' }
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @rdname eco.lmtree.mctree-summary
#' @aliases summary,eco.lmtree-method
#' @exportMethod summary

setMethod("summary", 	"eco.mlm", 	function(object) {
  
  anova.pres <- function(mod) {
    
    asteriscos <- function(x, y) {
      
      if((x < 0.05) && (x >= 0.01)) {
        y <- paste(y, "*")
      } else if((x < 0.01) && (x >= 0.001)) {
        y <- paste(y, "**")
      } else if(x < 0.001) {
        y <- paste(y, "***")
      } else if(x >= 0.05) {
        y <- paste(y, "ns")
      }
    }
    
    
    final <- as.data.frame(matrix(ncol = length(object@MLM), 
                                  nrow = ncol(object@DF2) + 1))
    rownames(final) <- c("Corrected model", colnames(object@DF2))
    colnames(final) <- colnames(object@DF1)
    pvalor <- rep(0, ncol(object@DF1))
    
    for(i in 1:length(object@MLM)) { 
      pvalor[i] <- object@ANOVA.MLM[[i]][, 5][1]   
    }
    pvalor <- simplify2array(pvalor)
    pvalor <- p.adjust(pvalor, method = "fdr")
    
    
    for(i in 1:length(object@MLM)) {
      
      an <- object@ANOVA.MLM[[i]]
      an <- as.data.frame(an)
      mm <- object@SUMMARY.MLM[[i]]
      if(is.null(mm$fstatistic)) {
        final[1, i] <- "ns"
        
      } else {
        
        final[1, i] <- paste("F", "(", mm$fstatistic[2],
                             ", ", mm$fstatistic[3], ")","",
                             "=", round(mm$fstatistic[1],
                                        3), sep = "")
        final[1, i] <- asteriscos(pvalor[i], final[1, i])
        
      }
      
      for(j in 1:(nrow(an))) {
        
        solapa <- match(rownames(an), rownames(final))
        if(is.na(solapa[j]) == FALSE) {
          final[solapa[j], i] <- paste("F", "(", an[, 1][j],
                                       ", ", an[, 1][nrow(an)], ")",
                                       "= ", round(an[, 4][j], 3),
                                       sep= "")
          final[solapa[j], i] <- asteriscos(an[, 5][j], 
                                            final[solapa[j], i])
        }
      }
    }
    final[is.na(final) == TRUE] = "-"
    final
  }
  
  cat("\n", "ANOVA's F statistics, degrees of freedom and P-VALUES", "\n\n")
  print(anova.pres(object@MLM))
  cat("ns: non significative, *P<0.05, **P<0.01, ***P<0.001, ",
      "- null model")
  
})


#-------------------------------------------------------------------#
#' @rdname eco.lmtree.mctree-summary
#' @aliases summary,eco.mctree-method
#' @exportMethod summary

setMethod("summary", 	"eco.mctree", 	function(object) {
  
  
  count <-0
  for(i in seq(along = object@TREES))
  {
    where.lev <- object@TREES[[i]]@where
    where.lev <- as.factor(where.lev)
    plot.tre<- max(levels(where.lev[[i]])) > 1
    if(plot.tre) {
      count <- count +1
    }
    if(plot.tre) {
      plot(object@TREES[[i]], main = names(object@TREES)[i])
    }
  }
  if(count == 0) {
    cat("\n", "There are not trees with significant splits to plot", "\n\n")
  }
  
})


#################################################

#-------------------------------------------------------------------#
#' show eco.IBD
#' @keywords internal
#' @rdname eco.IBD-methods
#' @aliases show,eco.IBD-method
#' @export


setMethod("show", "eco.IBD", function(object) {
  
  randtest <- (object@TEST)[1]
  metodo <- (object@METHOD)[2]
  
  cat("\n", 
      "##################","\n",
      " Kinship analysis", "\n", 
      "##################","\n\n",
      paste(" >","Method:", metodo), "\n",
      paste(" >","Number of simulations:", object@NSIM), "\n")
      if(randtest == "permutation") {
      cat(paste("  >","P adjust method:", object@PADJUST), "\n\n")
      } 
      if(metodo == "local") {
      cat(paste(" >", "Conditional:", (object@TEST)[2]), "\n\n")
        } 
      cat(
      "  Available information: ", "\n",
      paste(" >", paste(aue.access("SP", "name of this object"), ":", sep = ""), "SP analysis"), "\n", 
      paste(" >", paste(aue.access("OUT", "name of this object"), ":", sep = ""), "table with results"),"\n")
      cat("\n Results table in", aue.access("OUT",  "name of this object"), "\n")
      .printaccess()
      cat("\n")
})

#-------------------------------------------------------------------#


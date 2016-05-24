step <- function(model, ddf="Satterthwaite", type=3, alpha.random = 0.1, 
                 alpha.fixed = 0.05, reduce.fixed = TRUE, reduce.random = TRUE, 
                 fixed.calc=TRUE ,lsmeans.calc=TRUE, difflsmeans.calc=TRUE, 
                 test.effs=NULL, keep.effs = NULL,...)
{  
  if(!inherits(model, "lmerMod"))
    stop("The model is not linear mixed effects model")
  
  ddfs <- c("Satterthwaite", "Kenward-Roger")
  ind.ddf <- pmatch(tolower(ddf), tolower(ddfs))
  if(is.na(ind.ddf))  
    stop('Parameter ddf is wrongly specified')  
  else
    ddf <- ddfs[ind.ddf]
  
  result <- totalAnovaRandLsmeans(model = model, ddf = ddf , type = type,  
                                  alpha.random = alpha.random, 
                                  alpha.fixed = alpha.fixed,
                                  reduce.fixed = reduce.fixed, 
                                  reduce.random = reduce.random,
                                  fixed.calc = fixed.calc, 
                                  lsmeans.calc = lsmeans.calc,
                                  difflsmeans.calc = difflsmeans.calc, 
                                  isTotal = TRUE, 
                                  isTtest = FALSE, test.effs = test.effs,
                                  keep.effs = keep.effs)
  class(result) <- "step"
  result
}



print.step <- function(x, ...)
{
  
  if(!is.null(x$rand.table))
  {
    cat("\nRandom effects:\n") 
    x$rand.table[,"p.value"] <- format.pval(x$rand.table[,"p.value"], digits=4, 
                                            eps=1e-7)
    x$rand.table[,"Chi.sq"] <- round(x$rand.table[,"Chi.sq"], 2)
    print(x$rand.table)     
  } 
  if(is.null(x$anova.table)){
    
  }else{
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
          printCoefmat(x$lsmeans.table, dig.tst=3 ,
                       tst.ind=c(1:(which(colnames(x$lsmeans.table)=="Estimate")-1),
                                 which(colnames(x$lsmeans.table)=="DF")), 
                       digits=3 , P.values = TRUE, has.Pvalue=TRUE)
        }
        if(!is.null(x$diffs.lsmeans.table))
        {
          cat("\n Differences of LSMEANS:\n")
          printCoefmat(x$diffs.lsmeans.table, dig.tst=1  ,
                       tst.ind=c(1:(which(colnames(x$diffs.lsmeans.table)==
                                            "Estimate")-1),
                                 which(colnames(x$diffs.lsmeans.table)=="DF")),
                       digits=3 , P.values = TRUE, has.Pvalue=TRUE)
        }
        
      }    
    }
    else
      print(x$anova.table)
  }
  
  cat("\nFinal model:\n")
  print(x$model@call) 
}


plot.step <- function(x, main = NULL, cex = 1.4, 
                      which.plot = c("LSMEANS", "DIFF of LSMEANS"),
                      effs = NULL, mult = TRUE, ...)
{
  if(!is.null(x$lsmeans.table) && nrow(x$lsmeans.table)>0 && ("LSMEANS" %in% which.plot)){
    if(length(which.plot) == 1 && which.plot == "LSMEANS")
      return(plotLSMEANS(x$lsmeans.table, x$response, "LSMEANS", main = main, cex = cex,
                         effs = effs, mult = mult))
  }         
  if(!is.null(x$diffs.lsmeans.table) && nrow(x$diffs.lsmeans.table)>0 
     && ("DIFF of LSMEANS" %in% which.plot))
    plotLSMEANS(x$diffs.lsmeans.table, x$response, "DIFF of LSMEANS", 
                main = main, cex = cex, effs = effs, mult = mult)
}


lmer <- function(formula, data = NULL, REML = TRUE,
                 control = lmerControl(), start = NULL, verbose = 0L,
                 subset, weights, na.action, offset, contrasts = NULL,
                 devFunOnly = FALSE, ...)
{
  mc <- match.call()
  mc[[1]] <- quote(lme4::lmer)
  model <- eval.parent(mc)
  if(inherits(model, "merMod"))
    model <- as(model,"merModLmerTest")    
  return(model)
}




setMethod("anova", signature(object="merModLmerTest"),
          function(object, ..., ddf="Satterthwaite", type=3)  
          {
            
            mCall <- match.call(expand.dots = TRUE)
            dots <- list(...)
            modp <- if (length(dots))
              sapply(dots, is, "merModLmerTest") | sapply(dots, is, "merMod") | 
              sapply(dots, is, "lm") else logical(0)
            if (any(modp)) {
              return(callNextMethod())
            }
            else
            {
              cnm <- callNextMethod()
              if(!is.null(ddf) &&  ddf=="lme4") 
                return(cnm)              
{
                  table <- cnm 
                  
                  #errors in specifying the parameters
                  ddfs <- c("Satterthwaite", "Kenward-Roger")
                  ind.ddf <- pmatch(tolower(ddf), tolower(ddfs))
                  if(is.na(ind.ddf))  
                    stop('Parameter ddf is wrongly specified')  
                  else
                    ddf <- ddfs[ind.ddf]
                  
                  an.table <- tryCatch({totalAnovaRandLsmeans(model=object, 
                                                              ddf=ddf, 
                                                              type=type,
                                                              isAnova=TRUE, 
                                                              reduce.random=FALSE,
                                                              reduce.fixed=FALSE)$anova.table}
                                       , error = function(e) { NULL })
                  if(!is.null(an.table))
                  {
                    table <- an.table
                    
                    attr(table, "heading") <- 
                      paste("Analysis of Variance Table of type", as.roman(type) ,
                            " with ", ddf, 
                            "\napproximation for degrees of freedom")
                  }
                  else
                    message("anova from lme4 is returned\nsome computational error has occurred in lmerTest")
                  
                  
                  
                  class(table) <- c("anova", "data.frame")
                  return(table)
                }  

            }

          })

setMethod("summary", signature(object = "merModLmerTest"),
          function(object, ddf="Satterthwaite", ...)
          {
                      
            if(!is.null(ddf) && ddf=="lme4"){
              if(class(object) == "merModLmerTest")
                return(summary(as(object, "lmerMod")))
              #return(cl)
            }else{
              ## commented callNextMethod
              ## since it produces warning, summary cannot have multiple arguments
              ##cl <- callNextMethod()
              if(class(object) == "merModLmerTest")
                cl <- summary(as(object, "lmerMod"))
              #errors in specifying the parameters
              ddfs <- c("Satterthwaite", "Kenward-Roger")
              ind.ddf <- pmatch(tolower(ddf), tolower(ddfs))
              if(is.na(ind.ddf))  
                stop('Parameter ddf is wrongly specified')  
              else
                ddf <- ddfs[ind.ddf]
              
              tsum <- tryCatch( {totalAnovaRandLsmeans(model=object, 
                                                       ddf=ddf, 
                                                       isTtest=TRUE)$ttest}, 
                                error = function(e) { NULL })
              if(is.null(tsum)){
                message("summary from lme4 is returned\nsome computational error has occurred in lmerTest")
                return(cl)
              }
              coefs.satt <- cbind(cl$coefficients[,1:2, drop = FALSE], tsum$df, 
                                  tsum$tvalue, tsum$tpvalue)               
              cl$coefficients <- coefs.satt
              colnames(cl$coefficients)[3:5] <- c("df","t value","Pr(>|t|)")              
            }   
            
            cl$methTitle <- paste(cl$methTitle,  "\nt-tests use ", ddf, 
                                  "approximations to degrees of freedom")
            return(cl)
          }
          
)

#randTAB.default<-function(model, data, ...)
rand <- function(model, ...)
{
  if(!inherits(model, "lmerMod"))
    stop("The model is not linear mixed effects model")
  result <- totalAnovaRandLsmeans(model=model, isRand=TRUE, reduce.random=FALSE)  
  res <- list(rand.table=result$rand.table, isCorr = result$corr.intsl)
  class(res) <- "rand"
  res
}

print.rand <- function(x, ...)
{
  
  cat("Analysis of Random effects Table:\n")
  if(!is.null(x$rand.table))
    printCoefmat(x$rand.table, digits=3 , dig.tst=1  ,
                 tst.ind=c(which(colnames(x$rand.table)=="Chi.DF"),
                           which(colnames(x$rand.table)=="elim.num")), 
                 P.values=TRUE, has.Pvalue=TRUE)        
}





lsmeans <- function(model, test.effs=NULL, ...)
{
  if(!inherits(model, "lmerMod"))
    stop("The model is not linear mixed effects model")
  result <- totalAnovaRandLsmeans(model = model, ddf = "Satterthwaite", 
                                  isLSMEANS = TRUE, test.effs = test.effs, 
                                  reduce.random = FALSE, reduce.fixed = FALSE)  
  res <- list(lsmeans.table=result$lsmeans.table, response=result$response)
  class(res) <- "lsmeans"
  res 
}

print.lsmeans <- function(x, ...)
{
  
  cat("Least Squares Means table:\n")
  printCoefmat(data.matrix(x$lsmeans.table), dig.tst=1, 
               tst.ind=c(1:(which(colnames(x$lsmeans.table)=="Estimate")-1),
                         which(colnames(x$lsmeans.table)=="DF")), digits=3 , 
               P.values=TRUE, has.Pvalue=TRUE)       
}

plot.lsmeans <- function(x, main = NULL, cex = 1.4, effs = NULL, mult = TRUE, ...)
{
  
  #plots for LSMEANS
  if(!is.null(x$lsmeans.table) && nrow(x$lsmeans.table)>0)
    plotLSMEANS(x$lsmeans.table, x$response, "LSMEANS", main = main, cex = cex,
                effs = effs,  mult = mult)     
}

difflsmeans <- function(model, test.effs=NULL, ...)
{
  if(!inherits(model, "lmerMod"))
    stop("The model is not linear mixed effects model")
  result <- totalAnovaRandLsmeans(model = model, ddf = "Satterthwaite", 
                                  isDiffLSMEANS = TRUE, test.effs = test.effs, 
                                  reduce.random = FALSE, reduce.fixed = FALSE)  
  res <- list(diffs.lsmeans.table=result$diffs.lsmeans.table, 
              response=result$response)
  class(res) <- "difflsmeans"
  res 
}

print.difflsmeans <- function(x, ...)
{
  
  cat("Differences of LSMEANS:\n")
  printCoefmat(data.matrix(x$diffs.lsmeans.table), dig.tst=1, 
               tst.ind=c(1:(which(colnames(x$diffs.lsmeans.table)=="Estimate")-1),
                         which(colnames(x$diffs.lsmeans.table)=="DF")), digits=3 ,
               P.values=TRUE, has.Pvalue=TRUE)
  
}

plot.difflsmeans <- function(x, main = NULL, cex = 1.4, effs = NULL, 
                             mult = TRUE, ...)
{
  
  #plots for DIFF of LSMEANS
  if(!is.null(x$diffs.lsmeans.table) && nrow(x$diffs.lsmeans.table)>0)
    plotLSMEANS(x$diffs.lsmeans.table, x$response, "DIFF of LSMEANS", 
                main = main, cex = cex, effs = effs, mult = mult)   
}
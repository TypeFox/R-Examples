stemPlot.default <-
function(x, grpvar = NULL, varname = NULL, grpvarname = NULL, ...)
{     
    if (!is.numeric(x)) stop("variable must be numeric")
    
    if (is.null(varname))
       varname <- deparse(substitute(x))
       
    comCases <- complete.cases(x)
    nmiss <- length(x) - sum(comCases)
    if (nmiss > 0)
      cat("\n ", nmiss, "observation(s) removed due to missing values.\n")
              
     cat("\n***Stem and Leaf plot for ", varname, "***\n")
    
     if (!is.null(grpvar))
       {        
       
         if (is.null(grpvarname)) 
          grpvarname <- deparse(substitute(grpvar))  
         
         if (!is.factor(grpvar)) 
          {grpvar <- as.factor(grpvar)
          } else grpvar <- droplevels(grpvar)
         
         grp.levels <- levels(grpvar)    
       
       cat("   Grouped by levels of ", grpvarname, "\n")          
        for (i in grp.levels)
         {
          cat("\n   " ,i,"\n :")
          stem(x[grpvar == i],...)
    
         } #end for
         } else  { stem(x, ...)}
     
  
  invisible() 
}

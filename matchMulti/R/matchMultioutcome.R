matchMultioutcome <- function(obj, out.name = NULL, schl_id_name = NULL, treat.name = NULL, end.1=-1000, end.2=1000){
	
	
   pval.c <- pval_func(obj, out.name = out.name, schl_id_name = schl_id_name,  treat.name = treat.name, wt=FALSE)
   pval.p <- pval_func(obj, out.name = out.name, schl_id_name = schl_id_name,  treat.name = treat.name, wt=TRUE)
          
   ci1 <- uniroot(ci_func, c(end.1, end.2), obj=obj, out.name = out.name,
          schl_id_name = schl_id_name, treat.name = treat.name, alternative="less", alpha=.025)$root;
          
   ci2 <- uniroot(ci_func, c(end.1, end.2), obj=obj, out.name = "mathach",
          schl_id_name = "school",  treat.name = "sector", alternative="great", alpha=.025)$root; 
          
   pe <- uniroot(pe_func, c(end.1, end.2), obj=obj, out.name = "mathach",
          schl_id_name = "school", 
          treat.name = "sector")$root;   
          
   ci.lo <- min(ci1, ci2)
   ci.up <- max(ci1, ci2) 
   ci <- c(ci.lo, ci.up)          
     
	cat("Point Estimate is: ", pe, "\n")
	cat("95 confidence interval", ci,"\n")
	cat("One.sided p-value is: ", pval.c, "\n")
	res <- list(pval.c = pval.c, pval.p = pval.p, ci1=ci1, ci2=ci2, p.est=pe)
	return(res)
}



     
     
     


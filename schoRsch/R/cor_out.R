# -----------------------------------------------
# Assembles summary line for correlation analyses
# Requires output object from cor.test as input
# argument.
# -----------------------------------------------
cor_out <- function(coroutput, stats = FALSE, print = TRUE) {
  
  # ---------------------------------------------
  # (1) Assemble summary table
  # ---------------------------------------------

  if (coroutput$method == "Pearson's product-moment correlation") 
  {
    outtable <- data.frame(
    coefficient=format(round(coroutput$estimate,2),nsmall = 2),
    n = coroutput$parameter+2,
    p=format(round(coroutput$p.value,3),nsmall = 3)
    )
    
    
    if (stats == TRUE)  # wenn Statistik mit ausgegeben werden soll
    {
       statout = paste(", t(" , coroutput$parameter , ") = " , 
                       format(round(coroutput$statistic,2),nsmall = 2),sep="")  
    } else
    {
      statout = ""
    }
  } else
  {
    outtable <- data.frame(
    coefficient=format(round(coroutput$estimate,2),nsmall = 2),
    p=format(round(coroutput$p.value,3),nsmall= 3)
    )
    
    if (stats == TRUE)  # wenn Statistik mit ausgegeben werden soll
    {
      if (coroutput$method == "Kendall's rank correlation tau")
      {
        statout = paste(", z = " , format(round(coroutput$statistic,2),nsmall=2),sep="")  
      } else if (coroutput$method == "Spearman's rank correlation rho")
      {
        statout = paste(", S = " , format(round(coroutput$statistic,2),nsmall= 2),sep="")  
      }
    } else
    {
      statout = ""
    }
  }
  
  # ---------------------------------------------
  # (2) Format output table
  # ---------------------------------------------
  ####
  
  
  if (coroutput$method == "Pearson's product-moment correlation") 
  {  
     pcorr <- paste(", p = ", outtable$p, sep="")
     pcorr <- gsub("p = 1.000","p > .999", pcorr, fixed=TRUE)
     pcorr <- gsub("p = 0.000","p < .001", pcorr, fixed=TRUE)
     pcorr <- gsub("p = 0","p = ", pcorr, fixed=TRUE)
    
  outtext <- data.frame(
  Text=paste("r(",outtable$n,") = ", outtable$coefficient,statout, pcorr,sep=""));
  }  else if (coroutput$method == "Kendall's rank correlation tau")
  {
     pcorr <- paste(", p = ", outtable$p, sep="")
     pcorr <- gsub("p = 1.000","p > .999", pcorr, fixed=TRUE)
     pcorr <- gsub("p = 0.000","p < .001", pcorr, fixed=TRUE)
     pcorr <- gsub("p = 0","p = ", pcorr, fixed=TRUE)
    
     outtext <- data.frame(
     Text=paste("tau = " , outtable$coefficient,statout, pcorr, sep="")); 
  }  else if (coroutput$method == "Spearman's rank correlation rho")
  {
    
     pcorr <- paste(", p = ", outtable$p, sep="")
     pcorr <- gsub("p = 1.000","p > .999", pcorr, fixed=TRUE)
     pcorr <- gsub("p = 0.000","p < .001", pcorr, fixed=TRUE)
     pcorr <- gsub("p = 0","p = ", pcorr, fixed=TRUE)
     
     outtext <- data.frame(
     Text=paste("rho = " , outtable$coefficient, statout,pcorr,sep=""));
  }
  
  if (print==TRUE) {
    print(outtext);
  } else {
    outtext;  
  }
}

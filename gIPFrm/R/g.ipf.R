g.ipf <-
function(ModelMatrix, ObsTable, tol, estimand, adjustment)
{
  nc <- ncol(ModelMatrix);
  nr <- nrow(ModelMatrix);
  if(sum(ObsTable == 0)>0) 
     warning("Some of the observed values are zero; the model parameters may not be identifiable");
  if(length(ObsTable) != nc) 
     stop("Dimensions of the model matrix and the data vector do not match");
  est <- (estimand == "probabilities") + (estimand == "intensities");
  if( est != 1) stop("The estimand is not specified correctly.")
  adj <- (adjustment == "grid") + (adjustment == "bisection") + 
         (adjustment == "none");
  if( adj != 1) stop("The adjustment method is not specified correctly.")

  mleTable <- NULL;
  F <- NULL;
  iter <- 0;
  if (estimand == "intensities")
  {
     F <- ipf.gamma(ModelMatrix, ObsTable, 1, tol, "intensities");
     mleTable <- round(F$fitted.values,2);
     model.parameters <- round(F$model.parameters,2);
     }  else   {      
     F <- ipf.gamma(ModelMatrix, ObsTable, 1, tol, "probabilities")
     p <- F$fitted.values; 
     model.parameters <- F$parameters;
     p.sum <- sum(p); 
     obs.sum <- sum(ObsTable);
     Total <- p.sum/obs.sum;
     gamma.hat <- 1;
     if (abs(Total - 1) > tol)
     { 
           b <- suff.stat(ModelMatrix, ObsTable/obs.sum);
           
           gamma.left <- 1/(sum(b));
           gamma.right <- min(1/b); 
           if( adjustment == "grid")
           {
               UpdateGamma <-  grid.update(ModelMatrix, ObsTable, tol)
                                
           }
           else
           {
               UpdateGamma <- bisection.update(ModelMatrix, ObsTable, tol) 
           }

     
          
             
           F <- UpdateGamma$model.tilde;
           gamma.hat <- UpdateGamma$gamma.tilde;  
         }
        model.parameters <- F$model.parameters; 
        mleTable <- round(F$fitted.values,4);
      }
      chisqv <-  sum( (ObsTable-mleTable)^2/mleTable);
      LLratio <- 0;
         
      for(i in 1:nc)
      {
         if(ObsTable[i] > 0) 
           {LLratio <- LLratio + ObsTable[i]*log(ObsTable[i]/mleTable[i])};
      }
      LLratio <- 2*LLratio;

      df <- nc - qr(ModelMatrix)$rank;
      pv1 <- 1-pchisq(chisqv,df); 
      pv2 <- 1-pchisq(LLratio,df);

      if(estimand == "probabilities")
      {
         result <- list(model.matrix = ModelMatrix,
                     observed.data = ObsTable,
                     fitted.values = mleTable,
                     adjustment.for.subsets = round(gamma.hat,4), 
                     model.parameters = round(model.parameters,4), 
                     degrees.of.freedom = df,
                     chisq.statistic = round(chisqv,2),
                     p.value.chisq = round(pv1,2),
                     log.likelihood.ratio.statistic = round(LLratio,2),
                     p.value.log.likelihood.ratio = round(pv2,2))
      }
      else
      {
           if (sum(mleTable)==sum(ObsTable))
           {
              result <- list(model.matrix = ModelMatrix,
                        observed.data = ObsTable,
                        fitted.values = mleTable,
                        estimated.total = sum(mleTable),
                        adjustment.for.total = 1, 
                        model.parameters = round(model.parameters,4), 
                        degrees.of.freedom = df,
                        chisq.statistic = round(chisqv,2),
                        p.value.chisq = round(pv1,2),
                        log.likelihood.ratio.statistic = round(LLratio,2),
                        p.value.log.likelihood.ratio = round(pv2,2))
           }
           else
           {  
              warning("No p-values are produced. The distributions of the chisq statistic and the log likelihood ratio statistic are unknown");
              result <- list(model.matrix = ModelMatrix,
                        observed.data = ObsTable,                        fitted.values = mleTable,
                        estimated.total = sum(mleTable),
                        adjustment.for.total = round(sum(mleTable)/sum(ObsTable),4), 
                        model.parameters = round(model.parameters,4), 
                        degrees.of.freedom = df,
                        chisq.statistic = round(chisqv,2),
                        log.likelihood.ratio.statistic = round(LLratio,2))
            }
      }
   return(result)
}


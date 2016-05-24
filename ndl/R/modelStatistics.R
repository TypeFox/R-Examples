modelStatistics <- function(observed, predicted, frequency=NA, p.values, n.data, n.predictors, outcomes=levels(as.factor(observed)), 
    p.normalize=TRUE, cross.tabulation=TRUE, p.zero.correction=1/(NROW(p.values)*NCOL(p.values))^2)
{ if(p.zero.correction==0) warning("Loglikelihood and related statistics may be inestimable, if P=0 for any observed outcome.");
  N <- length(observed);
#  if(any(!(unique(observed) %in% colnames(p.values))))
#    p.values <- cbind(p.values, matrix(0,NROW(p.values),1,dimnames=list(NULL,setdiff(unique(observed),colnames(p.values)))))
  if(p.normalize)
    p.values <- p.values/apply(p.values,1,sum)
  p.outcomes <- colnames(p.values);
  d <- sapply(1:N, function(i)
    { p <- p.values[i,which(observed[i]==p.outcomes)]
      if(p==0) p = p.zero.correction
      if(all(is.na(frequency)))
        return(log(p))
      else
        return(as.vector(frequency)[i]*log(p))
    });
  loglikelihood.model <- sum(d);
  deviance.model <- -2 * loglikelihood.model;
  if(all(is.na(frequency)))
    n.outcomes <- sapply(outcomes, function(o) length(which(observed==o)))
  else
    n.outcomes <- sapply(outcomes, function(o) sum(as.vector(frequency)[which(observed==o)]))
  loglikelihood.null <- sum(sapply(outcomes, function(o) n.outcomes[o]*log(n.outcomes[o]/n.data)));
  deviance.null <- -2 * loglikelihood.null;
  R2.likelihood <- 1 - deviance.model/deviance.null;
#  R2.nagelkerke <- (1-(exp(loglikelihood.null)/exp(loglikelihood.model))^(2/n.data))/(1-(exp(loglikelihood.null)^(2/n.data)));
  R2.nagelkerke <- (1-exp(-2*(loglikelihood.model-loglikelihood.null)/n.data))/(1-exp(2*loglikelihood.null/n.data))
  AIC.model <- 2 * n.predictors - 2 * loglikelihood.model
  BIC.model <- n.predictors * log(n.data) - 2 * loglikelihood.model

  # calculate index of concordance C if response variable is binary
  if (length(outcomes) == 2)
     {
       if(all(is.na(frequency)))
         { predval = sort(unique(predicted))[1]
           binvec = as.numeric(as.character(observed)==predval)
           C.statistics = Hmisc::somers2(p.values[,predval], binvec)
         }
       else
         { predval = sort(unique(predicted))[1]
           binvec = rep(as.numeric(as.character(observed)==predval), frequency)
           C.statistics = Hmisc::somers2(rep(p.values[,predval],frequency), binvec)
         }

       statistics <- list(loglikelihood.null = loglikelihood.null, loglikelihood.model = loglikelihood.model, deviance.null = deviance.null, deviance.model = deviance.model, R2.likelihood = R2.likelihood, R2.nagelkerke = R2.nagelkerke, AIC.model = AIC.model, BIC.model = BIC.model, C=C.statistics[["C"]]);

     }
  else 

     statistics <- list(loglikelihood.null = loglikelihood.null, loglikelihood.model = loglikelihood.model, deviance.null = deviance.null, deviance.model = deviance.model, R2.likelihood = R2.likelihood, R2.nagelkerke = R2.nagelkerke, AIC.model = AIC.model, BIC.model = BIC.model)

  # Statistics based on crosstabulated of observed vs. predicted outcomes

  if(cross.tabulation)
    if(all(is.na(frequency)))
      { crosstable <- table(factor(observed, levels=outcomes), factor(predicted, levels=outcomes));
        statistics <- c(statistics, list(crosstable = crosstable), crosstableStatistics(crosstable))
      }
    else
      { crosstable <- matrix(0, length(outcomes), length(outcomes), dimnames=list(outcomes,outcomes));
        for(i in 1:N)
           crosstable[observed[i],predicted[i]] = crosstable[observed[i],predicted[i]] + as.vector(frequency)[i];
        statistics <- c(statistics, list(crosstable = crosstable), crosstableStatistics(crosstable))
      }

  return(statistics);

}

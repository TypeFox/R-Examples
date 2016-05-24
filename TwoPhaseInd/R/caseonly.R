

caseonly <- function(data,treatment,BaselineMarker,extra=NULL,fraction=0.5) {
  
  ### data: a data matrix with cases only
  ### BaselineMarker: baseline biomarker 
  ### extra: the extra covariates that may be adjusted for in case-only data
  ### fraction: the randomization fraction of active treatment assignment 
  
  if (!is.data.frame(data)) {
    stop("Argument data must be a data.frame object.")
  }
  else {
    colNames <- colnames(data)
    if (!(treatment %in% colNames)) {
      stop("Treatment variable was not found in the data.")
    }else {
      if (any(levels(factor(data[, treatment])) != c("0", 
                                                     "1"))) {
        warning("Treatment variable must be either 0 or 1 only.")
      }
    }
    if (!(BaselineMarker %in% colNames)) {
      stop("BaselineMarker variable was not found in the data.")
    }
    if (!is.null(extra)) {
      if (!any(extra %in% colNames)) {
        extraNotFound <- paste(extra[!(extra %in% colNames)], 
                               sep = "", collapse = ", ")
        stop(paste("Extra variable(s) was not found in the data:", 
                   extraNotFound))
      }
      tmp <- remove_rarevariants(data[, extra])
      if (any(tmp)) {
        idx <- tmp == TRUE
        toremove = NULL
        for (i in 1:length(idx)) {
          if (idx[i]) 
            toremove <- c(toremove, extra[idx[i]])
        }
        warnings(paste0(paste(toremove, sep = ", "), 
                        " were removed due to rare vairant"))
        extra <- extra[!idx]
        if (length(extra) == 0) 
          extra <- NULL
      }
    }  
  }
  #remove missing data
  dat <- remove_missingdata(data)$data
  if(!is.null(extra)){
    glmFormula <- as.formula(paste(treatment, "~",
                                   paste(BaselineMarker,
                                         ## extra variables
                                         paste(extra, collapse=" + "),
                                         sep=" + "
                                        )
                                  )
                            )
  }else{
    glmFormula <- as.formula(paste(treatment, "~", BaselineMarker,sep= " "))
  }
  
  fit <- glm(glmFormula,data=dat,family=binomial,offset=rep(log(fraction/(1-fraction)),nrow(dat)))
             
  ## RETURN THE RESULTS 
  ## first row: treatment effect in baselineMarker =0
  ## second row: treatment+BaselineMarker interaction
  tmpResult=summary(fit)$coefficients
  tmpResult=tmpResult[,-3]
  colnames(tmpResult)=c("beta","stder","pVal")
  tmpResult=tmpResult[1:2,]
  rownames(tmpResult)=c("treatment effect when baselineMarker=0","treatment+baselineMarker interaction")
  return(tmpResult)
             
             
}
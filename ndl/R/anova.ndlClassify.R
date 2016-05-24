# Comparison of lists of models or individual model predictors for single model

anova.ndlClassify <- function(object, ...,  statistic = "deviance", test = "Chisq")
{   if(statistic!="deviance")
      stop(paste(c("ANOVA for Naive Discriminative Learning not implemented for statistic: ",statistic),collapse=""))
    dotargs <- list(...)
    named <- if (is.null(names(dotargs)))
        rep(FALSE, length(dotargs))
    else (names(dotargs) != "")
    if (any(named))
        warning("the following arguments to 'anova.ndlClassify' are invalid and dropped: ",
            paste(deparse(dotargs[named]), collapse = ", "))
    dotargs <- dotargs[!named]
    is.ndlClassify <- unlist(lapply(dotargs, function(x) inherits(x,
        "ndlClassify")))
    dotargs <- dotargs[is.ndlClassify]
    c(list(object), dotargs);
    if (length(dotargs) > 0)
        return(anova.ndlClassifylist(c(list(object), dotargs), statistic=statistic, test=test))
   else 
       { formula <- object$formula;
         data <- object$data;
         frequency <- object$frequency
         cues <- attr(terms(formula),"term.labels") # N.B. if interactions are determined in the formula, such terms are also extracted
         n.cues = length(cues)
         outcome <- as.character(formula)[2]

         formula.model <- vector("character",n.cues+1); statistic.model <- statistic.difference <- significance <- vector("numeric",n.cues+1); df <- df.difference <- vector("integer",n.cues+1);
         formula.model[1] <- "NULL"; statistic.difference <- df.difference <- significance <- NA;

         for(i.cue in 1:n.cues)
            { formula.model[i.cue+1] <- paste(c(outcome, paste(cues[1:i.cue], collapse=" + ")), collapse=" ~ ")
              summary.model <- ndlStatistics(ndlClassify(as.formula(formula.model[i.cue+1]), data, frequency=frequency))
              statistic.model[i.cue+1] <- summary.model[["deviance.model"]] # Another statistic could be substituted here
              df[i.cue+1] <- summary.model[["df.model"]];
              if(i.cue==1)
                { statistic.model[1] <- summary.model[["deviance.null"]]; df[1] <- summary.model[["df.null"]]; }
               statistic.difference[i.cue+1] = statistic.model[i.cue]-statistic.model[i.cue+1]; df.difference[i.cue+1] = df[i.cue]-df[i.cue+1];
               significance[i.cue+1] = 1-pchisq(statistic.difference[i.cue+1], df.difference[i.cue+1])
            }
         effect.cues <- data.frame(df, statistic.model, df.difference, statistic.difference, significance)
         if(statistic=="deviance")
           { colnames(effect.cues) <- c("Resid. Df","Resid. Dev.","Df","Deviance","P(>|Chi|)");
             rownames(effect.cues) <- formula.model;
             title <- paste("Analysis of Deviance Table", "\n\nModel: ndlClassify", 
                      "\n\nResponse: ", outcome, "\n\nTerms added sequentially (first to last)\n\n", sep = "")
           }

         structure(effect.cues, heading = title, class = c("anova", "data.frame"))

       }
}

anova.ndlClassifylist <- function(object, ..., statistic = "deviance", test = "Chisq") 
{
    if(statistic!="deviance")
      stop(paste(c("ANOVA for Naive Discriminative Learning not implemented for statistic: ",statistic),collapse=""))

    responses <- as.character(lapply(object, function(x) {
        deparse(formula(x)[[2]])
    }))
    sameresp <- responses == responses[1]
    if (!all(sameresp)) {
        object <- object[sameresp]
        warning("models with response ", deparse(responses[!sameresp]),
            " removed because response differs from model 1")
    }
    ns <- sapply(object, function(x) NROW(x$data))
    if (any(ns != ns[1]))
        stop("models were not all fitted to the same size of dataset")
    nmodels <- length(object);
    if (nmodels == 1)
        return(anova.ndlClassify(object[[1]], statistic=statistic, test=test))
    df <- unlist(as.numeric(lapply(object, function(x) ndlStatistics(x)$df.model)))
    deviance.model <- as.numeric(lapply(object, function(x) ndlStatistics(x)$deviance.model)) # another statistic could be substituted here
    formula.model <- unlist(lapply(object, function(x) { f <- as.character(x$formula); return(paste(c(f[2]," ",f[1]," ",f[3]),collapse=""))}))
    df.difference <- c(NA,df[1:(nmodels-1)]-df[2:nmodels])
    deviance.difference <- c(NA,deviance.model[1:(nmodels-1)]-deviance.model[2:nmodels])
    significance <- c(NA, sapply(1:(nmodels-1), function(i) 1-pchisq(deviance.difference[i+1], df.difference[i+1])))

    effects.models <- data.frame(df, deviance.model, df.difference, deviance.difference, significance);
    colnames(effects.models) <- c("Resid. Df","Resid. Dev.","Df","Deviance","P(>|Chi|)");
    rownames(effects.models) <- formula.model;

    title <- "Analysis of Deviance Table\n"
    topnote <- paste("Model ", format(1:nmodels), ": ", formula.model, sep = "", collapse = "\n")

    structure(effects.models, heading=c(title, topnote), class =c("anova", "data.frame"))

}

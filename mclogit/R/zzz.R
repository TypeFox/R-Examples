.onLoad <- function(lib,pkg){
  if(requireNamespace("memisc",quietly = TRUE)){
    memisc::setSummaryTemplate(
      mclogit = c(
        "Likelihood-ratio" = "($LR:f1#)",
        #p             = "($p:#)",
        "Log-likelihood" = "($logLik:f1#)",
        Deviance      = "($deviance:f1#)",
        AIC           = "($AIC:f1#)",
        BIC           = "($BIC:f1#)",
        N             = "($N:d)"
      ),
      mclogitRandeff = c(
        #"Likelihood-ratio" = "($LR:f1#)",
        #p             = "($p:#)",
        #"Log-likelihood" = "($logLik:f1#)",
        Deviance      = "($deviance:f1#)",
        #AIC           = "($AIC:f1#)",
        #BIC           = "($BIC:f1#)",
        N             = "($N:d)"
      )
    )  
  }
}


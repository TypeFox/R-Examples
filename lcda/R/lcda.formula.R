###########################################################################
####                         lcda.formula                              ####
####  ===============================================================  ####
####  - estimation of one Latent-Class-Model per class                 ####
####  - determination of model selection criteria                      ####
###########################################################################



lcda.formula <- function(
                         formula,             # formula
                         data,                # data
                         ...
                        )
{

  modelf <- model.frame(formula, data=data)
  grouping <- modelf[,1]
  x <- modelf[,-1]

  lcda(x, grouping, ...)
}



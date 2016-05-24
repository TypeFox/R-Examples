###########################################################################
####                         cclcda.formula                            ####
####  ===============================================================  ####
####  - estimation of one Latent-Class-Model                           ####
####  - determination of model selection criteria                      ####
###########################################################################



cclcda.formula <- function(
                         formula,             # formula
                         data,                # data
                         ...
                        )
{

  modelf <- model.frame(formula, data=data)
  grouping <- modelf[,1]
  x <- modelf[,-1]

  cclcda(x, grouping, ...)
}



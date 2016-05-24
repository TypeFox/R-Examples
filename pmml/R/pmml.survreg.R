# PMML: Predictive Model Markup Language
#
# Copyright (c) 2009-2013, some parts by Togaware Pty Ltd and other by Zementis, Inc. 
#
# This file is part of the PMML package for R.
#
# The PMML package is free software: you can redistribute it and/or 
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 2 of 
# the License, or (at your option) any later version.
#
# The PMML package is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. Please see the
# GNU General Public License for details (http://www.gnu.org/licenses/).
######################################################################################

########################################################################
# Linear Model PMML exporter
#
# Implemented: 091015 Graham Williams based on pmml.lm

.pmml.survreg <- function(model,
                         model.name="Survival_Regression_Model",
                         app.name="Rattle/PMML",
                         description="Survival Regression Model",
                         copyright=NULL,
                         transforms=NULL,
			 unknownValue=NULL,
                         ...)
{
  if (! inherits(model, "survreg")) stop("Not a legitimate survreg object")
  
  # Collect the required information.

  # For a regression, all variables will have been used except those
  # with a NA coefficient indicating singularities. We mark
  # singularities as inactive shortly.

  terms <- attributes(model$terms)
  
  field <- NULL
  field$name <- names(terms$dataClasses)
  field$class <- terms$dataClasses
  
  # 091020 The target field is actually "risk" of the event, and it is
  # not an actual supplied variable. Notice that in pmml.lm we get the
  # target from field$name[1], which in our case here is "Surv(time,
  # status)". We could get the risk score as status - since risk is
  # the probability of the event occuring in comparison to the
  # population. But let's introduce a new variable, "risk" as numeric,
  # as the predicted variable.
  
  # field$name[1] <- sub(')', '', sub('Surv\\([^,]*, *', '', field$name[1]))
  field$name[1] <- 'risk'
  field$class[1] <- "numeric"
  names(field$class)[1] <- field$name[1]

  # 090103 Support transforms if available.
  
  orig.names <- field$name
  orig.class <- field$class

  number.of.fields <- length(field$name)

  target <- field$name[1]

  # 091020 TODO Do we want to do this for coxph? 090501 Identify those
  # who are singularities. For numerics, this is easy since the names
  # are just the variable names. For categorics this gets tricky
  # because the names include the levels. So we need to keep in
  # inactive the actual variable name, if all coefficients for that
  # variable are NAs.
  
  inactive <- names(which(is.na(coef(model))))
  active <- names(which(!is.na(coef(model))))

  # These are the actual variable names.
  
  tmp <- sapply(sapply(field$name, grep, inactive), length)
  inactive.vars <- names(tmp[tmp>0])
  tmp <- sapply(sapply(field$name, grep, active), length)
  active.vars <- names(tmp[tmp>0])
  
  # Now remove any which have any non-NA levels. This final list is
  # passed on as the definitive list of nonactive variables

  inactive <- setdiff(inactive.vars, active.vars)

  # 091020 Do we need to modify this for coxph?
  
  for (i in 1:number.of.fields)
    if (field$class[[field$name[i]]] == "factor")
      # 081004 gjw Test if the data is available in the model, as it
      # would be for a glm (but not an lm), since if the target
      # variable is categoric then the levels are not recorded in
      # xlevels for the target variable, so we will need to get the
      # levels from the data itself.
      if (is.null(model$data))
        field$levels[[field$name[i]]] <- model$xlevels[[field$name[i]]]
      else
        field$levels[[field$name[i]]] <- levels(model$data[[field$name[i]]])

  # PMML

  pmml <- .pmmlRootNode("4.2")

  # PMML -> Header

  pmml <- append.XMLNode(pmml, .pmmlHeader(description, copyright, app.name))

  # PMML -> DataDictionary

  pmml <- append.XMLNode(pmml, .pmmlDataDictionary(field,transformed=transforms))

  # PMML -> RegressionModel

  model.type <- "survreg"
  
  the.model <- xmlNode("RegressionModel",
                       attrs=c(modelName=model.name,
                         functionName="regression",
                         algorithmName="survreg",
                         targetFieldName=target))

  # PMML -> RegressionModel -> MiningSchema

  the.model <- append.XMLNode(the.model, .pmmlMiningSchema(field, target, transformed=transforms, unknownValue=unknownValue))

  # PMML -> TreeModel -> LocalTransforms

# Wen: sep 2013 : obsolete? We need to create surv-reg model to test this piece of code ...
 # if (.supportTransformExport(transforms))
 #   the.model <- append.XMLNode(the.model, .gen.transforms(transforms))

  # test of Zementis xform functions
  if(!is.null(transforms))
  {
    the.model <- append.XMLNode(the.model, .pmmlLocalTransformations(field, transforms, NULL))
  }
  
  # PMML -> RegressionModel -> RegressionTable

  coeff <- coefficients(model)
  coeffnames <- names(coeff)
  means <- model$means

  intercept <- coeff[['(Intercept)']]
  
  regTable <- xmlNode("RegressionTable",
                      attrs=c(intercept=intercept))
  
  # 080620 gjw The PMML spec (at least the Zementis validator)
  # requires NumericPredictors first and then
  # CategoricalPredictors. Simplest approach is to loop twice!!
  # Hopefully, this is not a significant computational expense.

  # For the coxph regression, we need to record the means. This is
  # then used to subtract from supplied value before multiplying by
  # the coefficient in the regression formula.

  for (i in 1:length(orig.names))
  {
    name <- orig.names[[i]]
    if (name == target) next
    klass <- orig.class[[name]]
    if (klass == 'numeric')
    {
      predictorNode <- xmlNode("NumericPredictor",
                               attrs=c(name=name,
                                 exponent="1",
                                 coefficient=as.numeric(coeff[which(coeffnames==name)]),
                                 mean=as.numeric(means[which(coeffnames==name)])))
      regTable <- append.XMLNode(regTable, predictorNode)
    }
  }

  for (i in 1:length(orig.names))
  {
    name <- orig.names[[i]]
    if (name == target) next
    klass <- orig.class[[name]]
    if (klass == 'factor')
    {
      levs <- model$xlevels[[name]]
      # 081019 gjw Add in a zero coefficient for the base level. In
      # this way, we communicate through the PMML which level is the
      # base. Can be useful in then comparing with the full list of
      # levels available for this variable and determining levels that
      # are just missing from the training. Note that xlevels does not
      # include any levels that were not modelled (i.e., missing
      # levels from the training data). We do this by iterating over
      # all the modelled levels (levs, i.e., all values in xlevels)
      # instead of all but the first level (levs[-1], i.e., the base
      # level). When we have the first level, we simply note the
      # coefficient as 0. 090306 This was updated to remove the
      # assumption that the first level has a 0 coefficient. This is
      # not the case in simple lm models (e.g., exampe(lm);
      # pmml(lm.D90)).
      for (l in levs)
      {
        tmp <- paste(name, l, sep='')
        # 090306 Change this test from one that assumes a 0
        # coefficient for the first level, to one that has a 0
        # coefficient for any missing level.
        coefficient <- ifelse(!length(which(coeffnames == tmp)), 0.00,
                              as.numeric(coeff[which(coeffnames == tmp)]))
        predictorNode <- xmlNode("CategoricalPredictor",
                                 attrs=c(name=name,
                                   value=l, coefficient=coefficient))
        regTable <- append.XMLNode(regTable, predictorNode)
      }
    }
  }
  
  the.model <- append.XMLNode(the.model, regTable)
  
  # Add to the top level structure.
  
  pmml <- append.XMLNode(pmml, the.model)
  
  return(pmml)
}

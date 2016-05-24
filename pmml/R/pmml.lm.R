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
# Implemented: 070528 rguha@indiana.edu based on Graham's template for
# handling rpart trees.
#
# Modified: 080201 by Zementis, Inc. (info@zementis.com) to add the
# capability to export binary logistic regression models using glm.
#
# Modified: 090103 by Graham Williams to add transforms framework.

pmml.lm <- function(model,
                    model.name="Linear_Regression_Model",
                    app.name="Rattle/PMML",
                    description="Linear Regression Model",
                    copyright=NULL,
                    transforms=NULL,
		    unknownValue=NULL,
                    dataset=NULL,
                    weights=NULL,
                    ...)
{
  if (! inherits(model, "lm")) stop("Not a legitimate lm object")

  # Collect the required information.

  # For a regression, all variables will have been used except those
  # with a NA coefficient indicating singularities. We mark
  # singularities as inactive shortly.

  terms <- attributes(model$terms)
  
  field <- NULL

  field$name <- names(terms$dataClasses)
  # 101009 Check for a "(weights)" data class and remove it. This
  # arises in the glm models when a Weight variable is used. Not sure
  # why the glm model records this here.
  weights <- which(field$name == "(weights)")
  if (length(weights)) field$name <- field$name[-weights]
  orig.names <- field$name

  field$class <- terms$dataClasses
  if (length(weights)) field$class <- field$class[-weights]
  orig.class <- field$class

  # 090103 Support transforms if available.
# Wen ... commented to see if it break anything ... august 2013
#  if (.supportTransformExport(transforms))
#  {
#    field <- .unifyTransforms(field, transforms)
 #   transforms <- .activateDependTransforms(transforms)
#  }
  number.of.fields <- length(field$name)

  target <- field$name[1]

  # 090501 Identify those who are singularities. For numerics, this is
  # easy since the names are just the variable names. For categorics
  # this gets tricky because the names include the levels. So we need
  # to keep in inactive the actual variable name, if all coefficients
  # for that variable are NAs.
  
  inactive <- names(which(is.na(coef(model))))
  active <- names(which(!is.na(coef(model))))

  # 110113 For the following grep, with certain Japanes characters we
  # see the string including a "[" in some encoding and causes the
  # grep to fail. We can't do the usual Encoding<- "UTF-8" trick since
  # the characters already look like UTF-8. But using enc2utf8 works -
  # does it hurt doing it always, or just when we have Japanese? Needs
  # testing.

  field$name <- enc2utf8(field$name)
  
  # These are the actual variable names. 110113 Should be using grepl
  # rather than grep and then tmp>0!
  
  tmp <- sapply(sapply(field$name, grep, inactive), length)
  inactive.vars <- names(tmp[tmp>0])
  tmp <- sapply(sapply(field$name, grep, active), length)
  active.vars <- names(tmp[tmp>0])
  
  # Now remove any which have any non-NA levels. This final list is
  # passed on as the definitive list of nonactive variables

  inactive <- setdiff(inactive.vars, active.vars)

  for (i in 1:number.of.fields)
  {
    # We don't need to bother with ylevels since lm doesn't do
    # factor predictions.
      
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
  }

  # PMML

  pmml <- .pmmlRootNode("4.2")

  # PMML -> Header

  pmml <- append.XMLNode(pmml, .pmmlHeader(description, copyright, app.name))

  # PMML -> DataDictionary

  pmml <- append.XMLNode(pmml, .pmmlDataDictionary(field, weights=weights,transformed=transforms))

  # PMML -> RegressionModel

  # Added by Zementis so that code can also export binary logistic
  # regression glm models built with binomial(logit). 090303 gjw This
  # looks dangerous, assuming the third argument is the model
  # type. For now, go with it, but set a default model type in case
  # the call has less than two arguments. A general lm model has data
  # as the third part of the call, thus we need to accept that as a
  # genuine model and not an unknown model type! For now, default
  # to generating lm PMML.

  if (model$call[[1]] == "lm")
    model.type <- "lm"
  else if (model$call[[1]] == "glm")
    model.type <- model$family$family
  else
    model.type <- "unknown"
  
  if (model.type == "binomial")
  {
    the.model <- xmlNode("RegressionModel",
                         attrs=c(modelName=model.name,
                           # 100915 Wen-Ching Lin of Zementis noted
                           # this was regression but should be
                           # classification.
                           functionName="classification",
                           algorithmName="glm",
                           normalizationMethod="softmax")) 
  }
  else if (model.type == "poisson")
  {
    the.model <- xmlNode("RegressionModel",
                         attrs=c(modelName=model.name,
                           functionName="regression",
                             algorithmName="glm",
                           normalizationMethod="exp")) 
  }
  else if (model.type == "gaussian")
  {
    the.model <- xmlNode("RegressionModel",
                         attrs=c(modelName=model.name,
                           functionName="regression",
                           algorithmName="glm")) 
  }		
  else if (model.type == "lm")
  {
    # The original code for linear regression models.
    the.model <- xmlNode("RegressionModel",
                         attrs=c(modelName=model.name,
                           functionName="regression",
                           algorithmName="least squares"))
  }
  else 
    stop("pmml.lm: Not a supported family object: ", model.type)

  # PMML -> RegressionModel -> MiningSchema

  the.model <- append.XMLNode(the.model, .pmmlMiningSchema(field, target, transformed=transforms, unknownValue=unknownValue))

  # PMML -> TreeModel -> Output

  the.model <- append.XMLNode(the.model, .pmmlOutput(field, target))

  #----------------------------------------------------
  # PMML -> TreeModel -> LocalTransforms

  # test of Zementis xform functions
  if(!is.null(transforms))
  {
    the.model <- append.XMLNode(the.model, .pmmlLocalTransformations(field, transforms, NULL))
  }

  # PMML -> RegressionModel -> RegressionTable

  coeff <- coefficients(model)

  # 100519 From Wen of Zementis For singularities the coefficient is
  # NA. From DMG the specification says:
  #
  # <xs:attribute name="coefficient" type="REAL-NUMBER" use="required" />
  #
  # So replace NAs with 0. The effect should be the same.
  
  coeff[is.na(coeff)] <- 0 

  coeffnames <- names(coeff)

  # 090306 Handle the case where the intercept is not in the
  # coefficients, and hence is 0?
  
  if (coeffnames[[1]] == "(Intercept)")
    intercept <- as.numeric(coeff[[1]])
  else
    intercept <- 0
  
  # Added by Graham Williams so that code identifies a targetCategory
  # for binary logistic regression glm models built with
  # binomial(logit) and 091009 adds an extra RegressionTable
  # (regTable2) to explicitly refer to the second target category, as
  # recommended in PMML 4.0 specs. and having an intercept of 0.0.

  regTable2 <- NULL
  
  if (model.type == "binomial")
  {
    # 090117 Identify the two possible values for the target variable,
    # and select the second as the target. Extend the PMML specs so I
    # can add the other value as well, since I need that when
    # generating C code to return a class rather than a probability.

    values <- sort(unique(model$data[[target]]))
    alternative.value <- as.character(values[1])
    target.value <- as.character(values[2])
    regTable <- xmlNode("RegressionTable",
                        attrs=c(targetCategory=target.value,
                          intercept=intercept))
    regTable2 <- xmlNode("RegressionTable",
                         attrs=c(targetCategory=alternative.value,
                           intercept="0.0"))
  }
  else
  {
    regTable <- xmlNode("RegressionTable",
                        attrs=c(intercept=intercept))
  }
  
  # 080620 gjw The PMML spec (at least the Zementis validator)
  # requires NumericPredictors first and then
  # CategoricalPredictors. Simplest approach is to loop twice!!
  # Hopefully, this is not a significant computational expense.

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
                                 coefficient=as.numeric(coeff[which(coeffnames==name)])))
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
                                   value=.markupSpecials(l), coefficient=coefficient))
        regTable <- append.XMLNode(regTable, predictorNode)
      }
    }
  }
  
  the.model <- append.XMLNode(the.model, regTable)
  if (! is.null(regTable2)) the.model <- append.XMLNode(the.model, regTable2)
  
  # Add to the top level structure.
  
  pmml <- append.XMLNode(pmml, the.model)
  
  return(pmml)
}

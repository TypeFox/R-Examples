

# Set global chunk options for knitr documentation
opts_chunk$set(fig.width=6.5, fig.height=4.5, size="small", tidy=FALSE)

# Generate R code from the sweave/knitr documentation using the purl command as shown below
purl(input="inst/doc/Causata-vignette.rnw", output="inst/doc/Causata-vignette.R", documentation=0)


library(Causata)
library(pROC)

# Set a random seed so random numbers are repeated
set.seed(87352)

# Load example data from the Causata package
data(df.causata)



# Set the dependent variable
dvname <- "has.responded.mobile.logoff_next.hour_466"
dv <- rep(0, nrow(df.causata))
dv[df.causata[[dvname]] == "true"] <- 1

# Remove unwanted variables from the data frame
df.causata[[dvname]] <- NULL



# Loop to remove variables with no variation
for (col in colnames(df.causata)){
  if (length(unique(df.causata[[col]])) == 1){
    # Single valued variable, remove it
    df.causata[[col]] <- NULL
  }
}



# Create a CausataData object
causataData <- CausataData(df.causata, dependent.variable=dv)



# Replace missing values
causataData <- CleanNaFromContinuous(causataData)
causataData <- CleanNaFromFactor(causataData)



# Merge levels in factors with # levels exceeding a threshold
causataData <- MergeLevels(causataData, max.levels=15)



# Replace outliers in authentication count variables
for (col in grep("^online.average.authentications.per.month", 
                 names(causataData$df), value=TRUE)){
  causataData <- ReplaceOutliers(causataData, col, upperLimit=200)
}



# Extract a set of independent variable names
indep.vars <- colnames(causataData$df)
indep.vars <- indep.vars[!(indep.vars == "dependent.variable") ]
# Build a formula string
formula.string <- paste("dependent.variable", "~",  
                        paste(indep.vars, collapse=" + "))
formula.object <- formula(formula.string)
# Build a model matrix
x.matrix <- model.matrix(formula.object, data=causataData$df)



# Split into training and testing data
idx.split <- sample(1:nrow(x.matrix), round(0.75 * nrow(x.matrix)))
x.matrix.train <- x.matrix[ idx.split, ]
x.matrix.test  <- x.matrix[-idx.split, ]
y.train <- causataData$df$dependent.variable[ idx.split ]
y.test  <- causataData$df$dependent.variable[-idx.split ]



# Build model, select alpha value near 1 
# since we expect most coefficients to be zero
cv.glmnet.obj <- cv.glmnet(x.matrix.train, y.train, alpha=0.8, 
                           family=c("binomial"))
plot(cv.glmnet.obj)



# Use the model to predict responses for training / testing data
predicted.train <- predict(cv.glmnet.obj, newx=x.matrix.train, 
                           type="response", s="lambda.min")
predicted.test  <- predict(cv.glmnet.obj, newx=x.matrix.test,  
                           type="response", s="lambda.min")



# Compute area under the ROC curve using the pROC package
roc.train <- roc(y.train, predicted.train)
roc.test  <- roc(y.test,  predicted.test )
cat("Training / testing area under ROC curve: ", 
    roc.train$auc, ", ", roc.test$auc, "\n")
plot(roc.train, 
     main=sprintf("ROC plot for glmnet model training data.  AUC=%6.4f", 
     roc.train$auc))



# Prepare to import the model into Causata
model.def <- ModelDefinition.Glmnet(
  cv.glmnet.obj, 
  causataData, 
  formula.object, 
  cv.glmnet.obj$lambda.min )
variable.def <- VariableDefinition(
  name = "score-test-model", 
  display.name = "Score: Test Model",
  description = "A logistic regression model trained with sampledata.",
  author = "Test User" )

# Generate a string of PMML representing the model and 
# preprocessing transformations.
# This step is not required to upload a model, 
# it's shown for illustration purposes only
pmml.string <- ToPmml(model.def, variable.def)

# Upload model to Causata.
# The parameters below are for illustration only.
causata.config <- CausataConfig(
  config.server.host = "123.456.789.10",
  config.server.port = 8002,
  username = "testuser",
  password = "1234abcd")

# The command below is commented since it requires a 
# Causata server connection.
#result <- UploadModel(causata.config, model.def, variable.def)



# extract nonzero coefficients
coefs.all <- as.matrix(coef(cv.glmnet.obj, s="lambda.min"))
idx <- as.vector(abs(coefs.all) > 0)
coefs.nonzero <- as.matrix(coefs.all[idx])
rownames(coefs.nonzero) <- rownames(coefs.all)[idx]



print(coefs.nonzero)



## ----echo = FALSE--------------------------------------------------------
library(mldr)

## ----install, eval = FALSE-----------------------------------------------
#  install.packages("mldr")

## ----mld_load, eval = FALSE----------------------------------------------
#  corel5k <- mldr("corel5k")

## ----mld_load_k, eval = FALSE--------------------------------------------
#  corel5k <- mldr("corel5k", label_amount = 374)

## ----mld_load_meka, eval = FALSE-----------------------------------------
#  imdb <- mldr("imdb", use_xml = FALSE)

## ----custom_mld, eval = FALSE--------------------------------------------
#  df <- data.frame(matrix(rnorm(1000), ncol = 10))
#  df$Label1 <- c(sample(c(0,1), 100, replace = TRUE))
#  df$Label2 <- c(sample(c(0,1), 100, replace = TRUE))
#  mymldr <- mldr_from_dataframe(df, labelIndices = c(11, 12), name = "testMLDR")

## ----summary-------------------------------------------------------------
summary(birds)

## ----measures------------------------------------------------------------
emotions$measures$num.attributes
genbase$measures$scumble

## ----labels--------------------------------------------------------------
birds$labels

## ----plot, fig.width=4, fig.height=4-------------------------------------
plot(emotions, type = "LH")

## ----plot_layout, fig.width=8, fig.height=11-----------------------------
layout(matrix(c(1,1,6,1,1,3,5,5,3,2,4,7), 4, 3, byrow = TRUE))

plot(emotions, type = "LC")
plot(emotions, type = "LH")
plot(emotions, type = "LB")
plot(emotions, type = "CH")
plot(emotions, type = "LSB")
plot(emotions, type = "AT")
plot(emotions, type = "LSH")

## ----plot_lc_custom, fig.width=4, fig.height=4---------------------------
plot(genbase, labelIndices = genbase$labels$index[1:11])

## ----transforms----------------------------------------------------------
emotionsbr <- mldr_transform(emotions, type = "BR")

emotionslp <- mldr_transform(emotions, type = "LP", emotions$labels$index[1:4])

## ----transf_eval, eval = FALSE-------------------------------------------
#  emo_lp <- mldr_transform(emotions, "LP")
#  library(RWeka)
#  classifier <- IBk(classLabel ~ ., data = emo_lp, control = Weka_control(K = 10))
#  evaluate_Weka_classifier(classifier, numFolds = 5)

## ----filter--------------------------------------------------------------
emotions$measures$num.instances
emotions[emotions$dataset$.SCUMBLE > 0.01]$measures$num.instances

## ----decoupling, eval = FALSE--------------------------------------------
#  mldbase <- mld[.SCUMBLE <= mld$measures$scumble]
#  
#  # Samples with coocurrence of highly imbalanced labels
#  mldhigh <- mld[.SCUMBLE > mld$measures$scumble]
#  majIndexes <- mld$labels[mld$labels$IRLbl < mld$measures$meanIR,"index"]
#  
#  # Deactivate majority labels
#  mldhigh$dataset[, majIndexes] <- 0
#  joined <- mldbase + mldhigh # Join the instances without changes with the filtered ones

## ----equals--------------------------------------------------------------
emotions[1:10] == emotions[20:30]
emotions == birds

## ----evaluation----------------------------------------------------------
# Get the true labels in emotions
predictions <- as.matrix(emotions$dataset[,emotions$labels$index])
# and introduce some noise
predictions[sample(1:593, 100),sample(1:6, 100, replace = TRUE)] <- sample(0:1, 100,
  replace = TRUE)
# then evaluate the predictive performance
res <- mldr_evaluate(emotions, predictions)
str(res)

## ----eval_plot, eval = FALSE---------------------------------------------
#  plot(res$ROC, main = "ROC curve for emotions") # Plot ROC curve

## ----gui, eval = FALSE---------------------------------------------------
#  mldrGUI()


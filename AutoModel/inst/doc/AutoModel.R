## ---- fig.width=6, fig.height=5------------------------------------------
library(AutoModel)
run_model("mpg", c("disp", "hp"), c("cyl", "wt"), c("drat", "qsec"), dataset=mtcars)

## ------------------------------------------------------------------------
formulas <- create_formula_objects("mpg", c("disp", "hp"), c("cyl", "wt"), c("drat", "qsec"))
models <- create_model_objects(formulas, dataset = mtcars)
models

## ---- fig.width=6, fig.height=5------------------------------------------
run_model("mpg", c("disp", "hp"), c("cyl", "wt"), c("drat", "qsec"), dataset=mtcars, assumptions.check = F)

## ---- fig.width=6, fig.height=5------------------------------------------
model_object <- run_model("mpg", c("disp", "hp"), c("cyl", "wt"), c("drat", "qsec"), dataset=mtcars)
model_object$Summary
model_object$Coefficients
model_object$Checks

## ------------------------------------------------------------------------
model_object$Summary$DeltaR2
model_object$Coefficients$estimate
model_object$Checks$Correlation.Matrix

## ----xtable, results='asis'----------------------------------------------
library(xtable)
sum_table <- xtable(model_object$SummaryDF)
coef_table <- xtable(model_object$CoefficientsDF)
print(sum_table, type="html")
print(coef_table, type="html")

## ------------------------------------------------------------------------
model <- run_model("am", c("disp", "hp"), c("cyl", "wt"), dataset = mtcars, type="binomial")
model$Summary
model$Coefficients
model$Class_Table

## ---- results='asis'-----------------------------------------------------
sum_table <- xtable(model$Summary)
coef_table <- xtable(model$Coefficients)
class_table <- xtable(model$Class_Table)
print(sum_table, type="html")
print(coef_table, type="html", digits = 4)
print(class_table, type="html")


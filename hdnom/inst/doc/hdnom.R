## ------------------------------------------------------------------------
library("hdnom")

## ------------------------------------------------------------------------
data("smart")
x = as.matrix(smart[, -c(1, 2)])
time = smart$TEVENT
event = smart$EVENT

library("survival")
y = Surv(time, event)

## ---- eval = FALSE-------------------------------------------------------
#  # Enable parallel parameter tuning
#  suppressMessages(library("doParallel"))
#  registerDoParallel(detectCores())
#  
#  aenetfit = hdcox.aenet(x, y, nfolds = 10, rule = "lambda.1se",
#                         seed = c(5, 7), parallel = TRUE)
#  names(aenetfit)

## ---- echo = FALSE-------------------------------------------------------
aenetfit = readRDS("aenetfit.rds")

## ------------------------------------------------------------------------
fit    = aenetfit$aenet_model
alpha  = aenetfit$aenet_best_alpha
lambda = aenetfit$aenet_best_lambda
adapen = aenetfit$pen_factor

## ---- fig.width = 8, fig.height = 8, out.width = 600, out.height = 600----
suppressMessages(library("rms"))
x.df = as.data.frame(x)
dd = datadist(x.df)
options(datadist = "dd")

nom = hdnom.nomogram(fit, model.type = "aenet", x, time, event, x.df,
                     lambda = lambda, pred.at = 365 * 2,
                     funlabel = "2-Year Overall Survival Probability")
plot(nom)

## ------------------------------------------------------------------------
val.int = hdnom.validate(x, time, event, model.type = "aenet",
                         alpha = alpha, lambda = lambda, pen.factor = adapen,
                         method = "bootstrap", boot.times = 10,
                         tauc.type = "UNO", tauc.time = seq(1, 5, 0.5) * 365,
                         seed = 42, trace = FALSE)
val.int
summary(val.int)

## ---- fig.width = 8, fig.height = 8, out.width = 600, out.height = 600----
plot(val.int)

## ---- fig.width = 8, fig.height = 8, out.width = 600, out.height = 600----
x_new = as.matrix(smart[, -c(1, 2)])[1001:2000, ]
time_new = smart$TEVENT[1001:2000]
event_new = smart$EVENT[1001:2000]

# External validation with time-dependent AUC
val.ext =
  hdnom.external.validate(aenetfit, x, time, event,
                          x_new, time_new, event_new,
                          tauc.type = "UNO",
                          tauc.time = seq(0.25, 2, 0.25) * 365)

val.ext
summary(val.ext)
plot(val.ext)

## ------------------------------------------------------------------------
cal.int = hdnom.calibrate(x, time, event, model.type = "aenet",
                          alpha = alpha, lambda = lambda, pen.factor = adapen,
                          method = "bootstrap", boot.times = 10,
                          pred.at = 365 * 5, ngroup = 3,
                          seed = 42, trace = FALSE)
cal.int
summary(cal.int)

## ---- fig.width = 8, fig.height = 8, out.width = 600, out.height = 600----
plot(cal.int, xlim = c(0.5, 1), ylim = c(0.5, 1))

## ---- fig.width = 8, fig.height = 8, out.width = 600, out.height = 600----
cal.ext =
  hdnom.external.calibrate(aenetfit, x, time, event,
                           x_new, time_new, event_new,
                           pred.at = 365 * 5, ngroup = 3)

cal.ext
summary(cal.ext)
plot(cal.ext, xlim = c(0.5, 1), ylim = c(0.5, 1))

## ---- fig.width = 8, fig.height = 8, out.width = 600, out.height = 600----
hdnom.kmplot(cal.int, group.name = c('High risk', 'Medium risk', 'Low risk'),
             time.at = 1:6 * 365)

hdnom.kmplot(cal.ext, group.name = c('High risk', 'Medium risk', 'Low risk'),
             time.at = 1:6 * 365)

## ------------------------------------------------------------------------
cal.int.logrank = hdnom.logrank(cal.int)
cal.int.logrank
cal.int.logrank$pval
cal.ext.logrank = hdnom.logrank(cal.ext)
cal.ext.logrank
cal.ext.logrank$pval

## ---- fig.width = 8, fig.height = 6.4, out.width = 600, out.height = 480----
cmp.val =
  hdnom.compare.validate(x, time, event,
                         model.type = c("lasso", "alasso"),
                         method = "cv", nfolds = 5, tauc.type = "UNO",
                         tauc.time = seq(0.25, 2, 0.25) * 365,
                         seed = 42, trace = FALSE)

cmp.val
summary(cmp.val)
plot(cmp.val)
plot(cmp.val, interval = TRUE)

## ---- fig.width = 8, fig.height = 6.4, out.width = 600, out.height = 480----
cmp.cal =
  hdnom.compare.calibrate(x, time, event,
                          model.type = c("lasso", "alasso"),
                          method = "cv", nfolds = 5,
                          pred.at = 365 * 9, ngroup = 5,
                          seed = 42, trace = FALSE)

cmp.cal
summary(cmp.cal)
plot(cmp.cal, xlim = c(0.3, 1), ylim = c(0.3, 1))

## ------------------------------------------------------------------------
predict(aenetfit, x, y, newx = x[101:105, ], pred.at = 1:10 * 365)


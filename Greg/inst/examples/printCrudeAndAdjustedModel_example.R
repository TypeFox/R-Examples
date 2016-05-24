# simulated data to use 
set.seed(10)
ds <- data.frame(
  ftime = rexp(200),
  fstatus = sample(0:1,200,replace=TRUE),
	Variable1 = runif(200),
  Variable2 = runif(200),
  Variable3 = runif(200),
  Variable4 = factor(sample(LETTERS[1:4], size=200, replace=TRUE)))

library(rms)
dd <- datadist(ds)
options(datadist="dd")

fit <- cph(Surv(ftime, fstatus) ~ Variable1 + Variable3 + Variable2 +  Variable4,
           data=ds, x=TRUE, y=TRUE)
printCrudeAndAdjustedModel(fit, order = c("Variable[12]", "Variable3"))
printCrudeAndAdjustedModel(fit, 
                           order=c("Variable3", "Variable4"),
                           add_references = TRUE, 
                           desc_column=TRUE)

# Now to a missing example
n <- 500
ds <- data.frame(
  x1 = factor(sample(LETTERS[1:4], size = n, replace = TRUE)),
  x2 = rnorm(n, mean = 3, 2),
  x3 = factor(sample(letters[1:3], size = n, replace = TRUE)))

ds$Missing_var1 <- factor(sample(letters[1:4], size=n, replace=TRUE))
ds$Missing_var2 <- factor(sample(letters[1:4], size=n, replace=TRUE))
ds$y <- rnorm(nrow(ds)) +
  (as.numeric(ds$x1)-1) * 1 +
  (as.numeric(ds$Missing_var1)-1)*1 + 
  (as.numeric(ds$Missing_var2)-1)*.5

# Create a messy missing variable
non_random_missing <- sample(which(ds$Missing_var1 %in% c("b", "d")), 
                             size = 150, replace=FALSE)
# Restrict the non-random number on the x2 variables
non_random_missing <- non_random_missing[non_random_missing %in%
                                           which(ds$x2 > mean(ds$x2)*1.5) &
                                           non_random_missing %in%
                                           which(ds$x2 > mean(ds$y))]
ds$Missing_var1[non_random_missing] <- NA

# Simple missing variable
ds$Missing_var2[sample(1:nrow(ds), size=50)] <- NA

# Setup the rms environment
ddist <- datadist(ds)
options(datadist = "ddist")

impute_formula <- 
  as.formula(paste("~",
                   paste(colnames(ds),
                         collapse="+")))

imp_ds <- aregImpute(impute_formula, data = ds, n.impute = 10)

fmult <- fit.mult.impute(y ~ x1 + x2 + x3 + 
                           Missing_var1 + Missing_var2, 
                         fitter = ols, xtrans = imp_ds, data = ds)

printCrudeAndAdjustedModel(fmult, 
                           impute_args = list(variance.inflation=TRUE,
                                              coef_change=list(type="diff",
                                                               digits=3)))


# Use some labels to prettify the output
# fro the mtcars dataset
data("mtcars")

label(mtcars$mpg) <- "Gas"
units(mtcars$mpg) <- "Miles/(US) gallon"

label(mtcars$wt) <- "Weight"
units(mtcars$wt) <- "10^3 kg" # not sure the unit is correct 

mtcars$am <- factor(mtcars$am, levels=0:1, labels=c("Automatic", "Manual"))
label(mtcars$am) <- "Transmission"

mtcars$gear <- factor(mtcars$gear)
label(mtcars$gear) <- "Gears"

# Make up some data for making it slightly more interesting
mtcars$col <- factor(sample(c("red", "black", "silver"), size=NROW(mtcars), replace=TRUE))
label(mtcars$col) <- "Car color"

require(splines)
fit_mtcar <- lm(mpg ~ wt + gear + col, data=mtcars)
printCrudeAndAdjustedModel(fit_mtcar, 
                           add_references=TRUE,
                           ctable=TRUE, 
                           desc_column = TRUE,
                           digits=1,
                           desc_args = caDescribeOpts(digits = 1,
                                                      colnames = c("Avg.")))

printCrudeAndAdjustedModel(fit_mtcar, 
                           add_references=TRUE,
                           desc_column=TRUE,
                           order=c("Interc", "gear"))

# Alterntive print - just an example, doesn't make sense to skip reference
printCrudeAndAdjustedModel(fit_mtcar, 
                           order=c("col", "gear"), 
                           groups=c("Color", "Gears"),
                           add_references=c("Black", NA),
                           ctable=TRUE)

# Now we can also combine models into one table using rbind()
mpg_model <- printCrudeAndAdjustedModel(lm(mpg ~ wt + gear + col, data=mtcars), 
                                    add_references=TRUE,
                                    ctable=TRUE, 
                                    desc_column = TRUE,
                                    digits=1,
                                    desc_args = caDescribeOpts(digits = 1,
                                                               colnames = c("Avg.")))
wt_model <- printCrudeAndAdjustedModel(lm(wt ~ mpg + gear + col, data=mtcars), 
                                    add_references=TRUE,
                                    ctable=TRUE, 
                                    desc_column = TRUE,
                                    digits=1,
                                    desc_args = caDescribeOpts(digits = 1,
                                                               colnames = c("Avg.")))

library(magrittr)
rbind(Miles = mpg_model, Weight = wt_model) %>% 
  htmlTable(caption="Combining models together with a table spanner element separating each model")

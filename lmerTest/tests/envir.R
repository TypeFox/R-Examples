require(lmerTest)
load(system.file("testdata","cltlike.RData", package="lmerTest"))
form <- "liking ~ pig.type + pig.type:landro + pig.type:lskatol +
 (1|pig) + (1|consumer)"
variables <- c("sequence", "gender")
get_models <- function(variables) {
  lapply(variables, function(v) {
    form2 <- as.formula(paste(form, v, sep=" + "))    
    lmer(form2, data=cltlike)
  })
}
get_models2 <- function(variables, data1) {
  lapply(variables, function(v) {
    form2 <- as.formula(paste(form, v, sep=" + "))
    print(form2)
    print(class(form2))
    environment(data1) <- environment(form2) ## in 2.0-30  was added otherwise fail anova in updateModel
    lmerTest::lmer(form2, data=data1)
  })
}

mlist <- get_models(variables)
ls1 <- lapply(mlist, anova)
stopifnot("Pr(>F)" %in% colnames(ls1[[1]]))



mlist <- get_models2(variables, data=cltlike)
ls2 <- lapply(mlist, anova, type = 1)
stopifnot("Pr(>F)" %in% colnames(ls2[[1]]))

# 
# update_call <- function (object, formula., ...) {
#   call <- object@call
#   
#   # Use update.formula to deal with formulas like . ~ .
#   if (!missing(formula.)) {
#     call$formula <- update.formula(formula(object), formula.)
#   }
#   
#   modify_call(call, dots(...))
# }
# update_model <- function(object, formula., ...) {
#   call <- update_call(object, formula., ...)
#   eval(call, parent.frame())
# }
# update_model(mod, formula = . ~ . + cyl)




########################################################################
# PV-t-test AUX FUNCTIONS
########################################################################

t.test.aux <- function(formula, 
                       design){
  # convert formula args
  tar.char <- as.character(formula)[2]
  grp.char <- as.character(formula)[3]
  tar.form <-as.formula(paste("~", tar.char))
  grp.form <-as.formula(paste("~", grp.char))
  # test input formula
  # formulate test cond (has to be true)
  # left side of formula
  test.cond.d1 <- 
    is.numeric(design$variable[[tar.char]])
  # right side of formula
  test.cond.d2 <- 
    ifelse(length(all.vars(grp.form))==0,
           grp.char %in% c("0","1"),
           dim(table(design$variable[[grp.char]]))==2)
  #   t.cond <- c(dim(table(design$variable[grp.char]))==2, 
  #               grp.char %in% c("0","1")
  #               )
  test.cond <- sum(test.cond.d1, test.cond.d2)
  if(sum(test.cond) != 2){
    stop("Wrong formula structure or group definition!")}
  # case distinction 
  # case 1 (formula type: outcome~0 or outcome~1)
  if(grp.char %in% c("0","1")){
    # calc descriptives (& t-test input)
    des.res.c1 <- data.frame(svymean(tar.form, 
                                     design, 
                                     na.rm=TRUE))
    rownames(des.res.c1) <- ""
    # calc N and sum of w
    N.dummy <- 
      complete.cases(design$variable[all.vars(formula)])
    number.of.cases <- sum(N.dummy)
    sum.of.weights <- sum(design$pweights[N.dummy==TRUE])
    N <- data.frame(number.of.cases, sum.of.weights)
    rownames(N) <- ""
    f.res <- list(DESC=des.res.c1, N=N, TEST=des.res.c1)
    return(f.res)
  }
  # else (formula type: outcome~group)
  else
  {
    # calc descriptives
    des.res.c2 <- as.data.frame(svyby(tar.form,
                                      grp.form,
                                      design,
                                      svymean, 
                                      na.rm=TRUE))
    # calc t-test input
    mod <- summary(svyglm(formula, design))
    mod.ests <- mod$coefficients[2,1:2]
    # calc N and sum of w
    number.of.cases <- table(design$variables[[grp.char]])
    sum.of.weights <- svytable(grp.form, design)
    N <- rbind(number.of.cases, sum.of.weights)
    f.res <- list(DESC=des.res.c2, N=N, TEST=mod.ests)
    return(f.res)
  }
}
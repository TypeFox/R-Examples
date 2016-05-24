######################################################################
# Model Profile
# Note: at the moment there is a bug in the tables package
#       (it fails when only one factor is used)
#       Duncan Mudorch is working on this. See email on Jan 2, 2014 (Lars)
######################################################################


modelProfile <- function(formula, 
                         data, 
                         groups = 10,
                         group_label = c("I", "D"), 
                         digits_numeric = 1,
                         digits_factor = 4,
                         exclude_na = FALSE, 
                         LaTex = FALSE) {
  
  ### Perform preliminary checks
  if (!inherits(formula, "formula"))
    stop("uplift: Method is only for formula objects.")
  
  if (!groups %in% c(5, 10, 20))
    stop("uplift: groups must be either 5, 10 or 20. Aborting...")
  
  if (!tolower(group_label) %in% c("i", "d"))
    stop("uplift: group_label must be either 'I' or 'D'. Aborting...")
  
  ### Extract formula elements
  mf <- match.call(expand.dots = FALSE)
  args <- match(c("formula", "data"), 
                names(mf), 0L)
  mf <- mf[c(1L, args)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.pass
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  
  data_class <- attributes(mt)$dataClasses[-1] # exclude response
  if (!all(unique(data_class) %in% c("numeric", "factor", "ordered")))
      stop("uplift: variable types in formula must be either numeric, integer, factor or ordered. Aborting...")
  num_vars <- which(data_class ==  "numeric")
  fac_vars <- which(data_class %in% c("factor", "ordered"))
  num_var_names <- attributes(mt)$term.labels[num_vars]
  fac_var_names <- attributes(mt)$term.labels[fac_vars]
  resp_var_name <- names(attributes(mt)$dataClasses[1])
  
  attr(mt, "intercept") <- 0 
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm)) 
      names(Y) <- nm
  }
  
  if (!is.numeric(Y)) 
    stop("uplift: the LHS of the model formula must be a numeric vector. Aborting...")
  
  if (group_label == "I")
  rank.Y <- rank(Y) else rank.Y <- rank(-Y)
  
  Group <- cut(rank.Y, breaks = quantile(rank.Y, probs = seq(0, 1, 1/groups)),
                labels = 1:groups, include.lowest = TRUE)
  
  ### Create data frame to evaluate tabular formula
  dframe <- data.frame(mf, Group)
  
  if (exclude_na) dframe_out <-"na.omit(dframe)" else dframe_out <-"dframe"
  
  ### Build tabular formula
  if (length(num_vars) != 0L) t1 <- paste("+", paste(num_var_names, collapse = " + ")) else 
                              t1 <- ""
  if (length(fac_vars) != 0L) {t2 <- paste("+ Format(digits=", digits_factor, ") * ((", 
                                          paste("Factor(", fac_var_names, ")", sep = "",  collapse = " + "),
                                          ") * ",
                                          "(Pctn. = Percent('col')))", sep = "")} else {
                               t2 <- ""}  
  
  tab.form <- paste("tabular((n=1) + Format(digits=", digits_numeric, ")", " * ((", 
                    resp_var_name,  
                    t1, ") * (Avg. = mean))", t2,
                     " ~ Justify(c) * (Group + 1),", 
                    " data = ", dframe_out, ")", sep ="")
  
  if (LaTex) tab.form <- paste("latex(", tab.form, ")", sep ="")
  
  ### Evaluate formula 
  res <- eval(parse(text=tab.form))
  
  return(res)    

}





  
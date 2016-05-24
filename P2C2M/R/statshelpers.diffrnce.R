statshelpers.diffrnce <-
function (post_dist, post_pred_dist) {
  ####################################
  # Function "statshelpers.diffrnce" #
  ####################################
  # Descr:    calculates the difference between post_distirical and post_pred_distulated
  # Deps:     -
  # I/p:      post_dist
  #           post_pred_dist

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> statshelpers.diffrnce", 
        fg="red"), sep="")
  }

  # Note: The difference is "post_distirical-post_pred_distulated". Since "post_distirical" is worse 
  # than "post_pred_distulated" whenever it does not conform to the coalescent model (i.e.
  # has larger values), significant differences will be more positive than 
  # non-significant differences.
  diff = post_dist - post_pred_dist
  # TFL converts from type "list" to type "double"; is important, because
  # is.infinite and is.nan can only work on type "double"
  diff = as.matrix(diff)

  # Removing diff. values that are infinite ("Inf")
  diff = ifelse(is.infinite(diff), NA, diff)
  # Removing diff. values that are "NaN" ("not a number", i.e. 0/0)
  diff = ifelse(is.nan(diff), NA, diff)

  return(diff)
}

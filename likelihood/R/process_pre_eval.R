#####################################################
# process_pre_eval
# Performs function calls on which the main functions
# are dependent.
#
# Author:  Lora Murphy, Cary Institute of Ecosystem Studies
# murphyl@caryinstitute.org
#####################################################
process_pre_eval<-function(pre_evals, eval_results, ...) {

  if (length(pre_evals) > 0) {
    for(i in 1:length(pre_evals)) {
      if (!is.null(pre_evals[[i]]$pre_eval)) {
        eval_results = process_pre_eval(pre_evals[[i]]$pre_eval, eval_results, ...)
      }
      if (is.null(list(...)) == FALSE && any(names(formals(pre_evals[[i]]$fun))=="...")) {
        for (j in 1:length(list(...))) {
          pre_evals[[i]]$call[[names(list(...))[[j]]]]<-list(...)[[j]]
        }
      } 
      eval_results[[pre_evals[[i]]$parname]]<-do.call(pre_evals[[i]]$fun, pre_evals[[i]]$call)
    }
  }
  eval_results
}
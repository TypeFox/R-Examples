`eqmccLoop` <-
function(data, outcome = "", neg.out = FALSE, conditions = "", 
      relation = "suf", n.cut = 1, icp = 1, ica = 1, 
      explain = "1", include = "", row.dom = FALSE, min.dis = TRUE, 
      omit = c(), dir.exp = "", details = FALSE, show.cases = FALSE, 
      inf.test = "", use.tilde = FALSE, use.letters = FALSE, ...) {


    check.object <- verify.mqca(data, outcome, conditions)
    
    conditions <- check.object$conditions
    
    data <- data[, unique(conditions, check.object$outcome)]
    
    eqmcc.list <- solution.list <- vector(mode = "list", length = length(outcome))
    names(eqmcc.list) <- names(solution.list) <- outcome
    
    for (i in seq(length(outcome))) {
        
        conditions <- names(data)[- which(names(data) == check.object$outcome[i])]
        
        # to do: try manipulating the match.call() object or something, instead of this
        eqmcc.list[[i]] <- eqmcc(
            data, outcome = outcome[i], neg.out=neg.out, conditions = conditions,
            relation = relation, n.cut=n.cut, icp = icp, ica = ica,
            explain = explain, include = include, row.dom = row.dom, min.dis = min.dis,
            omit = omit, dir.exp = dir.exp, details = details, show.cases = show.cases,
            inf.test = inf.test, use.tilde = use.tilde, use.letters = use.letters, ... = ...)
    }
    
    return(structure(eqmcc.list, class="mqca"))
    
}

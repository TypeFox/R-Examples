loadModule("CPG", TRUE)

evalqOnLoad({
##
## show-methods for reference objects
setMethod("show", "Rcpp_CTRL", function(object){
    cat("Control parameters used in optimization:\n\n")
    cat(paste("Maximum iterations:\t", object$params$maxiters,"\n"))
    cat(paste("Absolute tolerance:\t", object$params$abstol,"\n"))
    cat(paste("Relative tolerance:\t", object$params$reltol,"\n"))
    cat(paste("Feasible tolerance:\t", object$params$feastol,"\n"))
    cat(paste("Tracing progress:\t", object$params$trace,"\n"))
})
setMethod("show", signature = "Rcpp_DQP", function(object){
    title <- paste("* Definition of Quadratic Program *")
    row <- paste(rep("*", nchar(title)), collapse = "")
    cat("\n")
    cat(row, "\n")
    cat(paste(title, "\n"))
    cat(row, "\n")
    cat("\n")
    cat(paste("Count of variables in objective:", ncol(object$P), "\n"))
    cat(paste("Count of equality constraints:", nrow(object$A), "\n"))
    countcc <- object$cList$K
    cat(paste("Count of cone constraints:", countcc, "\n"))
    cc <- object$cList$cone
    cat("These consist of:\n")
    cat(paste("Constraints w.r.t. the nonnegative orthant:", max(0, sum(cc %in% "NNOC")), "\n"))
    cat(paste("Constraints w.r.t. the second-order cone:", max(0, sum(cc %in% "SOCC")), "\n"))
    cat(paste("Constraints w.r.t. the semidefinite cone:", max(0, sum(cc %in% "PSDC")), "\n"))
    cat("\n")
})
setMethod("show", signature = "Rcpp_DLP", function(object){
    title <- paste("* Definition of Linear Program *")
    row <- paste(rep("*", nchar(title)), collapse = "")
    cat("\n")
    cat(row, "\n")
    cat(paste(title, "\n"))
    cat(row, "\n")
    cat("\n")
    cat(paste("Count of variables in objective:", length(object$q), "\n"))
    cat(paste("Count of equality constraints:", nrow(object$A), "\n"))
    countcc <- object$cList$K
    cat(paste("Count of cone constraints:", countcc, "\n"))
    cc <- object$cList$cone
    cat("These consist of:\n")
    cat(paste("Constraints w.r.t. the nonnegative orthant:", max(0, sum(cc %in% "NNOC")), "\n"))
    cat(paste("Constraints w.r.t. the second-order cone:", max(0, sum(cc %in% "SOCC")), "\n"))
    cat(paste("Constraints w.r.t. the semidefinite cone:", max(0, sum(cc %in% "PSDC")), "\n"))
    cat("\n")
})
setMethod("show", signature = "Rcpp_DNL", function(object){
    title <- paste("* Definition of Linear Program *")
    row <- paste(rep("*", nchar(title)), collapse = "")
    cat("\n")
    cat(row, "\n")
    cat(paste(title, "\n"))
    cat(row, "\n")
    cat("\n")
    cat(paste("Count of variables in objective:", length(object$q), "\n"))
    cat(paste("Count of equality constraints:", nrow(object$A), "\n"))
    countcc <- object$cList$K
    cat(paste("Count of constraints:", countcc, "\n"))
    cc <- object$cList$cone
    cat("These consist of:\n")
    cat(paste("Constraints w.r.t. non-linearities:", object$cList$dims[1, 1], "\n"))
    cat(paste("Constraints w.r.t. the nonnegative orthant:", max(0, sum(cc %in% "NNOC")), "\n"))
    cat(paste("Constraints w.r.t. the second-order cone:", max(0, sum(cc %in% "SOCC")), "\n"))
    cat(paste("Constraints w.r.t. the semidefinite cone:", max(0, sum(cc %in% "PSDC")), "\n"))
    cat("\n")
})
setMethod("show", signature = "Rcpp_DCP", function(object){
    title <- paste("* Definition of Convex Program *")
    row <- paste(rep("*", nchar(title)), collapse = "")
    cat("\n")
    cat(row, "\n")
    cat(paste(title, "\n"))
    cat(row, "\n")
    cat("\n")
    cat(paste("Count of variables in objective:", nrow(object$x0) - 1L, "\n"))
    cat(paste("Count of equality constraints:", nrow(object$A), "\n"))
    if(object$cList$dims[1, 1] > 1){
        countcc <- object$cList$K
    } else {
        countcc <- 0
    }
    cat(paste("Count of constraints:", countcc, "\n"))
    cc <- object$cList$cone
    cat("These consist of:\n")
    cat(paste("Constraints w.r.t. non-linearities:", object$cList$dims[1, 1] - 1L, "\n"))
    cat(paste("Constraints w.r.t. the nonnegative orthant:", max(0, sum(cc %in% "NNOC")), "\n"))
    cat(paste("Constraints w.r.t. the second-order cone:", max(0, sum(cc %in% "SOCC")), "\n"))
    cat(paste("Constraints w.r.t. the semidefinite cone:", max(0, sum(cc %in% "PSDC")), "\n"))
    cat("\n")
})
setMethod("show", signature = "Rcpp_CPS", function(object){
    title <- "* Solution of Convex Program *"
    row <- paste(rep("*", nchar(title)), collapse = "")
    cat("\n")
    cat(row, "\n")
    cat(title, "\n")
    cat(row, "\n")
    cat("\n")
    state <- object$state
    cat(paste("Value of primal objective:", signif(state["pobj"]), "\n"))
    if(!is.na(state["dobj"])){
        cat(paste("Value of dual objective:", signif(state["dobj"]), "\n"))
    }
    if(!is.na(state["dgap"])){
        cat(paste("Value of duality gap:", signif(state["dgap"]), "\n"))
    }
    if(!is.na(state["rdgap"])){
        cat(paste("Value of relative duality gap:", signif(state["rdgap"]), "\n"))
     }
    if(!is.na(state["certp"])){
        cat(paste("Certificate of primal infeasibility:", signif(state["certp"]), "\n"))
    }
    if(!is.na(state["certd"])){
        cat(paste("Certificate of dual infeasibility:", signif(state["certd"]), "\n"))
    }
    if(!is.na(state["pslack"])){
        cat(paste("Value of smallest primal slack:", signif(state["pslack"]), "\n"))
    }
    if(!is.na(state["dslack"])){
        cat(paste("Value of smallest dual slack:", signif(state["dslack"]), "\n"))
    }
    cat(paste("Status of solution:", object$status, "\n"))
    cat(paste("Count of iterations:", object$niter, "\n\n"))
    cat("Solutions are contained in 'PDV'.\n")
    cat("Use 'getx()', 'gety()', 'gets()' and 'getz()', respectively.\n")
})
## cps-methods
setMethod("cps", signature = c("Rcpp_DLP", "Rcpp_CTRL"), function(cpd, ctrl){
    cpd$cps(ctrl)
})
## cps-methods
setMethod("cps", signature = c("Rcpp_DNL", "Rcpp_CTRL"), function(cpd, ctrl){
    cpd$cps(ctrl)
})
## cps-methods
setMethod("cps", signature = c("Rcpp_DQP", "Rcpp_CTRL"), function(cpd, ctrl){
    cpd$cps(ctrl)
})
## cps-methods
setMethod("cps", signature = c("Rcpp_DCP", "Rcpp_CTRL"), function(cpd, ctrl){
    cpd$cps(ctrl)
})
## gets-methods
setMethod("gets", signature = "Rcpp_PDV", function(object){
    object$s
})
setMethod("gets", signature = "Rcpp_CPS", function(object){
    pdv <- object$pdv
    sidx <- object$sidx
    if(nrow(sidx) > 1){
        sidx <- sidx + 1
        ans <- list()
        length(ans) <- nrow(sidx)
        for(i in 1:nrow(sidx)){
            ans[[i]] <- pdv$s[sidx[i, 1]:sidx[i, 2], 1]
        }
    } else {
       ans <- gets(pdv)
    }
    ans
})
## getz-methods
setMethod("getz", signature = "Rcpp_PDV", function(object){
    object$z
})
setMethod("getz", signature = "Rcpp_CPS", function(object){
    pdv <- object$pdv
    sidx <- object$sidx
    if(nrow(sidx) > 1){
        sidx <- sidx + 1
        ans <- list()
        length(ans) <- nrow(sidx)
        for(i in 1:nrow(sidx)){
            ans[[i]] <- pdv$z[sidx[i, 1]:sidx[i, 2], 1]
        }
    } else {
       ans <- getz(pdv)
    }
    ans
})
## getx-methods
setMethod("getx", signature = "Rcpp_PDV", function(object){
    object$x
})
setMethod("getx", signature = "Rcpp_CPS", function(object){
    pdv <- object$pdv
    getx(pdv)
})
## gety-methods
setMethod("gety", signature = "Rcpp_PDV", function(object){
    object$y
})
setMethod("gety", signature = "Rcpp_CPS", function(object){
    pdv <- object$pdv
    gety(pdv)
})
## other get-methods for Rcpp_CPS
setMethod("getstatus", signature = "Rcpp_CPS", function(object){
    object$status
})
setMethod("getstate", signature = "Rcpp_CPS", function(object){
    object$state
})
setMethod("getniter", signature = "Rcpp_CPS", function(object){
    object$niter
})
setMethod("getparams", signature = "Rcpp_CTRL", function(object){
    object$params
})

})

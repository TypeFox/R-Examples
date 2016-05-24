"print.relimplm" <-
function (x, ..., show.coeffs = ifelse(any(c("lmg", "pmvd") %in% x@type) & is.null(x@always), TRUE, FALSE)) 
{
    if (!(is(x, "relimplm"))) 
        stop("x must be the output from function calc.relimp")
    p <- length(slot(x, "namen")) - 1
    grouped <- FALSE
    if (length(slot(x, "groupdocu"))==2) {
        g <- length(slot(x,"groupdocu")[[1]])
        grouped <- TRUE
        p <- sum(list2vec(lapply(slot(x,"groupdocu")[[2]],length)))
      }
    cat("Response variable:", slot(x, "namen")[1], "\n")
    cat("Total response variance:", x@var.y, "\n")
    if (length(x@nobs)>0) cat("Analysis based on", x@nobs, "observations", "\n")

    cat("\n")
    cat(p + length(slot(x, "always")), "Regressors:", "\n", sep = " ")

    if (!is.null(x@always)) 
    {
      cat("Proportion of variance explained: ", round(100 * 
          x@R2, 2), "%", "\n", "\n", sep = "")
      if (length(slot(x, "always")) == 1) {
         cat("One Regressor always included in model:", "\n", 
             paste(slot(x, "alwaysnam"), collapse = " "), "\n")
         cat(round(100 * (x@R2 - x@R2.decomp), 2), "%", "of variance explained by this regressor", 
          "\n", sep = " ") 
      } else {
         cat(length(slot(x, "always")), "Regressors always included in model:", "\n", 
             paste(slot(x, "alwaysnam"), collapse = " "), "\n")     
         cat(round(100 * (x@R2 - x@R2.decomp), 2), "%", "of variance explained by these", 
          "\n", sep = " ")
      }

    if (grouped) {
        cat("\n","Some regressors combined in groups:", "\n")
        cat("        Group ", x@groupdocu[[1]][[1]], ":", x@groupdocu[[2]][[1]], "\n")
        for (j in 2:g) {
          if (length(slot(x,"groupdocu")[[2]][[j]]) > 1)
          cat("        Group ", x@groupdocu[[1]][[j]], ":", x@groupdocu[[2]][[j]], "\n")
        }
       cat("\n","Relative importance of", g, "(groups of) regressors assessed:", "\n",
           paste(slot(x, "groupdocu")[[1]], collapse = " "), "\n", "\n") 
       cat(round(100 * x@R2.decomp, 2), "%", "of variance decomposed among these", "\n",
           sep = " ")
      }
 
    else {
       cat("\n", "Relative importance of", p, "regressors assessed:", "\n",
           paste(slot(x, "namen")[2:(p + 1)], collapse = " "), "\n") 
       cat(round(100 * x@R2.decomp, 2), "% of variance decomposed among these", "\n",
           sep = " ")
         cat("\n")
      }
    } 

    if (is.null(x@always)) {
     if (grouped) {
        cat("Some regressors combined in groups:", "\n")
        cat("        Group ", x@groupdocu[[1]][[1]], ":", x@groupdocu[[2]][[1]], "\n")
        for (j in 2:g) {
        if (length(slot(x,"groupdocu")[[2]][[j]]) > 1)
        cat("        Group ", x@groupdocu[[1]][[j]], ":", x@groupdocu[[2]][[j]], "\n")
        }
        cat("\n", "Relative importance of", g, "(groups of) regressors assessed:", "\n",
           paste(slot(x, "groupdocu")[[1]], collapse = " "), "\n", "\n") 
        cat("Proportion of variance explained by model: ", round(100 * 
           x@R2, 2), "%", "\n", sep = "")
      }
      else {
       cat(paste(slot(x, "namen")[2:(p + 1)], 
          collapse = " "), "\n")
       cat("Proportion of variance explained by model: ", round(100 * 
          x@R2, 2), "%", "\n", sep = "")
      }
    }
    if (x@rela) cat("Metrics are normalized to sum to 100% (rela=TRUE).", "\n")
    else cat("Metrics are not normalized (rela=FALSE).", "\n")
    cat("\n")

    cat("Relative importance metrics:", "\n")
    cat("\n")
    type <- slot(x, "type")
    if (!grouped) print(matrix(cbind(x@lmg, x@pmvd, x@last, x@first, x@betasq, 
        x@pratt, x@genizi, x@car), p, length(type), dimnames = list(x@namen[2:(p + 
        1)], type)))
    else print(matrix(cbind(x@lmg, x@pmvd, x@last, x@first), 
        g, length(type), dimnames = list(x@groupdocu[[1]], type)))

    if (show.coeffs & !is.null(x@ave.coeffs)) {
      cat("\n")
      cat("Average coefficients for different model sizes:", "\n")
      cat("\n")
      print(x@ave.coeffs)
    }
}


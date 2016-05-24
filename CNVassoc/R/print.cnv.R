print.cnv <-
function (x, digits = 4, ...)
{
    cat("\n")
    if (!is.null(attr(x, "mixture"))) {
        cat("Inferred copy number variant by a quantitative signal\n")
        if (is.null(attr(x, "batches")))
          cat("   Method:", attr(x, "mixture")$method, " \n\n")
        else 
          cat("   Method:", attr(x, "mixture")[[1]]$method, " \n\n")        
    }
    else {
        cat("-. Copy number variant\n   Input data: called probabilities\n")
    }
    cat("-. Number of individuals:", length(x), "\n")
    cat("-. Copies", paste(attr(x, "num.copies"), collapse = ", "),
        "\n")
    if (!is.null(attr(x, "means"))){
        means <- attr(x, "means")
        if (is.vector(means))
          cat("-. Estimated means:", paste(round(attr(x, "means"), digits = digits), collapse = ", "), "\n")
        else{
          colnames(means)<-paste("CNV",attr(x, "num.copies"))
          cat("-. Estimated means:\n")
          print(round(means, digits = digits))
        }
    }
    if (!is.null(attr(x, "sds"))){
        sds <- attr(x, "sds")
        if (is.vector(sds))
          cat("-. Estimated variances:", paste(round(attr(x, "sds")^2, digits = digits), collapse = ", "), "\n")
        else{
          colnames(sds)<-paste("CNV",attr(x, "num.copies"))
          cat("-. Estimated variances:\n")
          print(round(sds^2, digits = digits))
        }
    }
    probs <- colMeans(attr(x, "probabilities"))
    cat("-. Estimated proportions:", paste(round(probs, digits = digits), collapse = ", "), "\n")
    if (!is.null(attr(x, "mixture"))){
      if (is.null(attr(x, "batches"))){
        if (!is.null(attr(x, "mixture")$P))
            cat("-. Goodness-of-fit test: p-value=", attr(x, "mixture")$P,"\n")
      }else{
        if (!is.null(attr(x, "mixture")[[1]]$P)){  
          cat("-. Goodness-of-fit test:\n")
          p.values <- cbind(unlist(lapply(attr(x, "mixture"), function(mixture.i) mixture.i$P)))
          rownames(p.values)<-paste("Batch",unique(attr(x, "batches")))
          colnames(p.values)<-"p-value"
          print(p.values)
        }
      }
    }
    if (length(attr(x, "copynum.range")) > 1) {
        if (!is.null(attr(x, "mixture")))
            cat("\n\n-. Note: number of classes has been selected using the best BIC\n")
    }
}


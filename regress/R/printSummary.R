summary.regress <- function(object, ...) object

print.regress <- function(x, digits=3, fixed.effects=T, ...)
  {
      cat("Likelihood kernel: K = ")
      if(length(x$kernel) == 1){
          cat(max(sign(x$kernel), 0))
      } else if(!is.null(x$Kcolnames)) cat(x$Kcolnames, sep="+")

      cat("\n\nMaximized log likelihood with kernel K is ",round(x$llik,digits),"\n",sep=" ")
      indent.lin <- max(nchar(dimnames(x$beta)[[1]]))
      indent.var <- max(nchar(x$Vnames))
      indent <- max(indent.lin,indent.var)

      extra.space <- ""
      space.var <- extra.space
      for(i in 0:(indent-indent.var)) space.var <- paste(space.var," ",sep="")
      space.lin <- extra.space
      for(i in 0:(indent-indent.lin)) space.lin <- paste(space.lin," ",sep="")

      coefficients <- cbind(x$beta,x$beta.se)
      dimnames(coefficients)[[2]] <- c("Estimate","Std. Error")
      coefficients <- round(coefficients,digits)
      if(fixed.effects) {
          cat("\nLinear Coefficients:\n")
          row.names(coefficients) <- paste(space.lin,dimnames(x$beta)[[1]],sep="")
          print(coefficients)
          cat("\n")
      } else {
          cat("\nLinear Coefficients: not shown\n\n")
      }

      ## New version of regress automatically converts to the linear
      ## scale - as if pos was a vector of zeroes

      var.coefficients <- cbind(x$sigma,sqrt(diag(as.matrix(x$sigma.cov))))
      row.names(var.coefficients) <- paste(space.var,x$Vnames,sep="")
      dimnames(var.coefficients)[[2]] <- c("Estimate","Std. Error")
      var.coefficients <- round(var.coefficients,digits)
      cat("Variance Coefficients:\n")
      print(var.coefficients)
      cat("\n")
  }


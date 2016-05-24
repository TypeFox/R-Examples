print.sensitivity.0d <- function(x, ...) {
  labs <- attr(x, "parameters")

  if(all(c('s0','s1') %in% names(labs))) {
    cat("Empty Principle Stratum: ",
        paste("S(", labs$z0, ") = ",labs$s0, ", S(", labs$z1, ") = ",
              labs$s1, sep=''),
        "\n")
  }
  
  cat("ACE:\t", paste("E(Y(", labs$z1, ") - Y(", labs$z0,") | S(", labs$z0,
                      ") = S(", labs$z1, ") = ",labs$selected,")", sep=''),
      "\n")

  ACE <- x$ACE
  print(ACE)

  ci.dim <- which(names(dimnames(x$ACE.ci)) == "ci.method")
  if(length(ci.dim) == 0)
    return(invisible(NULL))

  ci.method <- dimnames(x$ACE.ci)[[ci.dim]]
  ci.slice <- slice.index(x$ACE.ci, ci.dim)
  newdim <- dim(x$ACE.ci)[-ci.dim]
  newdimnames <- dimnames(x$ACE.ci)[-ci.dim]

  cat("\nACE confidence interval:\n")
  if("analytic" %in% ci.method) {
    cat("By analytic method\n")
    print(array(x$ACE.ci[ci.slice == which(ci.method == "analytic")],
                dim=newdim, dimnames=newdimnames))
  }

  if(all(c("bootstrap", "analytic") %in% ci.method))
    cat("\n")
  
  if("bootstrap" %in% ci.method) {
    cat("By bootstrap method, N = ", attr(x, 'N.boot'), "\n",sep='')
    print(array(x$ACE.ci[ci.slice == which(ci.method == "bootstrap")],
                dim=newdim, dimnames=newdimnames))
  }

  invisible(NULL)
}

print.sensitivity.1d <- function(x, ...) {
  labs <- attr(x, "parameters")

  if(all(c('s0','s1') %in% names(labs))) {
    cat("Empty Principle Stratum: ",
        paste("S(", labs$z0, ") = ",labs$s0, ", S(", labs$z1, ") = ",
              labs$s1, sep=''),
        "\n")
  }
  
  cat("SCE:\t", paste("P(T(", labs$z0, ") =< t | S(", labs$z0,
                      ") = S(", labs$z1, ") = ",labs$selected,") - ",
                      "P(T(", labs$z1, ") =< t | S(", labs$z0,
                      ") = S(", labs$z1, ") = ",labs$selected,")", sep=''),
      "\n")

  SCE <- x$SCE
  print(SCE)

  ci.dim <- which(names(dimnames(x$SCE.ci)) == "ci.method")
  if(length(ci.dim) == 0)
    return(invisible(NULL))
  
  ci.method <- dimnames(x$SCE.ci)[[ci.dim]]
  ci.slice <- slice.index(x$SCE.ci, ci.dim)
  newdim <- dim(x$SCE.ci)[-ci.dim]
  newdimnames <- dimnames(x$SCE.ci)[-ci.dim]

  cat("\nSCE confidence interval:\n")
  if("analytic" %in% ci.method) {
    cat("By analytic method\n")
    print(array(x$SCE.ci[ci.slice == which(ci.method == "analytic")],
                dim=newdim, dimnames=newdimnames))
  }

  if(all(c("bootstrap", "analytic") %in% ci.method))
    cat("\n")
  
  if("bootstrap" %in% ci.method) {
    cat("By bootstrap method, N = ", attr(x, 'N.boot'), "\n",sep='')
    print(array(x$SCE.ci[ci.slice == which(ci.method == "bootstrap")],
                dim=newdim, dimnames=newdimnames))
  }

  invisible(NULL)
}

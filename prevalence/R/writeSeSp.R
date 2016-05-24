writeSeSp <-
function(f, x){
  out <- character()

  if (x$d[1] == "pert"){
    pertK <- ifelse(is.null(x$p[4]), 4, x$p[4])
    pertM <- x$d[2]

    betaPERT <-
      betaPERT(a = x$p[1], m = x$p[2], b = x$p[3], k = pertK, method = pertM)

    if (pertM  == "branscum"){
      distSpec <-
        paste0("~ dbeta(", betaPERT$alpha, ", ", betaPERT$beta, ")")

    } else {
      fp <- paste0(gsub("\\]", "", gsub("\\[", "", f)), "p")
      dif <- betaPERT$b - betaPERT$a
      distSpec1 <-
        paste0("<- ", fp, " * ", dif, " + ", betaPERT$a)
      distSpec2 <-
        paste0("~ dbeta(", betaPERT$alpha, ", ", betaPERT$beta, ")")
      out <- c(out, paste(f, distSpec1))
      out <- c(out, paste(fp, distSpec2))
    }
  }

  if (x$d[1] == "fixed"){
    distSpec <- paste0("<- ", x$p)
    out <- c(out, paste(f, distSpec))
  }

  if (x$d[1] == "uniform"){
    distSpec <- paste0("~ dunif(", x$p[1], ", ", x$p[2], ")")
    out <- c(out, paste(f, distSpec))
  }

  if (x$d[1] == "beta"){
    distSpec <- paste0("~ dbeta(", x$p[1], ", ", x$p[2], ")")
    out <- c(out, paste(f, distSpec))
  }

  return(out)
}
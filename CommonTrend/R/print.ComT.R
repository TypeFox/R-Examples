"print.ComT" <-
function(x,  ...){
  title <- paste(x$method)
  cat("According to ", title,"\n", sep="")
  cat("\n")
  cat("Othogonal Complement of Beta:", "\n")
  print(x$othog.beta)
  cat("\n")
  cat("Othogonal Complement of Alpha:", "\n")
  print(x$othog.alpha)
  cat("\n")
  cat("Loading Vector:", "\n")
  print(x$loading.vector)
  }
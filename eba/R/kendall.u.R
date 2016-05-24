kendall.u <- function(M, cont.correct = FALSE){
  # Kendall's (1940, 1962) coefficient of agreement
  # Assumptions: equal number of obs. per pair
  #              one obs. per judge and pair
  # Maximum agreement: u = 1, the smaller u the less agreement
  # Chi2 test, H_0: agreement is by chance
  # m observers (judges), n stimuli
  # last mod: 17/OCT/2006 fw

  ifelse(cont.correct, c <- 1, c <- 0)
  Sigma <- sum( c(choose(M[upper.tri(M)],2), choose(M[lower.tri(M)],2)) )
  u <- 2*Sigma / ( choose(m<-M[1,2]+M[2,1],2) * choose(n<-dim(M)[1],2) ) - 1
  min.u <- ifelse(m %% 2, -1/m, -1/(m-1))
  chi <- 4/(m-2) * (Sigma-c - 1/2*choose(n,2) * choose(m,2) * (m-3)/(m-2))
  df <- round( choose(n,2) * (m*(m-1))/(m-2)^2 )
  out <- list(u=u, min.u=min.u, chi=chi, df=df, p=1-pchisq(chi, df),
              cont.correct=cont.correct)
  class(out) <- "kendall.u"
  out
}


print.kendall.u <- function(x, digits = max(3,getOption("digits")-4), ...){
  cat("\nKendall's u coefficient of agreement\n\n")
  cat("u = ",x$u, ", minimum u = ",x$min.u, "\nchi2 = ",x$chi, ", df = ",x$df,
      ", p-value = ",x$p, "\n",sep="")
  cat("alternative hypothesis: between-judges agreement is not by chance\n")
  cat("correction for continuity has ", ifelse(x$cont.correct, "", "not "),
      "been applied\n", sep="")
  cat("\n")
  invisible(x)
}

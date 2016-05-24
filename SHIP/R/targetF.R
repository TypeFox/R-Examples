targetF <-
function(x,genegroups) {
  covm <- cov(x)
  corm <- cor(x)
  var.only <- covm/corm
  diag(corm) <- 0
  cora  <- sum(corm)/sum(corm!=0)
  T <- cora*var.only
  diag(T)<- diag(covm)
  T
}


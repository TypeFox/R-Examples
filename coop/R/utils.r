check_use <- function(use)
{
  match.arg(tolower(use), c("everything", "all.obs", "complete.obs"))
}

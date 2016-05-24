#' Get git stamp (commit and branch) for a repository
#'
#' @param repo Git repo directory. If unspecified, then the current working
#' directory is used.
#'
#' @description
#'  The function returns the latest git commit and branch for the repo specified
#'  in \code{repo} or the current working directory if unspecified.
#'  Git is needed for the command to run.
#'  The functions makes it possible to include the latest git commit and branch
#'  in a run to be able to know exactly which code where used.
#'
#' @return character vector with latest commit and branch
#'
#' @export
#'
git_stamp <- function(repo = getwd()){
    start_wd <- getwd()
    setwd(repo)
    git_commit <- system(paste0("git rev-parse --verify HEAD"), intern = TRUE)
    git_branch <- system(paste0("git rev-parse --abbrev-ref HEAD"),
    intern = TRUE)
    git_message <- system(paste0("git log --format=%B -n 1 ", git_commit),
                        intern = TRUE)[1]
    setwd(start_wd)
    c("commit" = git_commit, "branch" = git_branch, "message" = git_message)
}

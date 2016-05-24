##    archivist package for R
##		archivist.github package for R
##
#' @title Clone Github Repository
#'
#' @description
#' \code{cloneGitHubRepo} is a wrapper around \code{git clone} and clones GitHub Repository
#' into the \code{repoDir} directory.
#' 
#' @details
#' To learn more about  \code{Archivist Integration With GitHub} visit \link{agithub}.
#' @param repoURL The remote repository to clone.
#' @param repoDir Local directory to clone to. If \code{NULL}, by default, creates a local directory,
#' which corresponds to the name after last \code{/} in \code{repoURL}.
#' @param default Sets cloned Repository as default Local and GitHub Repository. 
#' If \code{default = TRUE} then \code{repoDir} (last piece of \code{repoURL}) is set as default Local Repository 
#'  and for GitHub repository also the \code{user} from  \code{repoURL} is set as default GitHub user).
#' @param ... Further parameters passed to \link[git2r]{clone}.
#' 
#' @author 
#' Marcin Kosinski, \email{m.p.kosinski@@gmail.com}
#'
#' 
#' @examples 
#' \dontrun{
#' 
#' cloneGitHubRepo("https://github.com/MarcinKosinski/Museum")
#' cloneGitHubRepo("https://github.com/MarcinKosinski/Museum-Extra")
#' 
#' 
#' # empty Github Repository creation
#' 
#' library(httr)
#' myapp <- oauth_app("github",
#'                    key = app_key,
#'                    secret = app_secret)
#' github_token <- oauth2.0_token(oauth_endpoints("github"),
#'                                myapp,
#'                                scope = "public_repo")
#' # setting options                              
#' aoptions("github_token", github_token)
#' aoptions("name", user_name)
#' aoptions("password", user_password)
#' 
#' createEmptyGithubRepo("archive-test4")
#' setRemotebRepo(aoptions("name"), "archive-test4")
#' ## artifact's archiving
#' example <- 1:100
#' 
#' # archiving
#' archive(example) -> md5hash_path
#' 
#' ## proof that artifact is really archived
#' showRemoteRepo() # uses options from setGithubRepo
#' # let's remove przyklad
#' rm(example)
#' # and load it back from md5hash_path
#' aread(md5hash_path)
#' 
#' 
#' # clone example
#' unlink("archive-test", recursive = TRUE)
#' cloneGitHubRepo('https://github.com/MarcinKosinski/archive-test')
#' setRemoteRepo(aoptions("name"), "archive-test")
#' # equivalent is cloneGitHubRepo('https://github.com/MarcinKosinski/archive-test', default = TRUE)
#' # check if default is set with
#' # aoptions('repoDir'); aoptions('repo'); aoptions('user')
#' data(iris)
#' archive(iris)
#' showRemoteRepo()
#' 
#' 
#' }
#' @family archivist.github
#' @rdname cloneGitHubRepo
#' @export
cloneGitHubRepo <- function(repoURL, repoDir = NULL, default = FALSE, ...){
	
	stopifnot((is.character(repoDir) & length(repoDir) == 1) | is.null(repoDir))
	stopifnot( is.logical( default ), length( default ) == 1 )
	
	if (is.null(repoDir)) {
		repoDir <-tail(strsplit(repoURL,
														"/")[[1]],1)
	}
	
	if (!file.exists(repoDir)) {
		dir.create(repoDir)
	}
	git2r::clone(repoURL, repoDir, ...) -> repo2return
	
	
	if (default) {
		archivist::aoptions('repoDir', repoDir)
		archivist::aoptions('user', tail(strsplit(repoURL,
																			 "/")[[1]],2)[1])
		archivist::aoptions('repo',tail(strsplit(repoURL,
																			 "/")[[1]],1))
	}
	return(repo2return)
}

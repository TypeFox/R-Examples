##    archivist package for R
##    archivist.github package for R
##
#' @title Delete the Existing GitHub Repository
#'
#' @description
#' \code{deleteGitHubRepo} can delete whole GitHub-Repository or only archivist-like \link{Repository}
#' stored on a GitHub-Repository (then it requires more parameters to be specified).
#' 
#' @details
#' To learn more about  \code{Archivist Integration With GitHub} visit \link{agithub}.
#'  
#' @param deleteRoot A logical value that specifies if the repository root directory
#' should be deleted for Local Repository or for GitHub whether to delete whole GitHub-Repository.
#' @param unset A logical. If deleted \code{repoDir/repo} was set to be default Local/GitHub Repository
#' and \code{unset} is TRUE, then \code{repoDir/repo} is unset as a default Local/GitHub Repository (\code{aoptions('repoDir/repo', NULL, T)}).
#' @param repo While working with a Github repository. A character denoting GitHub repository name to be deleted.
#' @param github_token While working with a Github repository. An OAuth GitHub Token created with the \link{oauth2.0_token} function. To delete GitHub Repository you
#' need to have \code{delete_repo} scope set - see examples. See \link{archivist-github-integration}.
#' Can be set globally with \code{aoptions("github_token", github_token)}.
#' @param user While working with a Github repository. A character denoting GitHub user name. Can be set globally with \code{aoptions("user", user)}.
#'  See \link{archivist-github-integration}.
#' @param response A logical value. Should the GitHub API response be returned (only when \code{deleteRoot = TRUE}).
#' @param subdir Only when \code{deleteRoot = FALSE}. Subdirectory in which the archivist-Repository is stored on a GitHub Repository.
#' If it's in main directory, then specify to \code{NULL} (default).
#' @param password Only when \code{deleteRoot = FALSE}. While working with a Github repository. A character denoting GitHub user password. Can be set globally with \code{aoptions("password", password)}.
#' See \link{archivist-github-integration}.
#' 
#' @details
#' 
#' To delete GitHub Repository you
#' need to have \code{delete_repo} scope set - see examples.
#' 
#' 
#' @author 
#' Marcin Kosinski, \email{m.p.kosinski@@gmail.com}
#'
#' @examples
#' \dontrun{
#' 
#' ########################
#' #### GitHub version ####
#' ########################
#' 
#' library(httr)
#' myapp <- oauth_app("github",
#'                    key = app_key,
#'                    secret = app_secret)
#' github_token <- oauth2.0_token(oauth_endpoints("github"),
#'                                 myapp,
#'                                 scope = c("public_repo",
#'                                           "delete_repo"))
#' aoptions("github_token", github_token)
#' aoptions("user", user)
#' aoptions("password", password)
#' 
#' createGitHubRepo("Museum")
#' deleteGitHubRepo(repo = "Museum", deleteRoot = TRUE, response = TRUE)
#' 
#' }
#' 
#' @family archivist.github
#' @rdname deleteRepo
#' @export
deleteGitHubRepo <- function(repo,
														 github_token = aoptions("github_token"), 
														 user = aoptions("user"),
														 password = aoptions("password"),
														 unset = FALSE, 
														 deleteRoot = FALSE, 
														 subdir = NULL, 
														 response = aoptions("response")) {
	stopifnot(is.character(repo) & length(repo) ==1)
	stopifnot(is.character(user) & length(user)==1)
	
	if (deleteRoot) {
		httr::DELETE(url = file.path("https://api.github.com/repos",user,repo),
					 config = httr::config(token = github_token)
		) -> resp
	} else {
		tempfile() -> tmpDir
		# clone repo to tmpDir
		cloneGitHubRepo(file.path("https://github.com/", user, repo), 
										repoDir = tmpDir) -> clonedRepo
		# remove archivist-repository
		archivist::deleteLocalRepo(repoDir = file.path(tmpDir, 
																				ifelse(is.null(subdir), "", subdir)
		))
		# add changes to git 
		git2r::add(clonedRepo, c("backpack.db", "gallery/"))
		# message
		delete_commit <- git2r::commit(clonedRepo, "archivist Repository deletion.")
		# GitHub authorization
		# to perform pull and push operations
		cred <- git2r::cred_user_pass(user,
																	password)
		# push archivist-like Repository deletion to GitHub
		git2r::push(clonedRepo,
				 refspec = "refs/heads/master",
				 credentials = cred)
		
	}
	
	if (unset) {
		aoptions('repo', NULL,T)
	}
	
	
	if (response & deleteRoot){
		return(resp)
	}
	
	
}




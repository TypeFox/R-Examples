##    archivist package for R
##    archivist.github package for R
##
#' @title Create an Empty Repository on GitHub
#'
#' @description
#' 
#' \code{createGitHubRepo} is a GitHub version of \link{createLocalRepo} and creates a new GitHub repository 
#' with an empty \pkg{archivist}-like \link{Repository}. It also creates a Local \code{Repository} which is git-synchronized with
#' new GitHub repository. 
#' 
#' @details
#' To learn more about  \code{Archivist Integration With GitHub} visit \link{agithub}.
#' 
#' 
#' @details
#' At least one Repository must be initialized before using other functions from the \pkg{archivist.github} package. 
#' While working in groups, it is highly recommended to create a Repository on a shared Dropbox/GitHub folder.
#' 
#' All artifacts which are desired to be archived are going to be saved in the local Repository, which is an SQLite 
#' database stored in a file named \code{backpack}. 
#' After calling \code{saveToLocalRepo} function, each artifact will be archived in a \code{md5hash.rda} file. 
#' This file will be saved in a folder (under \code{repoDir} directory) named 
#' \code{gallery}. For every artifact, \code{md5hash} is a unique string of length 32 that is produced by
#' \link[digest]{digest} function, which uses a cryptographical MD5 hash algorithm.
#' 
#' To learn more about artifacts visit \link[archivist]{archivist-package}.
#' 
#' Created \code{backpack} database is a useful and fundamental tool for remembering artifact's 
#' \code{name}, \code{class}, \code{archiving date} etc. (the so called \link{Tags})
#' or for keeping artifact's \code{md5hash}.
#' 
#' Besides the \code{backpack} database, \code{gallery} folder is created in which all 
#' artifacts will be archived.
#' 
#' After every \code{saveToLocalRepo} call the database is refreshed. As a result, the artifact is available 
#' immediately in \code{backpack.db} database for other collaborators.
#' 
#' @param repoDir A character that specifies the directory for the Repository which is to be made. While working with GitHub Repository, this will
#' be the directory of the synchronized Local Repository, in which the new Local Repository will be created (is \code{NULL} then is the same as \code{repo}).
#' 
#' @param default If \code{default = TRUE} then \code{repoDir} (\code{repo}) is set as default local repository. Also the \code{user} is set as default GitHub user.
#' 
#' @param ... further arguments passed to \link{createLocalRepo} such as \code{force}.
#' 
#' @param repo While working with a Github repository. A character denoting new GitHub repository name. White spaces will be substitued with a dash.
#' @param github_token While working with a Github repository. An OAuth GitHub Token created with the \link{oauth2.0_token} function. See \link{archivist-github-integration}.
#' Can be set globally with \code{aoptions("github_token", github_token)}.
#' @param repoDescription While working with a Github repository. A character specifing the new GitHub repository description.
#' @param user While working with a Github repository. A character denoting GitHub user name. Can be set globally with \code{aoptions("user", user)}.
#'  See \link{archivist-github-integration}.
#' @param password While working with a Github repository. A character denoting GitHub user password. Can be set globally with \code{aoptions("password", password)}.
#' See \link{archivist-github-integration}.
#' @param readmeDescription While working with a Github repository. A character of the content of \code{README.md} file. By default a description of \link{Repository}.
#' Can be set globally with \code{aoptions("readmeDescription", readmeDescription)}. In order to omit 
#' \code{README.md} file set \code{aoptions("readmeDescription", NULL)}.
#' @param response A logical value. Should the GitHub API response be returned.
#' 
#' @author 
#' Marcin Kosinski, \email{m.p.kosinski@@gmail.com}
#'
#' @examples
#' \dontrun{
#' ## GitHub version
#' 
#' library(httr)
#' myapp <- oauth_app("github",
#'                    key = app_key,
#'                    secret = app_secret)
#' github_token <- oauth2.0_token(oauth_endpoints("github"),
#'                                 myapp,
#'                                 scope = "public_repo")
#' aoptions("github_token", github_token)
#' aoptions("user", user)
#' aoptions("password", password)
#' 
#' createGitHubRepo("Museum")
#' createGitHubRepo("Museum-Extras", response = TRUE)
#' createGitHubRepo("Gallery", readmeDescription = NULL)
#' createGitHubRepo("Landfill", 
#' repoDescription = "My models and stuff") 
#' createGitHubRepo("MuseumYYYY", repoDir = "Museum_YY")
#'         
#'         
#'         
#' # empty GitHub Repository creation
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
#' aoptions("user", user)
#' aoptions("password", password)
#' 
#' createGitHubRepo("archive-test4", default = TRUE)
#' ## artifact's archiving
#' przyklad <- 1:100
#' 
#' # archiving
#' archive(przyklad) -> md5hash_path
#' 
#' ## proof that artifact is really archived
#' showRemoteRepo() # uses options from setGithubRepo
#' # let's remove przyklad
#' rm(przyklad)
#' # and load it back from md5hash_path
#' aread(md5hash_path)
#' 
#' 
#' # clone example
#' unlink("archive-test", recursive = TRUE)
#' cloneGithubRepo('https://github.com/MarcinKosinski/archive-test')
#' setRemoteRepo(aoptions("user"), "archive-test")
#' data(iris)
#' archive(iris)
#' showRemoteRepo()
#' 
#' }
#' @family archivist.github
#' @rdname createEmptyRepo
#' @export
createGitHubRepo <- function(repo,
														 github_token = aoptions("github_token"), 
														 user = aoptions("user"),
														 repoDir = NULL,
														 #user.email = aoptions("user.email"),
														 password = aoptions("password"),
														 repoDescription = aoptions("repoDescription"),
														 readmeDescription = aoptions("readmeDescription"),
														 response = aoptions("response"),
														 default = FALSE,
														 ...){
	stopifnot(is.character(repo) & length(repo) ==1)
	stopifnot((is.character(repoDir) & length(repoDir) ==1) | (is.null(repoDir)))
	stopifnot(is.character(repoDescription) & length(repoDescription) ==1)
	#stopifnot(any(class(github_token) %in% "Token"))
	stopifnot(is.character(user) & length(user)==1)
	#stopifnot(is.character(user.email) & length(user.email)==1)
	stopifnot(is.character(password) & length(password)==1)
	stopifnot((is.character(readmeDescription) & length(readmeDescription)==1) |
							is.null(readmeDescription))
	stopifnot(is.logical(response) & length(response) ==1)
	
	
	stopifnot( is.logical( default ), length( default ) == 1 )
	repo <- gsub(pattern = " ", "-", repo)
	
	shortPath <- FALSE
	if(is.null(repoDir)) {
		shortPath <- TRUE
		repoDir <- repo
	}
	
	# httr imports are in archivist-package.R file
	# creating an empty GitHub Repository
	httr::POST(url = "https://api.github.com/user/repos",
			 encode = "json",
			 body = list(
			 	name = jsonlite::unbox(repo),
			 	description = jsonlite::unbox(repoDescription)
			 ),
			 config = httr::config(token = github_token)
	) -> resp
	
	
	# git2r imports are in the archivist-package.R
	#path <- repoDir
	dir.create(repoDir)
	
	if (!shortPath){
		dir.create(file.path(repoDir, repo))
		repoDir_path <- file.path(repoDir, repo)
	} else {
		repoDir_path <- repo
	}
	
	# initialize local git repository
	# git init
	repoDir_git2r <- git2r::init(repoDir_path)
	
	## Create and configure a user
	# git config - added to Note section
	#git2r::config(repo, ...) # if about to use, the add to archivist-package.R
	
	# archivist-like Repository creation
	archivist::createLocalRepo(repoDir = repoDir_path, ...)
	file.create(file.path(repoDir_path, "gallery", ".gitkeep"))
	# git add
	if (!is.null(readmeDescription)){
		file.create(file.path(repoDir_path, "README.md"))
		writeLines(aoptions("readmeDescription"), file.path(repoDir_path, "README.md"))
		git2r::add(repoDir_git2r, c("backpack.db", "gallery/", "README.md"))
	} else {
		git2r::add(repoDir_git2r, c("backpack.db", "gallery/"))
	}
	
	# git commit
	new_commit <- git2r::commit(repoDir_git2r, "archivist Repository creation.")
	
	# association of the local and GitHub git repository
	# git add remote
	git2r::remote_add(repoDir_git2r,
						 #"upstream2",
						 'origin',
						 file.path("https://github.com", user, paste0(repo, ".git")))
	
	# GitHub authorization
	# to perform pull and push operations
	cred <- git2r::cred_user_pass(user,
																password)
	
	# push archivist-like Repository to GitHub repository
	git2r::push(repoDir_git2r,
			 #name = "upstream2",
			 refspec = "refs/heads/master",
			 credentials = cred)
	
	if (response){
		return(resp)
	}
	
	if (default) {
		archivist::aoptions('repoDir',repoDir_path)
		archivist::aoptions('repo', repo)
		archivist::aoptions('user', user)
	}
	
}



checkDirectory <- function( directory, create = FALSE ){
	# check if global repository was specified by setLocalRepo
	if ( is.null(directory) ){
		
		directory <- aoptions("repoDir")
	}
	# check whether it is second call of checkDirectory 
	# (e.g CreatEmptyRepo + default = TRUE)
	#   if ( grepl("/$", x = directory , perl=TRUE) ){
	#     directory <- gsub(pattern = ".$", replacement = "",
	#                       x = directory, perl = TRUE)
	#   }
	# check property of directory
	if ( !create ){
		# check whether repository exists
		if ( !dir.exists( directory ) ){
			stop( paste0( "There is no such repository as ", directory ) )
		}
		# check if repository is proper (has backpack.db and gallery)
		if ( !all( c("backpack.db", "gallery") %in% list.files(directory) ) ){
			stop( paste0( directory, " is not a proper repository. There is neither backpack.db nor gallery." ) )
		}
	}
	# check if repoDir has "/" at the end and add it if not
	#   if ( !grepl("/$", x = directory , perl=TRUE) ){
	#     directory <- paste0(  directory, "/"  )
	#   }
	return( directory )
}

# checkDirectory2 <- function( directory ){
#   check if repoDir has "/" at the end and add it if not
#   if ( !grepl("/$", x = directory , perl=TRUE) ){
#     directory <- paste0(  directory, "/"  )
#   }
#   return( directory )
# }

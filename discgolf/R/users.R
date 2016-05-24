#' Work with users
#'
#' @name users
#' @param username A user name
#' @param type A type of user, one of active, staff, new, suspended, blocked, or suspect
#' @param name a name
#' @param email an email address
#' @param password a password
#' @param user_id a user id
#' @param new_username a username
#' @template args
#' @section users_list:
#' note that there is no paging, so if you have more than 100 users, you only
#' get the first 100. :sad panda:
#' @examples \dontrun{
#' # list a user
#' user('sckott')
#' user('cboettig')
#'
#' # list users
#' users_list('staff')
#' users_list('new')
#'
#' # create a user
#' (x <- user_create("jane doe", "jane@doe.com", "jane_doe", "afafasfdasdf"))
#'
#' # activate a user
#' user_activate(x$user_id)
#'
#' # upate email address
#' user_update_email('jane_doe', 'jane2@doe.com')
#'
#' # upate user name
#' user_update_username('jane_doe', 'jane_doe2')
#'
#' # delete a user
#' user_delete(x$user_id)
#' }

#' @export
#' @rdname users
user <- function(username, url = NULL, key = NULL, user = NULL, ...) {
  args <- dc(list(api_key = check_key(key), api_username = check_user(user)))
  disc_GET(check_url(url), sprintf("users/%s.json", username), args, ...)
}

#' @export
#' @rdname users
users_list <- function(type, url=NULL, key=NULL, user=NULL, ...){
  args <- dc(list(api_key = check_key(key), api_username = check_user(user)))
  disc_GET(check_url(url), sprintf("admin/users/list/%s.json", type), args, ...)
}

#' @export
#' @rdname users
user_create <- function(name, email, username, password, url=NULL, key=NULL, user=NULL, ...){
  args <- dc(list(api_key = check_key(key), api_username = check_user(user),
                  name = name, email = email, username = username, password = password))
  disc_POST(check_url(url), "users", args, ...)
}

#' @export
#' @rdname users
user_activate <- function(user_id, url=NULL, key=NULL, user=NULL, ...){
  args <- dc(list(api_key = check_key(key), api_username = check_user(user)))
  disc_PUT(check_url(url), sprintf("admin/users/%s/activate", user_id), args, ...)
}

#' @export
#' @rdname users
user_delete <- function(user_id, url=NULL, key=NULL, user=NULL, ...){
  args <- dc(list(api_key = check_key(key), api_username = check_user(user)))
  disc_DELETE(check_url(url), sprintf("admin/users/%s", user_id), args, ...)
}

#' @export
#' @rdname users
user_update_email <- function(username, email, url=NULL, key=NULL, user=NULL, ...){
  args <- dc(list(api_key = check_key(key), api_username = check_user(user),
                  email = email))
  disc_PUT(check_url(url), sprintf("users/%s/preferences/email", username), args, ...)
}

#' @export
#' @rdname users
user_update_username <- function(username, new_username, url=NULL, key=NULL, user=NULL, ...){
  args <- dc(list(api_key = check_key(key), api_username = check_user(user),
                  new_username = new_username))
  disc_PUT(check_url(url), sprintf("users/%s/preferences/username", username), args, ...)
}

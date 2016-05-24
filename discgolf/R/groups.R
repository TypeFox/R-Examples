#' Work with groups
#'
#' @export
#' @param name (character) A group name. required
#' @param id (numeric) A group id. required
#' @template args
#' @examples \dontrun{
#' # all groups
#' groups()
#'
#' # specific group by name
#' group_members("admins")
#' group_members("moderators")
#' group_members("trust_level_3")
#'
#' # create a group
#' (res <- group_create("group_testing"))
#'
#' # delete a group
#' group_delete(res$basic_group$id)
#' }
groups <- function(url = NULL, key = NULL, user = NULL, ...){
  args <- dc(list(api_key = check_key(key), api_username = check_user(user)))
  disc_GET(check_url(url), "admin/groups.json", args, ...)
}

#' @export
#' @rdname groups
group_members <- function(name, url = NULL, key = NULL, user = NULL, ...){
  args <- dc(list(api_key = check_key(key), api_username = check_user(user)))
  disc_GET(check_url(url), sprintf("groups/%s/members.json", name), args, ...)
}

#' @export
#' @rdname groups
group_create <- function(name, url = NULL, key = NULL, user = NULL, ...) {
  args <- dc(list(name = name, api_key = check_key(key), api_username = check_user(user)))
  disc_POST(check_url(url), "admin/groups", args, ...)
}

#' @export
#' @rdname groups
group_delete <- function(id, url = NULL, key = NULL, user = NULL, ...) {
  args <- dc(list(api_key = check_key(key), api_username = check_user(user)))
  disc_DELETE(check_url(url), sprintf("admin/groups/%s.json", id), args, ...)
}

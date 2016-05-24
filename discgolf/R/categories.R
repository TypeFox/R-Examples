#' Work with categories
#'
#' @name categories
#' @param category Category name. required
#' @param color A color by name or hex string. optional
#' @param text_color A color by name or hex string. optional
#' @param description Description of the category. optional
#' @param permissions Permissions - a list with the group name and permission_type
#' which is an integer: 1 = Full, 2 = Create Post, 3 = Read Only. optional
#' @param parent_category x. optional
#' @template args
#' @details Apprently there's no ability to delete categories via the API.
#' @examples \dontrun{
#' # all categories
#' categories()
#'
#' # a specfic category
#' category("questions")
#' category("packages")
#'
#' # latest topics for a category
#' category_latest_topics("packages")
#'
#' # top topics for a category
#' category_top_topics("packages")
#'
#' # new topics for a category
#' category_new_topics("packages")
#'
#' # create a category
#' category_create("stuff", "F7941D", "FFFFFF", "My new category")
#' }
categories <- function(url = NULL, key = NULL, user = NULL, ...){
  args <- dc(list(api_key = check_key(key), api_username = check_user(user)))
  disc_GET(check_url(url), "categories.json", args, ...)
}

#' @export
#' @rdname categories
category <- function(category, url = NULL, key = NULL, user = NULL, ...){
  args <- dc(list(api_key = check_key(key), api_username = check_user(user)))
  disc_GET(check_url(url), sprintf("c/%s.json", category), args, ...)
}

#' @export
#' @rdname categories
category_latest_topics <- function(category, url = NULL, key = NULL, user = NULL, ...){
  args <- dc(list(api_key = check_key(key), api_username = check_user(user)))
  disc_GET(check_url(url), sprintf("c/%s/l/latest.json", category), args, ...)
}

#' @export
#' @rdname categories
category_top_topics <- function(category, url = NULL, key = NULL, user = NULL, ...){
  args <- dc(list(api_key = check_key(key), api_username = check_user(user)))
  disc_GET(check_url(url), sprintf("c/%s/l/top.json", category), args, ...)
}

#' @export
#' @rdname categories
category_new_topics <- function(category, url = NULL, key = NULL, user = NULL, ...){
  args <- dc(list(api_key = check_key(key), api_username = check_user(user)))
  disc_GET(check_url(url), sprintf("c/%s/l/new.json", category), args, ...)
}

#' @export
#' @rdname categories
category_create <- function(category, color, text_color, description = NULL,
                            permissions = NULL, parent_category = NULL,
                            url = NULL, key = NULL, user = NULL, ...){

  args <- dc(list(name = category, color = color, text_color = text_color,
                  description = description, permissions = permissions,
                  parent_category_id = parent_category,
                  api_key = check_key(key), api_username = check_user(user)))
  disc_POST(check_url(url), "categories", args, ...)
}

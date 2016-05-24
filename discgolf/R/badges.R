#' Work with badges
#'
#' @export
#' @param username A current user name. required
#' @param name (character) A badge name to be created. required
#' @param badge_type_id (integer) A badge type id. one of 1 (gold), 2 (silver),
#' or 3 (bronze). required
#' @param description (character) A description(character) . optional
#' @param allow_title (logical) Allow a title. optional
#' @param multiple_grant (character) multiple grant. optional
#' @param icon (character) icon. optional
#' @param listable (logical) listable. optional
#' @param target_posts (character) . optional
#' @param query (character) query string. optional
#' @param enabled (logical) enabled or not. optional
#' @param auto_revoke (character) auto revoke. optional
#' @param badge_grouping_id (character) badge grouping id. optional
#' @param trigger (character) trigger. optional
#' @param show_posts (logical) show posts. optional
#' @param image (character) an image. optional
#' @param long_description (character) long description. optional
#' @template args
#' @examples \dontrun{
#' # all badges
#' badges()
#'
#' # specific group by name
#' badges_user("sckott")
#'
#' # create a group
#' (res <- badge_create(name = "funbadge", badge_type_id = 3))
#' }
badges <- function(url = NULL, key = NULL, user = NULL, ...){
  args <- dc(list(api_key = check_key(key), api_username = check_user(user)))
  disc_GET(check_url(url), "admin/badges.json", args, ...)
}

#' @export
#' @rdname badges
badges_user <- function(username, url = NULL, key = NULL, user = NULL, ...){
  args <- dc(list(api_key = check_key(key), api_username = check_user(user)))
  disc_GET(check_url(url), sprintf("users/%s/activity/badges.json", username), args, ...)
}

#' @export
#' @rdname badges
badge_create <- function(name, badge_type_id, description = NULL, allow_title = NULL,
                         multiple_grant = NULL, icon = NULL, listable = NULL,
                         target_posts = NULL, query = NULL, enabled = NULL,
                         auto_revoke = NULL, badge_grouping_id = NULL,
                         trigger = NULL, show_posts = NULL, image = NULL,
                         long_description = NULL, url = NULL, key = NULL, user = NULL, ...) {

  args <- dc(list(name = name, badge_type_id = badge_type_id, description = description,
                  allow_title = allow_title, multiple_grant = multiple_grant, icon = icon,
                  listable = listable, target_posts = target_posts, query = query,
                  enabled = enabled, auto_revoke = auto_revoke, badge_grouping_id = badge_grouping_id,
                  trigger = trigger, show_posts = show_posts, image = image,
                  long_description = long_description,
                  api_key = check_key(key), api_username = check_user(user)))
  disc_POST(check_url(url), "admin/badges.json", args, ...)
}

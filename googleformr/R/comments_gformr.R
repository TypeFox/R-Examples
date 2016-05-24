#' Function to Collect Helful Comments about Googleformr
#'
#' Function to Collect Helful Comments about googleformr.
#'
#' @param post_content the option of the function returned from gformr
#' @param active an option to make the function active or not
#' @export
#' @include gformr.R make_url.R get_form_id.R get_form_entry_ids.R get_form.R
#' @return invisible output of httr::POST
#' @examples
#' \dontrun{
#' comments_gformr(post_content = "It would useful to be able to distinguish
#' between lots of comments from one useR vs. many useRs")
#' }

comments_gformr <- gformr(form = "https://docs.google.com/forms/d/1doF8dvWKqbTBS_hQIj_QpYtGc1_fl8GdIN-X08hJ644/viewform"
                      , custom_reply  = "Thanks for supporting googleformr!")

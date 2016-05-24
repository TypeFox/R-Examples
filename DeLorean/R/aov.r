#' Perform an analysis of variance to select genes for the DeLorean model.
#'
#' @param dl DeLorean object.
#'
#' @export
#'
aov.dl <- function(dl) {
  dl$aov <- melt.expr(dl) %>%
    dplyr::left_join(dl$cell.meta) %>%
    dplyr::group_by(gene) %>%
    dplyr::do(broom::tidy(aov(x ~ capture, .))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(term == 'capture') %>%
    dplyr::arrange(p.value)
  dl
}

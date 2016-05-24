#' Estimate the cell sizes according to Anders & Huber
#' Differential expression analysis for sequence count data
#'
#' @param expr.l Melted expression values.
#'
anders.huber.cell.sizes <- function(expr.l) {
    (
        expr.l
        %>% group_by(gene)
        %>% dplyr::summarise(mu=mean(x))
        %>% left_join(expr.l)
        %>% group_by(cell)
        %>% dplyr::summarise(size=median(x-mu))
    )
}


#' Estimate the cell sizes per capture. Only uses genes that are expressed
#' in more than half the cells.
#'
#' @param expr.l Melted expression values.
#'
estimate.capture.cell.sizes <- function(expr.l) (
    expr.l
    %>% group_by(gene)
    %>% dplyr::summarise(prop.expr=mean(x>0))
    %>% filter(prop.expr > .5)
    %>% left_join(expr.l)
    %>% group_by(capture)
    %>% do(anders.huber.cell.sizes(.))
)


#' Estimate the cell sizes and adjust the expression by cell size.
#'
#' @param dl de.lorean object.
#'
#' @export
#'
adjust.by.cell.sizes <- function(dl) within(dl, {
    expr.l <- melt.expr(dl)
    cell.sizes <- estimate.capture.cell.sizes(
        expr.l
        %>% left_join(cell.meta %>% dplyr::select(cell, capture)))
    expr.before.adj <- expr
    expr <- cast.expr(
        expr.l
        %>% left_join(cell.sizes)
        %>% mutate(x=x-size))
})

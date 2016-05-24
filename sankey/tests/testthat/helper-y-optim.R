
y_optim_data <- function() {

  nodes <- data.frame(
    stringsAsFactors = FALSE,
    ID = paste0("n", 1:9),
    x = c(1,1,1,1, 2,2,2,2,2),
    y = c(1,2,3,4, 1,2,3,4,5)
  )

  edges <- data.frame(
    stringsAsFactors = FALSE,
    from = c("n1", "n1", "n2", "n3", "n3", "n4", "n4"),
    to   = c("n5", "n9", "n6", "n6", "n7", "n8", "n9")
  )

  list(nodes = nodes, edges = edges)
}

y_optim_data2 <- function() {

  nodes <- data.frame(
    stringsAsFactors = FALSE,
    ID = paste0("n", 1:9),
    x = c(1,1,1,1, 2,2,2,2,2)
  )

  edges <- data.frame(
    stringsAsFactors = FALSE,
    from = c("n1", "n1", "n2", "n3", "n3", "n4", "n4"),
    to   = c("n5", "n9", "n6", "n6", "n7", "n8", "n9")
  )

  nodes$size <- optimize_sizes(nodes, edges)

  list(nodes = nodes, edges = edges)

}


g1 <- function() {
  nodes <- data.frame(
    stringsAsFactors = FALSE,
    id = letters[1:5]
  )

  edges <- data.frame(
    stringsAsFactors = FALSE,
    from = c("a", "b", "c", "a", "b", "e"),
    to   = c("b", "a", "c", "e", "d", "a")
  )

  graph(nodes, edges)
}

g2 <- function() {
  nodes <- data.frame(
    stringsAsFactors = FALSE,
    id = letters[1:5]
  )

  edges <- data.frame(
    stringsAsFactors = FALSE,
    from = character(),
    to   = character()
  )

  graph(nodes, edges)
}

g3 <- function() {
  nodes <- data.frame(
    stringsAsFactors = FALSE,
    id = character()
  )

  edges <- data.frame(
    stringsAsFactors = FALSE,
    from = character(),
    to   = character()
  )

  graph(nodes, edges)
}

g4 <- function() {
  nodes <- data.frame(
    stringsAsFactors = FALSE,
    id = letters[1:5]
  )

  edges <- data.frame(
    stringsAsFactors = FALSE,
    from = c("a", "b", "a", "a", "b", "e"),
    to   = c("b", "a", "e", "e", "d", "a")
  )

  graph(nodes, edges)
}

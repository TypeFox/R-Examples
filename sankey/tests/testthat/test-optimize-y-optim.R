
context("Optimize y coords")

test_that("switch_nodes works", {

  nodes <- data.frame(
    stringsAsFactors = FALSE,
    id = c("foo", "bar", "foobar"),
    x = c(1,1,1),
    y = c(1,2,3)
  )

  nodes2_truth <- data.frame(
    stringsAsFactors = FALSE,
    id = c("foo", "bar", "foobar"),
    x = c(1,1,1),
    y = c(3,2,1)
  )

  nodes2 <- switch_nodes(nodes, x = 1, y1 = 1, y2 = 3)

  expect_equal(nodes2, nodes2_truth)
})


test_that("crossing_edge", {

  data <- y_optim_data()

  expect_true(crossing_edge(data$nodes, data$edges, 2, 3))
  expect_true(crossing_edge(data$nodes, data$edges, 2, 4))
  expect_true(crossing_edge(data$nodes, data$edges, 2, 5))
  expect_true(crossing_edge(data$nodes, data$edges, 2, 6))

  expect_false(crossing_edge(data$nodes, data$edges, 1, 2))
  expect_false(crossing_edge(data$nodes, data$edges, 1, 3))
  expect_false(crossing_edge(data$nodes, data$edges, 1, 4))
  expect_false(crossing_edge(data$nodes, data$edges, 1, 5))
  expect_false(crossing_edge(data$nodes, data$edges, 1, 6))
  expect_false(crossing_edge(data$nodes, data$edges, 1, 7))

  expect_false(crossing_edge(data$nodes, data$edges, 2, 1))
  expect_false(crossing_edge(data$nodes, data$edges, 2, 7))

  expect_false(crossing_edge(data$nodes, data$edges, 3, 1))
  expect_false(crossing_edge(data$nodes, data$edges, 3, 4))
  expect_false(crossing_edge(data$nodes, data$edges, 3, 5))
  expect_false(crossing_edge(data$nodes, data$edges, 3, 6))
  expect_false(crossing_edge(data$nodes, data$edges, 3, 7))

  expect_false(crossing_edge(data$nodes, data$edges, 4, 1))
  expect_false(crossing_edge(data$nodes, data$edges, 4, 3))
  expect_false(crossing_edge(data$nodes, data$edges, 4, 5))
  expect_false(crossing_edge(data$nodes, data$edges, 4, 6))
  expect_false(crossing_edge(data$nodes, data$edges, 4, 7))

  expect_false(crossing_edge(data$nodes, data$edges, 5, 1))
  expect_false(crossing_edge(data$nodes, data$edges, 5, 3))
  expect_false(crossing_edge(data$nodes, data$edges, 5, 4))
  expect_false(crossing_edge(data$nodes, data$edges, 5, 6))
  expect_false(crossing_edge(data$nodes, data$edges, 5, 7))

  expect_false(crossing_edge(data$nodes, data$edges, 6, 1))
  expect_false(crossing_edge(data$nodes, data$edges, 6, 3))
  expect_false(crossing_edge(data$nodes, data$edges, 6, 4))
  expect_false(crossing_edge(data$nodes, data$edges, 6, 5))
  expect_false(crossing_edge(data$nodes, data$edges, 6, 7))

  expect_false(crossing_edge(data$nodes, data$edges, 7, 1))
  expect_false(crossing_edge(data$nodes, data$edges, 7, 2))
  expect_false(crossing_edge(data$nodes, data$edges, 7, 3))
  expect_false(crossing_edge(data$nodes, data$edges, 7, 4))
  expect_false(crossing_edge(data$nodes, data$edges, 7, 5))
  expect_false(crossing_edge(data$nodes, data$edges, 7, 6))

})

test_that("crossing_edges", {

  data <- y_optim_data()

  expect_equal(crossing_edges(data$nodes, data$edges, c(), 3), 0)
  expect_equal(crossing_edges(data$nodes, data$edges, 2, c()), 0)
  expect_equal(crossing_edges(data$nodes, data$edges, c(), c()), 0)

  expect_equal(crossing_edges(data$nodes, data$edges, 2, 1), 0)
  expect_equal(crossing_edges(data$nodes, data$edges, 2, 3), 1)
  expect_equal(crossing_edges(data$nodes, data$edges, 2, 3:6), 4)

  expect_equal(crossing_edges(data$nodes, data$edges, c(1,2,7), 3:6), 4)
})

test_that("eval_node_pair", {

  data <- y_optim_data()

  expect_equal(eval_node_pair(data$nodes, data$edges, 1, 1, 2), 1)
  expect_equal(eval_node_pair(data$nodes, data$edges, 1, 1, 3), 2)
  expect_equal(eval_node_pair(data$nodes, data$edges, 1, 1, 4), 1)
  expect_equal(eval_node_pair(data$nodes, data$edges, 1, 2, 3), 0)
  expect_equal(eval_node_pair(data$nodes, data$edges, 1, 2, 4), 0)
  expect_equal(eval_node_pair(data$nodes, data$edges, 1, 3, 4), 0)

  expect_equal(eval_node_pair(data$nodes, data$edges, 2, 1, 2), 0)
  expect_equal(eval_node_pair(data$nodes, data$edges, 2, 1, 3), 0)
  expect_equal(eval_node_pair(data$nodes, data$edges, 2, 1, 4), 0)
  expect_equal(eval_node_pair(data$nodes, data$edges, 2, 1, 5), 0)
  expect_equal(eval_node_pair(data$nodes, data$edges, 2, 2, 3), 0)
  expect_equal(eval_node_pair(data$nodes, data$edges, 2, 2, 4), 0)
  expect_equal(eval_node_pair(data$nodes, data$edges, 2, 2, 5), 2)
  expect_equal(eval_node_pair(data$nodes, data$edges, 2, 3, 4), 0)
  expect_equal(eval_node_pair(data$nodes, data$edges, 2, 3, 5), 1)
  expect_equal(eval_node_pair(data$nodes, data$edges, 2, 4, 5), 1)
})

test_that("bring_up", {

  data <- y_optim_data()

  nodes1 <- bring_up(data$nodes, data$edges, 1, 4)
  nodes2 <- bring_up(data$nodes, data$edges, 1, 3)
  nodes3 <- bring_up(data$nodes, data$edges, 1, 2)

  nodesx <- data$nodes
  nodesx[1,]$y <- 2
  nodesx[2,]$y <- 1

  expect_equal(nodesx, nodes1)
  expect_equal(nodesx, nodes2)
  expect_equal(nodesx, nodes3)
})

test_that("bubble", {

  data <- y_optim_data()

  nodes1 <- bubble(data$nodes, data$edges, 1)
  nodesx1 <- data$nodes
  nodesx1$y <- c(3,1,2,4, 1:5)
  expect_equal(nodes1, nodesx1)

  nodes2 <- bubble(data$nodes, data$edges, 2)
  nodesx2 <- data$nodes
  nodesx2$y <- c(1:4, 1,2,3,5,4)
  expect_equal(nodes2, nodesx2)
})

test_that("optimize_y_optim", {

  data <- y_optim_data2()

  res <- optimize_y_optim(data$nodes, data$edges, gravity = "top")

})

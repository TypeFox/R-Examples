
context("BFS")

test_that("BFS", {

  G <- graph(list(
    "a" = "b",
    "b" = character(),
    "c" = "d",
    "d" = character(),
    "e" = "f",
    "f" = character()
  ))

  expect_equal(sort(bfs(G, "a")), c("a", "b"))
  expect_equal(sort(bfs(G, "c")), c("c", "d"))
  expect_equal(sort(bfs(G, c("a", "c"))), c("a", "b", "c", "d"))

  G <- graph(list(
    "7"  = c("11", "8"),
    "5"  = "11",
    "3"  = c("8", "10"),
    "11" = c("2", "9", "10"),
    "8"  = "9",
    "2"  = character(),
    "9"  = character(),
    "10" = character()
  ))

  expect_equal(bfs(G, character()), character())
  expect_equal(bfs(G, "10"), "10")
  expect_equal(bfs(G, "7"), c("7", "11", "8", "2", "9", "10"))
  expect_equal(
    sort(bfs(G, c("7", "5"))),
    c("10", "11", "2", "5", "7", "8", "9")
  )
})

test_that("BFS bug fixed", {
  bridges <- graph(list(
    "Altstadt-Loebenicht" = c(
      "Kneiphof",
      "Kneiphof",
      "Lomse"
    ),
    "Kneiphof" = c(
      "Altstadt-Loebenicht",
      "Altstadt-Loebenicht",
      "Vorstadt-Haberberg",
      "Vorstadt-Haberberg",
      "Lomse"
    ),
    "Vorstadt-Haberberg" = c(
      "Kneiphof",
      "Kneiphof",
      "Lomse"
    ),
    "Lomse" = c(
      "Altstadt-Loebenicht",
      "Kneiphof",
      "Vorstadt-Haberberg"
    )
  ))

  expect_equal(
    bfs(bridges),
    c("Altstadt-Loebenicht", "Kneiphof", "Lomse", "Vorstadt-Haberberg")
  )
})

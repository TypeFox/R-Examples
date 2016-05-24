##
## Ensure the number of evaluations conforms to the formula given in the
## documentation.
##

context("reg-n_evals")

## Trivial test function
f <- function(x) {
  c(sum(x), sum(x - 5))
}

eval_counter <- function(f) {
  force(f)
  N <- 0
  function(...) {
    N <<- N + 1
    f(...)
  }
}

for (p in c(4, 8, 96, 100, 120)) {
  for (g in c(1:10, 100, 400)) {
    tn <- sprintf("Number of evaluations for popsize=%i generations=%i", p, g)
    test_that(tn, {
      counting_f <- eval_counter(f)
      r <- nsga2(counting_f, idim=2, odim=2,
                 lower.bounds=c(-2, -2), upper.bounds=c(2, 2),
                 popsize=p, generations=g )
      expect_equal(environment(counting_f)$N, p * (g + 1))
    })
  }
}

for (p in c(4, 8, 96, 100, 120)) {
  for (g in c(1, 2, 5, 100)) {
    tn <- sprintf("Number of evaluations for popsize=%i generations=1:%i", p, g)
    test_that(tn, {
      counting_f <- eval_counter(f)
      generations <- 1:g
      r <- nsga2(counting_f, idim=2, odim=2,
                 lower.bounds=c(-2, -2), upper.bounds=c(2, 2),
                 popsize=p, generations=generations)
      expect_equal(environment(counting_f)$N, p * (max(generations) + 1))
    })
  }
}

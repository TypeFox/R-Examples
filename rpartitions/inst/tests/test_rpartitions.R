test_that("P() returns the number of partitions of Q with k or less parts", {
  D = list()
  expect_equal(P(D, 0, 0, F, F)[[2]], 1)
  expect_equal(P(D, 0, 123, F, F)[[2]], 1)
  expect_equal(P(D, 123, 123, F, F)[[2]], 2552338241)
  expect_equal(P(D, 123, 0, F, F)[[2]], 0)
  expect_equal(P(D, 12345, 123, T, F)[[2]], 488259162924433580696194373878466788895554319556195978121822180221785381227453675217103501271281020550492)
})  

test_that("NrParts() returns the number of partitions for a given total (q) and number of parts (n)", {
  expect_equal(NrParts(0), 1)
  expect_equal(NrParts(0, 1), 1)
  expect_equal(NrParts(100, use_c=F), 190569292)
  expect_equal(NrParts(100, use_c=T), 190569292)
  expect_equal(NrParts(100, 10, use_c=F), 2977866)
  expect_equal(NrParts(100, 10, use_c=T), 2977866)
})

test_that("conjugate() returns a conjugated vector of integers", {  
  expect_equal(conjugate(1), 1)
  expect_equal(conjugate(10), rep(1,10))
  expect_equal(conjugate(3:1), 3:1)
  expect_equal(conjugate(c(12,9,8,8,7,5,4,2,2,1)), 
                         c(10,9,7,7,6,5,5,4,2,1,1,1))
})

test_that("divide_and_conquer can discover the entire feasible set.", {
  q = 20
  n = 5
  answer = 84 # there 84 partitions of 20 having 5 parts
  sample_size = 1000
  methods = c('divide_and_conquer','multiplicity','top_down','bottom_up')
  D = list()
  for (method in methods) {
    partitions = rand_partitions(q, n, sample_size, method, D, zeros=F)
    feasibleset = unique(apply(partitions, 2, list))
    expect_equal(84, length(feasibleset))
  }  
})

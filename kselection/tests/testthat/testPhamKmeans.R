##
##  kselection tests with other kmeans function
##
##  Created by Daniel Rodriguez Perez on 8/1/2015.
##
##  Copyright (c) 2015 Daniel Rodriguez Perez.
##
##  This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <http://www.gnu.org/licenses/>
## 

context("Tests for kselection with other kmeans functions")

test_that("evaluate with amap", {
  skip_on_cran()
  
  if (!requireNamespace('amap')) {
    skip('No amap package')
  }
  
  set.seed(1000)
  x <- matrix(c(rnorm(100, 2, .1), rnorm(100, 3, .1),
                rnorm(100, -2, .1), rnorm(100, 1, .1),
                rnorm(100, 1, .1), rnorm(100, -3, .1),
                rnorm(100, -1, .1), rnorm(100, -2, .1)), 400, 2)
  k <- kselection(x, fun_cluster = amap::Kmeans, nstart = 10)
  
  expect_that(num_clusters(x), is_null())
  expect_that(num_clusters_all(x), is_null())
  
  expect_that(class(k), equals('Kselection'))
  expect_that(k$k, equals(4))
  expect_that(num_clusters(k), equals(4))
  
  valid_clusters <- which(get_f_k(k) < k$k_threshold)
  expect_that(num_clusters_all(k), equals(valid_clusters))
  
  valid_clusters <- which(get_f_k(k) < 1)
  k$k_threshold  <- 1
  expect_that(num_clusters_all(k), equals(valid_clusters))
  
  valid_clusters <- which(get_f_k(k) < 0.1)
  k$k_threshold  <- 0.1
  expect_that(num_clusters_all(k), equals(valid_clusters))
})

test_that("evaluate with FactoClass", {
  skip_on_cran()
  
  if (!requireNamespace('FactoClass')) {
    skip('No FactoClass package')
  }
  
  set.seed(1000)
  x <- matrix(c(rnorm(100, 2, .1), rnorm(100, 3, .1),
                rnorm(100, -2, .1), rnorm(100, 1, .1),
                rnorm(100, 1, .1), rnorm(100, -3, .1),
                rnorm(100, -1, .1), rnorm(100, -2, .1)), 400, 2)
  k <- kselection(x, fun_cluster = FactoClass::kmeansW, nstart = 10)
  
  expect_that(num_clusters(x), is_null())
  expect_that(num_clusters_all(x), is_null())
  
  expect_that(class(k), equals('Kselection'))
  expect_that(k$k, equals(4))
  expect_that(num_clusters(k), equals(4))
  
  valid_clusters <- which(get_f_k(k) < k$k_threshold)
  expect_that(num_clusters_all(k), equals(valid_clusters))
  
  valid_clusters <- which(get_f_k(k) < 1)
  k$k_threshold  <- 1
  expect_that(num_clusters_all(k), equals(valid_clusters))
  
  valid_clusters <- which(get_f_k(k) < 0.1)
  k$k_threshold  <- 0.1
  expect_that(num_clusters_all(k), equals(valid_clusters))
})

test_that("evaluate with LICORS", {
  skip_on_cran()
  
  if (!requireNamespace('LICORS')) {
    skip('No LICORS package')
  }
  
  set.seed(1000)
  x <- matrix(c(rnorm(100, 2, .1), rnorm(100, 3, .1),
                rnorm(100, -2, .1), rnorm(100, 1, .1),
                rnorm(100, 1, .1), rnorm(100, -3, .1),
                rnorm(100, -1, .1), rnorm(100, -2, .1)), 400, 2)
  k <- kselection(x, fun_cluster = LICORS::kmeanspp, nstart = 10)
  
  expect_that(num_clusters(x), is_null())
  expect_that(num_clusters_all(x), is_null())
  
  expect_that(class(k), equals('Kselection'))
  expect_that(k$k, equals(4))
  expect_that(num_clusters(k), equals(4))
  
  valid_clusters <- which(get_f_k(k) < k$k_threshold)
  expect_that(num_clusters_all(k), equals(valid_clusters))
  
  valid_clusters <- which(get_f_k(k) < 1)
  k$k_threshold  <- 1
  expect_that(num_clusters_all(k), equals(valid_clusters))
  
  valid_clusters <- which(get_f_k(k) < 0.1)
  k$k_threshold  <- 0.1
  expect_that(num_clusters_all(k), equals(valid_clusters))
})

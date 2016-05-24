test_that("SRD returns correct results",
          {
            cov.matrix.1 <- cov(matrix(rnorm(30*10), 30, 10))
            cov.matrix.2 <- cov(matrix(rnorm(30*10), 30, 10))
            srd.output <- SRD(cov.matrix.1, cov.matrix.2)
            expect_that(srd.output, is_a("SRD"))
            expect_that(names(srd.output), equals(c("output", "pc1", "model", "cormat")))
            expect_that(names(srd.output[[3]]), equals(c("quantiles", "interval", "code")))
            expect_that(dim(srd.output[[1]]), equals(c(10, 6)))
            expect_that(length(srd.output[[2]]), equals(10))
            expect_that(length(srd.output[[3]]), equals(3))
            expect_that(dim(srd.output[[4]]), equals(c(10, 10)))
            expect_that(srd.output[[1]], is_a("matrix"))
            expect_that(srd.output[[2]], is_a("numeric"))
            expect_that(srd.output[[3]], is_a("list"))
            expect_that(srd.output[[4]], is_a("matrix"))
            srd.scores.bool <- sapply(srd.output[[1]][,1], function(x) isTRUE(x > -1 & x < 1))
            srd.lower.bool <- sapply(srd.output[[1]][,2], function(x) isTRUE(x > -1 & x < 1))
            srd.upper.bool <- sapply(srd.output[[1]][,3], function(x) isTRUE(x > -1 & x < 1))
            srd.sd.bool <- sapply(srd.output[[1]][,4], function(x) isTRUE(x > 0 & x < 1))
            expect_that(sum(srd.scores.bool), equals(length(srd.scores.bool)))
            expect_that(sum(srd.lower.bool), equals(length(srd.scores.bool)))
            expect_that(sum(srd.upper.bool), equals(length(srd.scores.bool)))
            expect_that(mean(srd.output[[1]][,5]), equals(0))
            expect_that(mean(srd.output[[1]][,6]), equals(0))
            expect_that(mean(srd.output[[2]]), equals(0))
            expect_that(sum(diag(srd.output[[4]])), equals(10))
            expect_that(as.numeric(srd.output[[2]] < srd.output[[3]][[2]]), equals(srd.output[[3]][[3]]))
          }
)

test_that("SRD returns correct results on lists",
          {
            m.list <- RandomMatrix(10, 4)
            set.seed(42)
            srd.all <- SRD(m.list)
            set.seed(42)
            one.and.two  <- srd.all[1, 2][[1]]
            one.and.two.default <- SRD(m.list[[1]], m.list[[2]])
            expect_that(one.and.two, is_equivalent_to(one.and.two.default))
            expect_that(srd.all, is_a('matrix'))
            expect_that(isSymmetric(srd.all),is_true())
          }
)

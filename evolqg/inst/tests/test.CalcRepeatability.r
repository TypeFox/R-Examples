test_that("CalcRepeatability returns correct values",
          {
            num.ind = length(iris[,1])
            ID = as.factor(rep(1:num.ind, 2))
            ind.data = rbind(iris[,1:4], iris[,1:4]+array(rnorm(num.ind*4, 0, 0.1), dim(iris[,1:4])))
            models.list <- apply (ind.data, 2, function (vec){return (lm (vec ~ ID))})
            models.list <- lapply (models.list, anova)
            rep.itself <- function (lm.model){
                msq <- lm.model$'Mean Sq' ## 1 entre, 2 dentro
                s2a <- (msq[1] - msq[2])/2
                out <- s2a / (s2a + msq[2])
                return (out)
            }
            out <- sapply (models.list, rep.itself)
            expect_that(CalcRepeatability(ID, ind.data),equals(out))
          }
)

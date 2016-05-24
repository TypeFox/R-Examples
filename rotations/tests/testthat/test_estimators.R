rs<-rvmises(20)
Rs<-genR(rs)
Qs<-as.Q4(Rs)


context("Estimators")

#Do the estimator functions return correct objects?
expect_that(mean(Rs),is_a("SO3"))
expect_that(mean(Rs,type='geometric'),is_a("SO3"))
expect_that(median(Rs),is_a("SO3"))
expect_that(median(Rs,type='geometric'),is_a("SO3"))

expect_that(mean(Qs),is_a("Q4"))
expect_that(mean(Qs,type='geometric'),is_a("Q4"))
expect_that(median(Qs),is_a("Q4"))
expect_that(median(Qs,type='geometric'),is_a("Q4"))


#Estimator should be the same regardless of parametrization
expect_equal(mean(Qs),as.Q4(mean(Rs)))
expect_equal(mean(Qs,type='geometric'),as.Q4(mean(Rs,type='geometric')))
expect_equal(median(Qs),as.Q4(median(Rs)))
expect_equal(median(Qs,type='geometric'),as.Q4(median(Rs,type='geometric')))

expect_equal(as.SO3(matrix(as.SO3(mean(Qs)),3,3)),mean(Rs))
expect_equal(as.SO3(matrix(as.SO3(mean(Qs,type='geometric')),3,3)),mean(Rs,type='geometric'))
expect_equal(as.SO3(matrix(as.SO3(median(Qs)),3,3)),median(Rs))
expect_equal(as.SO3(matrix(as.SO3(median(Qs,type='geometric')),3,3)),median(Rs,type='geometric'))


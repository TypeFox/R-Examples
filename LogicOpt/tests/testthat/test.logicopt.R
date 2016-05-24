test_that("Logicopt basic functionality", {
   data(l.small)
   esp1 <- logicopt(l.small,4,3)
   expect_equal(length(esp1), 2)
   expect_equal(length(esp1[[2]]),1)
   expect_equal(ncol(esp1[[1]]),7)
   expect_equal(nrow(esp1[[1]]),esp1[[2]])

   esp2 <- logicopt(l.small,4,3,find_dc=TRUE)
   expect_equal(length(esp2), 2)
   expect_equal(length(esp2[[2]]),1)
   expect_equal(ncol(esp2[[1]]),7)
   expect_equal(nrow(esp2[[1]]),esp2[[2]])

   esp3 <- logicopt(l.small,4,3,find_dc=TRUE, mode="qm", exact_cover=TRUE)
   expect_equal(length(esp3), 2)
   expect_equal(length(esp3[[2]]),1)
   expect_equal(ncol(esp3[[1]]),7)
   expect_equal(nrow(esp3[[1]]),esp3[[2]])

   esp4 <- logicopt(l.small,4,3,find_dc=TRUE, mode="primes", exact_cover=TRUE)
   expect_equal(length(esp4), 2)
   expect_equal(length(esp4[[2]]),3)
   expect_equal(ncol(esp4[[1]]),7)
   expect_equal(nrow(esp4[[1]]),sum(esp4[[2]]))

   esp5 <- logicopt(l.small,4,3,find_dc=TRUE, mode="multi-min", exact_cover=TRUE)
   expect_equal(length(esp5), 2)
   expect_equal(ncol(esp5[[1]]),7)
   expect_equal(nrow(esp5[[1]]),sum(esp5[[2]]))

   esp6 <- logicopt(l.small,4,3,find_dc=TRUE, mode="multi-full", exact_cover=TRUE)
   expect_equal(length(esp6), 2)
   expect_equal(ncol(esp6[[1]]),7)
   expect_equal(nrow(esp6[[1]]),sum(esp6[[2]]))

   # sanity check interactive runs
   esp7 <- logicopt(l.small,4,3)
   expect_equal(esp1,esp7)
})

test_that("Logicopt MV tests", {
   data(l.partybans.1)
   values = num_input_values(l.partybans.1,5)
   mv1 <- logicopt(l.partybans.1,5,1,find_dc=TRUE,mode="qm",exact_cover=TRUE)
   expect_equal(length(mv1), 2)
   expect_equal(length(mv1[[2]]),1)
   expect_equal(ncol(mv1[[1]]),6)
   expect_equal(nrow(mv1[[1]]),mv1[[2]])

   mv2 <- logicopt(l.partybans.1,5,1,find_dc=TRUE,mode="qm",
     exact_cover=TRUE, input_sizes=values)
   expect_equal(mv1,mv2)

   # round trip!
   mv3 <- logicopt(mv2[[1]],5,1,find_dc=FALSE,mode="qm",
     exact_cover=TRUE,input_sizes=values)
   expect_equal(mv2[[2]],mv3[[2]])
})


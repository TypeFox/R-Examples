context("wrong_fc tests")

test_that("Regular fiscal codes", {
            ## right
            expect_false(wrong_fc("QWEASD34D12H221M"))
            ## wrong
            expect_true(wrong_fc("qWeASd34D12h 221M   "))
            expect_true(wrong_fc("QWEASD34D12H221X"))
            expect_true(wrong_fc("QWEASD34D12H221m"))
          })

test_that("Temporary fiscal codes", {
            ## right
            expect_false(wrong_fc("12312312312"))
            ## wrong
            expect_true(wrong_fc(" 12312312312 "))
            expect_true(wrong_fc("12312312315"))
          })

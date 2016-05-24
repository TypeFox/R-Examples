#############################################################################
# test-out_all.R
# 
# Testing input/output functions
#############################################################################

test_that("push and pull",
{
    ppo <- pubprint()

    push(ppo) <- t.test(1:10)
    push(ppo) <- t.test(1:100, 2:101)
    push(ppo, add = TRUE) <- 12.758690765

    push(ppo, item = "i1") <- t.test(1:30)
    push(ppo, item = "i2") <- t.test(1:30, 2:31)

    expect_identical(pull(ppo),
                     "(\\ensuremath{M\\ifmmode_{x}\\else\\textsubscript{x}\\fi=5.50, t[9]=5.74, p\\ifmmode<\\else\\textless\\fi.001})")
    expect_identical(pull(ppo),
                     "(\\ensuremath{M\\ifmmode_{x}\\else\\textsubscript{x}\\fi=50.50, M\\ifmmode_{y}\\else\\textsubscript{y}\\fi=51.50, t[198]=-0.24, p=.808, d=12.76})")


    expect_identical(pull(ppo, item = "i1"),
                     "(\\ensuremath{M\\ifmmode_{x}\\else\\textsubscript{x}\\fi=15.50, t[29]=9.64, p\\ifmmode<\\else\\textless\\fi.001})")
    expect_identical(pull(ppo, item = "i2"),
                     "(\\ensuremath{M\\ifmmode_{x}\\else\\textsubscript{x}\\fi=15.50, M\\ifmmode_{y}\\else\\textsubscript{y}\\fi=16.50, t[58]=-0.44, p=.662})")
    expect_identical(pull(ppo, item = "i2", mmode = FALSE),
                     "(M\\ifmmode_{x}\\else\\textsubscript{x}\\fi=15.50, M\\ifmmode_{y}\\else\\textsubscript{y}\\fi=16.50, t[58]=-0.44, p=.662)")
    expect_identical(pull(ppo, item = "i2", mmode = FALSE, concat = FALSE, separator = NULL),
                     c("M\\ifmmode_{x}\\else\\textsubscript{x}\\fi=15.50", 
                       "M\\ifmmode_{y}\\else\\textsubscript{y}\\fi=16.50", 
                       "t(58)=-0.44", 
                       "p=.662"))
})

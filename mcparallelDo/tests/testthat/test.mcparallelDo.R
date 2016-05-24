context("explicitly retrieved tests")
test_that("basic calls work",{
    mcparallelDo({2+2}, targetValue = "output4")
    mcparallelDo({4+4}, targetValue = "output8")
    mcparallelDo({Sys.sleep(5);8+8}, targetValue = "outputSleepy16")
    Sys.sleep(1)
    if (.Platform$OS.type == "unix") {
      #we only expect this code to do anything undex unix-alikes
      expect_equal(sum(mcparallelDoCheck()),2)
    }
    expect_equal(output4,4)
    expect_equal(output8,8)
    Sys.sleep(5)
    if (.Platform$OS.type == "unix") {
      #we only expect this code to do anything undex unix-alikes
      expect_equal(sum(mcparallelDoCheck()),1)
    }
    expect_equal(outputSleepy16,16)
})
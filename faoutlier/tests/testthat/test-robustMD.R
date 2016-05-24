context('robustMD')

test_that('robustMD run', {    
    set.seed(1)
    output <- robustMD(holzinger)
    expect_equal(as.numeric(output$mah[1:3]), c(6.679447, 4.261030, 15.056849),
                 tolerance = 1e-5)
    expect_is(output, 'robmah')
    expect_is(plot(output), 'trellis')
    expect_is(plot(output, type = 'qqplot'), 'trellis')
})
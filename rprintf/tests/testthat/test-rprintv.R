context("rprintv: variable-based formatting")

test_that("rprintv works as expected", {
    expect_equivalent(rprintv("$name:s", name = "a"), "a")
    expect_equivalent(rprintv("$num:d", num = 10), "10")
    expect_equivalent(rprintv("$pi:.2f", pi = pi), "3.14")
    expect_equivalent(rprintv("$pi:+.2f", pi = pi), "+3.14")
    expect_equivalent(rprintv("$pi: .2f", pi = pi), " 3.14")
    expect_equivalent(rprintv("$name1:s,$name", name = "a", name1 = "b"), "b,a")
    expect_equivalent(rprintv("$a,$b:s,$c:.2f,$b,$a", a = "a", b = "b", c = 1.4), "a,b,1.40,b,a")
    expect_equivalent(rprintv("$$name", name = "a"), "$name")
    
    p <- list(name = "Ken", age = 25)
    expect_equivalent(rprintn("{1} {2}", p), "Ken 25")
}) 

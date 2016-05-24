n.reps = 9


###############################################################
# Effective sample size functions
###############################################################

context("Sample size functions")


test_that("ess.weights throws proper errors", {
    w = numeric(0)
    expect_error(ess.weights(w    ))
    expect_error(ess.weights(w,"R"))
    expect_error(ess.weights(w,"C"))

    w = runif(4); w[2] = -w[2]
    expect_warning(ess.weights(w    ))
    expect_warning(ess.weights(w,"R"))
    expect_warning(ess.weights(w,"C"))
})

test_that("ess.weights works properly", {
    n = 4
    w = rep(1/n, n)
    expect_equal(ess.weights(w    ), n) 
    expect_equal(ess.weights(w,"R"), n) 
    expect_equal(ess.weights(w,"C"), n) 
})

test_that("ess.weights R matches C", {
    for (i in 1:n.reps) {
        w = runif(rpois(1,10)+1); w=w/sum(w)
        expect_equal(ess.weights(w,engine="R"),
                     ess.weights(w,engine="C"))
    }

})


test_that("ent.weights throws proper errors", {
    w = numeric(0)
    expect_error(ent.weights(w    ))
    expect_error(ent.weights(w,"R"))
    expect_error(ent.weights(w,"C"))

    w = runif(4); w[2] = -w[2]
    expect_warning(ent.weights(w    ))
    expect_warning(ent.weights(w,"R"))
    expect_warning(ent.weights(w,"C"))
})

test_that("ent.weights works properly", {
    n = 4
    w = rep(1/n, n)
    ent = -log2(1/n)
    expect_equal(ent.weights(w    ), ent) 
    expect_equal(ent.weights(w,"R"), ent) 
    expect_equal(ent.weights(w,"C"), ent) 
})

test_that("ent.weights R matches C", {
    for (i in 1:n.reps) {
        w = runif(rpois(1,10)+1); w=w/sum(w)
        expect_equal(ent.weights(w,engine="R"),
                     ent.weights(w,engine="C"))
    }

})


test_that("cov.weights throws proper errors", {
    w = numeric(0)
    expect_error(cov.weights(w    ))
    expect_error(cov.weights(w,"R"))
    expect_error(cov.weights(w,"C"))

    w = runif(4); w[2] = -w[2]
    expect_warning(cov.weights(w    ))
    expect_warning(cov.weights(w,"R"))
    expect_warning(cov.weights(w,"C"))
})

test_that("cov.weights works properly", {
    n = 4
    w = rep(1/n, n)
    answer = 0
    expect_equal(cov.weights(w    ), answer) 
    expect_equal(cov.weights(w,"R"), answer) 
    expect_equal(cov.weights(w,"C"), answer) 
})

test_that("cov.weights R matches C", {
    for (i in 1:n.reps) {
        w = runif(rpois(1,10)+1); w=w/sum(w)
        expect_equal(cov.weights(w,engine="R"),
                     cov.weights(w,engine="C"))
    }

})


###############################################################
# Resampling functions
###############################################################


context("Resampling functions")


test_that("stratified resampling throws errors", {
    w = runif(4); lw = log(w)
    expect_error(stratified.resample(lw))
    expect_error(stratified.resample( w,-1))
})

test_that("stratified resampling gets base case correct", {
    set.seed(1)
    w = runif(4); w=w/sum(w)
    n = 5
    id = c(1,3,3,4,4)
    set.seed(2); expect_equal(stratified.resample(w,n,   ), id)
    set.seed(2); expect_equal(stratified.resample(w,n,"R"), id)
    set.seed(2); expect_equal(stratified.resample(w,n,"C"), id)
})


test_that("stratified resampling: R matches C", {
    for (i in 1:n.reps) 
    {
        w = runif(rpois(1,10)+2); w=w/sum(w)
        n = rpois(1,10)+1
        seed = proc.time()
        set.seed(seed); mR = stratified.resample(w,n,"R")
        set.seed(seed); mC = stratified.resample(w,n,"C")
        expect_equal(mR,mC)
    }
})



########################## Multinomial ###############################
test_that("multinomial resampling throws errors", {
    w = runif(4); lw = log(w)
    expect_error(multinomial.resample(lw))
    expect_error(multinomial.resample( w,-1))
})

test_that("multinomial resampling gets base case correct", {
    set.seed(1)
    w = runif(4); w=w/sum(w)
    n = 5
    id = c(2,2,4,4,4)
    set.seed(2); expect_equal(multinomial.resample(w,n    ), id)
    set.seed(2); expect_equal(multinomial.resample(w,n,"R"), id)
    set.seed(2); expect_equal(multinomial.resample(w,n,"C"), id)
})


test_that("multinomial resampling: R matches C", {
    for (i in 1:n.reps) 
    {
        w = runif(rpois(1,10)+2); w=w/sum(w)
        n = rpois(1,10)+1
        seed = proc.time()
        set.seed(seed); mR = multinomial.resample(w,n,"R")
        set.seed(seed); mC = multinomial.resample(w,n,"C")
        expect_equal(mR,mC)
    }
})



test_that("systematic resampling throws errors", {
    w = runif(4); lw = log(w)
    expect_error(systematic.resample(lw))
    expect_error(systematic.resample( w,-1))
})


test_that("systematic resampling gets base case correct", {
    set.seed(1)
    w = runif(4); w=w/sum(w)
    n = 5
    id = c(1,2,3,4,4)
    set.seed(2); expect_equal(systematic.resample(w,n    ), id)
    set.seed(2); expect_equal(systematic.resample(w,n,"R"), id)
    set.seed(2); expect_equal(systematic.resample(w,n,"C"), id)
})


test_that("systematic resampling: R matches C", {
    for (i in 1:n.reps) 
    {
        w = runif(rpois(1,10)+2); w=w/sum(w)
        n = rpois(1,10)+1
        seed = proc.time()
        set.seed(seed); mR = systematic.resample(w,n,"R")
        set.seed(seed); mC = systematic.resample(w,n,"C")
        expect_equal(mR,mC)
    }
})



test_that("residual resampling throws errors", {
    w = runif(4); lw = log(w)
    expect_error(residual.resample(lw))
    expect_error(residual.resample( w,-1))
})

test_that("residual resampling gets base case correct", {
    set.seed(1)
    w = runif(4); w=w/sum(w)
    n = 5
    id = c(3,4,4,1,3)
    set.seed(2); expect_equal(residual.resample(w,n           ), id)
    set.seed(2); expect_equal(residual.resample(w,n,engine="R"), id)
    set.seed(2); expect_equal(residual.resample(w,n,engine="C"), id)
})

test_that("residual resampling: R matches C", {
    for (i in 1:n.reps) 
    {
        w = runif(rpois(1,10)+2); w=w/sum(w)
        n = rpois(1,10)+1
        seed = proc.time()
        set.seed(seed); mR = residual.resample(w,n,engine="R")
        set.seed(seed); mC = residual.resample(w,n,engine="C")
        expect_equal(mR,mC)
    }
})


context("Resample chooses proper method")

n=9
test_that("resample chooses stratified resampling", {
    w = runif(n); w=w/sum(w)
    seed = proc.time()
    set.seed(seed); m1 = resample(log(w), method="stratified")$indices
    set.seed(seed); m2 = stratified.resample(w,engine="R")
    expect_equal(m1,m2)
})

test_that("resample chooses multinomial resampling", {
    w = runif(n); w=w/sum(w)
    seed = proc.time()
    set.seed(seed); m1 = resample(log(w), method="multinomial")$indices
    set.seed(seed); m2 = multinomial.resample(w,engine="R")
    expect_equal(m1,m2)
})

test_that("resample chooses systematic resampling", {
    w = runif(n); w=w/sum(w)
    seed = proc.time()
    set.seed(seed); m1 = resample(log(w), method="systematic")$indices
    set.seed(seed); m2 = systematic.resample(w,engine="R")
    expect_equal(m1,m2)
})

test_that("resample chooses residual resampling", {
    w = runif(n); w=w/sum(w)
    seed = proc.time()
    set.seed(seed); m1 = resample(log(w), method="residual")$indices
    set.seed(seed); m2 = residual.resample(w,engine="R")
    expect_equal(m1,m2)
})


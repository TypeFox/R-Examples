# 
# Test some common-use cases of rollply
# 

library(plyr)
context('Common uses of rollply')



# Generate data
dat <- data.frame(time=seq.int(1000),
                  position=cumsum(rnorm(1000,0,10)))
rollav <- rollply(dat, ~ time, wdw.size = 10, 
                  summarise, position = mean(position))

test_that("1d moving averages works", { 
    expect_is(rollav, "data.frame")
})



# Generate three 2D random walks
dat <- ddply(data.frame(person=c('françois','nicolas','jacques')), ~ person, 
             summarise, 
             time=seq.int(1000),
             x=cumsum(rnorm(1000,0,1)),
             y=cumsum(rnorm(1000,0,1)))

# Smoothed trajectory over ten time-steps
rollav <- rollply(dat, ~ time | person, wdw.size=10, mesh.res=1000,
                  summarise, x=mean(x), y=mean(y))

test_that("2d moving averages with groups", { 
    expect_is(rollav, "data.frame")
})


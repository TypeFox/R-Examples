library(RUnit)
library(ROCR)

testFarg <- function() {
    ll <- list(arg1=c(1,2,3), arg2=c(4,5,6))
    print(str(.farg(ll, arg3=c(7,8,9)) ))
    checkEquals(.farg(ll, arg3=c(7,8,9)), list(arg1=c(1,2,3), arg2=c(4,5,6), arg3=c(7,8,9)))
    checkEquals(.farg(ll, arg1=c(1,4,3)), list(arg1=c(1,2,3), arg2=c(4,5,6)))
}

testGarg <- function() {
    ll <- list(arg1=list(1,2,3), arg2=list(4,5,6))
    checkEquals(.garg(ll, 'arg1'), 1)
    checkEquals(.garg(ll, 'arg1',2), 2)
    checkEquals(.garg(ll, 'arg2',3), 6)

    checkEquals(.garg(ll, 'arg3'), ll$arg3)
}

testSlice <- function() {
    ll <- list(arg1=list(c(1,2,3), c(2,3,4), c(3,4,5)),
               arg2=list('a', 'b', 'c'))
    checkEquals(.slice.run(ll, 1), list(arg1=c(1,2,3), arg2='a'))
    checkEquals(.slice.run(ll, 2), list(arg1=c(2,3,4), arg2='b'))
    checkEquals(.slice.run(ll, 3), list(arg1=c(3,4,5), arg2='c'))

    ll <- list(arg1=list(c(1,2,3), c(2,3,4), c(3,4,5)),
               arg2=c('a', 'b', 'c'))
    checkEquals(.slice.run(ll, 1), list(arg1=c(1,2,3), arg2=c('a', 'b', 'c')))
    checkEquals(.slice.run(ll, 2), list(arg1=c(2,3,4), arg2=c('a', 'b', 'c')))
    checkEquals(.slice.run(ll, 3), list(arg1=c(3,4,5), arg2=c('a', 'b', 'c')))
}

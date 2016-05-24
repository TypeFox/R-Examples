library(foreach)
m <- matrix(1:10/10, ncol=2, byrow=TRUE)
colnames(m) <- c('a', 'b')
d <- as.data.frame(m)
#       a    b
# [1,]  0.1  0.2
# [2,]  0.3  0.4
# [3,]  0.5  0.6
# [4,]  0.7  0.8
# [5,]  0.9  1.0
                 



test_that('lcut3 on matrix', {
    testedHedges <- c("ex", "si", "ve", "ml", "ro", "qr", "vr")
    hedgeNames   <- c("Ex", "Si", "Ve", "Ml", "Ro", "Qr", "Vr")

    smHedgeNames <- c("Ex", "Si", "Ve", "", "Ml", "Ro", "Qr", "Vr")
    meHedgeNames <- c("", "Ml", "Ro", "Qr", "Vr")
    biHedgeNames <- smHedgeNames

    attrs <- c(paste(smHedgeNames, 'Sm', sep=''),
               paste(meHedgeNames, 'Me', sep=''),
               paste(biHedgeNames, 'Bi', sep=''))

    res <- lcut3(m,
                context=c(0, 0.5, 1),
                hedges=testedHedges)
    expectedAttrs <- c(paste(attrs, '.a', sep=''),
                       paste(attrs, '.b', sep=''))

    expect_true(is.matrix(res))
    expect_equal(ncol(res), 42)
    expect_equal(nrow(res), 5)
    expect_equal(sort(colnames(res)), sort(expectedAttrs))
    expect_true(inherits(res, 'fsets'))
    expect_equivalent(vars(res), c(rep('a', 21), rep('b', 21)))
    expect_equal(sort(names(vars(res))), sort(expectedAttrs))
    expect_equal(sort(colnames(specs(res))), sort(expectedAttrs))
    expect_equal(sort(rownames(specs(res))), sort(expectedAttrs))

    s <- matrix(c(0,1,1,1,1,1,1,1, 0,0,0,0,0, 0,0,0,0,0,0,0,0,  #ExSm
                  0,0,1,1,1,1,1,1, 0,0,0,0,0, 0,0,0,0,0,0,0,0,  #SiSm
                  0,0,0,1,1,1,1,1, 0,0,0,0,0, 0,0,0,0,0,0,0,0,  #VeSm
                  0,0,0,0,1,1,1,1, 0,0,0,0,0, 0,0,0,0,0,0,0,0,  #Sm
                  0,0,0,0,0,1,1,1, 0,0,0,0,0, 0,0,0,0,0,0,0,0,  #MlSm
                  0,0,0,0,0,0,1,1, 0,0,0,0,0, 0,0,0,0,0,0,0,0,  #RoSm
                  0,0,0,0,0,0,0,1, 0,0,0,0,0, 0,0,0,0,0,0,0,0,  #QrSm
                  0,0,0,0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,0,0,0,  #VrSm

                  0,0,0,0,0,0,0,0, 0,1,1,1,1, 0,0,0,0,0,0,0,0,  #Me
                  0,0,0,0,0,0,0,0, 0,0,1,1,1, 0,0,0,0,0,0,0,0,  #MlMe
                  0,0,0,0,0,0,0,0, 0,0,0,1,1, 0,0,0,0,0,0,0,0,  #RoMe
                  0,0,0,0,0,0,0,0, 0,0,0,0,1, 0,0,0,0,0,0,0,0,  #QrMe
                  0,0,0,0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,0,0,0,  #VrMe

                  0,0,0,0,0,0,0,0, 0,0,0,0,0, 0,1,1,1,1,1,1,1,  #ExBi
                  0,0,0,0,0,0,0,0, 0,0,0,0,0, 0,0,1,1,1,1,1,1,  #SiBi
                  0,0,0,0,0,0,0,0, 0,0,0,0,0, 0,0,0,1,1,1,1,1,  #VeBi
                  0,0,0,0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,1,1,1,1,  #Bi
                  0,0,0,0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,1,1,1,  #MlBi
                  0,0,0,0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,0,1,1,  #RoBi
                  0,0,0,0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,0,0,1,  #QrBi
                  0,0,0,0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,0,0,0), #VrBi
                nrow=21,
                ncol=21)
    sfill <- matrix(0, nrow=21, ncol=21)
    s <- cbind(rbind(s, sfill), rbind(sfill, s))
    colnames(s) <- expectedAttrs
    rownames(s) <- expectedAttrs
    s <- s[colnames(res), colnames(res)]
    expect_equal(specs(res), s)
    expect_true(is.fsets(res))
})

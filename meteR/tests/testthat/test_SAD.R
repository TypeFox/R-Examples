context('test SAD')
## values from Newman et al.

valid <- data.frame(species = c('TAROFF', 'UNKSP', 'AQUCOE', 'OSMOCC', 'VIGMUL', 'CASOCC',
                                'GERRIC', 'SENAMP', 'ERISPE', 'EUCENG', 'DRASPE', 'AGOGLA', 'ERIELA',
                                'CASSUL', 'SOLMUL', 'LATLEU', 'ERIFOR', 'DELBAR', 'POLDOU', 'SENCRA',
                                'NOCMON', 'IPOAGG', 'VICAME', 'LIGPOR', 'POTGRA', 'FRAVES', 'HYMHOO',
                                'ANDSEP', 'LUPARG', 'HELQUI', 'VIONUT', 'BOEDRU'), 
                    abund = c(1, 1, 2, 2, 2, 6, 7, 7, 8, 8, 10, 11, 11, 12, 12, 13, 15, 17, 17, 19, 
                              20, 23, 33, 40, 47, 55, 55, 80, 86, 95, 101, 104), 
                    newmanPred = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 5, 6, 6, 8, 9, 11, 13, 15, 
                                   18, 21, 25, 31, 37, 46, 58, 74, 97, 138, 240),
                    newmanPredCDF = c(0.1998572, 0.1998572, 0.299091141, 0.299091141,
                                      0.299091141, 0.484756028, 0.512136717, 0.512136717, 0.535928274,
                                      0.535928274, 0.575698857, 0.592643463, 0.592643463, 0.608068043,
                                      0.608068043, 0.62220714, 0.647329142, 0.669093965, 0.669093965,
                                      0.68823307, 0.696985488, 0.720562864, 0.77917697, 0.808702375,
                                      0.832324099, 0.854213614, 0.854213614, 0.901096185, 0.909179025,
                                      0.919744831, 0.925915091, 0.928772959))

ourESF <- meteESF(valid$species, valid$abund)
ourSAD <- sad(ourESF)

test_that('meteESF returns correct class', {
  expect_is(ourESF, 'meteESF')
})

test_that('SAD returns correct class', {
  expect_is(ourSAD, 'meteDist')
})

test_that('rank SAD from empirical data works', {
  expect_true(all((sort(meteDist2Rank(ourSAD)) - valid$newmanPred)/valid$newmanPred < 0.05))
})

test_that('cdf SAD from epirical data works', {
  newmanCDF <- valid$newmanPredCDF[match(unique(valid$abund), valid$abund)]
  ourCDF <- ourSAD$p(valid$abund[match(unique(valid$abund), valid$abund)])
  
  expect_true(all(abs(ourCDF - newmanCDF) < 0.0005))
})



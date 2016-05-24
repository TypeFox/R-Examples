# This is file ../spam/tests/covmat.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     








options( echo=FALSE)
library( spam, warn.conflict=FALSE)  


cat("Results of the form '[1] TRUE' are from 'all.equal'\n\n")
h <- nearest.dist(100*1:10, 100*1:10+1:10, delta=10)
identical(cov.exp(1:10, 10), cov.exp(h, 10)@entries)
identical(cov.sph(1:10, 10), cov.exp(h, 10)@entries)
identical(cov.nug(1:10, 10), cov.exp(h, 10)@entries)
identical(cov.wu1(1:10, 10), cov.exp(h, 10)@entries)
identical(cov.wu1(1:10, 10), cov.exp(h, 10)@entries)
identical(cov.wu2(1:10, 10), cov.exp(h, 10)@entries)
identical(cov.wu3(1:10, 10), cov.exp(h, 10)@entries)
identical(cov.wend1(1:10, 10), cov.exp(h, 10)@entries)
identical(cov.wend2(1:10, 10), cov.exp(h, 10)@entries)
identical(cov.mat(1:10, 10), cov.exp(h, 10)@entries)

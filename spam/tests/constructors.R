# This is file ../spam/tests/constructors.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     


options( echo=FALSE)
library( spam, warn.conflict=FALSE)  

set.seed(14)
n <- 7
ln <- 20
A <- spam(0,n,n)
is <- sample(n,ln, replace=TRUE)
js <- sample(n,ln, replace=TRUE)

A[ unique( cbind(is,js)) ] <- 1:8


re <- A@rowpointers
rowpointers(A) <- re

# following will case error, thus the `try`
cat("===A series of errors caught by 'try':\n")
r <- re; r[1:2] <- rev(r[1:2]); try( rowpointers(A) <- r  )
r <- re; r[1] <- 0;             try( rowpointers(A) <- r  )
r <- re; r[n+1] <- 2;           try( rowpointers(A) <- r  )
r <- re; r[n+1] <- 20;          try( rowpointers(A) <- r  )
r <- c(rep(1,n),n+1);           try( rowpointers(A) <- r  )




ce <- A@colindices
colindices(A) <- ce

r <- ce; r[1:4] <- rev(r[1:4]); try( colindices(A) <- r  )
r <- ce; r[1] <- 0;             try( colindices(A) <- r  )
r <- ce; r[1] <- 20;            try( colindices(A) <- r  )

entries(A) <- A@entries
try( entries(A) <- as.logical(A@entries))
try( entries(A) <- c(r,1))
try( entries(A) <- r[-1])

try( dimension(A) <- c(1,2))

cat('===end of errors.\n')

options( echo=TRUE)

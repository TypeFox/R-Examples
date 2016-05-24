# Reordering eigenvalues
library(QZ, quiet = TRUE)

select <- c(TRUE, FALSE, FALSE, TRUE)
(ret <- qz(exAB1$A, exAB1$B, select = select))

### http://www.nag.com/lapack-ex/node119.html
select <- c(TRUE, FALSE, FALSE, TRUE)
(ret <- qz(exAB2$A, exAB2$B, select = select))
(ret <- ordqz(exAB2$A, exAB2$B, keyword = "ref"))
(ret <- ordqz(exAB2$A, exAB2$B, keyword = "cef"))

select <- c(TRUE, FALSE, FALSE, TRUE)
(ret <- qz(exA1$A, select = select))

### http://www.nag.com/lapack-ex/node89.html
select <- c(TRUE, FALSE, FALSE, TRUE)
(ret <- qz(exA2$A, select = select))
(ret <- ordqz(exA2$A, keyword = "lhp"))
(ret <- ordqz(exA2$A, keyword = "rhp"))
(ret <- ordqz(exA2$A, keyword = "ref"))
(ret <- ordqz(exA2$A, keyword = "cef"))


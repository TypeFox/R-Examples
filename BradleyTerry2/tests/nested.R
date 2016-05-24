## nested use of BTm
## (in response to Jing Hua Zhao's bug report)
library(BradleyTerry2)
myfun <- function(x) {
    c2b <- countsToBinomial(x)
    names(c2b) <- c("allele1", "allele2", "transmitted", "nontransmitted")
    btx <- BTm(cbind(transmitted, nontransmitted), allele1, allele2,
               ~allele, id = "allele", data = c2b)
}

x <- matrix(c(0,0, 0, 2, 0,0, 0, 0, 0, 0, 0, 0,
              0,0, 1, 3, 0,0, 0, 2, 3, 0, 0, 0,
              2,3,26,35, 7,0, 2,10,11, 3, 4, 1,
              2,3,22,26, 6,2, 4, 4,10, 2, 2, 0,
              0,1, 7,10, 2,0, 0, 2, 2, 1, 1, 0,
              0,0, 1, 4, 0,1, 0, 1, 0, 0, 0, 0,
              0,2, 5, 4, 1,1, 0, 0, 0, 2, 0, 0,
              0,0, 2, 6, 1,0, 2, 0, 2, 0, 0, 0,
              0,3, 6,19, 6,0, 0, 2, 5, 3, 0, 0,
              0,0, 3, 1, 1,0, 0, 0, 1, 0, 0, 0,
              0,0, 0, 2, 0,0, 0, 0, 0, 0, 0, 0,
              0,0, 1, 0, 0,0, 0, 0, 0, 0, 0, 0),nrow=12)
colnames(x) <- 1:12
rownames(x) <- 1:12

xx <- myfun(x)

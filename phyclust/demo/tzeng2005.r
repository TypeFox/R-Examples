library(phyclust, quiet = TRUE)

### Read data set of Crohn's disease (Hugot et al. 2001, table 2).
data.path <- tools:::file_path_as_absolute(
               system.file("./data/crohn.phy", package = "phyclust"))
my.snp <- read.phylip(data.path, code.type = "SNP")
X <- my.snp$org
case.control <- gsub("hap.*-", "", my.snp$seqname)

### Obtain grouping information for RD and FD.
ret <- haplo.post.prob(X, ploidy = 1)
X.id.FD <- ret$FD.id[apply(ret$FD.post, 1, which.max)]
X.id.RD <- ret$RD.id[apply(ret$RD.post, 1, which.max)]

### Transmited/non-transmited test for FD and RD.
DATA <- data.frame(CC = case.control, FD = X.id.FD, RD = X.id.RD)
(TNT.FD <- chisq.test(table(DATA[, c(1, 2)])))
(TNT.RD <- chisq.test(table(DATA[, c(1, 3)])))

### Chi-squared test for truncated distribution.
tmp <- as.data.frame(table(ret$haplo$hap1code))
tmp <- tmp[order(tmp$Freq, decreasing = TRUE)[1:ret$g.truncate],]
DATA.tr <- DATA[ret$haplo$hap1code %in% tmp$Var1,]
(TNT.tr <- chisq.test(table(DATA.tr[, c(1, 2)])))

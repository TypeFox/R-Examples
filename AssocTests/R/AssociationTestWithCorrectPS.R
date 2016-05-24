AssociationTestWithCorrectPS <- function(PSfile="file1.geno", ASSfile="file2.geno", PHEfile="file3.pheno",
                                             m.splits=20, miss.val=9, outfile="ass_test_result.txt")
{
    V <- pcoc(genoFile=PSfile, num.splits=m.splits, miss.val=miss.val)

    xStr <- readLines(con=ASSfile)
    n <- nchar(xStr[1])
    m <- length(xStr)

    y <- read.table(file=PHEfile)[,1]

    res <- rep(0,m)
    for (i in 1:m)
    {
        g <- Str2Num(xStr[i])
        is.na(g[g==miss.val]) <- TRUE

        temp <- glm(y~.+g, family=binomial(link="logit"),data=V)
        a.1 <- summary(temp)$coefficients
        t.1 <- dim(a.1)
        res[i] <- a.1[t.1[1], t.1[2]]
    }

    write(res, file=outfile, ncolumns=1)
}

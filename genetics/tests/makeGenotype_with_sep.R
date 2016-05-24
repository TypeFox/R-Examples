library(genetics)

set.seed(12345)

# Create a test data set where there are several genotypes in columns
# of the form "A/T".
test1 <- data.frame(Tmt=sample(c("Control","Trt1","Trt2"),20, replace=TRUE),
                G1=sample(c("A/T","T/T","T/A",NA),20, replace=TRUE),
                N1=rnorm(20),
                I1=sample(1:100,20,replace=TRUE),
                G2=paste(sample(c("134","138","140","142","146"),20,
                                replace=TRUE),
                         sample(c("134","138","140","142","146"),20,
                                replace=TRUE),
                         sep=" / "),
                G3=sample(c("A /T","T /T","T /A"),20, replace=TRUE),
                comment=sample(c("Possible Bad Data/Lab Error",""),20,
                               rep=TRUE)
                )
test1


# now automatically convert genotype columns
geno1 <- makeGenotypes(test1)
geno1



set.seed(12345)

# Create a test data set where there are several genotypes in columns
# of the form "A_T".
test1.b <- data.frame(Tmt=sample(c("Control","Trt1","Trt2"),20, replace=TRUE),
                G1=sample(c("A_T","T_T","T_A",NA),20, replace=TRUE),
                N1=rnorm(20),
                I1=sample(1:100,20,replace=TRUE),
                G2=paste(sample(c("134","138","140","142","146"),20,
                                replace=TRUE),
                         sample(c("134","138","140","142","146"),20,
                                replace=TRUE),
                         sep=" _ "),
                G3=sample(c("A _T","T _T","T _A"),20, replace=TRUE),
                comment=sample(c("Possible Bad Data/Lab Error",""),20,
                               rep=TRUE)
                )
# now automatically convert genotype columns
geno1.b <- makeGenotypes(test1.b, sep="_")

stopifnot(identical(geno1,geno1.b))


##' Test for the association between a biallelic SNP and a
##' quantitative trait using the maximum value of the three
##' nonparametric trend tests derived for the recessive, additive, and
##' dominant models. It is a robust procedure against the genetic
##' models.
##'
##' Under the null hypothesis of no association, the vector of the
##' three nonparametric tests under the recessive, additive, and
##' dominant models asymptotically follows a three-dimensional normal
##' distribution. The p-value can be calculated using the function
##' \link{pmvnorm} in the R package "\pkg{mvtnorm}".
##'
##' This test is different from the MAX3 test using in the function
##' \link{max3}. On one hand, the NMAX3 applies to the quantitative
##' traits association studies. However, the MAX3 is used in the
##' case-control association studies. On the other hand, the NMAX3 is
##' based on the nonparametric trend test. However, the MAX3 is based
##' on the Cochran-Armitage trend test.
##'
##' @title The NMAX3 based on the nonparametric trend test in a
##' quantitative trait association study
##' @param y a numeric vector of the observed quantitative trait
##' values in which the \emph{i}th element is the trait value of the
##' \emph{i}th subject.
##' @param g a numeric vector of the observed genotype values (\code{0}, \code{1},
##' or \code{2} denotes the number of risk alleles) in which the \emph{i}th
##' element is the genotype value of the \emph{i}th subject for a
##' biallelic SNP. \code{g} has the same length as \code{y}.
##' @return A list with class "\code{htest}" containing the following components:
##' \tabular{llll}{
##' \code{statistic} \tab \tab \tab \cr
##' \tab \tab \tab the observed value of the test statistic.\cr
##' \code{p.value} \tab \tab \tab \cr
##' \tab \tab \tab the p-value for the test.\cr
##' \code{alternative} \tab \tab \tab \cr
##' \tab \tab \tab a character string describing the alternative hypothesis.\cr
##' \code{method} \tab \tab \tab \cr
##' \tab \tab \tab a character string indicating the type of test performed.\cr
##' \code{data.name} \tab \tab \tab \cr
##' \tab \tab \tab a character string giving the names of the data.
##' }
##' @author Lin Wang, Wei Zhang, and Qizhai Li.
##' @references W Zhang and Q Li. Nonparametric Risk and Nonparametric
##' Odds in Quantitative Genetic Association Studies. \emph{Science
##' Reports (2nd revision)}. 2015.
##' @references B Freidlin, G Zheng, Z Li, and JL Gastwirth. Trend
##' Tests for Case-Control Studies of Genetic Markers: Power, Sample
##' Size and Robustness. \emph{Human Heredity}. 2002; 53:146-152.
##' @references WG Cochran. Some Methods for Strengthening the Common
##' Chi-Square Tests. \emph{Biometrics}. 1954; 10:417-451.
##' @references P Armitage. Tests for Linear Trends in Proportions and
##' Frequencies. \emph{Biometrics}. 1955; 11:375-386.
##' @examples
##' g <- rbinom(1500, 2, 0.3)
##' y <- 0.5 + 0.25 * g + rgev(1500, 0, 0, 5)
##' nmax3(y, g)
##' @export
nmax3 <- function(y, g)
{
    Z.R <- npt(y, g, 0)[[1]]
    Z.A <- npt(y, g, 0.5)[[1]]
    Z.D <- npt(y, g, 1)[[1]]

    T.max <- max(abs(Z.R), abs(Z.A), abs(Z.D))
    corr.mat <- CorrMatNRTest(y, g)

### calculate the p-value
    pval <- 1-mvtnorm::pmvnorm(lower=-c(T.max,T.max,T.max),upper=c(T.max,T.max,T.max),sigma=corr.mat)[1]

    a <- deparse(substitute(y))
    b <- deparse(substitute(g))
    structure( 
    list(statistic = c(NMAX3 = T.max), 
        p.value = pval, 
        alternative = "the phenotype is significantly associated with the genotype", 
        method = "Nonparametric MAX3 test", 
        data.name = paste(a, "and", b, sep=" ")
        ), 
    .Names=c("statistic", "p.value", "alternative", "method", "data.name"), 
    class="htest"
    )
}

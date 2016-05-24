##' Conduct the single-marker test in an association study to test for
##' the association between the genotype at a biallelic marker and a
##' trait.
##'
##' Single-marker analysis is a core in many gene-based or
##' pathway-based procedures, such as the truncated p-value
##' combination and the minimal p-value.
##' @title Single-marker test
##' @param y a numeric vector of the observed trait values in which
##' the \emph{i}th element is for the \emph{i}th subject. The elements
##' could be discrete (\code{0} or \code{1}) or continuous. The missing value is
##' represented by \code{NA}.
##' @param g a numeric vector of the observed genotype values (\code{0}, \code{1},
##' or \code{2} denotes the number of risk alleles) in which the \emph{i}th
##' element is for the \emph{i}th subject. The missing value is
##' represented by \code{NA}. \code{g} has the same length as \code{y}.
##' @param covariates an optional data frame, list or environment
##' containing the covariates used in the model. The default is \code{NULL},
##' that is, there are no covariates.
##' @param min.count a critical value to decide which method is used
##' to calculate the p-value when the trait is discrete and \code{covariates
##' = NULL}. If the minimum number of the elements given a specific
##' trait value and a specific genotype value is less than
##' \code{min.count}, the Fisher exact test is adopted; otherwise, the
##' Wald test is adopted.
##' @param missing.rate the highest missing rate of the genotype
##' values that this function can tolerate.
##' @param y.continuous logical. If \code{TRUE}, \code{y} is continuous;
##' otherwise, \code{y} is discrete.
##' @return \code{smt} returns a list with class "\code{htest}".
##'
##' If y is continuous, the list contains the following components:
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
##' \tab \tab \tab a character string giving the names of the data. \cr
##' \code{sample.size} \tab \tab \tab \cr
##' \tab \tab \tab a vector giving the numbers of the subjects with the genotypes \code{0}, \code{1}, and \code{2} (\code{n0}, \code{n1}, and \code{n2}, respectively).
##' }
##' If y is discrete, the list contains the following components:
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
##' \tab \tab \tab a character string giving the names of the data. \cr
##' \code{sample.size} \tab \tab \tab \cr
##' \tab \tab \tab a vector giving the number of the subjects with the trait value \code{1} and the genotype \code{0} (\code{r0}), the number of the subjects with the trait value \code{1} and the genotype \code{1} (\code{r1}), the number of the subjects with the trait value \code{1} and the genotype \code{2} (\code{r2}), the number of the subjects with the trait value \code{0} and the genotype \code{0} (\code{s0}), the number of the subjects with the trait value \code{0} and the genotype \code{1} (\code{s1}), and the number of the subjects with the trait value \code{0} and the genotype \code{2} (\code{s2}).\cr
##' \code{bad.obs} \tab \tab \tab \cr
##' \tab \tab \tab a vector giving the number of the missing genotype values with the trait value \code{1} (\code{r.miss}), the number of the missing genotype values with the trait value \code{0} (\code{s.miss}), and the total number of the missing genotype values (\code{n.miss}).
##' }
##' @author Lin Wang, Wei Zhang, and Qizhai Li.
##' @examples
##' y <- rep(c(0, 1), 25)
##' g <- sample(c(0, 1, 2), 50, replace = TRUE)
##' smt(y, g, covariates = NULL, min.count=5,
##'         missing.rate=0.20, y.continuous = FALSE)
##' @export
smt <- function(y, g, covariates=NULL, min.count=5, missing.rate=0.20, y.continuous=FALSE)
{
    # g is the genotype vector taking 0, 1 and 2
    # y is the outcome
    a <- deparse(substitute(y))
    b <- deparse(substitute(g))
        
    dex <- !is.na(g) & !is.na(y)
    g.valid <- g[dex]
    y.valid <- y[dex]

    # y is continuous
    if (y.continuous)
    {
        n0 <- sum(g.valid==0)
        n1 <- sum(g.valid==1)
        n2 <- sum(g.valid==2)

        if(is.null(covariates))
        {
            temp <- lm(y~g)
        }
        else
        {
            temp <- lm(y~.+g, data=covariates)
        }
        a.1 <- summary(temp)$coefficients
        t.1 <- dim(a.1)
        p.value <- a.1[t.1[1], t.1[2]]

        res <- structure( 
              list(p.value = p.value, 
                  alternative = "the phenotype is significantly associated with the genotype", 
                  method = "Single-marker test", 
                  data.name = paste(a, "and", b, sep=" "),
                  sample.size = c(n0=n0, n1=n1, n2=n2)
                  ), 
              .Names=c("p.value", "alternative", "method", "data.name", "sample.size"), 
              class="htest"
              )
        return(res)
    }
    else
    {
        r0 <- sum(g.valid==0 & y.valid==1)
        r1 <- sum(g.valid==1 & y.valid==1)
        r2 <- sum(g.valid==2 & y.valid==1)

        s0 <- sum(g.valid==0 & y.valid==0)
        s1 <- sum(g.valid==1 & y.valid==0)
        s2 <- sum(g.valid==2 & y.valid==0)

        addr <- !is.na(y)
        u <- g[addr]
        v <- y[addr]
        r.miss <- sum(is.na(u[v==1]))
        s.miss <- sum(is.na(u[v==0]))

        n.miss <- sum(is.na(g))

        # missing rate
        m.r <- sum(is.na(g))/length(g)
        if (m.r>=missing.rate)
        {
            res <- structure( 
                  list(p.value = -9999, 
                      alternative = "the phenotype is significantly associated with the genotype", 
                      method = "Single-marker test", 
                      data.name = paste(a, "and", b, sep=" "),
                      sample.size = c(r0=r0, r1=r1, r2=r2, s0=s0, s1=s1, s2=s2),
                      bad.obs = c(r.miss=r.miss, s.miss=s.miss, n.miss=n.miss)
                      ), 
                  .Names=c("p.value", "alternative", "method", "data.name", "sample.size", "bad.obs"), 
                  class="htest"
                  )
            return(res)
        }
        else
        {
            if (min(CalExpect(matrix(c(r0,r1,r2,s0,s1,s2),nrow=2,byrow=TRUE)))>=min.count)
            {
                if (is.null(covariates))
                {
                    temp <- glm(y~g, family=binomial(link="logit"))
                }
                else
                {
                    temp <- glm(y~.+g, family=binomial(link="logit"), data=covariates)
                }
                a.1 <- summary(temp)$coefficients
                t.1 <- dim(a.1)
                p.value <- a.1[t.1[1], t.1[2]]
            }
            else
            {
                if (s0 >= s2)
                {
                    H <- matrix(c(r0,r1+r2,s0,s1+s2),ncol=2, byrow=TRUE)
                    if (min(CalExpect(H))<min.count)
                    {
                        p.value <- fisher.test(H)$p.value
                    }
                    else
                    {
                        g[g==2] <- 1
                        if (is.null(covariates))
                        {
                            temp <- glm(y~g, family=binomial(link="logit"))
                        }
                        else
                        {
                            temp <- glm(y~.+g, family=binomial(link="logit"), data=covariates)
                        }
                        a.1 <- summary(temp)$coefficients
                        t.1 <- dim(a.1)
                        p.value <- a.1[t.1[1], t.1[2]]
                    }
                }
                else
                {
                    H <- matrix(c(r0+r1,r2,s0+s1,s2),ncol=2,byrow=TRUE)
                    if (min(CalExpect(H))<min.count)
                    {
                        p.value <- fisher.test(H)$p.value
                    }
                    else
                    {
                        g[g==1] <- 0
                        g[g==2] <- 1
                        if (is.null(covariates))
                        {
                            temp <- glm(y~g, family=binomial(link="logit"))
                        }
                        else
                        {
                            temp <- glm(y~.+g, family=binomial(link="logit"), data=covariates)
                        }
                        a.1 <- summary(temp)$coefficients
                        t.1 <- dim(a.1)
                        p.value <- a.1[t.1[1], t.1[2]]
                    }
                }
            }

            res <- structure( 
                  list(p.value = p.value, 
                      alternative = "the phenotype is significantly associated with the genotype", 
                      method = "Single-marker test", 
                      data.name = paste(a, "and", b, sep=" "),
                      sample.size = c(r0=r0, r1=r1, r2=r2, s0=s0, s1=s1, s2=s2),
                      bad.obs = c(r.miss=r.miss, s.miss=s.miss, n.miss=n.miss)
                      ), 
                  .Names=c("p.value", "alternative", "method", "data.name", "sample.size", "bad.obs"), 
                  class="htest"
                  )
            return(res)
        }
    }
}

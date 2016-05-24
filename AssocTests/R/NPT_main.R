##' Test for the association between a genetic variant and a
##' non-normal distributed quantitative trait based on the
##' nonparametric risk under a specific genetic model.
##'
##' For a non-normal distributed quantitative trait, three genetic
##' models (recessive, additive and dominant) used commonly are
##' defined in terms of the nonparametric risk (NR). The recessive,
##' additive, and dominant models can be classified based on the
##' nonparametric risks. More specifically, the recessive, additive,
##' and dominant models refer to NR20 > NR10 = 1/2, NR12 = NR10>1/2, and
##' NR10 = NR20 > 1/2, respectively, where NR10 and NR20 are the
##' nonparametric risks of the groups with the genotypes 1 and 2
##' relative to the group with the genotype 0, respectively, and NR12
##' is the nonparametric risk of the group with the genotype 2
##' relative to the group with the genotype 1.
##'
##' \code{varphi} can be \code{0}, \code{0.5}, or \code{1} for the recessive, additive, or
##' dominant model, respectively. When \code{varphi} is \code{0}, the test is
##' constructed under the recessive model by pooling together the
##' subjects with the genotypes 0 and 1. Similarly, when \code{varphi}
##' is \code{1}, the test is constructed under the dominant model by pooling
##' together the subjects with the genotypes 1 and 2. When
##' \code{varphi} is \code{0.5}, the test is based on the weighted sum of
##' NR10 and NR12.
##'
##' @title The nonparametric trend test based on the nonparametric
##' risk under a given genetic model
##' @param y a numeric vector of the observed quantitative trait
##' values in which the \emph{i}th element corresponds to the trait
##' value of the \emph{i}th subject.
##' @param g a numeric vector of the observed genotype values (\code{0}, \code{1},
##' or \code{2} denotes the number of risk alleles) in which the \emph{i}th
##' element is the genotype value of the \emph{i}th subject for a
##' biallelic SNP. \code{g} has the same length as \code{y}.
##' @param varphi a numeric value which represents the genetic
##' model. It should be \code{0}, \code{0.5}, or \code{1}, which indicates that the
##' calculation is performed under the recessive, additive, or
##' dominant model, respectively. The default is \code{0.5}.
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
##' @examples
##' g <- rbinom(1500, 2, 0.3)
##' y <- 0.5 + 0.25 * g + rgev(1500, 0, 0, 5)
##' npt(y, g, varphi = 0.5)
##' @export
npt <- function(y, g, varphi)
{

      if(varphi==0.5)
      {
             n0 <- sum(g==0)
             n1 <- sum(g==1)
             n2 <- sum(g==2)
             n  <- n0 + n1 + n2

             Y0 <- y[g==0]
             Y1 <- y[g==1]
             Y2 <- y[g==2]

        # genotype 1 relative to genotype 0
             d0  <- 1:n0
             d1  <- (n0+1):(n0+n1)
             RK1  <- rank(c(Y0,Y1))
             rr1 <- ( sum(RK1[d1]) - n1*(n1+1)/2 ) / (n0*n1)

        #genotype 2 relative to genotype 0
             d2  <- (n0+1):(n0+n2)
             RK2  <- rank(c(Y0,Y2))
             rr2 <- ( sum(RK2[d2])- n2*(n2+1)/2 ) / (n0*n2)

       # genotype 2 relative to genotype 1
             d3 <- (n1+1):(n1+n2)
             RK12 <- rank(c(Y1,Y2))
             r12 <- (sum(RK12[d3])-n2*(n2+1)/2 ) /(n1*n2)

       # standard deviation
             temp1 <- rep(NA, n0)
             temp2 <- rep(NA, n0)
             for (i in 1:n0)
             {
                temp1[i] <- ( sum(Y1>Y0[i])/n1 - 1/2 )^2
                temp2[i] <- ( sum(Y2>Y0[i])/n2 - 1/2 )^2
             }
             temp3 <- rep(NA, n1)
             temp4 <- rep(NA, n2)
             for (j in 1:n1)
             {
                temp3[j] <- ( sum(Y0<Y1[j])/n0 - 1/2 )^2
             }
             for (j in 1:n2)
             {
                temp4[j] <- ( sum(Y0<Y2[j])/n0 - 1/2 )^2
             }
             temp5 <- rep(NA, n1)
             temp6 <- rep(NA, n1)
             for (j in 1:n1)
             {
                  temp5[j] <- sum(Y0<Y1[j])/n0 - 1/2
                  temp6[j] <- sum(Y2>Y1[j])/n2 - 1/2
             }
             temp7 <- temp6^2
             temp8 <- rep(NA, n2)
             for(k in 1:n2)
             {
                 temp8[k] <- (sum(Y1<Y2[k])/n1 - 1/2)^2
             }

       ##the variance of the test for the Additive model
            rr1.var <- 1/(n0*n1) * ( (n1-1)/n0 *sum(temp1) + (n0-1)/n1*sum(temp3) + 1/4 )
            rr2.var <- 1/(n0*n2) * ( (n2-1)/n0 *sum(temp2) + (n0-1)/n2*sum(temp4) + 1/4 )
            r12.var <- 1/(n1*n2) * ( (n2-1)/n1 *sum(temp7) + (n1-1)/n2*sum(temp8) + 1/4 )
            cov12   <- (1/n1)^2*sum(temp5*temp6)

        # estimate the final probability P(Y2>Y1)
            a1 <-  sqrt((n0+n1)/rr1.var)
            a2 <-  sqrt((n1+n2)/r12.var)
            w1 <-  a1/(a1+a2)
            w2 <-  a2/(a1+a2)
            eta1 <- w1*rr1 + w2*r12

            eta1.var <- w1^2*rr1.var+ w2^2*r12.var+  2*w1*w2*cov12
            z <- (eta1-0.5)/sqrt(eta1.var)
            p.val <- 2*(1-pnorm(abs(z)))

      }
      else
      {
         #### under the recessive model or dominant model
             g[g==1] <- 2*varphi
             n0 <- sum(g==0)
             n2 <- sum(g==2)
             Y0 <- y[g==0]
             Y2 <- y[g==2]
         # genotype 2 relative to genotype 0
             d2  <- (n0+1):(n0+n2)
             RK2  <- rank(c(Y0,Y2))
             rr2 <- ( sum(RK2[d2])- n2*(n2+1)/2 ) / (n0*n2)

         # standard deviation
             temp2 <- rep(NA, n0)
             for (i in 1:n0)
             {
                  temp2[i] <- ( sum(Y2>Y0[i])/n2 - 1/2 )^2
             }
             temp4 <- rep(NA, n2)
             for (j in 1:n2)
             {
                  temp4[j] <- ( sum(Y0<Y2[j])/n0 - 1/2 )^2
             }


             rr2.var <- 1/(n0*n2) * ( (n2-1)/n0 *sum(temp2) + (n0-1)/n2*sum(temp4) + 1/4 )
             rr2.sd <- sqrt(rr2.var)
             z <- (rr2-0.5)/rr2.sd
             p.val <- 2*(1-pnorm(abs(z)))

      }

    a <- deparse(substitute(y))
    b <- deparse(substitute(g))
    structure( 
    list(statistic = c(NPT = z), 
        p.value = p.val, 
        alternative = "the phenotype is significantly associated with the genotype", 
        method = "Nonparametric trend test", 
        data.name = paste(a, "and", b, sep=" ")
        ), 
    .Names=c("statistic", "p.value", "alternative", "method", "data.name"), 
    class="htest"
    )
}

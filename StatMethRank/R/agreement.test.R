#' Test for Agreement of Ranking Data Among Groups
#'
#' This function performs a test of agreement among groups.
#'
#' @param data a data frame of the frequencies of all possible rankings 
#' given by different groups 
#' @param method whether the test is based on Spearman metric or Kendall metric
#' @return a list of test statistics
#' @export
#' @author Li Qinglong <liqinglong0830@@163.com>
#' @examples
#' data(Sutton)
#' agreement.test(Sutton, method = "spearman")
#' agreement.test(Sutton, method = "kendall")
#' @references Intergroup Diversity and Concordance for Ranking Data: An Approach via 
#' Metrics for Permuatations, Paul D. Feigin and Mayer Alvo
agreement.test <- function(data, method = c("spearman", "kendall"))
{
    
    facItems = dim(data)[1] # number of n!
    nMax = 15
    nItems = sum(facItems / cumprod(1:nMax) >= 1) # solve n from n! 
    nGroups = dim(data)[2] - nItems # number of populations

    n = colSums(data[, nItems + seq(nGroups)]) # number of observations in each group
    N = sum(n) # number of total observations
    lambda = n / N # Proportion in each group
    # relative frequency of occurrence of each ranking
    p = mapply("/", data[nItems + seq(nGroups)], n)

    method = match.arg(method)
    if (method == "spearman")
    {
        T = t(data[, 1:nItems]) - (nItems + 1) / 2
        C = nItems * (nItems^2 - 1) / 12
    }
    else if (method == "kendall")
    {
        Tkj = matrix(rep(c(t(data[, 1:nItems])), rep(0:(nItems - 1), facItems)), nrow = nItems)
        Tki = matrix(rep(c(t(data[, 1:nItems])), rep((nItems - 1):0, facItems)), nrow = nItems)
        T = sign(Tkj - Tki)
        C = nItems * (nItems - 1) / 2
    }

    # Relative frquency vector 
    f = p %*% lambda
    # Diversity apportionment
    H_total = C - t(f) %*% t(T) %*% T %*% f
    H_within = C - sum(colSums((T %*% p)^2) * lambda)
    H_between = H_total - H_within
    # Relative Diversity
    alpha = H_within / H_total
    # Relative Similarity
    pho = (C - H_total) / (C - H_within)


    if (nGroups == 2)
    {
        Sigma_separate = matrix(0, nrow = facItems, ncol = facItems)
        Sigma_pooled = matrix(0, nrow = facItems, ncol = facItems)
        for (i in 1:2)
        {
            Sigma_separate = Sigma_separate + (diag(p[, i]) - p[, i] %*% t(p[, i])) / (n[i] - 1)
            Sigma_pooled = Sigma_pooled + (diag(p[, i])- p[, i] %*% t(p[, i])) * n[i]
        }
        Sigma_separate = Sigma_separate * N
        Sigma_pooled = Sigma_pooled * N / (N-2) * sum(1 / n)

        # Separate tests    
        D_separate = ginv(T %*% Sigma_separate %*% t(T))
        # Degree of freedom of chi-square test statistic
        df = qr(T %*% Sigma_separate %*% t(T))$rank
        # chi-square test statistic
        chisq_separate = N * t(p[, 1] - p[, 2]) %*% t(T) %*% D_separate %*% T %*% (p[, 1] - p[, 2])
        # F separate test statistic(df, min(n1,n2)-df)
        F_separate = (min(n) - df) / ((min(n) - 1) * df) * chisq_separate
        
        p_value = 1 - pchisq(chisq_separate, df)
        chisq_separate = list(stat = chisq_separate, df = df, p_value = p_value)
        # Degree of freedom of F test statistic
        df = c(df, min(n) - df)
        p_value = 1 - pf(F_separate, df[1], df[2])
        F_separate = list(stat = F_separate, df = df, p_value = p_value)

        # Pooled tests
        # D_pooled = MASS::ginv(T %*% Sigma_pooled %*% t(T))
		D_pooled = ginv(T %*% Sigma_pooled %*% t(T))
        # Degree of freedom of chi-square test statistics
        df = qr(T %*% Sigma_pooled %*% t(T))$rank
        # chi-square test statistic
        chisq_pooled = N * t(p[, 1] - p[, 2]) %*% t(T) %*% D_pooled %*% T %*% (p[, 1] - p[, 2])
        # F pooled test statistic(df, N-1-df)
        F_pooled = (N - df - 1) / ((N - 2) * df) * chisq_pooled    

        p_value = 1 - pchisq(chisq_pooled, df)
        chisq_pooled = list(stat = chisq_pooled, df = df, p_value = p_value)
        df = c(df, N - df - 1)
        p_value = 1 - pf(F_pooled, df[1], df[2])
        F_pooled = list(stat = F_pooled, df = df, p_value = p_value)

        res = list( Within_Diversity = H_within,
                    Between_Diversity = H_between,
                    Total_Diversity = H_total,
                    Relative_Diversity = alpha,
                    Relative_Similarity = pho,
                    Separate_Chi_square = chisq_separate,
                    F_separate = F_separate,
                    Pooled_Chi_square = chisq_pooled,
                    F_pooled = F_pooled
                )
    }    
    else
    {
        res = list( Within_Diversity = H_within,
                    Between_Diversity = H_between,
                    Total_Diversity = H_total,
                    Relative_Diversity = alpha,
                    Relative_Similarity = pho
                )

    }
    return(res)
}
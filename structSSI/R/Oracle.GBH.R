
  # This function performs the Group Benjamini-Hochberg procedure at level alpha.
  #
  # Input: 1) unadjp - [vector numerics between 0 and 1] - The unadjusted p-values that have been computed by testing
  # multiple hypotheses tests. This is the same as what one would input into the
  # multtest function mt.rawp2adjp(), for example.
  # 2) group.index - [same as other groups vectors] - A vector that contains the group labeling for each hypothesis where
  # the index of each hypothesis is the same as the index used for the unadjusted
  # p-values above. For example, if the first 3 hypotheses were in group 1 and
  # the last 2 were in group 2, then we would input (1, 1, 1, 2, 2) as the ``groups''
  # vector.
  # 3) pi.groups - [vector numerics between 0 and 1, length = # groups] - This is the
  # vector of proportion of true null hypotheses that corresponds
  # to each group. In the adaptive procedure, these proportions are estimated from
  # the data. For example, if the hypotheses were divided into two groups, and
  # we found that the proportion of true null hypotheses in the first group was 0.9
  # and the proportion of true null hypotheses in the second group was 0.6, then we
  # would input pi.groups = c(0.9, 0.6). This is the aspect of the Group-Benjamini
  # Hochberg procedure that gives it more information to operate on than the
  # standard BH procedure: We are able to direct attention to those hypotheses
  # that are in groups where the proportion of true null hypotheses is estimated
  # to be low.
  #
  # Output: 1) GBH.results - A list including (a) the sorted adjusted p-values and (b) the
  # indices of the original hypotheses /uadjusted p-values that these adjusted p-values
  # correspond to (c) the hypotheses that were rejected, labeled by their original
  # unadj.p indexing, and (d) the hypotheses that were not rejected, also labeled
  # by their original unadj.p indexing.

Oracle.GBH <- function(unadj.p, group.index, pi.groups, alpha = 0.05){
    if(!all(group.index %in% names(pi.groups))) {
        stop('Names of pi.groups vector must match
              the elements of groups vector.')
    }
    if(length(unadj.p) != length(group.index)) {
        stop('Length of p values vector does not match
              group indexing vector.')
    }
    
    p.weighted <- unadj.p
    N <- length(unadj.p)

    names.sort <- sort(names(pi.groups))
    n_g <- table(group.index)
    pi0 <- 1/N * sum(n_g[names.sort] * pi.groups[names.sort])

    # The first part of the procedure involves weighting p-values. This is where the
    # known group structure information is being explicitly accounted for.
    pi.groups.match <- pi.groups[as.character(group.index)]
    p.weighted <- unadj.p * (pi.groups.match / (1 - pi.groups.match))

    # The second part of the procedure is exactly like Benjamini-Hochberg in that it
    # is a step-up procedure where we compared ordered p-values to some constant factor
    # times alpha, where the constant is determined by the position of the p-value in
    # the ordered list.

    sorting.weighted.p <- sort(p.weighted, index.return = TRUE)
    p.weighted <- sorting.weighted.p$x
    p.weighted.index <- sorting.weighted.p$ix

    if(pi0 < 1) {
        adjp.temp <- N * (1 - pi0) * p.weighted / 1 : N
    } else {
        adjp.temp <- p.weighted
    }
    
    adjp <- StepUp(adjp.temp)

    p.vals <- data.frame('unadjp' = unadj.p[p.weighted.index],
                         'adjp' = adjp,
                         'group' = group.index[p.weighted.index],
                         'adj.significance' = SignificanceStars(alpha, adjp))
    rownames(p.vals) <- names(unadj.p)[p.weighted.index]
    GBH.result <- new('GBH', p.vals = p.vals,
                      pi0 = pi.groups, adaptive = F, alpha = alpha)
    return(GBH.result)
}

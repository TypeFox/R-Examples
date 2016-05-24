######################################################################################################################

# Function: MultipleSequenceGatekeepingAdj.
# Argument: rawp, Raw p-value.
#           par, List of procedure parameters: vector of family (1 x m) Vector of component procedure labels ('BonferroniAdj.global' or 'HolmAdj.global' or 'HochbergAdj.global' or 'HommelAdj.global') (1 x nfam) Vector of truncation parameters for component procedures used in individual families (1 x nfam)

# Description: Computation of adjusted p-values for gatekeeping procedures based on the modified mixture methods (ref Dmitrienko et al. (2014))

MultipleSequenceGatekeepingAdj = function(rawp, par) {
  # Determine the function call, either to generate the p-value or to return description
  call = (par[[1]] == "Description")

  if (any(call == FALSE) | any(is.na(call))) {
    # Error check
    if (is.null(par[[2]]$family)) stop("Analysis model: Multiple sequence gatekeeping procedure: Hypothesis families must be specified.")
    if (is.null(par[[2]]$proc)) stop("Analysis model: Multiple sequence gatekeeping procedure: Procedures must be specified.")
    if (is.null(par[[2]]$gamma)) stop("Analysis model: Multiple sequence gatekeeping procedure: Gamma must be specified.")

    # Number of p-values
    nhyp = length(rawp)
    # Extract the vector of family (1 x m)
    family = par[[2]]$family
    # Number of families in the multiplicity problem
    nfam = length(family)
    # Number of null hypotheses per family
    nperfam = nhyp/nfam

    # Extract the vector of procedures (1 x m)
    proc = paste(unlist(par[[2]]$proc), ".global", sep = "")
    # Extract the vector of truncation parameters (1 x m)
    gamma = unlist(par[[2]]$gamma)

    # Simple error checks
    if (nhyp != length(unlist(family)))
      stop("Multiple-sequence gatekeeping adjustment: Length of the p-value vector must be equal to the number of hypothesis.")
    if (length(proc) != nfam)
      stop("Multiple-sequence gatekeeping adjustment: Length of the procedure vector must be equal to the number of families.") else {
        for (i in 1:nfam) {
          if (proc[i] %in% c("BonferroniAdj.global", "HolmAdj.global", "HochbergAdj.global", "HommelAdj.global") == FALSE)
            stop("Multiple-sequence gatekeeping adjustment: Only Bonferroni (BonferroniAdj), Holm (HolmAdj), Hochberg (HochbergAdj) and Hommel (HommelAdj) component procedures are supported.")
        }
      }
    if (length(gamma) != nfam)
      stop("Multiple-sequence gatekeeping adjustment: Length of the gamma vector must be equal to the number of families.") else {
        for (i in 1:nfam) {
          if (gamma[i] < 0 | gamma[i] > 1)
            stop("Multiple-sequence gatekeeping adjustment: Gamma must be between 0 (included) and 1 (included).") else if (proc[i] == "bonferroni.global" & gamma[i] != 0)
              stop("Multiple-sequence gatekeeping adjustment: Gamma must be set to 0 for the global Bonferroni procedure.")
        }
      }
    # Number of intersection hypotheses in the closed family
    nint = 2^nhyp - 1

    # Construct the intersection index sets (int_orig) before the logical restrictions are applied.  Each row is a vector of binary indicators (1 if the hypothesis is
    # included in the original index set and 0 otherwise)
    int_orig = matrix(0, nint, nhyp)
    for (i in 1:nhyp) {
      for (j in 0:(nint - 1)) {
        k = floor(j/2^(nhyp - i))
        if (k/2 == floor(k/2))
          int_orig[j + 1, i] = 1
      }
    }

    # Construct the intersection index sets (int_rest) and family index sets (fam_rest) after the logical restrictions are applied.  Each row is a vector of binary
    # indicators (1 if the hypothesis is included in the restricted index set and 0 otherwise)
    int_rest = int_orig
    fam_rest = matrix(1, nint, nhyp)
    for (i in 1:nint) {
      for (j in 1:(nfam - 1)) {
        for (k in 1:nperfam) {
          # Index of the current null hypothesis in Family j
          m = (j - 1) * nperfam + k
          # If this null hypothesis is included in the intersection hypothesis all dependent null hypotheses must be removed from the intersection hypothesis
          if (int_orig[i, m] == 1) {
            for (l in 1:(nfam - j)) {
              int_rest[i, m + l * nperfam] = 0
              fam_rest[i, m + l * nperfam] = 0
            }
          }
        }
      }
    }

    # Number of null hypotheses from each family included in each intersection before the logical restrictions are applied
    korig = matrix(0, nint, nfam)

    # Number of null hypotheses from each family included in the current intersection after the logical restrictions are applied
    krest = matrix(0, nint, nfam)

    # Number of null hypotheses from each family after the logical restrictions are applied
    nrest = matrix(0, nint, nfam)

    # Compute korig, krest and nrest
    for (j in 1:nfam) {
      # Index vector in the current family
      # index = which(family == j)
      index = family[[j]]
      korig[, j] = apply(as.matrix(int_orig[, index]), 1, sum)
      krest[, j] = apply(as.matrix(int_rest[, index]), 1, sum)
      nrest[, j] = apply(as.matrix(fam_rest[, index]), 1, sum)
    }

    # Vector of intersection p-values
    pint = rep(1, nint)

    # Matrix of component p-values within each intersection
    pcomp = matrix(0, nint, nfam)

    # Matrix of family weights within each intersection
    c = matrix(0, nint, nfam)

    # P-value for each hypothesis within each intersection
    p = matrix(0, nint, nhyp)

    # Compute the intersection p-value for each intersection hypothesis
    for (i in 1:nint) {

      # Compute component p-values
      for (j in 1:nfam) {
        # Consider non-empty restricted index sets
        if (krest[i, j] > 0) {
          # Restricted index set in the current family
          int = int_rest[i, family[[j]]]
          # Set of p-values in the current family
          pv = rawp[family[[j]]]
          # Select raw p-values included in the restricted index set
          pselected = pv[int == 1]
          # Total number of hypotheses used in the computation of the component p-value
          tot = nrest[i, j]
          pcomp[i, j] = do.call(proc[j], list(pselected, tot, gamma[j]))
        } else if (krest[i, j] == 0)
          pcomp[i, j] = 1
      }

      # Compute family weights
      c[i, 1] = 1
      for (j in 2:nfam) {
        c[i, j] = c[i, j - 1] * (1 - errorfrac(krest[i, j - 1], nrest[i, j - 1], gamma[j - 1]))
      }

      # Compute the intersection p-value for the current intersection hypothesis
      pint[i] = pmin(1, min(pcomp[i, ]/c[i, ]))
      # Compute the p-value for each hypothesis within the current intersection
      p[i, ] = int_orig[i, ] * pint[i]

    }

    # Compute adjusted p-values
    adjustedp = apply(p, 2, max)
    result = adjustedp
  }
  else if (call == TRUE) {
    family = par[[2]]$family
    nfam = length(family)
    proc = unlist(par[[2]]$proc)
    gamma = unlist(par[[2]]$gamma)
    test.id=unlist(par[[3]])
    proc.par = data.frame(nrow = nfam, ncol = 4)
    for (i in 1:nfam){
      proc.par[i,1] = i
      proc.par[i,2] = paste0("{",paste(test.id[family[[i]]], collapse = ", "),"}")
      proc.par[i,3] = proc[i]
      proc.par[i,4] = gamma[i]
    }
    colnames(proc.par) = c("Family", "Hypotheses", "Component procedure", "Truncation parameter")
    result=list(list("Multiple-sequence gatekeeping"),list(proc.par))
  }



  return(result)
}
# End of MultipleSequenceGatekeepingAdj
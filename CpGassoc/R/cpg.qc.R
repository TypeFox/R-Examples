cpg.qc <-
function (beta.orig, siga, sigb, pval, p.cutoff = 0.001, cpg.miss = NULL,
    sample.miss = NULL, constant100 = FALSE, sig.return = FALSE)
{
    avg_sig = colMeans(siga + sigb, na.rm = TRUE)
    experimentwide_median_sig = median(avg_sig)
    flag_to_remove = (avg_sig < 0.5 * experimentwide_median_sig |
        avg_sig < 2000)
    gc()
    print(paste("Removed", sum(flag_to_remove), "samples with low signal"))
    siga = siga[, !flag_to_remove]
    sigb = sigb[, !flag_to_remove]
    const <- ifelse(constant100, 100, 0)
    beta.new = sigb/(siga + sigb + const)
    if (!sig.return) {
        rm(siga, sigb)
    }
    gc()
    beta.orig = beta.orig[, !flag_to_remove]
    beta.new <- as.matrix(beta.new)
    beta.new[is.na(beta.orig)] = NA
    pval = pval[, !flag_to_remove]
    beta.new[pval > p.cutoff] = NA
    if (sig.return) {
        siga <- as.matrix(siga)
        sigb <- as.matrix(sigb)
        siga[which(is.na(beta.orig) | pval > p.cutoff)] = NA
        sigb[which(is.na(beta.orig) | pval > p.cutoff)] = NA
    }
    rm(beta.orig, pval)
    gc()
    if (!is.null(sample.miss) | !is.null(cpg.miss)) {
        missing.location <- which(is.na(beta.new))
        cpg.missing <- missing.location%%nrow(beta.new)
        cpg.missing[which(cpg.missing == 0)] <- nrow(beta.new)
        cpg.missing.table <- table(cpg.missing)/ncol(beta.new)
        if (!is.null(cpg.miss)) {
            remove.cpg <- which(cpg.missing.table > cpg.miss)
            remove.cpg <- as.numeric(names(remove.cpg))
            if (length(remove.cpg) > 0) {
                beta.new <- beta.new[-remove.cpg, ]
                if (sig.return) {
                  siga <- siga[-remove.cpg, ]
                  sigb <- sigb[-remove.cpg, ]
                }
            }
            print(paste("Removed", length(remove.cpg), "CpG sites with missing data for >",
                cpg.miss, "of samples"))
        }
        if (!is.null(sample.miss)) {
            missing.location <- which(is.na(beta.new))
            sample.missing <- (missing.location%/%nrow(beta.new)) + 1
            sample.missing.table <- table(sample.missing)/nrow(beta.new)
            remove.samp <- which(sample.missing.table > sample.miss)
            remove.samp <- as.numeric(names(remove.samp))
            if (length(remove.samp) > 0) {
                beta.new <- beta.new[, -remove.samp]
                if (sig.return) {
                  siga <- siga[, -remove.samp]
                  sigb <- sigb[, -remove.samp]
                }
            }
            print(paste("Removed", length(remove.samp), "samples with missing data for >",
                sample.miss, "of CpG sites"))
        }
    }
    if (!sig.return) {
        beta.new
    }
    else {
        list(beta.new, siga, sigb)
    }
}

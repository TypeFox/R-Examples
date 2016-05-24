CommT.anova <-
function (in_df) {
  # Descr:  calculate the anova
  # I/p:    in_df = data frame
  # O/p:    a table

    out_list = list()
    for (i in sprintf("%02d", 1:length(unique(in_df[,'gene_id'])))) {
        handle = in_df
        target = paste("gene0", i, sep="")
        handle[which(handle[,2]==target),4] = target
        handle[which(handle[,2]!=target),4] = "blocked"
        aov_results = summary(aov(KF_dist ~ grouping_var + Error(gene_id), data = handle))
        out_list[[target]] = aov_results$'Error: gene'[[1]][1,5]
    }

    out_m = t(as.data.frame(out_list))
    out_m = round(out_m, digits=4)
    colnames(out_m) = "Pr(>F)"

    return(out_m)
}

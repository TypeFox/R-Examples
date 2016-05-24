CommT.viz <-
function (in_df, title_str="a_project_name_here", alpha=0.05, legend_text, legend_pos) {
  # Descr:  visualize the tree distances
  # Deps:   ggplot2::ggplot
  # I/p:    in_df = input dataframe
  #         title_str = string with title name
  #         alpha = significance level
  #         legend_text = input legend
  #         legend_pos = legend position
  # O/p:    a plot

  # 0. Parse legend position information
    annot_x_pos = legend_pos$annot_x_pos
    annot_y_pos = legend_pos$annot_y_pos
    xlim_thres = legend_pos$xlim_thres_pos

  # 1. Define colors
    color_specs = CommT.plotcolors(n=2)

  # 2. Assess significance
    in_df[,ncol(in_df)+1] = "insign"
    label_sign = paste("gene", sprintf("%03d", which(legend_text < alpha)), sep="")
    in_df[which(in_df[,2]==label_sign),ncol(in_df)] = "signif"
    colnames(in_df)[ncol(in_df)] = "significance"

  # Special to avoid error 'no visible binding for global variable' during compilation
  # For details, see: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    KF_dist = gene_id = significance = NULL

  # 3. Generate plot
    plot_handle = ggplot2::ggplot(data=in_df) +
    geom_density(aes(x=KF_dist, group=gene_id, color=factor(significance), line=2)) +
    xlim(0, xlim_thres) +
    #ylim(0, 25) +
    theme_bw() +
    scale_colour_manual(name="significance", values=rev(color_specs)) +
    ggtitle(paste(title_str, ", alpha=", alpha, "\n", sep="")) +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5)))

  # 4. Add annotations
  #    Note: x-position easy to define, because KF distances between 0 and 1
    n_entries = length(legend_text)
    plot_handle = plot_handle + annotate("text", x=annot_x_pos, y=annot_y_pos, label=paste("gene", colnames(legend_text), sep="     "), color=color_specs[2], fontface="bold")
    for (i in 1:n_entries) {
        if (legend_text[i] < alpha) {
            plot_handle = plot_handle + annotate("text", x=annot_x_pos, y=annot_y_pos-(annot_y_pos*i/n_entries), label=paste(rownames(legend_text)[i], sprintf("%03f", legend_text[i]), sep="   "), color=color_specs[1], size=4)
        }
        if (legend_text[i] >= alpha) {
            plot_handle = plot_handle + annotate("text", x=annot_x_pos, y=annot_y_pos-(annot_y_pos*i/n_entries), label=paste(rownames(legend_text)[i], sprintf("%03f", legend_text[i]), sep="   "), color=color_specs[2], size=4)
        }
    }
  # 3. Return plot
    return(plot_handle)
}

CommT.legendpos <-
function (in_data) {
  # Descr:  generate location coordinates for legend
  # I/p:    in_data
  # O/p:    a list
  
    annot_x_list = annot_y_list = xlim_thres_list = c()
    for (i in split(in_data, in_data$gene_id)) {
        xlim_thres_list = c(annot_x_list, quantile(i[,"KF_dist"], probs = c(0.9999)))
        annot_y_list = c(annot_y_list, max(table(cut(i[,"KF_dist"], breaks=100))))
    }
    xlim_thres_pos = quantile(xlim_thres_list, probs = c(0.95))
    annot_x_pos = quantile(xlim_thres_list, probs = c(0.80))[[1]]
    annot_y_pos = max(annot_y_list)
    
    out_l = list("xlim_thres_pos"=xlim_thres_pos, "annot_x_pos"=annot_x_pos, "annot_y_pos"=annot_y_pos)

    return(out_l)
}

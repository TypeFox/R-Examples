plot.vip <- function(x, ncomp = 1, ...) {
  if(x$val.method == "none" | x$val.method == "loo") {
    df <- x$VIP[x$VIP$ncomp %in% ncomp, ]
    print(with(df, ggplot(df, aes(x = reorder(variables, -VIP, mean), y = VIP)) + 
            theme_bw() + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            geom_point(aes(col = variables), size = 5) + 
            ggtitle("Variable Importance in the Projection") + 
            facet_wrap(~ncomp) + 
            xlab("Variable") + 
            theme(legend.position = "none") + 
            geom_hline(yintercept = 1, col = "blue") + 
            theme(strip.text.x = element_text(size = 15, face = "bold")) + 
            theme(plot.title = element_text(size = 20)) + 
            theme(axis.title.x = element_text(size = 20)) +
            theme(axis.title.y = element_text(size = 20, angle = 90)) + 
            theme(axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, face = "bold")) + 
            theme(axis.text.y = element_text(size = 15, angle = 0, face = "bold"))))
  } else {
    df <- x$VIP.s[x$VIP.s$ncomp %in% ncomp, ]
    print(with(df, ggplot(df, aes(reorder(variables, -abs(values), mean), values)) + 
            theme_bw() + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
            geom_line(aes(col = variables)) + 
            geom_hline(yintercept = 1) + 
            ylab("VIP") + 
            xlab("Variable") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            ggtitle("Variable Importance in the Projection") + 
            facet_wrap(~ncomp) + 
            theme(strip.text.x = element_text(size = 15, face = "bold")) + 
            theme(legend.position = "none") + 
            theme(plot.title = element_text(size = 20)) + 
            theme(axis.title.x = element_text(size = 20)) +
            theme(axis.title.y = element_text(size = 20, angle = 90)) + 
            theme(axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, face = "bold")) + 
            theme(axis.text.y = element_text(size = 15, angle = 0, face = "bold"))))
  }
}

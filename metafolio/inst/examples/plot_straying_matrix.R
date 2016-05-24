# plot an example straying matrix

library(colorspace)
M <- generate_straying_matrix(10, 0.02, 0.1)
library(reshape2)
M <- melt(M)
shades <- rev(sequential_hcl(20))
cuts <- seq(min(M$value), max(M$value), length.out = 20)
M$colour <- shades[findInterval(M$value, cuts)]

#image(generate_straying_matrix(10, 0.01, 0.3), col = rev(sequential_hcl(11)))

plot_straying_matrix <- function(x, pal) {

par(mfrow = c(10, 10), cex = 0.7, mar = c(0,0,0,0), oma = c(4, 4, .5, .5))
for(i in 1:10) {
  for(j in 1:10) {
    #if(j < i){
      plot(1, 1, xlim = c(0, 1), ylim = c(0, 1), yaxs = "i", xaxs = "i", axes = FALSE, pch = 20, col = "grey20", type = "n")
    with(x[x$Var1 == i & x$Var2 == j, ], rect(0, 0, 1, 1, col = colour, border = NA))
    box(col = "grey50")
    #}else{
      #plot(1,1, type = "n", xlab = "", ylab = "", axes = FALSE)
    #}
    if(i == 10) mtext(j, side = 1, line = 1, cex = 0.8, col = pal[j])
    if(j == 1) mtext(i, side = 2, las = 1, line = 1, cex = 0.8, col = pal[i])
  }
}
    mtext("Receiving population", side = 1, line = 2.5, cex = 0.8, outer = TRUE)
    mtext("Source population", side = 2, line = 2.5, cex = 0.8, outer = TRUE, las = 0)
}

pdf_eps("stray-matrix", width = 4, height = 4, type = TYPE)
plot_straying_matrix(M, pal = col_pal)
dev.off()


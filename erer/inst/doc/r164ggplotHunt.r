# A. Run the program and generate the default graph
setwd("C:/aErer"); source("r072sunSJAF.r", echo = FALSE)
names(p1); class(p1); plot.maTrend

windows(width = 4, height = 3, pointsize = 9); bringToTop(stay = TRUE)
par(mai = c(0.7, 0.7, 0.1, 0.1), family = "serif"); plot(p1)

# B. Draw the graph by ggplot2
# B1. Data preparation: ggplot need a data frame 
pr <- p1$trend; class(pr); head(pr)
ya <- seq(from = 0.10, to = 0.45, by = 0.05)
yb <- substr(sprintf("%.2f", ya), start = 2, stop = 4)
mv <- colMeans(p1$q$w$x)[p1$nam.c]
             
# B2. Create a graph in ggplot
library(ggplot2)
f1 <- ggplot(data = pr, mapping = aes(x = HuntYrs)) +
  geom_line(aes(y = all.pr)) +
  geom_line(aes(y = Nonres_d1.pr), linetype = 2) + 
  geom_line(aes(y = Nonres_d0.pr), linetype = 3) +
  geom_vline(xintercept = mv, linetype = 4)

f2 <- f1 + 
  scale_x_continuous(name = "Hunting experience (Year)", 
    breaks = seq(from = 0, to = 70, by = 10)) +
  scale_y_continuous(name = "Prob(Insurance purchase = yes)",
    breaks = ya, labels = yb) +   
  theme_classic(base_size = 9, base_family = "serif") +  
  annotate(geom = "text", label = c("Nonresident", "All", "Resident"), 
    x = c(39, 50, 60), y = c(0.35, 0.22, 0.18), family = "serif", size= 3)
str(f2); class(f2)
  
# B3. Show the graph on a screen device
windows(width = 4, height = 3); f2

# B4. Save the graph on a file device
pdf(file = "fig_ggplotHunt.pdf", width = 4, height = 3)
f2
dev.off()
ggsave(filename = "fig_ggplotHunt.pdf", plot = f2, width = 4, height = 3)
# Title: Graph codes for Sun (2011 FPE)
library(apt); library(ggplot2); setwd('C:/aErer'); data(daVich)

# -------------------------------------------------------------------------
# A. Data for graphs: value, price, and t5$path
prVi <- daVich[, 1]; prCh <- daVich[, 2]
vaVi <- daVich[, 3]; vaCh <- daVich[, 4] 
(date <- as.Date(time(daVich), format = '%Y/%m/%d'))
(value <- data.frame(date, vaCh, vaVi))
(price <- data.frame(date, prVi, prCh))
(t5 <- ciTarThd(y=prVi, x=prCh, model = 'mtar', lag = 3, th.range = 0.15))

# -------------------------------------------------------------------------
# B. Traditonal graphics
# Figure 1 Import values from China and Vietnam
win.graph(width = 5, height = 2.8, pointsize = 9); bringToTop(stay = TRUE)
par(mai = c(0.4, 0.5, 0.1, 0.1), mgp = c(2, 1, 0), family = "serif")
plot(x = vaCh, lty = 1, lwd = 1, ylim = c(0, 60), xlab = '',
  ylab = 'Montly import value($ million)', axes = FALSE)
box(); axis(side = 1, at = 2002:2010) 
axis(side = 2, at = c(0, 20, 40, 60), las = 1)
lines(x = vaVi, lty = 4, lwd = 1)
legend(x = 2008.1, y = 59, legend = c('China', 'Vietnam'), 
  lty = c(1, 4), box.lty = 0)
fig1.base <- recordPlot()

# Figure 2 Import prices from China and Vietnam
win.graph(width = 5, height = 2.8, pointsize = 9)
par(mai = c(0.4, 0.5, 0.1, 0.1), mgp = c(2, 1, 0), family = "serif")
plot(x = prCh, lty = 1, type = 'l', lwd = 1, ylim = range(prCh, prVi),
  xlab = '', ylab = 'Monthly import price ($/piece)' )
lines(x = prVi, lty = 3, type = 'l', lwd = 1)
legend(x = 2008.5, y = 175, legend = c('China', 'Vietnam'), 
  lty = c(1, 3), box.lty = 0)

# Figure 3 Sum of dquared errors by threshold value from MTAR
win.graph(width = 5.1, height = 3.3, pointsize = 9)
par(mai = c(0.5, 0.5, 0.1, 0.1),  mgp = c(2.2, 1, 0), family = "serif")
plot(formula = path.sse ~ path.thr, data = t5$path, type = 'l',
     ylab = 'Sum of Squared Errors', xlab = 'Threshold value')

# -------------------------------------------------------------------------     
# C. ggplot for three figures
pp <- theme(axis.text   = element_text(size = 8, family = "serif")) +
      theme(axis.title  = element_text(size = 9, family = "serif")) +
      theme(legend.text = element_text(size = 9, family = "serif")) +
      theme(legend.position = c(0.85, 0.9) ) +
      theme(legend.key = element_rect(fill = 'white', color = NA)) +
      theme(legend.background = element_rect(fill = NA, color = NA))      

fig1 <- ggplot(data = value, aes(x = date)) +
  geom_line(aes(y = vaCh, linetype = 'China')) +
  geom_line(aes(y = vaVi, linetype = 'Vietnam')) + 
  scale_linetype_manual(name = '', values = c(1, 3)) +       
  scale_x_date(name = '', labels = as.character(2002:2010), breaks = 
    as.Date(paste(2002:2010, '-1-1', sep = ''), format = '%Y-%m-%d')) +
  scale_y_continuous(limits = c(0, 60), 
    name = 'Monthly import value ($ million)') + pp

fig2 <- ggplot(data = price, aes(x = date)) +
  geom_line(aes(y = prCh, linetype = 'China')) +
  geom_line(aes(y = prVi, linetype = 'Vietnam')) + 
  scale_linetype_manual(name = '', values = c(1, 3))+       
  scale_x_date(name = '', labels = as.character(2002:2010), breaks = 
    as.Date(paste(2002:2010, '-1-1', sep = ''), format = '%Y-%m-%d')) +
  scale_y_continuous(limits = c(98, 180), 
    name = 'Monthly import price ($/piece)') + pp

fig3 <- ggplot(data = t5$path) +
  geom_line(aes(x = path.thr, y = path.sse)) +
  labs(x = 'Threshold value', y = 'Sum of squared errors') +
  scale_y_continuous(limits = c(5000, 5700)) +
  scale_x_continuous(breaks = c(-10:7)) +
  theme(axis.text  = element_text(size = 8, family = "serif")) + 
  theme(axis.title = element_text(size = 9, family = "serif"))

# -------------------------------------------------------------------------
# D. Show on screen devices or save on file devices
pdf(file = 'OutBedFig1base.pdf', width = 5, height = 2.8, pointsize = 9)
replayPlot(fig1.base); dev.off()

windows(width = 5, height = 2.8); fig1 
windows(width = 5, height = 2.8); fig2 
windows(width = 5, height = 2.8); fig3 

ggsave(fig1, filename = 'OutBedFig1ggplot.pdf', width = 5, height = 2.8)
ggsave(fig2, filename = 'OutBedFig2ggplot.pdf', width = 5, height = 2.8)
ggsave(fig3, filename = 'OutBedFig3ggplot.pdf', width = 5, height = 2.8) 
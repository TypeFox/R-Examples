# A. Data preparation 
library(erer); data(daBedRaw, daBed)
tos <- window(daBedRaw[, "vWD"], start = c(2001, 1), end = c(2008, 12))
sha <- daBed[, c('sCN', 'sVN', 'sID')] * 100
bed <- ts.union(tos / 10 ^ 6, sha)
colnames(bed) <- c("TotExp", colnames(sha)); head(bed)
 
# B. Left side
windows(width = 5.5, height = 3, pointsize = 9); bringToTop(stay = TRUE)
par(mai = c(0.45, 0.45, 0.45, 0.4), omi = c(0.2, 0.2, 0.2, 0.2),
  mgp = c(2, 0.6, 0), family = "serif")
ts.plot(bed[, -1], lty = 1:3, ylim = range(bed[, -1]) + c(-2, 2), 
  xlab = '', ylab = 'Share (%)')
text(x = 2002.5, y = c(53, 15, 3), labels = 1:3)
(usr.left <- par("usr"))  # usr for left side
  
# C1. Right side: approach by par(new = TRUE)
par(new = TRUE)
plot(x = bed[, 'TotExp'], type = "l", lty = 1, col = "gray80", lwd = 2, 
  axes = FALSE, ann = FALSE, ylim = range(bed[, 1]) + c(-5, 5))   
axis(4); mtext("Total Expenditure ($ million)", side = 4, line = 2)
text(x = 2002.5, y = 73, labels = 4)
  
# C2. Right side: approach by resetting user coordinates
# old.usr <- par("usr")
# new.usr <- c(old.usr[1:2], range(bed[, 1]) + c(-5, 5))
# par(usr = new.usr)
# lines(bed[, 'TotExp'], lty = 1, col = "gray80", lwd = 2) 
# axis(4); mtext("Total Expenditure ($ million)", side = 4, line = 2)

# D. Regions and clipping
abline(v = 2004.6, lty = 2, xpd = NA)     # current whole device region
abline(v = 2004.8, lty = 2, xpd = TRUE)   # current figure region
abline(v = 2005.0, lty = 2, xpd = FALSE)  # current plot region; default

# smaller clipping area on the plotting region
(usr.right <- par("usr"))  # usr for right side
clip(x1 = 2004, x2 = 2008, y1 = 30, y2 = 70); par("usr")  # same usr  
abline(v = 2005.5, lty = 2)          # a very short vertical line
do.call("clip", as.list(usr.right))  # restore the original clipping
abline(v = 2005.7, lty = 2)

# E. Legend: two choices
cc <- c("black", "black", "black", "gray80")
legend(x = 2001.5, y = 157, ncol = 2,
  box.lty = 0, lty = c(1:3, 1), lwd = c(1, 1, 1, 2), col = cc, xpd = TRUE,
  legend = c("1. Share - China", "2. Share - Vietnam            ", 
    "3. Share - Indonesia", "4. Total expenditure"))
    
legend(x = "bottom", inset = -0.3, ncol = 4,
  box.lty = 0, lty = c(1:3, 1), lwd = c(1, 1, 1, 2), col = cc, xpd = TRUE,
  legend = c("1. China", "2. Vietnam", "3. Indonesia", "4. Expenditure"))
out <- recordPlot()  # save the graph on the screen device

# F. Save a PDF version
pdf(file = "C:/aErer/fig_region.pdf", width = 5.5, height = 3, 
  pointsize = 9, family = "serif")
replayPlot(out); dev.off()
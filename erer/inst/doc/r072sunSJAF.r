# Title: R program for Sun et al. (2007 SJAF)
# Date: January - May 2006

# -------------------------------------------------------------------------
# Brief contents
# 0. Libraries and global setting
# 1. Import raw data in csv format
# 2. Descriptive statistics
# 3. Logit regression and figures
# 4. Export results

# -------------------------------------------------------------------------
# 0. Libraries and global setting 
library(erer)  # functions: bsTab(), maBina(), maTrend()
wdNew <- 'C:/aErer'  # Set up working directory
setwd(wdNew); getwd(); dir()

# -------------------------------------------------------------------------
# 1. Import raw data in csv format
daInsNam <- read.table(file = 'RawDataIns1.csv', header = TRUE, sep = ',')
daIns <- read.table(file = 'RawDataIns2.csv', header = TRUE, sep = ',')
class(daInsNam); dim(daInsNam); print(daInsNam); class(daIns); dim(daIns) 
head(daIns); tail(daIns); daIns[1:3, 1:5]

# -------------------------------------------------------------------------
# 2. Descriptive statistics
(insMean <- round(x = apply(X = daIns, MARGIN = 2, FUN = mean), digits =2))
(insCorr <- round(x = cor(daIns),  digits = 3))
table.3 <- cbind(daInsNam, Mean = I(sprintf(fmt="%.2f", insMean)))[-14, ]
rownames(table.3) <- 1:nrow(table.3)
print(table.3, right = FALSE)

# -------------------------------------------------------------------------
# 3. Logit regression and figures
# 3.1 Logit regression
ra <- glm(formula = Y ~ Injury + HuntYrs + Nonres + Lspman + Lnong + 
                    Gender + Age + Race + Marital + Edu + Inc + TownPop,
          family = binomial(link = 'logit'), data = daIns, x = TRUE)
fm.fish <- Y ~ Injury + FishYrs + Nonres + Lspman + Lnong + 
  Gender + Age + Race + Marital + Edu + Inc + TownPop
rb <- update(object = ra, formula = fm.fish)
names(ra); summary(ra)
(ca <- data.frame(summary(ra)$coefficients))
(cb <- data.frame(summary(rb)$coefficients))

# 3.2 Marginal effect
(me <- maBina(w = ra))
(u1 <- bsTab(w = ra, need = '2T'))
(u2 <- bsTab(w = me$out, need = '2T'))
table.4 <- cbind(u1, u2)[, -4]
colnames(table.4) <- c('Variable', 'Coefficient', 't-ratio', 
  'Marginal effect', 't-ratio')
table.4

# 3.3 Figures: probability response curve
(p1 <- maTrend(q = me, nam.d = 'Nonres', nam.c = 'HuntYrs'))
(p2 <- maTrend(q = me, nam.d = 'Nonres', nam.c = 'Age'))
(p3 <- maTrend(q = me, nam.d = 'Nonres', nam.c = 'Inc'))

# Show one graph on screen device
windows(width = 4, height = 3, pointsize = 9)
bringToTop(stay = TRUE)
par(mai = c(0.7, 0.7, 0.1, 0.1), family = 'serif')
plot(p1)

# Save three graphs on file device
fname <- c('OutInsFig1a.png', 'OutInsFig1b.png', 'OutInsFig1c.png')
pname <- list(p1, p2, p3)
for (i in 1:3) {
  png(file = fname[i], width = 4, height = 3,
    units = 'in', pointsize = 9, res = 300)
  par(mai = c(0.7, 0.7, 0.1, 0.1), family = 'serif')
  plot(pname[[i]])
  dev.off()
}

# -------------------------------------------------------------------------
# 4. Export results 
write.table(x = table.3, file = 'OutInsTable3.csv', sep = ',')
write.table(x = table.4, file = 'OutInsTable4.csv', sep = ',')
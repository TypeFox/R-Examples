# Title: R Program for Wan et al. (2010 JAAE ); last revised Feb. 2010
library(RODBC); library(erer); library(xlsx)
options(width = 120); setwd("C:/aErer")

# -------------------------------------------------------------------------
# 1. Data import and transformation
# 1.1 Import raw data in Microsoft Excel format
dat <- odbcConnectExcel2007('RawDataAids.xlsx')
  sheet <- sqlTables(dat); sheet$TABLE_NAME
  impo <- sqlFetch(dat, "dataImport")
  expe <- sqlFetch(dat, "dataExp")
odbcClose(dat)
names(impo); names(expe); head(impo); tail(expe)

# 1.2 Expenditure data for Hausman test
ex <- ts(data = expe[, -c(1, 2)], start = c(1959, 1), end = c(2009, 12), 
  frequency = 12)
Exp <- window(ex, start = c(2001, 1), end = c(2008, 12), frequency = 12)
head(Exp); bsStat(Exp)

# 1.3 Raw import data, date selection, and transformation for AIDS
BedRaw <- ts(data = impo[, -c(1, 2)], start = c(1996, 1), 
  end = c(2008, 12), frequency = 12)
lab8 <- c("CN", "VN", "ID", "MY", "CA", "BR", "IT")
dumm <- list(dum1 = c(2003, 10, 2003, 10), dum2 = c(2004, 7, 2004, 7),
  dum3 = c(2005, 1, 2005, 1))
imp8 <- aiData(x = BedRaw, label = lab8, label.tot = "WD", 
  prefix.value = "v", prefix.quant = "q", 
  start = c(2001, 1), end = c(2008, 12), dummy = dumm)
imp5 <- update(imp8, label = c("CN", "VN", "ID", "MY")); names(imp5)
Bed  <- imp8$out; colnames(Bed)[18:20] <- c("dum1", "dum2", "dum3")

# 1.4 Three datasets saved in 'erer' library already
# Results in Wan (2010 JAAE ) can be reproduced with saved data directly.
data(daExp, daBedRaw, daBed); str(daExp); str(Exp) 

# -------------------------------------------------------------------------
# 2. Descriptive statistics (Table 1)
lab8 <- c("CN", "VN", "ID", "MY", "CA", "BR", "IT")
pig <- aiData(x = daBedRaw, label = lab8, label.tot = "WD", 
  prefix.value = "v", prefix.quant = "q", 
  start = c(2001, 1), end = c(2008, 12))
hog <- cbind(pig$share * 100, pig$price, pig$m / 10 ^ 6)
colnames(hog) <- c(paste("s", lab8, sep = ""), "sRW", 
                   paste("p", lab8, sep = ""), "pRW", "Expend")
dog <- bsStat(hog, two = TRUE, digits = 3)$fstat[, -6]
colnames(dog) <- c("Variable", "Mean", "St. Dev.", "Minimum", "Maximum")
dog[, -1] <- apply(X = dog[, -1], MARGIN = 2, 
  FUN = function(x) {sprintf(fmt="%.3f", x)})
(table.1 <- dog)

# -------------------------------------------------------------------------
# 3. Monthly expenditure and import shares by country (Figure 1)
tos <- window(daBedRaw[, "vWD"], start = c(2001, 1), end = c(2008, 12))
tot <- tos / 10 ^ 6
sha <- daBed[, c('sCN', 'sVN', 'sID', 'sMY', 'sCA', 'sBR', 'sIT')] * 100
y <- ts.union(tot, sha); colnames(y) <- c('TotExp', colnames(sha))
windows(width = 5.5, height = 5, family = 'serif', pointsize = 11)
plot(x = y, xlab = "", main = "", oma.multi = c(2.5, 0, 0.2, 0))

# ------------------------------------------------------------------------- 
# 4. Hausman test and revised data
# 4.1 Getting started with a static AIDS model 
sh <- paste("s", c(lab8, "RW"), sep = "")
pr <- paste("lnp", c(lab8, "RW"), sep = "")
du3 <- c("dum1", "dum2", "dum3"); du2 <- du3[2:3]
rSta <- aiStaFit(y = daBed, share = sh, price = pr, shift = du3, 
  expen = "rte", omit = "sRW", hom = TRUE, sym = TRUE)
summary(rSta)

# 4.2 Hausman test and new data
(dg <- daExp[, "dg"])
rHau <- aiStaHau(x = rSta, instr = dg, choice = FALSE)
names(rHau); colnames(rHau$daHau); colnames(rHau$daFit); rHau
two.exp <- rHau$daFit[, c("rte", "rte.fit")]; bsStat(two.exp, digits = 4)
plot(data.frame(two.exp)); abline(a = 0, b = 1)
daBedFit <- rHau$daFit

# -------------------------------------------------------------------------
# 5. Static and dynamic AIDS models
# 5.1 Diagnostics and coefficients (Table 2, 3, 4)
hSta  <- update(rSta, y = daBedFit, expen = "rte.fit")
hSta2 <- update(hSta, hom = FALSE, sym = F); lrtest(hSta2$est, hSta$est)
hSta3 <- update(hSta, hom = FALSE, sym = T); lrtest(hSta2$est, hSta3$est)
hSta4 <- update(hSta, hom = TRUE,  sym = F); lrtest(hSta2$est, hSta4$est)
hDyn  <- aiDynFit(hSta)
hDyn2 <- aiDynFit(hSta2); lrtest(hDyn2$est, hDyn$est)
hDyn3 <- aiDynFit(hSta3); lrtest(hDyn2$est, hDyn3$est)
hDyn4 <- aiDynFit(hSta4); lrtest(hDyn2$est, hDyn4$est)
(table.2 <- rbind(aiDiag(hSta), aiDiag(hDyn)))
(table.3 <- summary(hSta)); (table.4 <- summary(hDyn))

# 5.2 Own-price elasticities (Table 5)
es <- aiElas(hSta); ed <- aiElas(hDyn); esm <- edm <- NULL
for (i in 1:7) {
  esm <- c(esm, es$marsh[c(i * 2 - 1, i * 2), i + 1])
  edm <- c(edm, ed$marsh[c(i * 2 - 1, i * 2), i + 1])
}
MM <- cbind(es$expen[-c(15:16), ], esm, ed$expen[-c(15:16), 2], edm)
colnames(MM) <- c("Country", "LR.exp", "LR.Marsh", "SR.exp", "SR.Marsh")
(table.5 <- MM)

# 5.3 Cross-price elasticities (Table 6)
(table.6a <- es$hicks[-c(15:16), -9])
(table.6b <- ed$hicks[-c(15:16), -9])
for (j in 1:7) {
  table.6a[c(j * 2 - 1, j * 2), j + 1] <- "___"
  table.6b[c(j * 2 - 1, j * 2), j + 1] <- "___"
}
rown <- rbind(c("Long-run",  rep("", times = 7)), 
              c("Short-run", rep("", times = 7)))
colnames(rown) <- colnames(table.6a)
(table.6 <- rbind(rown[1, ], table.6a, rown[2, ], table.6b))

# 5.4 Alternative specifications
summary(uSta1 <- update(hSta, shift = du2)); aiElas(uSta1)
summary(uDyn1a <- aiDynFit(uSta1)); aiElas(uDyn1a)
summary(uDyn1b <- aiDynFit(uSta1, dum.dif = TRUE))

# -------------------------------------------------------------------------
# 6. Export five tables
# Table in csv format
(output <- listn(table.1, table.2, table.3, table.4, table.5, table.6))
write.list(z = output, file = "OutAidsTable.csv")

# Table in excel format
name <- paste("table", 1:6, sep = ".")
for (i in 1:length(name)) {
  write.xlsx(x = get(name[i]), file = "OutAidsTable.xlsx",
    sheetName = name[i], row.names = FALSE, append = as.logical(i - 1))
}
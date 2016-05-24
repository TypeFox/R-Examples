# 1. Prepare data
library(erer); data(daIns)
daInsH <- daIns[, -14]  # with HuntYrs
daInsF <- daIns[, -3]   # with FishYrs
fm.hunt <- Y ~ Injury + HuntYrs + Nonres + Lspman + Lnong + 
  Gender + Age + Race + Marital + Edu + Inc + TownPop
fm.hunt2 <- bsFormu(name.y = "Y", name.x = colnames(daIns)[-c(1, 14)])
lg <- binomial(link = "logit") 

# 2. Run regressions 
HA <- glm(formula = Y ~ Injury + HuntYrs + Nonres + Lspman + Lnong + 
                    Gender + Age + Race + Marital + Edu + Inc + TownPop,
          family = binomial(link = "logit"), data = daIns, x = TRUE)
HB <- glm(family = lg, data = daIns,  x = TRUE, formula = fm.hunt)
HC <- glm(family = lg, data = daInsH, x = TRUE, formula = Y ~ .)
HD <- glm(family = lg, data = daInsH, x = TRUE, formula = formula(daInsH))
HE <- glm(family = lg, data = daIns,  x = TRUE, formula = formula(daInsH))
FA <- glm(family = lg, data = daInsF, x = TRUE, formula = formula(daInsF))
FB <- update(object = HA, formula = formula(daInsF))

# 3. Understand outputs
class(HA)
head(x = methods(class = "glm"), n = 4)
names(HA); names(summary(HA))
(aic.value <- summary(HA)$aic)

# 4. Extract selected outputs
HH <- summary(HA)$coefficients; head(HH)
FF <- summary(FA)$coefficients

ta <- bsTab(w = HH, need = "2T")
tb <- bsTab(w = FF, need = "2T")
ab <- merge(x = ta, y = tb, by = "Variable", all = TRUE, sort = FALSE)

out <- ab[c(1, 2, 13, 14, 3:12), ]
out[1, 1] <- "Constant"
colnames(out) <- c("Variable", "A_estimate", "A_t.ratio",
  "B_estimate", "B_t.ratio")
out[is.na(out)] <- "___"; out
## ----pkgs, eval = TRUE, echo = TRUE, message = FALSE---------------------
library(popEpi)
library(Epi)
library(survival)

## ------------------------------------------------------------------------

data(sire)
## NOTE: recommended to use factor status variable
x <- Lexis(entry = list(FUT = 0, AGE = dg_age, CAL = get.yrs(dg_date)), 
           exit = list(CAL = get.yrs(ex_date)), 
           data = sire[sire$dg_date < sire$ex_date, ],
           exit.status = factor(status, levels = 0:2, 
                                labels = c("alive", "canD", "othD")), 
           merge = TRUE)

## pretend some are male
set.seed(1L)
x$sex <- rbinom(nrow(x), 1, 0.5)

## observed survival - explicit method
st <- survtab(Surv(time = FUT, event = lex.Xst) ~ sex, data = x, 
              surv.type = "surv.obs",
              breaks = list(FUT = seq(0, 5, 1/12)))

## observed survival - easy method (assumes lex.Xst in x is the status variable)
st <- survtab(FUT ~ sex, data = x, 
              surv.type = "surv.obs",
              breaks = list(FUT = seq(0, 5, 1/12)))

## printing gives the used settings and 
## estimates at the middle and end of the estimated
## curves; more information available using summary()
st


## ------------------------------------------------------------------------
plot(st, col = c("blue", "red"))

## ----popmort-------------------------------------------------------------
data(popmort)
pm <- data.frame(popmort)
names(pm) <- c("sex", "CAL", "AGE", "haz")
head(pm)

## ----survtab_e2----------------------------------------------------------
st.e2 <- survtab(Surv(time = FUT, event = lex.Xst) ~ sex, data = x, 
                 surv.type = "surv.rel", relsurv.method = "e2",
                 breaks = list(FUT = seq(0, 5, 1/12)),
                 pophaz = pm)

## ------------------------------------------------------------------------
plot(st.e2, y = "r.e2", col = c("blue", "red"))

## ----survtab_pp----------------------------------------------------------
st.pp <- survtab(Surv(time = FUT, event = lex.Xst) ~ sex, data = x, 
                 surv.type = "surv.rel", relsurv.method = "pp",
                 breaks = list(FUT = seq(0, 5, 1/12)),
                 pophaz = pm)

## ------------------------------------------------------------------------
plot(st.e2, y = "r.e2", col = c("blue", "red"), lty = 1)
lines(st.pp, y = "r.pp", col = c("blue", "red"), lty = 2)

## ----survtab_adjust------------------------------------------------------
## an age group variable
x$agegr <- cut(x$dg_age, c(0, 60, 70, 80, Inf), right = FALSE)

## using "internal weights" - see ?ICSS for international weights standards
w <- table(x$agegr)
w

w <- list(agegr = as.numeric(w))

## ----survtab_adjust_2----------------------------------------------------
st.as <- survtab(Surv(time = FUT, event = lex.Xst) ~ sex + adjust(agegr), 
                 data = x, weights = w,
                 surv.type = "surv.rel", relsurv.method = "e2",
                 breaks = list(FUT = seq(0, 5, 1/12)),
                 pophaz = pm)

## ------------------------------------------------------------------------
plot(st.as, y = "r.e2.as", col = c("blue", "red"))

## ----weights_examples, eval = TRUE---------------------------------------
list(sex = c(0.4, 0.6), agegr = c(0.2, 0.2, 0.4, 0.2))

wdf <- merge(0:1, 1:4)
names(wdf) <- c("sex", "agegr")
wdf$weights <- c(0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.1, 0.1)
wdf

## ----survtab_adjust_3----------------------------------------------------
st.as <- survtab(Surv(time = FUT, event = lex.Xst) ~ sex, 
                 adjust = "agegr",
                 data = x, weights = w,
                 surv.type = "surv.rel", relsurv.method = "e2",
                 breaks = list(FUT = seq(0, 5, 1/12)),
                 pophaz = pm)

## ----survtab_cause-------------------------------------------------------
st.ca <- survtab(Surv(time = FUT, event = lex.Xst) ~ 1, 
                 data = x, 
                 surv.type = "surv.cause",
                 breaks = list(FUT = seq(0, 5, 1/12)))

st.pp <- survtab(Surv(time = FUT, event = lex.Xst) ~ 1, data = x, 
                 surv.type = "surv.rel", relsurv.method = "pp",
                 breaks = list(FUT = seq(0, 5, 1/12)),
                 pophaz = pm)

plot(st.ca, y = "surv.obs.canD", col = "blue")
lines(st.pp, y = "r.pp", col = "red")

## ----survtab_cif---------------------------------------------------------
st.cif <- survtab(Surv(time = FUT, event = lex.Xst) ~ 1, 
                  data = x, 
                  surv.type = "cif.obs",
                  breaks = list(FUT = seq(0, 5, 1/12)))
plot(st.cif, y = "CIF_canD", conf.int = FALSE)
lines(st.cif, y = "CIF_othD", conf.int = FALSE, col = "red")

## ----survtab_relcif------------------------------------------------------
st.cir <- survtab(Surv(time = FUT, event = lex.Xst) ~ 1, 
                  data = x, 
                  surv.type = "cif.rel",
                  breaks = list(FUT = seq(0, 5, 1/12)),
                  pophaz = pm)
plot(st.cif, y = "CIF_canD", conf.int = FALSE, col = "blue")
lines(st.cir, y = "CIF.rel", conf.int = FALSE, col = "red")

## ------------------------------------------------------------------------
sire$sex <- rbinom(nrow(sire), size = 1, prob = 0.5)
ag <- lexpand(sire, birth = "bi_date", entry = "dg_date", exit = "ex_date",
              status = "status", breaks = list(fot = seq(0, 5, 1/12)), 
              aggre = list(sex, fot))
head(ag)

## ----survtab_ag_example1-------------------------------------------------
st <- survtab_ag(fot ~ sex, data = ag, surv.type = "surv.obs",
                 surv.method = "hazard",
                 d = c("from0to1", "from0to2"), pyrs = "pyrs")

## ----survtab_ag_example2-------------------------------------------------
st <- survtab_ag(fot ~ sex, data = ag, surv.type = "surv.obs",
                 surv.method = "lifetable",
                 d = c("from0to1", "from0to2"), n = "at.risk",
                 n.cens = "from0to0")

## ----survtab_ag_cause----------------------------------------------------
st.ca <- survtab_ag(fot ~ sex, data = ag, surv.type = "surv.cause",
                    surv.method = "hazard",
                    d = list(canD = from0to1, othD = from0to2), pyrs = "pyrs")
plot(st.ca, y = "surv.obs.canD", col = c("blue", "red"))

## ------------------------------------------------------------------------
ag <- lexpand(sire, birth = "bi_date", entry = "dg_date", exit = "ex_date",
              status = "status", breaks = list(fot = seq(0, 5, 1/12)), 
              pophaz = popmort, pp = TRUE,
              aggre = list(sex, fot))

st.pp <- survtab_ag(fot ~ sex, data = ag, surv.type = "surv.rel",
                    surv.method = "hazard", relsurv.method = "pp",
                    d = list(from0to1 + from0to2), pyrs = "pyrs",
                    d.pp = list(from0to1.pp + from0to2.pp),
                    d.pp.2 = list(from0to1.pp.2 + from0to2.pp.2),
                    pyrs.pp = "ptime.pp", d.exp.pp = "d.exp.pp")
plot(st.pp, y = "r.pp", col = c("blue", "red"))


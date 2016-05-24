## Replicate analysis in Leblang (2003)
data("leblang2003")

## NOTE: Convergence tolerance is set to 1e-4 to reduce testing runtime on
## CRAN; do not reduce tolerance in real applications!
m1 <- egame12(outcome ~
              capcont + lreserves + overval + creditgrow + USinterest + service
              + contagion + prioratt - 1 |
              1 |
              1 |
              unifgov + lexports + preelec + postelec + rightgov + realinterest
              + capcont + lreserves,
              data = leblang2003,
              link = "probit",
              type = "private",
              reltol = 1e-4)

summary(m1)

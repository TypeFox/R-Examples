## ----load.data.fsm-------------------------------------------------------
# Load and attach datafsm into your R session, making its functions available:
library(datafsm)

## ----load.printr, include=FALSE------------------------------------------
#library(printr)

## ----cite.datafsm, eval=TRUE, include=TRUE-------------------------------
citation("datafsm")

## ----fake.data, eval=TRUE, include=TRUE----------------------------------
cdata <- data.frame(period = rep(1:10, 1000),
                    outcome = NA,
                    my.decision1 = NA,
                    other.decision1 = NA)
#
# Prisoner's dilemma
#
pd_outcome <- function(player_1, player_2) {
  #
  # 1 = C
  # 2 = D
  #
  player_1  + 1
}

tit_for_tat <- function(last_round_self, last_round_opponent) {
    last_round_opponent
}

noisy_tit_for_tat <- function(last_round_self, last_round_opponent) {
  if (runif(1,0,1) <= 0.1) {
    sample(0:1,1)
  } else {
    last_round_opponent
  }
}

for (i in seq_along(cdata$period)) {
  if (cdata$period[i] == 1) {
    my.decision <- sample(0:1,1)
    other.decision <- sample(0:1,1)
  } else{
    my.last <- my.decision
    other.last <- other.decision
    my.decision <- tit_for_tat(my.last, other.last)
    other.decision <- noisy_tit_for_tat(other.last, my.last)
  }
  cdata[i,-1] <- c(pd_outcome(my.decision, other.decision), my.decision, other.decision)
}

## ----cdata.table, eval=TRUE, echo=FALSE, results='asis'------------------
knitr::kable(head(cdata, 11))

## ----evolve.model, eval=TRUE, include=TRUE, message=FALSE, warning=FALSE, results='hide'----
res <- evolve_model(cdata)

## ----plot.fsm, eval=TRUE, include=TRUE, fig.width=5, fig.height=4, fig.cap="Result of `plot()` method call on `ga_fsm` object."----
summary(res)
plot(res, action_label = c("C", "D"), transition_label = c('cc','cd','dc','dd'))

## ----plot.evolution, eval=TRUE, include=TRUE, fig.width=8, fig.height=6, fig.cap="Result of `plot()` method call on ga object, which is obtained by calling `estimation_details()` on `ga_fsm` object."----
suppressMessages(library(GA))
plot(estimation_details(res))

## ----get.help, eval=TRUE, include=TRUE, render.args = list(help = list(sections = list("usage", "arguments", "details", "references")))----
?evolve_model

## ----session.info, eval=TRUE, include=TRUE-------------------------------
sessionInfo()


## ----size="small"--------------------------------------------------------
library(coala)
list_simulators()

## ------------------------------------------------------------------------
model <- coal_model(10, 1) +
  feat_mutation(5, model = "IFS") +
  sumstat_nucleotide_div()
check_model(model)
model

## ------------------------------------------------------------------------
activate_ms(priority = 500)
model


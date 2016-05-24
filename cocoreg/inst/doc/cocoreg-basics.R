## ---- include = FALSE----------------------------------------------------
library(cocoreg)

## ---- fig.show='hold', fig.height = 5------------------------------------
library(cocoreg)
dc <- create_syn_data_toy()
ccr <- cocoreg(dc$data)
shared.by.all.df <- variation_shared_by(dc, 'all') #only on synthetic datasets

ggplot_dflst(dc$data, ncol = 1)
ggplot_dflst(ccr$data, ncol = 1)

## ---- fig.width = 7.1, fig.height = 7.5----------------------------------
ggplot_dclst(list(observed = dc$data, shared = shared.by.all.df, cocoreg = ccr$data))

## ---- fig.width = 7.1, fig.height = 7.5, message = FALSE-----------------
library(reshape) #importing from namespace does not work as expected
ggcompare_dclst(list(shared = shared.by.all.df, cocoreg = ccr$data))

## ---- fig.width = 7.1, fig.height = 7.5----------------------------------
ggplot_dclst(list(observed = dc$data, shared = shared.by.all.df, cocoreg = ccr$data),
              legendMode = 'all')

## ---- message = FALSE, fig.width = 7.1, fig.height = 5-------------------
ggplot_dflst(dc$data, ncol=1)
ggplot_df(dc$data[[1]])


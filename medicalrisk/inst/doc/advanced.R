## ------------------------------------------------------------------------
library(medicalrisk)
library(plyr)
data(vt_inp_sample)
cm_df <- generate_comorbidity_df(vt_inp_sample, icd9mapfn=icd9cm_charlson_quan)
cci_df <- generate_charlson_index_df(cm_df)
rsi_df <- ddply(vt_inp_sample, .(id), function(x) { icd9cm_sessler_rsi(x$icd9cm) } )
num_icd9_df <- count(vt_inp_sample, c('id'))
num_icd9_df <- rename(num_icd9_df, c("freq" = "num_icd9"))
wide_df <- merge(merge(merge(merge(
  rsi_df, cci_df), 
    cm_df), 
      unique(vt_inp_sample[,c('id','scu_days','drg','mdc')])),
        num_icd9_df)

## ---- results='asis', echo=F---------------------------------------------
library(knitr)
kable(head(wide_df[1:13], n=5), digits=3, table.attr='id="wide_df_table"')
kable(head(wide_df[c(1,14:ncol(wide_df))], n=5), digits=3, table.attr='id="wide_df_table"')

## ----hist, fig.width=10, fig.height=4, cache=F---------------------------
library(reshape2)
library(ggplot2)
# generate a 100 pt x 17 comorbidity table (1700 rows)
cm_melted <- melt(cm_df, id.vars=c('id'), variable.name='cm')
# get rid of all the false entries
cm_melted <- cm_melted[cm_melted$value,]
## count only flags that are true
ggplot(cm_melted, aes(cm, fill=cm)) + 
  geom_bar() + 
  scale_fill_discrete()

## ----hist_chrnlung, fig.width=6, fig.height=4, cache=T-------------------
# make a histogram dataframe for all the icd-9 codes
icd9cm_df <- count(vt_inp_sample, vars='icd9cm')
# create a charlson comorbidity map for all icd-9 codes
icd9cm_charlson_df <- icd9cm_charlson_quan(icd9cm_df$icd9cm)
# isolate just the chrnlung icd_9_cm codes
icd9cm_chrnlung <- row.names(icd9cm_charlson_df[icd9cm_charlson_df$chrnlung,])
# create a hist df
icd9cm_chrnlung_hist <- icd9cm_df[icd9cm_df$icd9cm %in% icd9cm_chrnlung,]
# plot it
ggplot(icd9cm_chrnlung_hist, aes(icd9cm, freq)) + 
  geom_bar(stat="identity")

## ----coincide_icd9-------------------------------------------------------
# create base dataset
pairs <- unique(
  vt_inp_sample[vt_inp_sample$icd9cm %in% icd9cm_chrnlung, c('id','icd9cm')])

# create coincidence matrix
t <- table(
    ddply(pairs, c('id'), function(x) { if (length(x$icd9cm) > 1) {
      data.frame(t(combn(as.character(x$icd9cm),2))) } })[c('X1','X2')])

## ---- results='asis', echo=F---------------------------------------------
kable(t, table.attr='id="chrnlung_coincidence_table"')

## ----coincide_cm---------------------------------------------------------
# create coincidence matrix
t <- table(
    ddply(cm_melted, c('id'), function(x) { if (length(x$cm) > 1) {
      data.frame(t(combn(as.character(x$cm),2))) } })[c('X1','X2')])
# sort it
t <- t[order(rownames(t)),order(colnames(t))]

## ---- results='asis', echo=F---------------------------------------------
kable(t, table.attr='id="cm_coincidence_table"')

## ----fig.width=8, fig.height=4, cache=T----------------------------------
m <- melt(t)
ggplot(m[m$value>0,], aes(X1,X2)) + stat_sum(aes(group=value))

## ----scatterplot, fig.width=10, fig.height=3, cache=T--------------------
library(grid)
library(gridExtra)
p.inhosp <- ggplot(wide_df, aes(rsi_inhosp, index)) + geom_point() + geom_smooth(method=lm) +
  scale_y_continuous(limits=c(-3,10))
p.30dpod <- ggplot(wide_df, aes(rsi_30dpod, index)) + geom_point() + geom_smooth(method=lm) +
  scale_y_continuous(limits=c(-3,10))
p.1yrpod <- ggplot(wide_df, aes(rsi_1yrpod, index)) + geom_point() + geom_smooth(method=lm) +
  scale_y_continuous(limits=c(-3,10))
grid.arrange(p.inhosp, p.30dpod, p.1yrpod, nrow=1)


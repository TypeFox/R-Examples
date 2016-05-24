## ------------------------------------------------------------------------
library(medicalrisk)
library(plyr)
data(vt_inp_sample)
x <- count(vt_inp_sample, c('id'))
cat("average count of ICD codes per patient is: ", mean(x$freq))
y <- count(vt_inp_sample, c('icd9cm'))

## ---- results='asis'-----------------------------------------------------
library(knitr)
kable(head(y[order(-y$freq),], n=5), row.names=F,
      caption='Top 5 most popular ICD-9-CM codes in this dataset')

## ---- results='asis'-----------------------------------------------------
kable(
  icd9cm_charlson_quan(c('D25000')))
kable(
  icd9cm_elixhauser_ahrq37(c('D25000')))
kable(
  icd9cm_rcri(c('D25000')))

## ---- results='asis'-----------------------------------------------------
cases <- vt_inp_sample[vt_inp_sample$id %in% 1:5, c('id','icd9cm')]
cases_with_cm <- merge(cases, icd9cm_charlson_quan(levels(cases$icd9cm)), 
   by.x="icd9cm", by.y="row.names", all.x=TRUE)
 
# generate crude comorbidity summary for each patient
kable(
  ddply(cases_with_cm, .(id), 
   function(x) { data.frame(lapply(x[,3:ncol(x)], any)) }),
  row.names=F)

## ---- results='asis'-----------------------------------------------------
kable(
  generate_comorbidity_df(cases, icd9mapfn=icd9cm_charlson_quan))

## ---- results='asis'-----------------------------------------------------
data(charlson_weights_orig)
kable(
  t(charlson_weights_orig))

## ---- results='asis'-----------------------------------------------------
data(charlson_weights)
kable(
  t(charlson_weights))

## ---- results='asis'-----------------------------------------------------
kable(
  generate_charlson_index_df(generate_comorbidity_df(cases)), row.names=F)

## ------------------------------------------------------------------------
ddply(cases, .(id), function(x) { icd9cm_sessler_rsi(x$icd9cm) } )


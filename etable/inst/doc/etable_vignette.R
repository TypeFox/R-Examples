### R code from vignette source 'etable_vignette.rnw'

###################################################
### code chunk number 1: etable_vignette.rnw:13-17
###################################################
#setwd('D:/WORK/Schulz/R/R_Package/EasyTable-package/Vignette')
options(width=80, continuation="   ")
library(xtable)
library(etable)


###################################################
### code chunk number 2: a
###################################################

set.seed(31415)
age    <- round(runif(10000, 20, 80))
dec    <- cut(age, c(20, 30, 40, 50, 60, 70, 80))
weight <- rnorm(10000, mean=80, sd=10)
height <- rnorm(10000, mean=1.6, sd=0.1)
bmi    <- weight/height^2
bmi_q  <- cut(bmi, quantile(bmi, c(0, 0.25, 0.5, 0.75, 1)))
sex    <- factor(as.factor(rbinom(10000, 1, 0.5)), labels=c('Men', 'Women'))
ethnic <- factor(as.factor(rbinom(10000, 1, 0.75)),labels=c('Other','Caucasian'))
stage  <- as.factor(rbinom(10000, 2, 0.3)+1)
disease<- factor(as.factor(rbinom(10000, 1, 0.1)), labels=c('no','yes'))
treat  <- factor(as.factor(rbinom(10000, 1, 0.2)), labels=c('no','yes'))
ws     <- abs(rnorm(10000))
d<-data.frame(sex,age,dec,ethnic,weight,height,bmi,bmi_q,stage,disease,treat,ws)

print(str(d))


###################################################
### code chunk number 3: etable_vignette.rnw:110-113
###################################################
Tab <- tabular.ade(x_vars='bmi', rows='sex', rnames='Sex', data=d, FUN=iqr_cell)
print(xtable(Tab, caption='Median (Q1/Q3) of BMI'), include.rownames=F,
include.colnames=F, caption.placement="top")


###################################################
### code chunk number 4: etable_vignette.rnw:125-129
###################################################
Tab<-tabular.ade(x_vars='bmi', rows='sex', rnames='Sex', cols='ethnic',
cnames='Ethnicity', data=d, FUN=iqr_cell)
print(xtable(Tab, caption='Median (Q1/Q3) of BMI'), include.rownames=F,
include.colnames=F, caption.placement="top")


###################################################
### code chunk number 5: etable_vignette.rnw:144-148
###################################################
Tab<-tabular.ade(x_vars='bmi', rows=c('sex', 'dec'), rnames=c('Sex', 'Decades'),
cols='ethnic', cnames='Ethnicity', data=d, FUN=iqr_cell)
print(xtable(Tab, caption='Median (Q1/Q3) of BMI'), include.rownames=F,
include.colnames=F, caption.placement="top")


###################################################
### code chunk number 6: etable_vignette.rnw:163-169
###################################################
Tab<-tabular.ade(x_vars='bmi',
rows=c('dec','bmi_q'),   rnames=c('Decades','BMI Quantiles'),
cols=c('sex', 'ethnic'), cnames=c('Sex', 'Ethnicity'), data=d, FUN=n_cell)

print(xtable(Tab, caption='N of Obs.'), include.rownames=F, include.colnames=F,
caption.placement="top")


###################################################
### code chunk number 7: etable_vignette.rnw:184-191
###################################################
Tab<-tabular.ade(x_vars='bmi', xname='BMI',
rows=c('sex','ethnic','disease','treat'),
rnames=c('Sex', 'Ethnicity', 'Disease', 'Treatment'),
data=d, FUN=quantile_cell, probs=0.95)

print(xtable(Tab, caption='95th quantile of BMI'), include.rownames=F,
include.colnames=F,caption.placement="top")


###################################################
### code chunk number 8: etable_vignette.rnw:287-297
###################################################
my_cell<- function(x, y, z, w, cell_ids, row_ids, col_ids, vnames, vars, n_min)
{
out<- format(mean(x[cell_ids], na.rm=TRUE), digits = 3)
return(out)
}

tab<-tabular.ade(x_vars='age', rows='sex', rnames='Sex', cols='dec',
cnames='Decades', data=d, FUN=my_cell)
print(xtable(tab, caption='Mean Age'), include.rownames=F,
include.colnames=F,caption.placement="top")


###################################################
### code chunk number 9: etable_vignette.rnw:324-339
###################################################
my_cell<- function(x, y, z, w, cell_ids, row_ids, col_ids, vnames, vars, n_min)
{
out<- NULL
tab<-table(x[cell_ids])
for(i in 1:length(tab)){
out<- paste(out, levels(x)[i],': ' ,tab[i], sep='')
if(i<length(tab)) out<- paste(out, ', ', sep='')
}
return(out)
}

tab<-tabular.ade(x_vars='sex', rows='dec', rnames='Decades',
cols='stage', cnames='Stage', data=d, FUN=my_cell)
print(xtable(tab, caption='Frequencies'), include.rownames=F,
include.colnames=F,caption.placement="top")


###################################################
### code chunk number 10: etable_vignette.rnw:379-405
###################################################
b_cell<- function(x, y, z, w, cell_ids, row_ids, col_ids, vnames, vars, n_min)
{
out<- NULL
if(length(unique(x))==2){
lv<-levels(x)
n <-sum(x[cell_ids]==lv[2])
N <-sum(table(x[cell_ids]))
out<-paste(levels(x)[2], ': ',format((n/N)*100, digits=3),'% (N:',n , ')',sep='')
}
if(!is.factor(x) & length(unique(x))> 2){
quant <- format(quantile(x[cell_ids], c(0.25, 0.5, 0.75), na.rm=TRUE), digits=3)
out<- paste(quant[1], ' (',quant[2],'/',quant[3],')', sep='')
}
if(is.factor(x) & length(unique(x))> 2){
lv<-levels(x)
n <-table(x[cell_ids])
N <-sum(table(x[cell_ids]))
out<- paste(lv, ': ', format((n/N)*100,  digits=3), '%', collapse=' | ', sep='')
}
return(out)
}

tab<-tabular.ade(x_vars=c('bmi','ethnic','stage'),xname=c('BMI','Ethnicity','Stages'),
cols='sex', cnames='Sex', data=d, FUN=b_cell)
print(xtable(tab, caption='Diverse variables'), include.rownames=F,
include.colnames=F,caption.placement="top")


###################################################
### code chunk number 11: etable_vignette.rnw:455-470
###################################################
t_test_cell<- function(x, y, z, w, cell_ids, row_ids, col_ids, vnames, vars, n_min)
{
v <- x[cell_ids]
group <- y[cell_ids]
test<-t.test(v[which(group==levels(group)[1])], v[which(group==levels(group)[2])])
mdiff<- format(diff(test$estimate), digits=3)
p<- base:::format.pval(test$p.value, digits=2, eps=0.0001)
out<- paste('Diff: ', mdiff, ', p-value: ', p, sep='')
return(out)
}

tab<-tabular.ade(x_vars='bmi', xname='BMI', y_vars='ethnic', yname='Ethnicity',
rows='dec', rnames='Decades', cols='sex', cnames='Sex', data=d, FUN=t_test_cell)
print(xtable(tab, caption='T test for BMI means between Ethnicity'),
include.rownames=F, include.colnames=F,caption.placement="top")


###################################################
### code chunk number 12: etable_vignette.rnw:483-490
###################################################
vars    <-c('age', 'weight', 'height', 'bmi')
vlabels <-c('Age', 'Weight', 'Height', 'BMI')

tab<-tabular.ade(x_vars=vars, xname=vlabels, y_vars=vars, yname=vlabels, data=d,
FUN=corr_p_cell, digits=2)
print(xtable(tab, caption='Pearson correlation'), include.rownames=F,
 include.colnames=F,caption.placement="top")


###################################################
### code chunk number 13: etable_vignette.rnw:503-510
###################################################
vars    <-c('age', 'weight', 'height', 'bmi')
vlabels <-c('Age', 'Weight', 'Height', 'BMI')

tab<-tabular.ade(x_vars=vars, xname=vlabels,
cols=c('sex','stage'), cnames=c('Sex','Stage'), data=d, FUN=quantile_cell)
print(xtable(tab, caption='Medians'), include.rownames=F,
include.colnames=F,caption.placement="top")


###################################################
### code chunk number 14: etable_vignette.rnw:554-559
###################################################

tab<-tabular.ade(x_vars='sex', rows=c('treat', 'ALL'), rnames=c('Treatment'),
cols=c('disease', 'ALL'), cnames=c('Disease'), data=d, FUN=n_cell, alllabel='both')
print(xtable(tab, caption='Contingency table'), include.rownames=F,
include.colnames=F,caption.placement="top")


###################################################
### code chunk number 15: etable_vignette.rnw:576-581
###################################################

tab<-tabular.ade(x_vars='sex', rows=c('sex', 'ALL', 'ethnic', 'stage'),
rnames=c('Sex','Ethnicity', 'Stage'), w='ws', data=d, FUN=n_cell, digits=1)
print(xtable(tab, caption='weighted N'), include.rownames=F,
include.colnames=F,caption.placement="top")


###################################################
### code chunk number 16: etable_vignette.rnw:608-619
###################################################
vars    <-c('age', 'weight', 'height', 'bmi')
vlabels <-c('Age', 'Weight', 'Height', 'BMI')

keywords  <-c('MIN', 'MAX', 'MEAN', 'SD', 'CV', 'SKEW',     'KURT')
keylabels <-c('Min', 'Max', 'Mean', 'SD', 'CV', 'Skewness', 'Kurtosis')


tab<-tabular.ade(x_vars=vars, xname=vlabels,
y_vars=keywords, yname=keylabels, data=d, FUN=stat_cell)
print(xtable(tab, caption='Various statistics'), include.rownames=F,
 include.colnames=F,caption.placement="top")


###################################################
### code chunk number 17: etable_vignette.rnw:635-642
###################################################
keywords  <-c('N', 'MIN', 'MAX', 'MEAN', 'SD')
keylabels <-c('N', 'Min', 'Max', 'Mean', 'SD')

tab<-tabular.ade(x_vars=vars, xname=vlabels, y_vars=keywords, yname=keylabels,
rows=c('sex','ALL','ethnic'), rnames=c('Sex','Ethnicity'), data=d, FUN=stat_cell)
print(xtable(tab, caption='Various statistics'), include.rownames=F,
 include.colnames=F,caption.placement="top")


###################################################
### code chunk number 18: etable_vignette.rnw:659-667
###################################################
keywords  <-c('N', 'MIN', 'MAX', 'MEAN', 'SD')
keylabels <-c('N', 'Min', 'Max', 'Mean', 'SD')


tab<-tabular.ade(x_vars=keywords, xname=keylabels, y_vars=vars, yname=vlabels,
rows=c('sex', 'ALL'), rnames=c('Sex'),data=d, FUN=stat_cell)
print(xtable(tab, caption='Various statistics'), include.rownames=F,
 include.colnames=F,caption.placement="top")


###################################################
### code chunk number 19: etable_vignette.rnw:696-707
###################################################
vars    <-c('age', 'weight', 'height', 'bmi')
vlabels <-c('Age', 'Weight', 'Height', 'BMI')

keywords  <-c('N', 'MEDIAN', 'IQR')
keylabels <-c('N', 'Median', 'IQR')

tab<-tabular.ade(x_vars=vars, xname=vlabels, y_vars=keywords, yname=keylabels,
rows=c('sex', 'ALL'), rnames=c('Sex'),cols=c('ethnic'),cnames=c('Ethnicity'),
w='ws',data=d,FUN=stat_cell)
print(xtable(tab, caption='Various statistics'), include.rownames=F,
include.colnames=F,caption.placement="top")



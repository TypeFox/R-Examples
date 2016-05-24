## ----lib, results="hide"-------------------------------------------------
## Loading needed libraries
library(funModeling)
data(heart_disease)

## ----df_status-----------------------------------------------------------
my_data_status=df_status(heart_disease)

## ----df_status3----------------------------------------------------------
# Removing variables with 60% of zero values
vars_to_remove=subset(my_data_status, my_data_status$p_zeros > 60)
vars_to_remove["variable"]

## Keeping all except vars_to_remove 
heart_disease_2=heart_disease[, !(names(heart_disease) %in% vars_to_remove[,"variable"])]


## ----df_status4----------------------------------------------------------
my_data_status[order(-my_data_status$p_zeros),]

## ----variable_importance1, results="hide", fig.height=4, fig.width=8-----
cross_gender=cross_plot(heart_disease, str_input="gender", str_target="has_heart_disease")

## ----variable_importance2, results="hide", fig.height=4, fig.width=12----
cross_plot(heart_disease, str_input="max_heart_rate", str_target="has_heart_disease")

## ----variable_importance3------------------------------------------------
heart_disease$oldpeak_2=equal_freq(var=heart_disease$oldpeak, n_bins = 3)
summary(heart_disease$oldpeak_2)

## ----variable_importance4, results="hide", fig.height=4, fig.width=8-----
cross_oldpeak_2=cross_plot(heart_disease, str_input="oldpeak_2", str_target="has_heart_disease", auto_binning = F)

## ----variable_importance5, results="hide", fig.height=4, fig.width=12----
heart_disease$max_heart_rate_2=equal_freq(var=heart_disease$max_heart_rate, n_bins = 10)
cross_plot(heart_disease, str_input="max_heart_rate_2", str_target="has_heart_disease")

## ----variable_importance6, results="hide", fig.height=4, fig.width=10----
heart_disease$max_heart_rate_3=equal_freq(var=heart_disease$max_heart_rate, n_bins = 5)
cross_plot(heart_disease, str_input="max_heart_rate_3", str_target="has_heart_disease")

## ----several_cross_plot1, eval=FALSE-------------------------------------
#  cross_plot(heart_disease, str_input="max_heart_rate_3", str_target="has_heart_disease", path_out="my_plots")

## ----several_cross_plot2, eval=FALSE-------------------------------------
#  vars_to_analyze=c("age", "oldpeak", "max_heart_rate")

## ----several_cross_plot3, eval=FALSE-------------------------------------
#  cross_plot(data=heart_disease, str_target="has_heart_disease", str_input=vars_to_analyze)

## ----model_perfomance1---------------------------------------------------
## Training and test data. Percentage of training cases default value=80%.
index_sample=get_sample(data=heart_disease, percentage_tr_rows=0.8)

## Generating the samples
data_tr=heart_disease[index_sample,] 
data_ts=heart_disease[-index_sample,]


## Creating the model only with training data
fit_glm=glm(has_heart_disease ~ age + oldpeak, data=data_tr, family = binomial)


## ----model_perfomance2,  fig.height=3, fig.width=4-----------------------
## Performance metrics for Training Data
model_performance(fit=fit_glm, data = data_tr, target_var = "has_heart_disease")

## Performance metrics for Test Data
model_performance(fit=fit_glm, data = data_ts, target_var = "has_heart_disease")


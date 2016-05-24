## ----load_modellingTools, results = "hide", warning = FALSE, message = FALSE----
library(modellingTools)

## ----install_otherpackages, results = "hide", warning = FALSE, message = FALSE----
# install.packages(c("dplyr","magrittr"))
library(dplyr)
library(magrittr)

## ----iris_tbl------------------------------------------------------------
d <- tbl_df(iris)
d[,1]
d[,c(1,2)]

## ----using_table---------------------------------------------------------
data(CO2)
CO2 %$% table(Plant)

## ----using_table2--------------------------------------------------------
CO2 %$% table(Plant,Type)

## ----basic_procfreq------------------------------------------------------
proc_freq(iris,"Species")

## ----basic_procfreq2-----------------------------------------------------
proc_freq(iris,"Sepal.Length")

## ----create_large_dframe-------------------------------------------------
library(foreach) # Check this out if you haven't already
# Try for 1,000 observations, for illustration purposes
d <- data_frame(v1 = times(1e04) %do% if (runif(1) < .1) NA else rnorm(1))

## ----procfreq_large_dframe-----------------------------------------------
proc_freq(d,"v1",bins = 3)

## ----mtcars_cor----------------------------------------------------------
x <- mtcars[,-1]
y <- mtcars[,1]
apply(x,2,function(a,b) cor(a,b),b = y)

## ----mtcars_gettopcorrs--------------------------------------------------
get_top_corrs(mtcars,"mpg")

## ----seatbelts_simplebin-------------------------------------------------
data(Seatbelts)
d <- tbl_df(as.data.frame(Seatbelts))
d

d_bin <- simple_bin(d,bins = 3)
d_bin

## ----traintest_simplebin, warnings = FALSE-------------------------------
d %<>% mutate(on_train = times(n()) %do% if (runif(1) < .7) 1 else 0)
d_train <- d %>% filter(on_train == 1) %>% select(-on_train)
d_test <- d %>% filter(on_train == 0) %>% select(-on_train)

d_split_binned <- simple_bin(d_train,test = d_test,bins = 3)
d_split_binned

## ------------------------------------------------------------------------
d_train_bin <- simple_bin(d_train,bins = 3)
train_cutpoints <- binned_data_cutpoints(d_train_bin)
train_cutpoints

d_test_bin <- simple_bin(d_test,bins = train_cutpoints)
d_test_bin

identical(d_test_bin,d_split_binned$test)

## ----vectorbin-----------------------------------------------------------
x <- rnorm(100)
vector_bin(x,bins = 3)
vector_bin(x,bins = 3,type = "width")
vector_bin(x,bins = c(-1,0,1))

## ----levels_x------------------------------------------------------------
x <- vector_bin(rnorm(100),bins = 3)
levels(x)

## ----vectorbin_x---------------------------------------------------------
get_vector_cutpoints(x)

## ----binneddatacutpoints-------------------------------------------------
binned_data_cutpoints(d_bin)

## ------------------------------------------------------------------------
# Print this in your own R console:
x <- create_model_matrix(d_train_bin)

# Prettier output:
create_model_matrix(d_train_bin,matrix_out = FALSE)


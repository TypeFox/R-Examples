## ----get-----------------------------------------------------------------

library(PP)

data(fourpl_df)

dim(fourpl_df)

head(fourpl_df)

diff_par <- attr(fourpl_df,"diffpar")
slope_par <- attr(fourpl_df,"slopes")


## ----check---------------------------------------------------------------

# extract items and transform the data.frame to matrix
itmat <- as.matrix(fourpl_df[,-(1:2)])


# are there any full scores?
fullsc <- apply(itmat,1,function(x) (sum(x,na.rm=TRUE)+sum(is.na(x))) == length(x))
any(fullsc)

# are there 0 scores?
nullsc <- apply(itmat,1,function(x) sum(x,na.rm=TRUE) == 0)
any(nullsc)

# are there missing values? how many and where?
nasc <- apply(itmat,1,function(x) sum(is.na(x)))
any(nasc > 0)
#which(nasc)

# is our data dichotomous as expected?
apply(itmat,2,function(x) table(x))

# are there any duplicates?
rdup <- duplicated(itmat)
sum(rdup)



## ----decide--------------------------------------------------------------
library(PP)

res1plmle <- PP_4pl(respm = itmat,thres = diff_par, slopes = slope_par, type = "mle")

summary(res1plmle)


## ----edit----------------------------------------------------------------

dafest <- data.frame(fourpl_df,res1plmle$resPP$resPP)

head(dafest,10)

## ----rerun---------------------------------------------------------------
library(PP)

res1plmle <- PP_4pl(respm = itmat,thres = diff_par, slopes = slope_par, type = "wle")

summary(res1plmle)


## ----setup, include = FALSE----------------------------------------------
library(qlcData)

## ----nominalData---------------------------------------------------------
data <- data.frame(
    size = c('large','very small','small','large','very small','tiny'),
    kind = c('young male','young female','old female','old male',NA,'young female'),
    row.names = paste('obs', 1:6, sep='')
    )
data

## ----merge---------------------------------------------------------------
# Specifying a recoding
recoding <- list(
  list(
    recodingOf = 'size',
    attribute = 'newSize',
    values = c('large','small'),
    link = c(1,2,2,2)
    )
  )
# Do the actual recoding and see the results
recode(data, recoding)

## ----mergesplit----------------------------------------------------------
# Specifying the recoding of three different new attributes
recoding <- list(
  list(
    recodingOf = 'size',
    attribute = 'size',
    values = c('large','small'),
    link = c(1,2,2,2)
    ),
  list(
    recodingOf = 'kind',
    attribute = 'gender',
    values = c('female','male'),
    link = c(1,2,1,2)
    ),
  list(
    recodingOf = 'kind',
    attribute = 'age',
    values = c('old','young'),
    link = c(1,1,2,2)
    )
  )
# Do the recoding and show it
newdata <- recode(data, recoding)
newdata

## ----mergeattributes-----------------------------------------------------
back_recoding <- list(
  list(
    recodingOf = c('size','age'),
    attribute = 'size+age',
    values = c('large+old','large+young','small+old','small+young'),
    link = c(1,3,0,2,4,0,0,0,0,0)
    )
  )
recode(newdata, back_recoding)

## ----expandgrid----------------------------------------------------------
expand.grid(c(levels(newdata$size),NA),c(levels(newdata$age),NA))

## ----yamloutput, echo = FALSE--------------------------------------------
cat(yaml::as.yaml(write.recoding(list(1,c(1,2)),data,yaml=F)))

## ----shortcuts-----------------------------------------------------------
short_recoding <- list(
  # same as first example at the start of this vignette
  # using abbreviations and a different order
  list(
    r = 'size',
    a = 'newSize',
    l = c(1,2,2,2),
    v = c('large','small')
    ),
  # same new attribute, but with automatically generated names
  list(
    r = 'size',
    l = c(1,2,2,2)
    ),
  # keep original attribute in column 2 of the data
  list(
    r = 2
    ),
  # add three times the first original attribute
  # senseless, but it illustrates the possibilities
  list(
    d = c(1,1,1)
    )
  )
recode(data, short_recoding)

## ----expand_shortcuts_notrun, eval = FALSE-------------------------------
#  read.recoding(short_recoding, file = yourFile , data =  data)

## ----expand_shortcuts, echo = FALSE--------------------------------------
meta <- list(
  title = NULL,
  author = NULL,
  date = format(Sys.time(),"%Y-%m-%d")
  )
recoding <- read.recoding(short_recoding, file = NULL , data =  data)
outfile <- c(meta, list(recoding = recoding))
cat(yaml::as.yaml(outfile))


## ----echo=FALSE,message=FALSE--------------------------------------------
knitr::opts_chunk$set(comment=NA,
               fig.align="center",
               results="markup")

## ------------------------------------------------------------------------
knit_print.codebook <-function(x,...) 
  knitr::asis_output(format_html(x,...))

knit_print.descriptions <-function(x,...) 
  knitr::asis_output(format_html(x,...))

knit_print.ftable <-function(x,options,...)
  knitr::asis_output(
    format_html(x,
                digits=if(length(options$ftable.digits))
                          options$ftable.digits
                       else 0,
                ...))
# We can now adjust the number of digits after the comma
# for each column e.g. by adding an `ftable.digits` option
# to an R chunk, as in ```{r,ftable=c(2,2,0)}

knit_print.mtable <-function(x,...)
  knitr::asis_output(format_html(x,...))

## ---- message=FALSE------------------------------------------------------
library(memisc)
options(digits=3)
nes1948.por <- unzip(system.file("anes/NES1948.ZIP",package="memisc"),
                     "NES1948.POR",exdir=tempfile())

## ------------------------------------------------------------------------
nes1948 <- spss.portable.file(nes1948.por)
print(nes1948)

## ------------------------------------------------------------------------
names(nes1948)

## ------------------------------------------------------------------------
description(nes1948)

## ---- eval=FALSE---------------------------------------------------------
#  codebook(nes1948)

## ------------------------------------------------------------------------
codebook(nes1948[1:5])

## ------------------------------------------------------------------------
vote.48 <- subset(nes1948,
              select=c(
                  v480018,
                  v480029,
                  v480030,
                  v480045,
                  v480046,
                  v480047,
                  v480048,
                  v480049,
                  v480050
                  ))

## ------------------------------------------------------------------------
str(vote.48)

## ------------------------------------------------------------------------
vote.48 <- rename(vote.48,
                  v480018 = "vote",
                  v480029 = "occupation.hh",
                  v480030 = "unionized.hh",
                  v480045 = "gender",
                  v480046 = "race",
                  v480047 = "age",
                  v480048 = "education",
                  v480049 = "total.income",
                  v480050 = "religious.pref"
        )

## ----eval=FALSE----------------------------------------------------------
#  vote.48 <- subset(nes1948,
#                    select=c(
#                      vote           = v480018,
#                      occupation.hh  = v480029,
#                      unionized.hh   = v480030,
#                      gender         = v480045,
#                      race           = v480046,
#                      age            = v480047,
#                      education      = v480048,
#                      total.income   = v480049,
#                      religious.pref = v480050
#                    ))

## ------------------------------------------------------------------------
codebook(vote.48)

## ------------------------------------------------------------------------
vote.48 <- within(vote.48,{
  vote3 <- recode(vote,
    1 -> "Truman",
    2 -> "Dewey",
    3:4 -> "Other"
    )
  occup4 <- recode(occupation.hh,
    10:20 -> "Upper white collar",
    30 -> "Other white collar",
    40:70 -> "Blue collar",
    80 -> "Farmer"
    )
  relig3 <- recode(religious.pref,
    1 -> "Protestant",
    2 -> "Catholic",
    3:5 -> "Other,none"
    )
   race2 <- recode(race,
    1 -> "White",
    2 -> "Black"
    )
  })

## ------------------------------------------------------------------------
ftable(xtabs(~vote3+occup4,data=vote.48))

## ---- ftable.digits=c(2,2,2,0)-------------------------------------------
gt1 <- genTable(percent(vote3)~occup4,data=vote.48)
## For knitr-ing, we use ```{r, ftable.digits=c(2,2,2,0)} here.
ftable(gt1,row.vars=2)

## ---- ftable.digits=c(2,2,2,0)-------------------------------------------
gt2 <- genTable(percent(vote3)~relig3,data=vote.48)
ftable(gt2,row.vars=2)

## ---- ftable.digits=c(2,2,2,0)-------------------------------------------
gt3 <- genTable(percent(vote3)~race2,data=vote.48)
ftable(gt3,row.vars=2)

## ---- ftable.digits=c(2,2,2,0)-------------------------------------------
gt4 <- genTable(percent(vote3)~total.income,data=vote.48)
ftable(gt4,row.vars=2)

## ---- ftable.digits=c(2,2,2)---------------------------------------------
## For knitr-ing, we use ```{r, ftable.digits=c(2,2,2)} here.
inc.tab <- genTable(percent(vote3,ci=TRUE)~total.income,data=vote.48)
ftable(inc.tab,row.vars=c(3,2))

## ---- ftable.digits=c(2,2,2)---------------------------------------------
occup.tab <- genTable(percent(vote3,ci=TRUE)~occup4,data=vote.48)
ftable(occup.tab,row.vars=c(3,2))

## ------------------------------------------------------------------------
vote.48 <- within(vote.48,{
  contrasts(occup4) <- contr("treatment",base = 3)
  contrasts(total.income) <- contr("treatment",base = 4)
  })

## ------------------------------------------------------------------------
model1 <- glm((vote3=="Truman")~occup4,data=vote.48,
              family="binomial")
model2 <- glm((vote3=="Truman")~total.income,data=vote.48,
              family="binomial")
model3 <- glm((vote3=="Truman")~occup4+total.income,data=vote.48,
              family="binomial")
model4 <- glm((vote3=="Truman")~relig3,data=vote.48,
              family="binomial")
model5 <- glm((vote3=="Truman")~occup4+relig3,data=vote.48,
              family="binomial")

## ------------------------------------------------------------------------
mtable(model1,model2,model3,summary.stats=c("Nagelkerke R-sq.","Deviance","AIC","N"))

## ------------------------------------------------------------------------
relabel(mtable(
            "Model 1"=model1,
            "Model 2"=model2,
            "Model 3"=model3,
            summary.stats=c("Nagelkerke R-sq.","Deviance","AIC","N")),
          UNDER="under",
          "AND OVER"="and over",
          occup4="Occup. class",
          total.income="Income",
          gsub=TRUE
          )

## ------------------------------------------------------------------------
relabel(mtable(
              "Model 1"=model1,
              "Model 4"=model4,
              "Model 5"=model5,
              summary.stats=c("Nagelkerke R-sq.","Deviance","AIC","N")),
            occup4="Occup. class",
            relig3="Religion",
            gsub=TRUE
            )

## ----echo=FALSE----------------------------------------------------------
rm(knit_print.codebook,
   knit_print.descriptions,
   knit_print.ftable,
   knit_print.mtable)


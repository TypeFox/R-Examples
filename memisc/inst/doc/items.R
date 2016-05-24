## ----echo=FALSE,message=FALSE--------------------------------------------
knitr::opts_chunk$set(comment=NA,
               fig.align="center",
               results="markup")

## ----echo=FALSE----------------------------------------------------------
set.seed(20)
vote.probs <- c(.361,.29,.23,1-.361-.29-.23)
decis.probs <- c(vote.probs*.651,1-.651)
all.probs <- c(.9*decis.probs,.1*rep(1/3,3))
voteint <- sample(c(1:4,9,97,98,99),
                  prob=all.probs,
                  replace=TRUE,
                  size=200)                
library(memisc)

## ------------------------------------------------------------------------
voteint

## ------------------------------------------------------------------------
# This is to be run *after* memisc has been loaded.
labels(voteint) <- c(Conservative       =  1,
                     Labour             =  2,
                     "Liberal Democrat" =  3, # We have whitespace in the label, 
                     "Other Party"      =  4, # so we need quotation marks
                     "Will not vote"    =  9,
                     "Don't know"       = 97,
                     "Answer refused"   = 98,
                     "Not applicable"   = 99)

## ------------------------------------------------------------------------
class(voteint)
str(voteint)
voteint

## ------------------------------------------------------------------------
labels(voteint)

## ------------------------------------------------------------------------
voteint <- relabel(voteint,
                   "Conservative"     = "Cons",
                   "Labour"           = "Lab",
                   "Liberal Democrat" = "LibDem",
                   "Other Party"      = "Other",
                   "Will not vote"    = "NoVote",
                   "Don't know"       = "DK",
                   "Answer refused"   = "Refused",
                   "Not applicable"   = "N.a.")

## ------------------------------------------------------------------------
labels(voteint)
voteint
str(voteint)

## ------------------------------------------------------------------------
missing.values(voteint) <- c(97,98,99)

## ------------------------------------------------------------------------
voteint

## ------------------------------------------------------------------------
missing.values(voteint) <- missing.values(voteint) + 9

## ------------------------------------------------------------------------
missing.values(voteint)

## ------------------------------------------------------------------------
as.numeric(voteint)[1:30]
as.factor(voteint)[1:30]

## ------------------------------------------------------------------------
missing.values(voteint) <- NULL
missing.values(voteint)
as.numeric(voteint)[1:30]

## ------------------------------------------------------------------------
valid.values(voteint) <- 1:4
valid.values(voteint)
missing.values(voteint)

## ------------------------------------------------------------------------
valid.range(voteint) <- c(1,9)
missing.values(voteint)

## ------------------------------------------------------------------------
description(voteint) <- "Vote intention"
description(voteint)

## ------------------------------------------------------------------------
wording(voteint) <- "Which party are you going to vote for in the general election next Tuesday?"
wording(voteint)
annotation(voteint)
annotation(voteint)["wording"]

## ------------------------------------------------------------------------
codebook(voteint)

## ------------------------------------------------------------------------
voteint1 <- voteint
voteint1[sample(length(voteint),size=20)] <- c(rep(5,13),rep(7,7))

## ------------------------------------------------------------------------
codebook(voteint1)

## ------------------------------------------------------------------------
wild.codes(voteint1)

## ------------------------------------------------------------------------
voteint2 <- voteint
labels(voteint2) <- NULL # This deletes all value labels
codebook(voteint2)

## ------------------------------------------------------------------------
measurement(voteint2) <- "interval"
codebook(voteint2)

## ----results='asis'------------------------------------------------------
show_html(codebook(voteint))

## ------------------------------------------------------------------------
Data <- data.set(
          vote = sample(c(1,2,3,4,8,9,97,99),
                        size=300,replace=TRUE),
          region = sample(c(rep(1,3),rep(2,2),3,99),
                          size=300,replace=TRUE),
          income = round(exp(rnorm(300,sd=.7))*2000)
          )

## ------------------------------------------------------------------------
Data

## ------------------------------------------------------------------------
options(show.max.obs=5)
Data
# Back to the default
options(show.max.obs=25)

## ----eval=FALSE----------------------------------------------------------
#  print(Data)

## ------------------------------------------------------------------------
Data <- within(Data,{
  description(vote) <- "Vote intention"
  description(region) <- "Region of residence"
  description(income) <- "Household income"
  wording(vote) <- "If a general election would take place next Tuesday,
                    the candidate of which party would you vote for?"
  wording(income) <- "All things taken into account, how much do all
                    household members earn in sum?"
  foreach(x=c(vote,region),{
    measurement(x) <- "nominal"
    })
  measurement(income) <- "ratio"
  labels(vote) <- c(
                    Conservatives         =  1,
                    Labour                =  2,
                    "Liberal Democrats"   =  3,
                    "Other"               =  4,
                    "Don't know"          =  8,
                    "Answer refused"      =  9,
                    "Not applicable"      = 97,
                    "Not asked in survey" = 99)
  labels(region) <- c(
                    England               =  1,
                    Scotland              =  2,
                    Wales                 =  3,
                    "Not applicable"      = 97,
                    "Not asked in survey" = 99)
  foreach(x=c(vote,region,income),{
    annotation(x)["Remark"] <- "This is not a real survey item, of course ..."
    })
  missing.values(vote) <- c(8,9,97,99)
  missing.values(region) <- c(97,99)

  # These to variables do not appear in the
  # the resulting data set, since they have the wrong length.
  junk1 <- 1:5
  junk2 <- matrix(5,4,4)
  
})

## ------------------------------------------------------------------------
Data

## ------------------------------------------------------------------------
EnglandData <- subset(Data,region == "England")
EnglandData

## ------------------------------------------------------------------------
codebook(Data)

## ----results='asis'------------------------------------------------------
show_html(codebook(Data))

## ------------------------------------------------------------------------
DataFr <- as.data.frame(Data)
## Looking a the data frame structure
str(DataFr)
## Looking at the first 25 observations
DataFr[1:25,]

## ------------------------------------------------------------------------
xtabs(~vote+region,data=DataFr)

## ------------------------------------------------------------------------
xtabs(~vote+region,data=Data)

## ------------------------------------------------------------------------
xtabs(~vote+region,data=within(Data, 
                               vote <- include.missings(vote)))

## ----results='asis'------------------------------------------------------
show_html(codebook(DataFr))

## ------------------------------------------------------------------------
load(system.file("gles/gles2013work.RData",package="memisc"))

## ----results='asis'------------------------------------------------------
with(gles2013work,
     show_html(codebook(bula)))

## ------------------------------------------------------------------------
gles2013work <- within(gles2013work,
                       east.west <- recode(bula,
                                          East = 1 <- c(3,4,8,13,14,16),
                                          West = 2 <- c(1,2,5:7,9:12,15)
                                          ))

## ------------------------------------------------------------------------
xtabs(~bula+east.west,data=gles2013work)

## ------------------------------------------------------------------------
x <- 1:10
xc <- cases(x <= 3, 
            x > 3 & x <= 7, 
            x > 7)
data.frame(x,xc)

## ------------------------------------------------------------------------
xn <- cases(1 <- x <= 3, 
            2 <- x > 3 & x <= 7, 
            3 <- x > 7)
data.frame(x,xn)

## ------------------------------------------------------------------------
gles2013work <- within(gles2013work,{

  candidate.vote <- cases(
              wave == 1 & intent.turnout == 6 -> postal.vote.candidate,
              wave == 1 & intent.turnout %in% 4:5 -> 900,
              wave == 1 & intent.turnout %in% 1:3 -> voteint.candidate,
              wave == 2 & turnout == 1 -> vote.candidate,
              wave == 2 & turnout == 2 -> 900
            )

  list.vote <- cases(
              wave == 1 & intent.turnout == 6 -> postal.vote.list,
              wave == 1 & intent.turnout %in% 4:5 -> 900,
              wave == 1 & intent.turnout %in% 1:3 -> voteint.list,
              wave == 2 & turnout ==1 -> vote.list,
              wave == 2 & turnout ==2 -> 900
            )
})

## ------------------------------------------------------------------------
gles2013work <- within(gles2013work,{
  candidate.vote <- recode(as.item(candidate.vote),
                      "CDU/CSU"   =  1 <- 1,
                      "SPD"       =  2 <- 4,
                      "FDP"       =  3 <- 5,
                      "Grüne"     =  4 <- 6,
                      "Linke"     =  5 <- 7,
                      "NPD"       =  6 <- 206,
                      "Piraten"   =  7 <- 215,
                      "AfD"       =  8 <- 322,
                      "Other"     = 10 <- 801,
                      "No Vote"   = 90 <- 900,
                      "WN"        = 98 <- -98,
                      "KA"        = 99 <- -99
                  )
  list.vote <- recode(as.item(list.vote),
                      "CDU/CSU"   =  1 <- 1,
                      "SPD"       =  2 <- 4,
                      "FDP"       =  3 <- 5,
                      "Grüne"     =  4 <- 6,
                      "Linke"     =  5 <- 7,
                      "NPD"       =  6 <- 206,
                      "Piraten"   =  7 <- 215,
                      "AfD"       =  8 <- 322,
                      "Other"     = 10 <- 801,
                      "No Vote"   = 90 <- 900,
                      "WN"        = 98 <- -98,
                      "KA"        = 99 <- -99
                  )
  
   missing.values(candidate.vote) <- 98:99
   missing.values(list.vote) <- 98:99
   measurement(candidate.vote) <- "nominal"
   measurement(list.vote) <- "nominal"
})

## ----width=120-----------------------------------------------------------
xtabs(~list.vote+east.west,data=gles2013work)
xtabs(~list.vote+candidate.vote,data=gles2013work)


## ----setup, include=FALSE------------------------------------------------
# library(printr)
# knitr::opts_chunk$set(cache=TRUE)
# knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_github("nikita-moor/ldatuning")

## ---- message=FALSE------------------------------------------------------
library("ldatuning")

## ---- message=FALSE------------------------------------------------------
library("topicmodels")
data("AssociatedPress", package="topicmodels")
dtm <- AssociatedPress[1:10, ]

## ---- eval=TRUE----------------------------------------------------------
result <- FindTopicsNumber(
  dtm,
  topics = seq(from = 2, to = 15, by = 1),
  metrics = c("Griffiths2004", "CaoJuan2009", "Arun2010", "Deveaud2014"),
  method = "Gibbs",
  control = list(seed = 77),
  mc.cores = 2L,
  verbose = TRUE
)

## ---- echo=FALSE---------------------------------------------------------
knitr::kable(result)

## ---- fig.width=6, fig.height=3, results="hide"--------------------------
FindTopicsNumber_plot(result)

## ---- fig.width=9, fig.height=5, echo=FALSE------------------------------
result <- read.csv(file = "files/APress.csv", header = TRUE)
FindTopicsNumber_plot(result[result$topics < 500, ])


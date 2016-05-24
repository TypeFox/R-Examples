## ----echo=FALSE,message=FALSE--------------------------------------------
knitr::opts_chunk$set(comment=NA,
               fig.align="center",
               results="markup")

## ------------------------------------------------------------------------
tab.Class.Age <- xtabs(Freq~Class+Age,data=Titanic)
tab.Survived.Class.Age <- xtabs(Freq~Survived+Class+Age,data=Titanic)
tab.Survived.Class.Sex <- xtabs(Freq~Survived+Class+Sex,data=Titanic)
tab.Survived.Class <- xtabs(Freq~Survived+Class,data=Titanic)
tab.Survived.Sex <- xtabs(Freq~Survived+Sex,data=Titanic)
tab.Survived.Age <- xtabs(Freq~Survived+Age,data=Titanic)
tab.Survived <- xtabs(Freq~Survived,data=Titanic)

## ----message=FALSE-------------------------------------------------------
library(memisc)

## ------------------------------------------------------------------------
(ftab.Survived.Age <- ftable(tab.Survived.Age))
(ftab.Survived.Sex <- ftable(tab.Survived.Sex))

## ------------------------------------------------------------------------
cbind(ftab.Survived.Age,
      ftab.Survived.Sex)

## ------------------------------------------------------------------------
cbind(ftab.Survived.Age,
      ftab.Survived.Sex,
      Total=tab.Survived)

## ----results='asis'------------------------------------------------------
show_html(ftab.Survived.Age)

## ----results='asis'------------------------------------------------------
show_html(ftab.Survived.Sex)

## ----results='asis'------------------------------------------------------
show_html(
  cbind(ftab.Survived.Age,
      ftab.Survived.Sex,
      Total=tab.Survived)
)

## ------------------------------------------------------------------------
knit_print.ftable_matrix <-function(x,options,...)
  knitr::asis_output(
    format_html(x,
                digits=if(length(options$ftable.digits))
                          options$ftable.digits
                       else 0,
                ...))

## ------------------------------------------------------------------------
cbind(ftab.Survived.Age,
      ftab.Survived.Sex,
      Total=tab.Survived)

## ------------------------------------------------------------------------
rm(knit_print.ftable_matrix)

## ----results='asis'------------------------------------------------------
show_html(
  cbind(ftab.Survived.Age,
      ftab.Survived.Sex,
      Total=tab.Survived),
  varinfront=FALSE
)

## ----results='asis'------------------------------------------------------
show_html(
  cbind(ftab.Survived.Age,
      ftab.Survived.Sex,
      Total=tab.Survived),
  varontop=FALSE
)

## ------------------------------------------------------------------------
ftab.Age.Survived <- ftable(tab.Survived.Age,col.vars=1)
ftab.Sex.Survived <- ftable(tab.Survived.Sex,col.vars=1)
ftab.Class.Survived <- ftable(tab.Survived.Class,col.vars=1)

rbind(
  ftab.Age.Survived,
  ftab.Sex.Survived,
  ftab.Class.Survived,
  Total=tab.Survived
)

## ----results='asis'------------------------------------------------------
show_html(
  rbind(
    ftab.Age.Survived,
    ftab.Sex.Survived,
    ftab.Class.Survived,
    Total=tab.Survived
  )
)

## ----results='asis'------------------------------------------------------
ptab.Survived.Age<-percentages(Survived~Age,data=Titanic)
ptab.Survived.Sex<-percentages(Survived~Sex,data=Titanic)
ptab.Survived.Class<-percentages(Survived~Class,data=Titanic)

fptab.Age.Survived <- ftable(ptab.Survived.Age,col.vars=1)
fptab.Sex.Survived <- ftable(ptab.Survived.Sex,col.vars=1)
fptab.Class.Survived <- ftable(ptab.Survived.Class,col.vars=1)

show_html(
  rbind(
    fptab.Age.Survived,
    fptab.Sex.Survived,
    fptab.Class.Survived
  ),
  digits=1
)

## ----results='asis'------------------------------------------------------
tab.Age <- xtabs(Freq~Age,data=Titanic)
tab.Sex <- xtabs(Freq~Sex,data=Titanic)
tab.Class <- xtabs(Freq~Class,data=Titanic)


show_html(
  rbind(
    cbind(fptab.Age.Survived,Total=tab.Age),
    cbind(fptab.Sex.Survived,Total=tab.Sex),
    cbind(fptab.Class.Survived,Total=tab.Class)
  ),
  digits=c(1,0) # One digit after dot for percentages 
                # no digits for total counts.
)

## ------------------------------------------------------------------------
toLatex(
  rbind(
    cbind(fptab.Age.Survived,Total=tab.Age),
    cbind(fptab.Sex.Survived,Total=tab.Sex),
    cbind(fptab.Class.Survived,Total=tab.Class)
  ),
  digits=c(1,0) # One digit after dot for percentages 
                # no digits for total counts.
)


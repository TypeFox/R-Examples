DATA_DIR <- "../../data"

## Anorexia data
anorexia <- read.csv("../extdata/anorexia.csv")
save(anorexia, file=file.path(DATA_DIR, "anorexia.rda"))

## anorexia.b <- anorexia[ anorexia$therapy == "b", ]
## weights <- anorexia.b$after - anorexia.b$before
UNdata <- transform(read.csv("../extdata/UN-data.csv", na=".."),
                    Nation = gsub("-", " ", Nation))
save(UNdata, file=file.path(DATA_DIR, "UNdata.rda"))

## student survey
student.survey <- transform(read.csv("../extdata/student-survey.csv"),
                            re = ordered(re, labels = c("never", "occasionally", "most weeks", "every week")),
                            pi = ordered(pi, labels = c("very liberal", "liberal", "slightly liberal",
                            "moderate", "slightly conservative", "conservative",
                            "very conservative")))
for (i in c("ve", "ab", "aa", "ld")) {
    student.survey[[i]] <- (ifelse(student.survey[[i]] == "u", NA,
                                   student.survey[[i]]) == "y")

}
save(student.survey, file=file.path(DATA_DIR, "student.survey.rda"))

## Crime data
crime2005 <- read.csv("../extdata/crime2005.csv")
save(crime2005, file=file.path(DATA_DIR, "crime2005.rda"))

## Statewide Crime
statewide.crime.2 <- read.csv("../extdata/statewide-crime-2.csv")
save(statewide.crime.2, file=file.path(DATA_DIR, "statewide.crime.2.rda"))

## OECD data
oecd.data <- read.csv("../extdata/oecd-data.csv")
save(oecd.data, file=file.path(DATA_DIR, "oecd.data.rda"))

## Mental Impairment
mental.impairment <- read.csv("../extdata/mental-impairment.csv")
names(mental.impairment) <- c("impair", "life", "ses")
save(mental.impairment, file=file.path(DATA_DIR, "mental.impairment.rda"))

##' birth-rates
birth.rates <- read.csv("../extdata/birth-rates-3rd.csv")
save(birth.rates, file=file.path(DATA_DIR, "birth.rates.rda"))

##' fertility-gdp
fertility.gdp <- read.csv("../extdata/fertility-gdp.csv")
save(fertility.gdp, file=file.path(DATA_DIR, "fertility.gdp.rda"))

##' fl-crime
fl.crime <- read.csv("../extdata/fl-crime.csv")
save(fl.crime, file=file.path(DATA_DIR, "fl.crime.rda"))

##' house-selling-price-2
house.selling.price.2 <- read.csv("../extdata/house-selling-price-2.csv")
save(house.selling.price.2, file=file.path(DATA_DIR, "house.selling.price.2.rda"))

##' house-selling-price
house.selling.price <- read.csv("../extdata/house-selling-price.csv")
save(house.selling.price, file=file.path(DATA_DIR, "house.selling.price.rda"))

##' inc.ed.race.13p1
inc.ed.race.13p1 <- read.csv("../extdata/inc-ed-race-13p1.csv")
comment(inc.ed.race.13p1$inc) <- "income (thousands of dollars)"
comment(inc.ed.race.13p1$educ) <- "number of years of education (where 12 = high school graduate, 16 = college graduate)"
comment(inc.ed.race.13p1$race) <- "racial-ethnic group (b=Black, h=Hispanic, w=White)"
comment(inc.ed.race.13p1$z1) <- "racial group (1=black, 0=white)"
comment(inc.ed.race.13p1$z2) <- "ethnic group (1=Hispanic, 0=non-Hispanic)"
save(inc.ed.race.13p1, file=file.path(DATA_DIR, "inc.ed.race.13p1.rda"))

##' income-credit
income.credit <- read.csv("../extdata/income-credit.csv")
comment(income.credit$Income) <- "annual income (euros)"
comment(income.credit$n) <- "number of subjects"
comment(income.credit$credit) <- "number possessing a credit card"
save(income.credit, file=file.path(DATA_DIR, "income.credit.rda"))

##' US.pop.size
us.pop.size <- read.csv("../extdata/us-pop-size.csv")
comment(us.pop.size$decade) <- "Decade (0=1890, 11=2000)"
comment(us.pop.size$population) <- "Population (millions)"
save(us.pop.size, file=file.path(DATA_DIR, "us.pop.size.rda"))

##' Zagat
zagat <- read.csv("../extdata/zagat.csv", header=TRUE, sep=",")
save(zagat, file=file.path(DATA_DIR, "zagat.rda"))

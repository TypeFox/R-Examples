context("ISO8601 Date conversion on Windows and Linux")

#Generate some dates
d <- as.Date("2001-01-01")
d2 <- as.Date(c("2001-01-01","2002-05-01"))

test_that("formatting date", expect_that(formatDate(d,"W%V-%Y"), equals("W01-2001")))

test_that("Formatting date vectors with ISO8601 and UK conventions",
          expect_that(formatDate(d2,"W%V-%G / W%W-%Y / %d-%m-%Y"),
                      equals(c("W01-2001 / W01-2001 / 01-01-2001","W18-2002 / W17-2002 / 01-05-2002"))))

test_that("Formatting date vectors with roman letters for quarters",
          expect_that(formatDate(d2,"%G\n%OQ"), equals(c("2001\nI","2002\nII"))))


#Some checks for the atChange
dates <- seq(as.Date("2007-01-01"),as.Date("2013-01-01"),by="1 week")
#Format with conversion string
x <- as.numeric(formatDate(dates,"%m"))
xm1 <- as.numeric(formatDate(dates[1]-7,"%m"))

#At change
test_that("atChange function works for %m",expect_that(
    atChange(x,xm1), equals(
    c(1L, 6L, 10L, 14L, 19L, 23L, 27L, 32L, 36L, 40L, 45L, 49L, 54L, 
      58L, 62L, 67L, 71L, 75L, 80L, 84L, 88L, 93L, 97L, 101L, 106L, 
      110L, 114L, 119L, 123L, 127L, 132L, 136L, 141L, 145L, 149L, 154L, 
      158L, 162L, 166L, 171L, 175L, 180L, 184L, 188L, 193L, 197L, 201L, 
      206L, 210L, 215L, 219L, 223L, 227L, 232L, 236L, 240L, 245L, 249L, 
      254L, 258L, 262L, 267L, 271L, 275L, 280L, 284L, 288L, 293L, 297L, 
      301L, 306L, 310L))))

#Test every second change function
test_that("at2ndChange function works for %m",expect_that(
    at2ndChange(x,xm1),equals(
        c(1L, 10L, 19L, 27L, 36L, 45L, 54L, 62L, 71L, 80L, 88L, 97L, 
          106L, 114L, 123L, 132L, 141L, 149L, 158L, 166L, 175L, 184L, 193L, 
          201L, 210L, 219L, 227L, 236L, 245L, 254L, 262L, 271L, 280L, 288L, 
          297L, 306L))))


#### Year formatting
x <- as.numeric(formatDate(dates,"%Y"))
xm1 <- as.numeric(formatDate(dates[1]-7,"%Y"))

test_that("atMedian function works for %Y",expect_that(
    atMedian(x,xm1),equals(c(26L, 79L, 131L, 183L, 235L, 287L))))

test_that("at2ndChange function works for %Y",expect_that(
    dates[at2ndChange(x,xm1)],equals(as.Date(c("2007-01-01","2009-01-05","2011-01-03")))))


#Does this look at expected (hard to check with testthat
data("rotaBB")
plot(rotaBB, xaxis.tickFreq=list("%Y"=atChange), xaxis.labelFreq=list("%Y"=at2ndChange),xaxis.labelFormat="%Y",xlab="time (months)")

#Test quarter formatting
test_that(formatDate(d2,"%Q"), equals(c("1","2")))

test_that(formatDate(d2,"%q"), equals(c("1","31")))

test_that(as.character(d2 - as.numeric(formatDate(d2,"%q")) + 1),
          equals(c("2001-01-01","2002-04-01")))


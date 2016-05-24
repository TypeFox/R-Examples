library(testthat)
library(LaF)

context("laf_to_ffdf")

test_that(
    "laf_to_ffdf",
    {
      lines <- c(
        " 1M 1.45Rotterdam ",
        " 2F12.00Amsterdam ",
        " 3  .22 Berlin    ",
        "  M22   Paris     ",
        " 4F12345London    ",
        " 5M     Copenhagen",
        " 6M-12.1          ",
        " 7F   -1Oslo      ")
      
      data <- data.frame(
        id=c(1,2,3,NA,4,5,6,7),
        gender=as.factor(c("M", "F", NA, "M", "F", "M", "M", "F")),
        x=c(1.45, 12, 0.22, 22, 12345, NA, -12.1, -1),
        city=c("Rotterdam", "Amsterdam", "Berlin", "Paris", 
               "London", "Copenhagen", "", "Oslo"),
        stringsAsFactors=FALSE
      )
      
      tmp.fwf <- tempfile("tmp.fwf")
      writeLines(lines, con=tmp.fwf, sep="\n")
        laf <- laf_open_fwf(filename=tmp.fwf, 
            column_types=c("integer", "categorical", "double", "string"),
            column_widths=c(2,1,5,10)
            )

        foo <- laf_to_ffdf(laf, nrow=2)
        
        # Volgende gaat mis omdat levels niet geappend worden; tenzij het 
        # bestand al een keer helemaal gelezen is en alle levels van begin af
        # aan bekend zijn. 
        # Dus dit gaat goed: 
        bar <- laf_to_ffdf(laf, nrows=1, columns=c(1, 2, 3))
        # maar dit niet
        laf2 <- laf_open_fwf(filename=tmp.fwf, 
            column_types=c("integer", "categorical", "double", "string"),
            column_widths=c(2,1,5,10)
            )
        bar2 <- laf_to_ffdf(laf2, nrows=1, columns=c(1, 2, 3))
        
        # Wat ook niet goed gaat is als er al gelezen is uit het bestand. Hij
        # begint dan met lezen vanaf laatste positie. In laf_to_ffdf ontbreekt
        # namelijk een begin(laf)
        
        bar3 <- laf_to_ffdf(laf2, nrows=1, columns=c(1, 2, 3))
})


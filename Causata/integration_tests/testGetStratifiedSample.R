library(testthat)

if (!nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) {
  library(causata)
  source('/tmp/causata_environment.R')
} else {
  source('/tmp/causata_environment.R')
  # causata.config.host <- "localhost"
  # causata.config.port <- 64030
}

source("utils.R")

test_that("Selecting one customer by customer attribute", 
  with.local.connection(function(conn) {
    causata.variables <- list("customer-id"="SELECT PROPERTY customer-id", "total-spend"="INCLUDES purchase")
    with.primary.variables(conn, causata.config, causata.variables, {      
      
      variables <- GetMetadata(conn)
      expect_that('customer-id' %in% variables$system.name, is_true())
      expect_that('total-spend' %in% variables$system.name, is_true())

      queryA <- Query() + WithVariables('customer-id', 'total-spend') + Limit(200)
      data <- GetCausataData(conn, queryA)
      
      newdata = GetStratifiedSample(conn, queryA, "total.spend_All.Past", causata.data=data)
      
      expect_that(dim(newdata$df)[1], equals(200))
      
      weights <- unique(newdata$df$weight)
      expect_that(length(weights), equals(2))
      expect_that(min(weights), equals(1))
      expect_that(max(weights) > 1, is_true())
      
      dv <- newdata$df$total.spend_All.Past
      frequency.table <- table(newdata$df$total.spend_All.Past)

      expect_that(frequency.table["0"], equals(100))
      expect_that(frequency.table["1"], equals(100))
    }) 
  })
)

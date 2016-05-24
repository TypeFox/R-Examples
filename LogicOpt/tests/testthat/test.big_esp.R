'test_esp_local' <- function(esp_file)
{
    ett <- logicopt(esp_file=esp_file,,mode="espresso")[[1]]
    qtt <- logicopt(esp_file=esp_file,,mode="qm", exact_cover=FALSE)[[1]]

    return(paste('ESPRESSO =',nrow(ett),'QM =',nrow(qtt)))
}

test_that("Logicopt big tests", {
   skip_on_cran()

   esp_file <- system.file("extdata/espresso/small.esp", package="LogicOpt") 
   tt1 <- logicopt(esp_file=esp_file)
   tt2 <- logicopt(tt1[[1]],4,3)
   expect_equal(tt1[2],tt2[2])

   esp_file <- system.file("extdata/espresso/cmu.esp", package="LogicOpt") 
   tt1 <- logicopt(esp_file=esp_file)
   tt2 <- logicopt(tt1[[1]],4,1,input_sizes=c(2,2,2,2))
   expect_equal(tt1[2],tt2[2])

   #esp_file <- system.file("extdata/espresso/pdc.esp", package="LogicOpt") 
   #esp_str  <- test_esp_local(esp_file)
   #expect_equal(esp_str,"ESPRESSO = 145 QM = 96")

   esp_file <- system.file("extdata/espresso/cmu.esp", package="LogicOpt") 
   esp_str  <- test_esp_local(esp_file)
   expect_equal(esp_str,"ESPRESSO = 3 QM = 3")

   esp_file <- system.file("extdata/espresso/ex1010.esp", package="LogicOpt") 
   tt1 <- logicopt(esp_file=esp_file)
   tt2 <- logicopt(tt1[[1]],10,10)
   expect_equal(tt1[2],tt2[2])
   #esp_str  <- test_esp_local(esp_file)
   #expect_equal(esp_str,"ESPRESSO = 284 QM = 256")

   esp_file <- system.file("extdata/espresso/health.esp", package="LogicOpt") 
   esp_str  <- test_esp_local(esp_file)
   expect_equal(esp_str,"ESPRESSO = 3 QM = 3")

   esp_file <- system.file("extdata/espresso/misex3.esp", package="LogicOpt") 
   esp_str  <- test_esp_local(esp_file)
   expect_equal(esp_str,"ESPRESSO = 690 QM = 656")

})


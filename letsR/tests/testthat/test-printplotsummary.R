context("Test for general methods")

# Path
data(PAM)

test_that("Summary and print PresenceAbsence object", {

print(PAM)
summarPAM <- summary(PAM)  
print.summary.PresenceAbsence(summarPAM) 

})


test_that("plot object", {
  
  plot(PAM)
  plot(PAM, name = "Phyllomedusa atelopoides")
  plot(PAM, world = FALSE)
  plot(PAM, col_rich = rainbow)
  
})

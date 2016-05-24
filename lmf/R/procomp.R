procomp <-
function(a, uage)
{
  #"a" is a data set with fecundity and survival of individuals with
  #known age
  cbind(uage, tapply(a$recruits, a$age, mean, na.rm = TRUE), tapply(a$survival,
                                                                    a$age, mean, na.rm = TRUE))
}

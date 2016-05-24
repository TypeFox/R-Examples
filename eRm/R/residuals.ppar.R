`residuals.ppar` <-
function(object,...)
# computes standardized residuals
# for object of class "ppar" (from person.parameter)
{
  result <- itemfit(object)$st.res
  result
}


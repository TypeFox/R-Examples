"anova.gamlist" <-
function(object, ..., test = c("none", "Chisq", "F", "Cp")){
  test=match.arg(test)
  class(object)="glmlist"
  anova(object, test = test)
}

print.raov <- function (x, digits = max(5, .Options$digits - 2), ...) {
  cat("\nRobust ANOVA Table\n")
  Table<-round(x$tab,digits)
  print(Table,...)
}

p.ci <- function(lower, upper){
  lower <- meta:::rmSpace(lower)
  upper <- meta:::rmSpace(upper)
  ##
  ifelse (lower=="NA" & upper=="NA",
          "",
          paste(" [", format(lower, justify="right"),
                "; ", format(upper, justify="right"), "]", sep=""))
}

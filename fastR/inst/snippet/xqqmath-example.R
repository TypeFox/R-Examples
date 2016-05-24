x = rnorm(100)
p <- qqmath(~x, main='QQ-plot using qqmath()')
q <- xqqmath(~x, main='QQ-plot using xqqmath()')

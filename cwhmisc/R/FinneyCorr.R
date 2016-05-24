FinneyCorr <- function(s,n) {
# Finney''s correction factor K in
#   x	=	e^{\ln x} \cdot K, \textrm{with} K = e^{s_{\ln x}^2/2} \cdot \left \{1-\frac{s^2}{4n}(s^2+2)+\frac{s^4}{96n^2}(3s^4+44s^2+84) \right \}.
  s2 <- s^2
  exp(s2/2)*(1-s2*(s2+2)/(4*n)+s2^2*(3*s2^2+44*s2+84)/(96*n^2))
}

FC.lm <- function(lmobj) {
  sumo <- summary.lm(lmobj,correlation=FALSE)
  FinneyCorr(sumo$sigma,length(sumo$residuals))
}

R2.lm <- function(lmobj) {
  summary.lm(lmobj,correlation=FALSE)$r.squared
}

s.lm <- function(lmobj) {summary.lm(lmobj,correlation=FALSE)$sigma}

summaryFs <- function(lmobj){
  print(summary(lmobj))
  prinE("FC.lm(lmobj)")
  prinE("s.lm (lmobj)")
}

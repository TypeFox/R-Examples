#' SchumakerIndInterval
#'
#' This creates quadratic coefficients for one interval of a domain.
#' This is an internal function that is called from the Schumaker function.
#' @param z This is the y value at edges of an interval.
#' @param s This is the slope at edges of an interval.
#' @param Smallt This is x values at the edge of an interval.
#' @return The location of the knot and quadratic coefficients for an interval.
SchumakerIndInterval = function(z,s,Smallt){
  # The SchumakerIndInterval function takes in each interval individually
  # and returns the location of the knot as well as the quadratic coefficients in each subinterval.

  # Judd (1998), page 232, Lemma 6.11.1 provides this if condition:
  if (sum(s)*(Smallt[2]-Smallt[1]) == 2*(z[2]-z[1])){
    tsi = Smallt[2]
  } else {
    # Judd (1998), page 233, Algorithm 6.3 along with equations 6.11.4 and 6.11.5 provide this whole section
    delta = (z[2] -z[1])/(Smallt[2]-Smallt[1])
    Condition = ((s[1]-delta)*(s[2]-delta) >= 0)
    Condition2 = abs(s[2]-delta) < abs(s[1]-delta)
    if (Condition){
      tsi = sum(Smallt)/2
    } else {
      if (Condition2){
        tsi = (Smallt[1] + (Smallt[2]-Smallt[1])*(s[2]-delta)/(s[2]-s[1]))
      } else {
        tsi = (Smallt[2] + (Smallt[2]-Smallt[1])*(s[1]-delta)/(s[2]-s[1]))
      }
    }
  }
  # Judd (1998), page 232, 3rd last equation of page.
  alpha = tsi-Smallt[1]
  beta = Smallt[2]-tsi
  # Judd (1998), page 232, 4th last equation of page.
  sbar = (2*(z[2]-z[1])-(alpha*s[1]+beta*s[2]))/(Smallt[2]-Smallt[1])
  # Judd (1998), page 232, 3rd equation of page. (C1, B1, A1)
  Coeffs1 = c((sbar-s[1])/(2*alpha), s[1], z[1])
  if (beta == 0){
    Coeffs2 = c(Coeffs1)
  } else {
    # Judd (1998), page 232, 4th equation of page. (C2, B2, A2)
    Coeffs2 = c((s[2]-sbar)/(2*beta), sbar , Coeffs1 %*% c(alpha^2, alpha, 1))
  }

  data.frame(tsi = c(tsi, tsi), C = c(Coeffs1[1], Coeffs2[1]), B =  c(Coeffs1[2], Coeffs2[2]), A = c(Coeffs1[3], Coeffs2[3]) )
}

StereoWeb <-
function(){
    my.lambda.sc = seq(from = 0, to = pi, by = pi / 36) - pi / 2
    for (j in seq(from = -1 * pi / 2 + pi / 18, to = pi / 2 - pi / 18, by = pi / 18)) {
        #phi = j
        x = (sqrt(2) / 2) * sqrt(2/(1 + cos(j) * cos(my.lambda.sc))) * cos(j) * sin(my.lambda.sc)
        y = (sqrt(2) / 2) * sqrt(2/(1 + cos(j) * cos(my.lambda.sc))) * sin(j)
        lines(x, y, lwd = .5, col = '#dddddd')
    }
    my.phi <- seq(from = -1 * pi / 2, to = pi / 2, by = pi / 36)
    for (j in seq(from = pi / 18, to = pi - pi / 18, by = pi / 18)) {
        x = (sqrt(2) / 2) * sqrt(2/(1 + cos(my.phi) * cos(j - pi / 2))) * cos(my.phi) * sin(j - pi / 2)
        y = (sqrt(2) / 2) * sqrt(2/(1 + cos(my.phi) * cos(j - pi / 2))) * sin(my.phi)
        lines(x, y, lwd = .5, col = '#dddddd')
    }
  lines(c(-0.025, 0.025), c(0, 0), lwd = .5)
  lines(c(0, 0), c(-0.025, 0.025), lwd = .5)
}

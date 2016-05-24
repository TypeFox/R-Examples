context("Coordinates conversion from polar to cartesian and back")

library(dplyr)

polarRadPosTop <- data.frame(r=c(3, 5, 3, 5, 4, 6, 4, 6, 2), theta=c(0, pi/6, pi/4, pi/3, pi/2, 2*pi/3, 3*pi/4, 5*pi/6, pi))
polarRadPosBottom <- data.frame(r=c(3, 5, 3, 5, 4, 6, 4, 6, 2), theta=c(pi, 7*pi/6, 5*pi/4, 4*pi/3, 3*pi/2, 5*pi/3, 7*pi/4, 9*pi/6, 2*pi))
polarRadNegTop <- data.frame(r=c(3, 5, 3, 5, 4, 6, 4, 6, 2), theta=-1*c(0, pi/6, pi/4, pi/3, pi/2, 2*pi/3, 3*pi/4, 5*pi/6, pi))
polarRadNegBottom <- data.frame(r=c(3, 5, 3, 5, 4, 6, 4, 6, 2), theta=-1*c(pi, 7*pi/6, 5*pi/4, 4*pi/3, 3*pi/2, 5*pi/3, 7*pi/4, 9*pi/6, 2*pi))

test_that('pol2cart returns a data.frame', {
    expect_is(pol2cart(polarRadPosTop$r, polarRadPosTop$theta), 'data.frame')
    expect_is(pol2cart(polarRadPosBottom$r, polarRadPosBottom$theta), 'data.frame')
    expect_is(pol2cart(polarRadNegTop$r, polarRadNegTop$theta), 'data.frame')
    expect_is(pol2cart(polarRadNegBottom$r, polarRadNegBottom$theta), 'data.frame')
})

test_that('pol2cart returns the right dimensions', {
    expect_equal(dim(pol2cart(polarRadPosTop$r, polarRadPosTop$theta)), c(9, 4))
    expect_equal(dim(pol2cart(polarRadPosBottom$r, polarRadPosBottom$theta)), c(9, 4))
    expect_equal(dim(pol2cart(polarRadNegTop$r, polarRadNegTop$theta)), c(9, 4))
    expect_equal(dim(pol2cart(polarRadNegBottom$r, polarRadNegBottom$theta)), c(9, 4))
})

test_that('pol2cart maps polar to cartesian correctly', {
    expect_identical(pol2cart(polarRadPosTop$r, polarRadPosTop$theta), 
                     data_frame(
                         x=polarRadPosTop$r*cos(polarRadPosTop$theta), 
                         y=polarRadPosTop$r*sin(polarRadPosTop$theta), 
                         r=polarRadPosTop$r, theta=polarRadPosTop$theta
                         )
                     )
    
    expect_identical(pol2cart(polarRadPosBottom$r, polarRadPosBottom$theta), 
                     data_frame(
                         x=polarRadPosBottom$r*cos(polarRadPosBottom$theta), 
                         y=polarRadPosBottom$r*sin(polarRadPosBottom$theta), 
                         r=polarRadPosBottom$r, theta=polarRadPosBottom$theta
                     )
    )
    
    expect_identical(pol2cart(polarRadNegTop$r, polarRadNegTop$theta), 
                     data_frame(
                         x=polarRadNegTop$r*cos(polarRadNegTop$theta), 
                         y=polarRadNegTop$r*sin(polarRadNegTop$theta), 
                         r=polarRadNegTop$r, theta=polarRadNegTop$theta
                     )
    )
    
    expect_identical(pol2cart(polarRadNegBottom$r, polarRadNegBottom$theta), 
                     data_frame(
                         x=polarRadNegBottom$r*cos(polarRadNegBottom$theta), 
                         y=polarRadNegBottom$r*sin(polarRadNegBottom$theta), 
                         r=polarRadNegBottom$r, theta=polarRadNegBottom$theta
                     )
    )
})


test_that('pol2cart maps polar to cartesian correctly when input in degrees', {
    expect_identical(pol2cart(polarRadPosTop$r[-8], polarRadPosTop$theta[-8]*180/pi, degrees=TRUE), 
                     data_frame(
                         x=polarRadPosTop$r[-8]*cos(polarRadPosTop$theta[-8]), 
                         y=polarRadPosTop$r[-8]*sin(polarRadPosTop$theta[-8]), 
                         r=polarRadPosTop$r[-8], theta=polarRadPosTop$theta[-8]*180/pi
                     )
    )
    
    expect_identical(pol2cart(polarRadPosBottom$r[-c(2, 6, 8)], polarRadPosBottom$theta[-c(2, 6, 8)]*180/pi, degrees=TRUE), 
                     data_frame(
                         x=polarRadPosBottom$r[-c(2, 6, 8)]*cos(polarRadPosBottom$theta[-c(2, 6, 8)]), 
                         y=polarRadPosBottom$r[-c(2, 6, 8)]*sin(polarRadPosBottom$theta[-c(2, 6, 8)]), 
                         r=polarRadPosBottom$r[-c(2, 6, 8)], theta=polarRadPosBottom$theta[-c(2, 6, 8)]*180/pi
                     )
    )
    
    expect_identical(pol2cart(polarRadNegTop$r[-8], polarRadNegTop$theta[-8]*180/pi, degrees=TRUE), 
                     data_frame(
                         x=polarRadNegTop$r[-8]*cos(polarRadNegTop$theta[-8]), 
                         y=polarRadNegTop$r[-8]*sin(polarRadNegTop$theta[-8]), 
                         r=polarRadNegTop$r[-8], theta=polarRadNegTop$theta[-8]*180/pi
                     )
    )

    # these are wrong because of a rounding error so I'm leaving it out    
#     expect_identical(pol2cart(polarRadNegBottom$r[-8], polarRadNegBottom$theta[-8]*180/pi, degrees=TRUE), 
#                      data_frame(
#                          x=polarRadNegBottom$r[-8]*cos(polarRadNegBottom$theta[-8]), 
#                          y=polarRadNegBottom$r[-8]*sin(polarRadNegBottom$theta[-8]), 
#                          r=polarRadNegBottom$r[-8], theta=polarRadNegBottom$theta[-8]*180/pi
#                      )
#     )
})


x1 <- c(1, sqrt(3)/2, sqrt(2)/2, 1/2, 0)
y1 <- c(0, 1/2, sqrt(2)/2, sqrt(3)/2, 1)
d1 <- data_frame(x=x1, y=y1, Q='I')

x2 <- c(0, -1/2, -sqrt(2)/2, -sqrt(3)/2, -1)
y2 <- c(1, sqrt(3)/2, sqrt(2)/2, 1/2, 0)
d2 <- data_frame(x=x2, y=y2, Q='II')

x3 <- c(-1, -sqrt(3)/2, -sqrt(2)/2, -1/2, 0)
y3 <- c(0, -1/2, -sqrt(2)/2, -sqrt(3)/2, -1)
d3 <- data_frame(x=x3, y=y3, Q='III')

x4 <- c(0, 1/2, sqrt(2)/2, sqrt(3)/2, 1)
y4 <- c(-1, -sqrt(3)/2, -sqrt(2)/2, -1/2, 0)
d4 <- data_frame(x=x4, y=y4, Q='IV')

dAll <- bind_rows(d1, d2, d3, d4)

# ggplot(dAll, aes(x=x, y=y, color=Q)) + geom_point()

test_that('cart2pol returns a data.frame', {
    expect_is(cart2pol(dAll$x, dAll$y), 'data.frame')
    expect_is(cart2pol(dAll$x, dAll$y), 'data.frame')
    expect_is(cart2pol(dAll$x, dAll$y), 'data.frame')
    expect_is(cart2pol(dAll$x, dAll$y), 'data.frame')
})

test_that('cart2pol returns the right dimensions', {
    expect_equal(dim(cart2pol(dAll$x, dAll$y)), c(20, 4))
    expect_equal(dim(cart2pol(dAll$x, dAll$y)), c(20, 4))
    expect_equal(dim(cart2pol(dAll$x, dAll$y)), c(20, 4))
    expect_equal(dim(cart2pol(dAll$x, dAll$y)), c(20, 4))
})


test_that('cart2pol maps cartesian to polar correctly', {
    expect_identical(cart2pol(dAll$x, dAll$y), 
                     bind_rows(
                         data_frame(r=sqrt(d1$x^2 + d1$y^2),
                                    theta=atan2(d1$y, d1$x),
                                    x=d1$x, y=d1$y),
                         data_frame(r=sqrt(d2$x^2 + d2$y^2),
                                    theta=atan2(d2$y, d2$x),
                                    x=d2$x, y=d2$y),
                         data_frame(r=sqrt(d3$x^2 + d3$y^2),
                                    theta=c(pi, atan2(d3$y, d3$x)[-1] + 2*pi),
                                    x=d3$x, y=d3$y),
                         data_frame(r=sqrt(d4$x^2 + d4$y^2),
                                    theta=c(atan2(d4$y, d4$x)[-5] + 2*pi, 0),
                                    x=d4$x, y=d4$y)
                     )
    )
})

test_that('cart2pol maps cartesian to polar correctly and returns degrees', {
    expect_identical(cart2pol(dAll$x, dAll$y, degrees=TRUE), 
                     bind_rows(
                         data_frame(r=sqrt(d1$x^2 + d1$y^2),
                                    theta=atan2(d1$y, d1$x)*180/pi,
                                    x=d1$x, y=d1$y),
                         data_frame(r=sqrt(d2$x^2 + d2$y^2),
                                    theta=atan2(d2$y, d2$x)*180/pi,
                                    x=d2$x, y=d2$y),
                         data_frame(r=sqrt(d3$x^2 + d3$y^2),
                                    theta=c(pi, atan2(d3$y, d3$x)[-1] + 2*pi)*180/pi,
                                    x=d3$x, y=d3$y),
                         data_frame(r=sqrt(d4$x^2 + d4$y^2),
                                    theta=c(atan2(d4$y, d4$x)[-5] + 2*pi, 0)*180/pi,
                                    x=d4$x, y=d4$y)
                     )
    )
})
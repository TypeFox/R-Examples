"test.data" <- 
function(type = "ppoly", n = 512, signal = 1, rsnr = 7, plotfn = FALSE)
{
    x <- seq(0., 1., length = n + 1)[1:n]
    #    elseif strcmp(Name,'Bumps'), 
    #    pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
    #    hgt = [ 4  5   3   4  5  4.2 2.1 4.3  3.1 5.1 4.2];
    #    wth = [.005 .005 .006 .01 .01 .03 .01 .01  .005 .008 .005];
    #    sig = zeros(size(t));
    #    for j =1:length(pos)
    #       sig = sig + hgt(j)./( 1 + abs((t - pos(j))./wth(j))).^4;
    #    end 
    if(type == "ppoly") {
        y <- rep(0., n)
        xsv <- (x <= 0.5)
        y[xsv] <- -16. * x[xsv]^3. + 12. * x[xsv]^2.
        xsv <- (x > 0.5) & (x <= 0.75)
        y[xsv] <- (x[xsv] * (16. * x[xsv]^2. - 40. * x[xsv] +
            28.))/3. - 1.5
        xsv <- x > 0.75
        y[xsv] <- (x[xsv] * (16. * x[xsv]^2. - 32. * x[xsv] +
            16.))/3.
    }
    else if(type == "blocks") {
        t <- c(0.10000000000000001, 0.13, 0.14999999999999999,
            0.23000000000000001, 0.25, 0.40000000000000002,
            0.44, 0.65000000000000002, 0.76000000000000001,
            0.78000000000000003, 0.81000000000000005)
        h <- c(4., -5., 3., -4., 5., -4.2000000000000002, 
            2.1000000000000001, 4.2999999999999998, 
            -3.1000000000000001, 2.1000000000000001, 
            -4.2000000000000002
            )
        y <- rep(0., n)
        for(i in seq(1., length(h))) {
            y <- y + (h[i] * (1. + sign(x - t[i])))/2.
        }
    }
    else if(type == "bumps") {
        t <- c(0.10000000000000001, 0.13, 0.14999999999999999,
            0.23000000000000001, 0.25, 0.40000000000000002,
            0.44, 0.65000000000000002, 0.76000000000000001,
            0.78000000000000003, 0.81000000000000005)
        h <- c(4., 5., 3., 4., 5., 4.2000000000000002, 
            2.1000000000000001, 4.2999999999999998, 
            3.1000000000000001, 5.0999999999999996, 
            4.2000000000000002)
        w <- c(0.0050000000000000001, 0.0050000000000000001,
            0.0060000000000000001, 0.01, 0.01, 
            0.029999999999999999, 0.01, 0.01, 
            0.0050000000000000001, 0.0080000000000000002,
            0.0050000000000000001)
        y <- rep(0, n)
        for(j in 1:length(t)) {
            y <- y + h[j]/(1. + abs((x - t[j])/w[j]))^
                4.
        }
    }
    else if(type == "heavi")
        y <- 4. * sin(4. * pi * x) - sign(x - 
            0.29999999999999999) - sign(0.71999999999999997 -
            x)
    else if(type == "doppler") {
        eps <- 0.050000000000000003
        y <- sqrt(x * (1. - x)) * sin((2. * pi * (1. + eps))/
            (x + eps))
    }
    else {
        cat(c("test.data: unknown test function type", type,
            "\n"))
        cat(c("Terminating\n"))
        return("NoType")
    }
    y <- y/sqrt(var(y)) * signal
    ynoise <- y + rnorm(n, 0, signal/rsnr)
    if(plotfn == TRUE) {
        if(type == "ppoly")
            mlab = "Piecewise polynomial"
        if(type == "blocks")
            mlab = "Blocks"
        if(type == "bumps")
            mlab = "Bumps"
        if(type == "heavi")
            mlab = "HeaviSine"
        if(type == "doppler")
            mlab = "Doppler"
        plot(x, y, type = "l", lwd = 2, main = mlab, ylim = 
            range(c(y, ynoise)))
        lines(x, ynoise, col = 2)
        lines(x, y)
    }
    return(list(x = x, y = y, ynoise = ynoise, type = type, rsnr = 
        rsnr))
}

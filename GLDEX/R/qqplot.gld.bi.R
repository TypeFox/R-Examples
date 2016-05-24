"qqplot.gld.bi" <-
function (data, fit, param1, param2, len = 10000, name = "", 
    corner = "topleft", type = "", range = c(0, 1), xlab = "", 
    main = "") 
{
    q <- seq(0, 1, length = len)
    minr<-min(range)
    maxr<-max(range)
 
    first.len <- len * fit[9]
    second.len <- len - first.len
    first.q <- seq(0, 1, length = first.len)
    second.q <- seq(0, 1, length = second.len)
    theo.quantile <- sort(c(qgl(first.q, fit[1], fit[2], fit[3], 
        fit[4], param1), qgl(second.q, fit[5], fit[6], fit[7], 
        fit[8], param2)))
    fitted.q <- (theo.quantile[!is.inf(theo.quantile)])


    qq<-fit[9]*pgl(fitted.q,fit[1:4],param=param1)+(1-fit[9])*pgl(fitted.q,
    fit[5:8],param=param2)

    fitted.q<-fitted.q[qq>=minr & qq<=maxr]
    qq<-qq[qq>=minr & qq<=maxr]

    data.q <- quantile(data, prob = qq)

    if (type == "") {
        if (is.na(xlab) == TRUE) {
            xlab = "Quantile"
        }
        if (is.na(main) == TRUE) {
            main = paste("Quantile plot for", name)
        }
        plot(qq, data.q, xlab = xlab, ylab = "Data values", main = main, 
            xlim = c(minr, maxr), ylim = c(min(fitted.q, data.q), max(fitted.q, 
                data.q)), col = 1, type = "l")
        lines(qq, fitted.q, col = 6)
        legend(corner, c("Data", "Fitted Distribution"), col = c(1, 
            6), lty = 1)
    }
    if (type == "str.qqplot") {
        if ((xlab) == "") {
            xlab = "Data"
        }
        if ((main) == "") {
            main = paste("Quantile plot for", name)
        }
        qqplot(data.q, fitted.q, xlab = xlab, 
            ylab = "Theoretical", main = main)
        abline(0, 1, col = 2)
    }
}



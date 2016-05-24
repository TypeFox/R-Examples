"qqplot.gld"<-
function (data, fit, param, len = 10000, name = "", corner = "topleft", 
    type = "", range = c(0, 1), xlab = "", main = "") 
{
    q <- seq(min(range), max(range), length = len)
    data.q <- quantile(data, prob = q)
    fitted.q <- qgl(q, fit[1], fit[2], fit[3], fit[4], param)
    if (type == "") {
        if (is.na(xlab) == TRUE) {
            xlab = "Quantile"
        }
        if (is.na(main) == TRUE) {
            main = paste("Quantile plot for", name)
        }
        plot(q, data.q, xlab = xlab, ylab = "Data values", main = main, 
            xlim = c(min(range), max(range)), 
            ylim = c(min(fitted.q[!is.inf(fitted.q)], 
                data.q), max(fitted.q[!is.inf(fitted.q)], data.q)), 
            col = 1, type = "l")
        lines(q[!is.inf(fitted.q)], fitted.q[!is.inf(fitted.q)], 
            col = 6)
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
        qqplot(data.q[!is.inf(fitted.q)], fitted.q[!is.inf(fitted.q)], 
            xlab = xlab, ylab = "Theoretical", main = main)
        abline(0, 1, col = 2)
    }
}


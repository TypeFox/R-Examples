## A graphical method
## migrated from S-Plus code  by Dr. Jiayang Sun.
eda <- function(x)
{
    x.na <- is.na(x)
    if(any(x.na)){
        x <- x[!x.na];
        warning("Missing value found")
    }
    par(mfrow = c(2., 3.))
    qqnorm(x)
    qqline(x)
    boxplot(x)
    title("Boxplot")
    hist(x, main = "Histogram")
    iqd <- summary(x)[5.] - summary(x)[2.]
    plot(density(x, width = 2. * iqd),
         main = "Density Plot",
         ylab = "Density", type = "l")
    ts.plot(x)
    title("Time Series Plot")
    acf(x)
    invisible()
}

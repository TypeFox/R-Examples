library(rpanel)
# plots a univariate normal distribution Initialize the gnuplot handles
h1 <- Gpinit()
h2 <- Gpinit()
# set gnuplot's additional search directories, to the extdata directory from Rgnuplot (default)
Gpsetloadpath(h1)

# change gnuplot's working directory to be the same as R's working directory (default)
Gpsetwd(h1)
# read gnuplot script
gpfile <- system.file("extdata/normal-r.gnu", package = "Rgnuplot")
s2 <- Gpfile2string(gpfile)
# set default values
d1m <- 0.5
d1d <- 0.5
Gpcmd(h1, "flag=1\nd1m=0.0\nd1d=0.5\n" %s% s2)
Gpcmd(h2, "flag=0\nd1m=0.0\nd1d=0.5\n" %s% s2)

b4close <- function() {
    # close gnuplot handles
    h1 <- Gpclose(h1)
    h2 <- Gpclose(h2)
    print("bye")
}

plot.normal <- function(panel) {
    with(panel, {
        s1 <- "d1m=" %s% d1m %s% "\nd1d=" %s% d1d
        Gpresetplot(h1)
        Gpcmd(h1, "flag=1\n" %s% s1 %s% "\n" %s% s2)
        Gpresetplot(h2)
        Gpcmd(h2, "flag=0\n" %s% s1 %s% "\n" %s% s2)
    })
    panel
}
# tktoplevel<- tcltk::tktoplevel()
rp.normal <- function() {
    # x11(width = 3, height = 3)
    pname <- rpanel::rp.control("Univariate normal distribution", d1m = 0, d1d = 0.5, dev.number = dev.cur())
    rpanel::rp.slider(pname, d1m, -5, 5, initval = 0, title = "u:", action = plot.normal)
    rpanel::rp.slider(pname, d1d, 0, 5, initval = 0.5, title = "var(x):", action = plot.normal)  #    , quitbutton = TRUE
    # rp.button(pname,'Close', action=b4close)
    rpanel::rp.button(pname, {
        h1 <- Gpclose(h1)
        h2 <- Gpclose(h2)
        break
    }, title = "Close", quitbutton = TRUE)
}

rp.normal()

# close gnuplot handles h1<-Gpclose(h1) h2<-Gpclose(h2) 

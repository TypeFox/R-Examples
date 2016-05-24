# Help page, example, demo
?plot          # initiate the help page
example(plot)  # run examples at the example section
demo()         # list all demos availabe by package
demo(persp)    # demo for one function

# Math typing
aa <- 1000
windows(width = 5.5, height = 3, pointsize = 10, family = "serif")
par(mai = c(0.4, 0.4, 0.1, 0.1))
plot(1:5, type = "n")  # an empty plot generated
text(x = 2, y = 4, labels = 
  expression(paste("This is ", alpha %*% beta, ".", sep = "")))
text(x = 2, y = 3, labels = expression(delta == aa))        # bad on aa
text(x = 2, y = 2, labels = bquote(expr = delta == .(aa)))  # good
text(x = 4, y = 3, labels = substitute(expr = aa))  # good for "aa" symbol
text(x = 4, y = 2, labels = aa)  # good for the value of "aa"
myplot <- recordPlot()
 
pdf(file = "c:/aErer/fig_math.pdf", width = 6, height = 3, pointsize = 10,
  family = "serif")   
replayPlot(myplot); dev.off()    

# Interactive functions
windows(); plot(1:10); locator(n = 3); identify(x = 3, y = 3)

# Build a graph with low-level functions
set.seed(12); dat <- rnorm(30) 
windows() 
plot(dat, type = 'n', axes = FALSE, ann = FALSE)  # an empty plot
box(); axis(side = 1); axis(side = 2); points(x = dat)
title(main = "Random numbers from the normal distribution")
title(xlab = "x", ylab = "y")
legend(x = 5, y = 2, legend = "random numbers", pch = 1)
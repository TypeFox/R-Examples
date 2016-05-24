plot1 <- xyplot(births~dayofyear,groups=dayofyear%%7,births78,
			auto.key=list(columns=3))
plot2 <- bwplot(births~factor(dayofyear%%7),groups=dayofyear%%7,births78)

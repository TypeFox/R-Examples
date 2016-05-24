create.KM.plot <- function (riskgroup = NULL, survtime = NULL, survstat = NULL, file.name = NULL, main.title = "", resolution = 100) {
	
	# make appropriate data structure for coxph
	all.data <- data.frame("riskgroup" = riskgroup, "survtime" = survtime, "survstat" =  survstat);

	coxmodel.obj <- summary(
		coxph(
			Surv(survtime, survstat) ~ riskgroup, 
			data = all.data
			)
		);
	survfit.obj <- survfit(
		Surv(survtime, survstat) ~ riskgroup, 
		data = all.data
		);

	# set the graphics driver
	current.type <- getOption("bitmapType");
	options(bitmapType = "cairo");

	png(
		filename = file.name,
		height = 5,
		width = 5,
		units = "in",
		res = resolution
		#,compression = "lzw"
		);
	par(
		mfrow = c(1,1),
		tcl = -0.1, 
		las = 1,
		font.lab = 2,
		font.axis = 2,
		mar = c(2.5, 2.3, 1, 0.2),
		mgp = c(1.4, 0.3, 0),
		oma = c(0, 0, 0, 0)
		);

	plot(survfit.obj, xlab = "Time (Years)", ylab = "Fraction of Cohort", col = c(1,2), lwd = 2);
	text(1, 0.1, paste("HR=", round(coxmodel.obj$conf.int[1,1], digits=2)), adj = c(0,0), font = 2);
	text(1, 0, paste("P=", format(coxmodel.obj$coef[1,5], scientific = TRUE, digits=3)), adj = c(0,0), font = 2);
	title(main.title, cex.main = 1, font.main = 2);

	dev.off();
	options(bitmapType = current.type);

	}
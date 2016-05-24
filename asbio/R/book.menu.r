book.menu <- function () 
{
    tclServiceMode(TRUE)
    tt <- tktoplevel()
    tkwm.title(tt, "Foundational and Applied Statistics for Biologists using R")
    topMenu <- tkmenu(tt)
    tkconfigure(tt, menu = topMenu, width = 800, height = 150)
  
#----------------------- Menu labels ---------------------#

# Book menu
	Ch1 <- tkmenu(topMenu, tearoff = FALSE)
	Ch2 <- tkmenu(topMenu, tearoff = FALSE)
	Ch3 <- tkmenu(topMenu, tearoff = FALSE)
	Ch4 <- tkmenu(topMenu, tearoff = FALSE)
	Ch5	<- tkmenu(topMenu, tearoff = FALSE)
	Ch6 <- tkmenu(topMenu, tearoff = FALSE)
	Ch7 <- tkmenu(topMenu, tearoff = FALSE)
	Ch8 <- tkmenu(topMenu, tearoff = FALSE)
	Ch9 <- tkmenu(topMenu, tearoff = FALSE)
	Ch10 <- tkmenu(topMenu, tearoff = FALSE)
	Ch11 <- tkmenu(topMenu, tearoff = FALSE)
	
# Biology menu
	Biology <- tkmenu(topMenu, tearoff = FALSE)

# Self test questions
     questions1 <- tkmenu(topMenu, tearoff = FALSE)
     questions2 <- tkmenu(topMenu, tearoff = FALSE)
	questions3 <- tkmenu(topMenu, tearoff = FALSE)
	questions4 <- tkmenu(topMenu, tearoff = FALSE)
	questions5 <- tkmenu(topMenu, tearoff = FALSE)
	questions6 <- tkmenu(topMenu, tearoff = FALSE)
	questions7 <- tkmenu(topMenu, tearoff = FALSE)
	questions8 <- tkmenu(topMenu, tearoff = FALSE)
	questions9 <- tkmenu(topMenu, tearoff = FALSE)
	questions10 <- tkmenu(topMenu, tearoff = FALSE)
	questions11 <- tkmenu(topMenu, tearoff = FALSE)
	
# Distributions
	dists1 <- tkmenu(topMenu, tearoff = FALSE) # Ch 3 (Ch 4, Ch 5, Ch 6)

#------------------- Menu contents -----------------------#

# Chapter 1
	tkadd(Ch1, "cascade", label = "Self-test questions", menu = questions1) #Names question menus
		tkadd(questions1, "command", label = "Logic", command = substitute(see.logic())) # Ch 1 question
	tkadd(Ch1, "command", label = "Accuracy and precision", command = substitute(see.accPrec.tck())) # Ch 1
	tkadd(topMenu, "cascade", label = "Chapter 1", menu = Ch1) # Names Ch 1

# Chapter 2
	tkadd(Ch2, "cascade", label = "Self-test questions", menu = questions2) #Names question menus
		tkadd(questions2, "command", label = "Probability", command = substitute(selftest.prob.tck1())) # Ch 2 question
	tkadd(Ch2, "command", label = "Bayesian analysis of discrete data", 
        command = substitute(print(Bayes.disc.tck()))) # Ch 2
	tkadd(Ch2, "command", label = "Coin flips", command = substitute(anm.coin.tck())) # Ch 2
	tkadd(Ch2, "command", label = "Die throws", command = substitute(anm.die.tck())) # Ch 2
	tkadd(Ch2, "command", label = "Venn diagrams", command = substitute(Venn.tck())) # Ch 2
	tkadd(topMenu, "cascade", label = "Chapter 2", menu = Ch2) # Names Ch 2

# Chapter 3
	tkadd(Ch3, "cascade", label = "Self-test questions", menu = questions3) #Names question menus
		tkadd(questions3, "command", label = "Pdfs", command = substitute(selftest.pdfs.tck1())) # Ch 3 question
	tkadd(Ch3, "command", label = "Continuous pdf conceptualization", command = substitute(see.pdf.conc.tck())) # Ch 3 
	tkadd(Ch3, "command", label = "Exponential power function", 
		command = substitute(see.exppower.tck())) # Ch 3
      tkadd(Ch3, "command", label = "Pdf depiction", command = substitute(see.pdfdriver.tck())) # Ch 3 (also in Ch 4, 5, 6)
	tkadd(Ch3, "cascade", label = "Pdf probability", menu = dists1) # Names dist1 (also in Ch 4, 5, 6) 
	tkadd(topMenu, "cascade", label = "Chapter 3", menu = Ch3) # Names Ch 3
	
# Chapter 4
	tkadd(Ch4, "cascade", label = "Self-test questions", menu = questions4) #Names question menus
		tkadd(questions4, "command", label = "Statistics and inference", command = substitute(selftest.stats.tck1())) # Ch 4 question
	tkadd(Ch4, "command", label = "Cognitive illusions", command = substitute(example(illusions))) # Ch 4
	tkadd(Ch4, "command", label = "Least squares", command = substitute(anm.ls.tck())) # Ch 4
	tkadd(Ch4, "command", label = "Log-likelihood", command = substitute(anm.loglik.tck())) # Ch 4
	tkadd(Ch4, "command", label = "M-estimation", command = substitute(see.M())) # Ch 4	
	tkadd(topMenu, "cascade", label = "Chapter 4", menu = Ch4) # Names Ch 4

    # Shaded Distributions
	tkadd(dists1, "command", label = "Normal", command = substitute(shade.norm.tck()))  
	tkadd(dists1, "command", label = "Binomial", command = substitute(shade.bin.tck()))
	tkadd(dists1, "command", label = "Chi-square", command = substitute(shade.chi.tck()))
	tkadd(dists1, "command", label = "F", command = substitute(shade.F.tck()))
	tkadd(dists1, "command", label = "Poisson", command = substitute(shade.poi.tck()))
	tkadd(dists1, "command", label = "t", command = substitute(shade.t.tck()))

# Chapter 5
	tkadd(Ch5, "cascade", label = "Self-test questions", menu = questions5) #Names question menus
		tkadd(questions5, "command", label = "Sampling distributions", 
			command = substitute(selftest.sampd.tck1())) # Ch 5 question
		tkadd(questions5, "command", label = "Confidence intervals", 
			command = substitute(selftest.conf.tck1())) # Ch 5 question
	tkadd(Ch5, "command", label = "Confidence intervals", command = substitute(anm.ci.tck())) # Ch 5
	tkadd(Ch5, "command", label = "MCMC simulation", command = substitute(anm.mc.bvn.tck())) # Ch 5
	tkadd(Ch5, "command", label = "Sampling distribution basics", 
        command = substitute(samp.dist.mech.tck())) # Ch 5 
	tkadd(Ch5, "command", label = "Sampling distributions", 
        command = substitute(samp.dist.method.tck())) # Ch 5
	tkadd(topMenu, "cascade", label = "Chapter 5", menu = Ch5) # Names Ch 5
	
# Chapter 6
	tkadd(Ch6, "cascade", label = "Self-test questions", menu = questions6) #Names question menus
	tkadd(questions6, "command", label = "Null hypothesis tests", 
			command = substitute(selftest.H0.tck1())) # Ch 6 question
	tkadd(Ch6, "command", label = "Power", command = substitute(see.power.tck())) # Ch 6
	tkadd(Ch6, "command", label = "t-test mechanics", command = substitute(see.ttest.tck())) # Ch 6
	tkadd(Ch6, "command", label = "Type I and II error", command = substitute(see.typeI_II())) # Ch 6
	tkadd(topMenu, "cascade", label = "Chapter 6", menu = Ch6) # Names Ch 6

# Chapter 7
	tkadd(Ch7, "cascade", label = "Self-test questions", menu = questions7) #Names question menus
	   tkadd(questions7, "command", label = "Experimental and sampling design", 
     command = substitute(selftest.se.tck1()))   # Ch 7 question
	tkadd(Ch7, "command", label = "Experimental designs", command = substitute(anm.ExpDesign.tck())) # Ch 7
	tkadd(Ch7, "command", label = "Regression (Add/delete points)", 
        command = substitute(see.adddel())) # Ch 7 (also in Ch 9)
	tkadd(Ch7, "command", label = "Regression (Move points)", 
        command = substitute(see.move())) # Ch 7 (also in Ch 9)
	tkadd(Ch7, "command", label = "Sampling designs", command = substitute(anm.samp.design.tck())) # Ch 7
	tkadd(topMenu, "cascade", label = "Chapter 7", menu = Ch7) # Names Ch 7

# Chapter 8
	tkadd(Ch8, "cascade", label = "Self-test questions", menu = questions8) #Names question menus
	tkadd(questions8, "command", label = "Correlation", 
			command = substitute(selftest.corr.tck1())) # Ch 8 question
	tkadd(Ch8, "command", label = "Effect of range on correlation", command = substitute(see.cor.range.tck())) # Ch 8
	tkadd(Ch8, "command", label = "Pearson correlation sampling distribution", command = substitute(see.r.dist.tck())) # Ch 8
	tkadd(topMenu, "cascade", label = "Chapter 8", menu = Ch8) # Names Ch 8

# Chapter 9
	tkadd(Ch9, "cascade", label = "Self-test questions", menu = questions9) #Names question menus
	tkadd(questions9, "command", label = "Approaches", 
			command = substitute(selftest.regapproaches.tck1())) # Ch 9 question
	tkadd(questions9, "command", label = "Diagnostics", 
			command = substitute(selftest.regdiag.tck1())) # Ch 9 question
	tkadd(Ch9, "command", label = "Linear models (Regression)", 
		command = substitute(see.lmr.tck())) # Ch 9
	tkadd(Ch9, "command", label = "Non-linear models", 
		command = substitute(see.nlm())) # Ch 9
	tkadd(Ch9, "command", label = "OLS estimation in regression", 
		command = substitute(anm.ls.reg.tck())) # Ch 9
	tkadd(Ch9, "command", label = "Regression (Add/delete points)", 
		command = substitute(see.adddel())) # Ch 9 (also in Ch 7)
	tkadd(Ch9, "command", label = "Regression (Move points)", 
		command = substitute(see.move())) # Ch 9 (also in Ch 7)
	tkadd(Ch9, "command", label = "Regression mechanics", 
		command = substitute(see.regression.tck())) # Ch 9
	tkadd(Ch9, "command", label = "ROC mechanics", 
		command = substitute(see.roc.tck())) # Ch 9
	tkadd(Ch9, "command", label = "Smoothers", 
		command = substitute(see.smooth.tck())) # Ch 9
	tkadd(topMenu, "cascade", label = "Chapter 9", menu = Ch9) # Names Ch 9
	
# Chapter 10
	tkadd(Ch10, "cascade", label = "Self-test questions", menu = questions10) #Names question menus
	tkadd(questions10, "command", label = "Simultaneous inference", 
			command = substitute(selftest.ANOVAsiminf.tck1())) # Ch 10 question
	tkadd(questions10, "command", label = "Mixed effect models", 
			command = substitute(selftest.ANOVAmixed.tck1()))
	tkadd(questions10, "command", label = "Type II and III SS", 
			command = substitute(selftest.typeIISS.tck1()))
	tkadd(Ch10, "command", label = "ANOVA mechanics", 
		command = substitute(see.anova.tck())) # Ch 10
	tkadd(Ch10, "command", label = "ANCOVA mechanics", 
		command = substitute(see.ancova.tck())) # Ch 10  
	tkadd(Ch10, "command", label = "Fixed effects in mixed effect model", 
		command = substitute(see.mixedII())) # Ch 10
	tkadd(Ch10, "command", label = "Linear models (ANOVA)", 
		command = substitute(see.lma.tck())) # Ch 10
	tkadd(Ch10, "command", label = "Type I, II, III sums of squares", 
		command = substitute(see.lmu.tck())) # Ch 10
	tkadd(Ch10, "command", label = "Random effects", 
		command = substitute(see.rEffect.tck())) # Ch 10
	tkadd(topMenu, "cascade", label = "Chapter 10", menu = Ch10) # Names Ch 10
	
# Chapter 11
	tkadd(Ch11, "cascade", label = "Self-test questions", menu = questions11) #Names question menus
	tkadd(Ch11, "command", label = "Multinomial distribution", 
		command = substitute(see.mnom.tck())) # Ch 11
	tkadd(Ch11, "command", label = "Simpson's paradox", 
		command = substitute(example(paik))) # Ch 11
	tkadd(topMenu, "cascade", label = "Chapter 11", menu = Ch11) # Names Ch 11

# Biology
  
	tkadd(Biology, "command", label = "Geometric growth", 
		command = substitute(anm.geo.growth.tck()))
	tkadd(Biology, "command", label = "Exponential growth", 
		command = substitute(anm.exp.growth.tck()))
	tkadd(Biology, "command", label = "Hardy Weinberg equilibrium", 
		command = substitute(see.HW.tck()))    
	tkadd(Biology, "command", label = "Logistic growth", 
		command = substitute(anm.log.growth.tck()))          
	tkadd(Biology, "command", label = "Lotka-Volterra competition", 
		command = substitute(anm.LVc.tck()))
	tkadd(Biology, "command", label = "Lotka-Volterra exploitation", 
		command = substitute(anm.LVe.tck()))
	tkadd(Biology, "command", label = "Matrix population models", 
		command = substitute(anm.TM.tck()))       
	tkadd(topMenu, "cascade", label = "Biology", menu = Biology)
}


`hglm2` <-
	function(meanmodel = NULL, data = NULL, family = gaussian(link = identity),
             rand.family = gaussian(link = identity), method = "EQL", 
             conv = 1e-6, maxit = 50, startval = NULL,
             X.disp = NULL, disp = NULL, link.disp = "log", 
             weights = NULL, fix.disp = NULL, offset = NULL, sparse = TRUE,
			 vcovmat = FALSE, calc.like = FALSE, RandC = NULL,
			 bigRR = FALSE, verbose = FALSE, ...) UseMethod("hglm2")


AmpSim <- function(cyc = 1:35, b.eff = -25, bl = 0.05, ampl = 1, 
		   Cq = 20, noise = FALSE, nnl = 0.025, 
		   nnl.method = "constant") {
#   tmp.warn <- getOption("warn")
#   options(warn = -1)
  if (min(cyc) < 1) 
    stop("The minimum cycle value must 1 or larger.")
  if (Cq < 1) 
    stop("The Cq value must larger or equal to 1.")
  if (nnl < 0 || nnl > 10) 
    stop("nnl must be within 0 and 10.")

# Define the model used to simulate the amplification curve
# based on a 5-parameter sigmoidal function

  fluo <- bl + (
		(ampl - bl) / 
		(1 + exp(b.eff * (log(cyc) - log(Cq))))
		)

# Decide if noise is added and how noise is added to the 
# simulated curve
# The noise is simulated by the rnorm function and used
# defined values

 if (noise) {
 	mean.noise <- mean(fluo) * nnl
  	sd.noise <- sd(fluo) * nnl
  if (nnl.method == "increase") {
	noise.sim <- sort(rnorm(length(fluo), 
	mean = mean.noise, 
	sd = sd.noise))
  }
  if (nnl.method == "decrease") {
	noise.sim <- sort(rnorm(length(fluo), 
			mean = mean.noise, 
			sd = sd.noise), decreasing = TRUE)
  }
  if (nnl.method == "constant") {
	noise.sim <- rnorm(length(fluo), 
			mean = mean.noise, 
			sd = sd.noise)
  }
    
# Add the noise the the simulated curve and combine the output
    fluo <- fluo + noise.sim
    res <- data.frame(cyc, fluo)
  } else (res <- data.frame(cyc, fluo))

  res
}

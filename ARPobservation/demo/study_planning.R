set.seed(123)
# set parameters
stream_length <- c(5,10)
sessions <- 4
phi <- 0.25
zeta <- 1 / 60
prev_ratio <- 0.4
replicates <- 5

# create parameter combinations

# function for simulating ABAB data
ABAB <- function(phi, zeta, prev_ratio, stream_length, sessions, replicates) {
  Condition <- rep(rep(c("Baseline","NCR"), times = 2), each = sessions)
  Phase <- rep(c("A1","B1","A2","B2"), each = sessions) 
  phi_vec <- rep(phi * rep(c(1, prev_ratio), times = 2), each = sessions)
  zeta_vec <- rep(zeta, 4 * sessions)
  BS <- r_behavior_stream(n = 4 * sessions * replicates, 
                          mu = phi_vec / zeta_vec, lambda = (1 - phi_vec) / zeta_vec, 
                          F_event = F_exp(), 
                          F_interim = F_exp(), 
                          stream_length = 60 * stream_length)
  CDR <- continuous_duration_recording(BS)
  MTS <- momentary_time_recording(BS, interval_length = 20)
  
  repnr <- rep(1:replicates, each = 4 * sessions)
  data.frame(repnr, session = 1:(4 * sessions), Phase, Condition, CDR, MTS)
}
# apply function to each combination of parameter values
library(plyr)
ABAB_sim <- mdply(data.frame(stream_length), ABAB, phi=phi, zeta = zeta, prev_ratio = prev_ratio, replicates = replicates, sessions = sessions)

# plot results
library(reshape2)
ABAB_sim <- melt(ABAB_sim, measure.vars = c("CDR","MTS"))
ABAB_sim <- within(ABAB_sim,{
  session_length <- ordered(paste(stream_length, "m"), levels = c("5 m","10 m"))
  system <- factor(variable, labels = c("Continuous recording", "20 s MTS"))
  replication <- paste("Replication", repnr)
})
library(ggplot2)
qplot(x = session, y = value, color = Condition, linetype = Condition, group = Phase,
      data = ABAB_sim,
      facets = replication ~ session_length + system, 
      geom = c("point","line"), 
      xlab = "Session", ylab = "% Duration") + 
  labs(color="Phase", linetype = "Phase") + theme_bw()

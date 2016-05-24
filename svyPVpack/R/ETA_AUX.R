########################################################################
# ETA AUX FUNCTIONS
########################################################################
eta.func <- function(x, y, svydes){
  withReplicates(svydes, function(w,data){
    mod <- lm(x~factor(y), weights = w)
    anova.res <- anova(mod)
    SS <- anova.res[,2]
    eta2 <- SS[1]/sum(SS)
    sqrt(eta2)})}
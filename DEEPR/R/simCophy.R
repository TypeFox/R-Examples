simCophy <-
function (group_alphas = c(1,1,1,1), cophy_number = 20, mean_event_count = 33) {
group_probs<-rdirichlet(cophy_number,group_alphas)
totals<-rpois(1,mean_event_count*cophy_number)
group_events<-matrix(rmultinom(1,totals,group_probs),ncol=4)
return(group_events)
}
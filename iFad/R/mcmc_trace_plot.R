mcmc_trace_plot <-
function(tau_g_chain,plot_file_name,index){

# library(coda)

  tau_g_trace<-1/sqrt(cbind(tau_g_chain[[1]],tau_g_chain[[2]])[,index])


  mh.draws <- mcmc(tau_g_trace)

  pdf(plot_file_name)
  traceplot(mh.draws)
  dev.off()

}


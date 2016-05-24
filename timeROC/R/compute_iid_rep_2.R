# {{{ Technical details
# t                       : time at which we compute the i.i.d decomposition (Influence function)
# n                       : sample size
# cause                   : the value that indicates the main event of interest
# F01t                    : the cumulative incidence function of the main event at time t
# weights                 : object weights, output of the main function, order by order(T) (by default with ipcw() function of package pec)
# T                       : vector of observed failure times, order by order(T)
# delta                   : vector of indicator of status (0 for censoring, 1 for type of event one, 2 for type of event two and so on...),order by order(T)
# marker                  : vector ofmarker values,order by order(T)
# times                   : vector of times for wich we compute the AUCs
#
## CAUTION : T,delta,marker,weights should be order by order(T)
#
# }}}
compute_iid_decomposition<-function(t,n,cause,F01t,St,weights,T,delta,marker,MatInt0TcidhatMCksurEff){ 
  # indicator vectors 
  Cases<-(T< t & delta==cause)
  Controls_1<-(T> t)
  Controls_2<-(T< t &  delta!=cause & delta!=0)
  if(sum(Controls_2)>0){
    compute_iid_decomposition_competing_risks(t,n,cause,F01t,St,weights,T,delta,marker,MatInt0TcidhatMCksurEff) 
  }else{
    compute_iid_decomposition_survival(t,n,cause,F01t,St,weights,T,delta,marker,MatInt0TcidhatMCksurEff)
  }  
}

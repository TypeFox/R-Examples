RelativeRisk.Equality <-
function(alpha,beta,or,k,pt,pc){
n.c<-(qnorm(1-alpha/2)+qnorm(1-beta))^2*(1/(k*pt*(1-pt))+1/(pc*(1-pc)))/(log(or))^2
n.c
}

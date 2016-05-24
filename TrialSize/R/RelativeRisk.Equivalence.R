RelativeRisk.Equivalence <-
function(alpha,beta,or,k,pt,pc,margin){
n.c<-(qnorm(1-alpha)+qnorm(1-beta/2))^2*(1/(k*pt*(1-pt))+1/(pc*(1-pc)))/(margin-abs(log(or)))^2
n.c
}

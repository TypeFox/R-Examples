Cox.NIS <-
function(alpha,beta,loghr,p1,p2,d,margin){

n=(qnorm(alpha)+qnorm(1-beta))^2/((loghr-margin)^2*p1*p2*d)

}

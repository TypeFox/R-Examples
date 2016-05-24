Cox.Equivalence <-
function(alpha,beta,loghr,p1,p2,d,margin){

n=(qnorm(alpha/2)+qnorm(1-beta))^2/((margin-abs(loghr))^2*p1*p2*d)

}

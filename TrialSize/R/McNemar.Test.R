McNemar.Test <-
function(alpha,beta,psai,paid){

n<-(qnorm(1-alpha/2)*(psai+1)+qnorm(1-beta)*sqrt((psai+1)^2-(psai-1)^2*paid))^2/((psai-1)^2*paid)

}

eva.hat<-function(x,a=0.5,b=0.5)
{
# 0<a<1, b<1
# if b<a then marginal is unimodal
# if a^2 < b < a then not star unimodal

d<-length(x)
eta<-sum(x^2)     #vektorin x pituuden nelio
normvakio<-((2*pi)^d*(a^(-d)-b))^(-1)
tulos<-normvakio*(exp(-a^2*eta/2)-b*exp(-eta/2))

return(tulos)
}

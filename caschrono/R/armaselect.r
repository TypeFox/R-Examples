armaselect = function(y, max.p= 15, max.q=15, nbmod = 10)
{
 matlag <- function(y, maxlag)
 {
 # entrée matrice colonne ou vecteur
 # sortie matrice de maxlag colonne des valeurs retardée
 # contient un triangle de manquants
 y = as.matrix(y) ; n= nrow(y)
 x = matrix(1,nrow= n, ncol=1)
 for(i in 1:maxlag )
 # i=1
 { x= cbind(x, Lag(y,i))}
 x = x[,-1]
 colnames(x) = paste("Lag_",as.character(1:maxlag),sep="")
 x
 }
# calcule une autorégression longue puis des régressions linéaires sur y passés et résidus
# pour chacune, BIC
# inter.p = c(p1,p2)
# inter.q = c(q1,q2) intervalles de valeurs p et q
n = length(y); pmaxi = floor(min(n-1, 10*log10(n)) )
# centrage de la série
yc = y-mean(y)
# residus d'une longue autorégression
z.tilde= ar(yc, aic = FALSE, order.max = pmaxi,   method= "yule-walker", demean=FALSE)$resid
# pmaxi résidus manquants au début
## remplacés par les pmaxi derniers résidus de la régression sur série retournée
#aa1 = ar(rev(yc),aic = FALSE, order.max = pmaxi,   method= "yule-walker", demean=FALSE)
#z.tilde[1:pmaxi] = rev(aa1$resid)[1:pmaxi]
## organisation en matrice
## y retardés
## calcul de la matrice des y retardés
## ici on   prédit à l'horizon pmaxi la série retardée pour avoir des non manquants dans y
## prédiction par date croissante
#pred.rev= rev(predict(aa1, n.ahead=max.p)$pred)
#yret = matlag(c(pred.rev,yc),max.p)[-(1:max.p),]
yret = matlag(yc,max.p)
# organisation des résidus
resret = matlag(z.tilde,max.q)
# boucle sur toutes les régressions p  :   0:max.p et q  :  0:max.q
# stockage data frame de 3 colonnes : p,q, bic
# ordres nuls
# p et q nuls
mm = lm(yc ~ 0)
sbc= nrow(resret)*log(var(mm$residuals))
# autre critère d'information
bic = AIC(mm, k = log(nrow(resret)))
resul = matrix(NA,nrow = (max.p+1)*(max.q+1), ncol=4)
colnames(resul) = c("p","q","bic","sbc")
ili=1 ; resul[ili,] = c(0,0,bic,sbc)
# q nul
for(ip in 1:max.p)
{
ili = ili+1
iq=0
mm = lm(yc ~ yret[,1:ip]-1)
bic = AIC(mm, k = log(nrow(resret)))
sbc= nrow(resret)*log(var(mm$residuals))+ ip*log(nrow(resret))
resul[ili,] = c(ip,0,bic,sbc)
}
#p nul
for(iq in 1:max.q)
{
ili = ili+1 ; ip=0
mm = lm(yc ~ resret[,1:iq]-1)
bic = AIC(mm, k = log(nrow(resret)))
sbc= nrow(resret)*log(var(mm$residuals))+ iq*log(nrow(resret))
resul[ili,] = c(0,iq,bic,sbc)
}
for(ip in 1:max.p)
{
for(iq in 1:max.q)
{
ili=ili+1
mm = lm(yc ~ yret[,1:ip]+ resret[,1:iq] -1)
bic = AIC(mm, k = log(nrow(resret)))
sbc= nrow(resret)*log(var(mm$residuals))+ (ip+iq)*log(nrow(resret))
resul[ili,] = c(ip,iq,bic,sbc)
}
}
ordre = order(resul[,4])
sbc_opt = resul[ordre,][1:nbmod,c(1,2,4)]
colnames(sbc_opt ) = c("p","q","sbc")
sbc_opt
}
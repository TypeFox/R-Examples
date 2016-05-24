set.seed(2)
library(mvtnorm)


par(mfrow=c(2,2))

Wetness<-function(T, Tmin, Topt, Tmax, Wmin, Wmax) {

	fT<-((Tmax-T)/(Tmax-Topt))*(((T-Tmin)/(Topt-Tmin))^((Topt-Tmin)/(Tmax-Topt)))
	W <- Wmin/fT
	W[W>Wmax | T<Tmin | T>Tmax]<-Wmax 
	return(W)

	}

par(mfrow=c(2,2))

############Parameter distribution###############
TAB<-read.table("MagareyParam.txt", header=T)

meanP<-mean(log(TAB))
SIG<-cov(log(TAB))

print(meanP)
print(cor(log(TAB)))

Num<-10000

Z<-rmvnorm(Num,meanP,SIG)
Zexp<-exp(Z)

Tmin_vec<-Zexp[,1]
Topt_vec<-Zexp[,3]
Tmax_vec<-Zexp[,2]
Wmin_vec<-Zexp[,4]
Wmax_vec<-Zexp[,5]

############Uncertainty analysis#################

#Graphique initial

W<-Wetness(seq(5,35, by=0.1),Tmin=mean(Tmin_vec), Topt=mean(Topt_vec), Tmax=mean(Tmax_vec), Wmin=mean(Wmin_vec), Wmax=mean(Wmax_vec))

plot(seq(5,35, by=0.1), W, xlab="Temperature (degC)", ylab="Wetness duration (h)", type="l", lwd=3, ylim=c(0,50) )
text(7,45,"A")

for (k in 1:Num) {
while(Tmin_vec[k]>Topt_vec[k] | Topt_vec[k]>Tmax_vec[k] | Wmin_vec[k]>Wmax_vec[k]) {
Z<-rmvnorm(1,meanP,SIG)
Zexp<-exp(Z)
Tmin_vec[k]<-Zexp[1]
Topt_vec[k]<-Zexp[3]
Tmax_vec[k]<-Zexp[2]
Wmin_vec[k]<-Zexp[4]
Wmax_vec[k]<-Zexp[5]

}
}

#plot(Wmin_vec, Wmax_vec, xlab="Wmin", ylab="Wmax")
#plot(Tmin_vec, Tmax_vec, xlab="Tmin", ylab="Tmax")

plot(c(0), c(0), pch=" ", xlab="Temperature (degC)", ylab="Wetness duration (h)", xlim=c(5, 35), ylim=c(0, 50))
text(30,45,"B")

T_vec<-seq(from=5, to=35, by=0.1)
W_mat<-matrix(nrow=Num, ncol=length(T_vec))

for (i in 1:Num) {

	W_mat[i,]<-Wetness(T_vec, Tmin_vec[i], Topt_vec[i], Tmax_vec[i], Wmin_vec[i], Wmax_vec[i])

	if (i<20) {lines(T_vec, W_mat[i,])}

}

med_vec<-apply(W_mat, 2, quantile, 0.5)
Q0.01_vec<-apply(W_mat, 2, quantile, 0.01)
Q0.1_vec<-apply(W_mat, 2, quantile, 0.1)
Q0.9_vec<-apply(W_mat, 2, quantile, 0.9)
Q0.99_vec<-apply(W_mat, 2, quantile, 0.99)

plot(c(0), c(0), pch=" ", xlab="Temperature (degC)", ylab="Wetness duration (h)", xlim=c(5, 35), ylim=c(0, 150))
text(7,140,"C")
lines(T_vec, med_vec, lwd=3)
lines(T_vec, Q0.9_vec, lty=2)
lines(T_vec, Q0.1_vec, lty=2)
lines(T_vec, Q0.99_vec, lty=9)
lines(T_vec, Q0.01_vec, lty=9)

hist(W_mat[,221],xlab=paste("Wetness duration (h) for T=", T_vec[221],"degC"), main=" ", xlim=c(0,90))
text(80,7000,"D")
summary(W_mat[,221])
quantile(W_mat[,221], 0.1)
quantile(W_mat[,221], 0.9)
sd(W_mat[,221])/mean(W_mat[,221])

T_vec[51]
summary(W_mat[,51])
quantile(W_mat[,51], 0.1)
quantile(W_mat[,51], 0.9)
sd(W_mat[,51])/mean(W_mat[,51])

#########Sensitivity analysis#############

Indices_mat<-matrix(nrow=5, ncol=length(T_vec))

for (i in 1:length(T_vec)) {

Indices_mat[1,i]<-cor(Tmin_vec,W_mat[,i])
Indices_mat[2,i]<-cor(Topt_vec,W_mat[,i])
Indices_mat[3,i]<-cor(Tmax_vec,W_mat[,i])
Indices_mat[4,i]<-cor(Wmin_vec,W_mat[,i])
Indices_mat[5,i]<-cor(Wmax_vec,W_mat[,i])
}

dev.new()
par(mfrow=c(1,2))
plot(T_vec,Indices_mat[1,], xlab="Temperature (degC)", ylab="Correlation", ylim=c(-0.6,0.6), type="l")
lines(T_vec,Indices_mat[2,], lty=4, lwd=3)
lines(T_vec, Indices_mat[3,], lwd=3)

lines(c(5,10),c(-0.3,-0.3))
text(14,-0.3, "Tmin")
lines(c(5,10),c(-0.4,-0.4), lty=4, lwd=3)
text(14,-0.4, "Topt")
lines(c(5,10),c(-0.5,-0.5), lwd=3)
text(14,-0.5, "Tmax")

plot(T_vec,Indices_mat[4,], xlab="Temperature (degC)", ylab="Correlation", ylim=c(0.2,1), type="l")
lines(T_vec,Indices_mat[5,], lwd=3)
lines(c(5,10),c(0.4,0.4))
text(16,0.4, "Wmin")
lines(c(5,10),c(0.3,0.3), lwd=3)
text(16,0.3, "Wmax")






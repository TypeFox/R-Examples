distTraitspace <-
function (result, obs = NULL, byrow = TRUE){
if(is.null(obs)& byrow) obs <- result$true.p$p
if(is.null(obs)& !byrow) obs <- result$true.p$p_dist

# Uniform scaling
if (byrow){
result <- result$predicted.p$P_level_1_level_3
scaled.result <- (apply(result,1,norm<-function(x){return (x/sum(x))}))
scaled.obs <- (apply(obs,1,norm<-function(x){return (x/sum(x))}))
}else{
result <- result$predicted.p$P_level_1_level_3_dist
scaled.result <- (apply(result,2,norm<-function(x){return (x/sum(x))}))
scaled.obs <- (apply(obs,2,norm<-function(x){return (x/sum(x))}))
}


#### true distance
# euclidean:
# Usual distance between the two vectors (2 norm aka L_2), sqrt(sum((x_i - y_i)^2)).
EUC.dist = sqrt(apply((scaled.result - scaled.obs)^2, 2 , sum))

# manhattan:
# Absolute distance between the two vectors (1 norm aka L_1).
MAN.dist = apply(abs(scaled.result - scaled.obs), 2 , sum)

# Hellinger Distance
# sqrt(P) - sqrt(Q)
sqrt.xy = sqrt(scaled.result) - sqrt(scaled.obs)
# 2-norm of sqrt(P) - sqrt(Q)
norm.2 = apply(abs(sqrt.xy)^2,2,sum)^0.5
# Uniform scaling
H.dist = (2)^-0.5*norm.2

# Kullback Leibler divergence
# D1The K-L distance of 'spec2' with respect to 'spec1' (i.e. D(spec1 || spec2))
# D2The K-L distance of 'spec1' with respect to 'spec2' (i.e. D(spec2 || spec1))
# DThe symmetric K-L distance (i.e. D = 0.5*(D1+D2))

# my own code of Kullback Leibler divergence
# ln(P(i)/Q(i))
log.ratio = log(scaled.result+0.0001,2)- log(scaled.obs+0.0001,2)  # non zero
# sum(ln(ratio)P(i))
KL.div.D1 = apply((log.ratio)*scaled.result, 2, sum)
# KL.div.D2 = apply((-log.ratio)*scaled.obs, 2, sum)
# KL.div.D = 0.5 * (KL.div.D1 + KL.div.D2)


# bhattacharya
# Bhattacharyya coefficient
B.ceof = apply(sqrt((scaled.result+0.0001) * (scaled.obs+0.0001)),2,sum)
# Bhattacharyya distance
B.dist = -log(B.ceof)

#### check
# print(round(sqrt(1-B.ceof),2) == round(H.dist,2))

#### return
result = list(EUC.dist = EUC.dist, MAN.dist = MAN.dist, H.dist = H.dist, KL.div.D1 = KL.div.D1, B.dist = B.dist) 

return(result)
}

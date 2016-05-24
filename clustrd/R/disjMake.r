disjMake <- function(obs) {

obs=apply(as.matrix(obs),2,as.numeric)

#in case of binary input
minobs = min(obs)
maxobs = max(obs)
  
if ((minobs==0) & (maxobs==1)) {
    obs = obs + 1
}

n = nrow(obs)
p = ncol(obs)
Q = p

mods = (apply(obs, 2, max))

wz = which(mods == 0)
if (length(wz) != 0) {
  mods = mods[-wz]
  Q = length(mods)
}

cum_mods = (cumsum(mods))
J = sum(mods)

dZ = matrix(0, nrow = n, ncol = J)
for (i in 1:n) {
  dZ[i, obs[i, 1]] = 1
}

for (j in 2:Q) {
  for (i in 1:n) {
    if (obs[i, j] != 0) {
      dZ[i, cum_mods[j - 1] + obs[i, j]] = 1
    } else {
      dZ[i, cum_mods[j - 1] + 1] = 0
    }
  }
}
if ((minobs==1) & (maxobs==0)) {
  dZ=abs(dZ-1)
}
out = list()
out$J = J
out$Q = Q
out$dZ = dZ
out
} 
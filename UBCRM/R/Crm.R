Crm <-
function(Dk, prior, target = 1/3, nptmax = 24, nmaxmtd = 6, nmaxdose = nptmax, sd = 1.34, approach = "bayes", model = "power", method = "fpost", nextlevel = "ntarget", upskipping = F, downskipping = F, lastdose = NA){#browser()
# Bayesian method
if (approach == "bayes") {
if (upskipping | downskipping){
if(!lastdose %in% Dk$dose){
stop("lastdose specified not available.")
}
}
mtd <- NA
# Method 1: estimation of the mean of a then use of psi:
if (method == "fpost"){
# Power model
if (model == "power") {
sgl <- prior
	norm<- integrate(function(a){Lp(a,data=Dk,sgl=sgl)*fp(a, sd)},-Inf,Inf)$value
	f_post<- function(a){a*Lp(a,data=Dk,sgl=sgl)*fp(a, sd)/norm}
	a_hat<- integrate(f_post,-Inf,Inf)$value
p_post<- psip(sgl,a_hat)
}

# Tangent model
else if (model == "tangent") {
sgl <- ait1(prior)
	norm<- integrate(function(a){Lt(a,data=Dk,sgl=sgl)*ft(a)},0,Inf)$value
	f_post<- function(a){a*Lt(a,data=Dk,sgl=sgl)*ft(a)/norm}
	a_hat<- integrate(f_post,0,Inf)$value
p_post<- psit(sgl,a_hat)
}

# Logisitic model
else if (model == "logistic") {
sgl <- ail1(prior)
	norm<- integrate(function(a){Ll(a,data=Dk,sgl=sgl)*ft(a)},0,Inf)$value
	f_post<- function(a){a*Ll(a,data=Dk,sgl=sgl)*ft(a)/norm}
	a_hat<- integrate(f_post,0,Inf)$value
p_post<- psil(sgl,a_hat)
}
else
stop("model specified not available.")
}

# Method 2a: estimation of the mean of the a posteriori probability to do a DLT
## Probability estimation without singleton calculations
else if (method == "ppostp") {
# Power model
if (model == "power") {
psi <- psip
sgl <- prior
	norm<- integrate(function(a){Lp(a,data=Dk,sgl=sgl)*fp(a, sd)},-Inf,Inf)$value
	f_post<- function(a){Lp(a,data=Dk,sgl=sgl)*fp(a, sd)/norm}
p_post<- sapply(sgl,FUN = function(s){integrate(function(a){psip(s,a)*f_post(a)},-Inf,Inf)$value})
}

# Tangent model
else if (model == "tangent") {
sgl <- ait1(prior)
	norm<- integrate(function(a){Lt(a,data=Dk,sgl=sgl)*ft(a)},0,Inf)$value
	f_post<- function(a){Lt(a,data=Dk,sgl=sgl)*ft(a)/norm}
p_post <- sapply(sgl,FUN = function(s){integrate(function(a){psit(s,a)*f_post(a)},0,Inf)$value})
}

# Logistic model
else if (model == "logistic") {
sgl <- ail1(prior)
	norm<- integrate(function(a){Ll(a,data=Dk,sgl=sgl)*fl(a)},0,Inf)$value
	f_post<- function(a){Ll(a,data=Dk,sgl=sgl)*fl(a)/norm}
p_post <- sapply(sgl,FUN = function(s){integrate(function(a){psil(s,a)*f_post(a)},0,Inf)$value})
}
else
stop("model specified not available.")
}

# Method 2a: estimation of the mean of the a posteriori probability to do a DLT
## Probability estimation with singleton calculations
else if (method == "pposts") {
# Power model
if (model == "power") {
sgl <- aip(prior, sd = sd)
	norm<- integrate(function(a){Lp(a,data=Dk,sgl=sgl)*fp(a, sd)},-Inf,Inf)$value
	f_post<- function(a){Lp(a,data=Dk,sgl=sgl)*fp(a, sd)/norm}
p_post<- sapply(sgl,FUN = function(s){integrate(function(a){psip(s,a)*f_post(a)},-Inf,Inf)$value})
}

# Tangent model
else if (model == "tangent") {
sgl <- ait2(prior)
# Method 2: estimation of the mean of the a posteriori probability to do a DLT
	norm<- integrate(function(a){Lt(a,data=Dk,sgl=sgl)*ft(a)},0,Inf)$value
	f_post<- function(a){Lt(a,data=Dk,sgl=sgl)*ft(a)/norm}
p_post <- sapply(sgl,FUN = function(s){integrate(function(a){psit(s,a)*f_post(a)},0,Inf)$value})
}

# Logistic model
else if (model == "logistic") {
sgl <- ail2(prior)
	norm<- integrate(function(a){Ll(a,data=Dk,sgl=sgl)*fl(a)},0,Inf)$value
	f_post<- function(a){Ll(a,data=Dk,sgl=sgl)*fl(a)/norm}
p_post <- sapply(sgl,FUN = function(s){integrate(function(a){psil(s,a)*f_post(a)},0,Inf)$value})
}
else
stop("model specified not available.")
}
else
stop("method specified not available.")


# Option nextdose
## Inferior to target
if (nextlevel == "utarget"){
nextdose <- ifelse(p_post[1] < target, Dk$dose[max(which(p_post < target))], NA)
}
## Near target with euclidian distance
else if (nextlevel == "ntarget"){
nextdose <- Dk$dose[which.min(abs(p_post - target))]
}
else
stop("nextlevel specified not available.")


# Dose escalation et de-escalation options
if (upskipping)
nextdose <- ifelse(nextdose > lastdose, lastdose + 1, nextdose)
if (downskipping)
nextdose <- ifelse(nextdose < lastdose, lastdose - 1, nextdose)
}

# Likelihood maximum
if (approach == "mle") {
mtd <- NA
npt<- Dk$npt
ndlt<- Dk$ndlt
if (sum(npt) == 0) {
nextdose <- 1
p_post <- prior
}
if (sum(npt) > 0) {
if (sum(ndlt) > 0){
sgl <- prior
a_hat <- optimize(function(a){Lp(a,Dk,sgl)}, c(-10, 10), tol = 1e-04, maximum = TRUE)$max
p_post <- psip(sgl,a_hat)

# Next dose
if (nextlevel == "utarget"){
nextdose <- ifelse(p_post[1] < target, Dk$dose[max(which(p_post < target))], NA)
}
else if (nextlevel == "ntarget"){
nextdose <- Dk$dose[which.min(abs(p_post - target))]
}
else
stop("nextlevel specified not available.")}

else {
if (!(lastdose %in% Dk$dose) | Dk$npt[lastdose] == 0)
stop("lastdose specified not available.")
nextdose <- lastdose + 1
p_post <- prior
}
}
}

# Stop criterion
if (sum(Dk$npt) >= nptmax | !(nextdose %in% Dk$dose) | max(Dk$npt) == nmaxdose){
mtd <- ifelse(nextdose == (length(Dk$dose) + 1), lastdose, nextdose)
nextdose <- NA
}
if (nextdose %in% Dk$dose)
if (Dk$npt[nextdose] == nmaxmtd) {
mtd <- nextdose
nextdose <- NA
}

list(nextdose = nextdose, mtd = mtd, prob = p_post)      
}

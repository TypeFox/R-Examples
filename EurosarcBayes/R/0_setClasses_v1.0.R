

trialDesign_binom_one=setClass("trialDesign_binom_one",slots=c(reviews="numeric", success = "numeric", failure="numeric", eta = "numeric", zeta = "numeric", alpha = "numeric", power = "numeric", exp.p0 = "numeric", exp.p1 = "numeric", p0 = "numeric", p1 = "numeric"))


trialDesign_binom_two=setClass("trialDesign_binom_two",slots=c(reviews="numeric",data="data.frame",cutpoints="data.frame",precision="numeric",decision="list",post.futility="numeric",post.efficacy="numeric", post.toxicity ="numeric",post.no_toxicity="numeric",graph="list"))


binom_two_bryantday=setClass("binom_two_bryantday",slots=c(optimal="data.frame",minmax="data.frame",all.fit="data.frame"))
binom_two_singlestage=setClass("binom_two_singlestage",slots=c(optimal="data.frame",output="data.frame"))


finetab <-
function(nom_cov, t_id, c_id) {
	Cat = nom_cov
	Units = rep(0, length(Cat))
	Units[t_id] = 1
	Units[c_id] = 2
	tab = table(Cat, Units)[, 2:3]
	colnames(tab) = c("T", "C")
	tab
}

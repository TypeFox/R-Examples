mod.lt <- function(child.value, child.mort=4, e0.target=NULL, adult.mort=NULL, sex="female", alpha=0){
	#data(MLTobs)
	class <- hmd.DA(x=child.value, sex=sex, child.mort=child.mort, adult.mort=adult.mort)
	class <- as.numeric(class$classification)	
	# If e0.target is provided (not null), then get alpha from that
	# If e0.target is NULL (default), use alpha  
	a.out <- if (is.null(e0.target)) alpha else alpha.e0(pattern=class, e0.target=e0.target, sex=sex)

	mx.out <- mortmod(pattern=class, alpha=a.out, sex=sex)
	lt.out <- lt.mx(nmx=exp(mx.out), sex=sex)
	return(structure(c(lt.out, list(alpha=a.out, sex=sex, family=class)), class='LifeTable'))
	}
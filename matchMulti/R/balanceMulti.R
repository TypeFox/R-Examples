balanceMulti <-
function(match.obj, student.cov = NULL, school.cov = NULL){
	treatment <- match.obj$treatment
	
	#if no school or student covariates are provided, compute balance on all
	if (is.null(student.cov)) student.cov = setdiff(colnames(match.obj$raw), c(treatment, match.obj$school.id))
	if(is.null(school.cov)) school.cov = setdiff(colnames(match.obj$raw), c(treatment, match.obj$school.id))
	
	#student balance
	student.bal <- balanceTable(match.obj$raw[c(student.cov,treatment)], match.obj$matched[c(student.cov, treatment)], treatment)
	
	#school balance
	schools.raw <- students2schools(match.obj$raw, c(school.cov, treatment), match.obj$school.id)
	schools.matched <- schools.raw[which(schools.raw[[match.obj$school.id]] %in% as.vector(match.obj$school.match)),]
	
	#schools.treat <- which(schools.raw[[treatment]] == 1)
	#schools.ctrl <- which(schools.raw[[treatment]] == 0)
	#match.clean <- match.obj$school.match$matches[,1,drop = FALSE]
	#match.clean <- match.clean[which(!is.na(match.clean)),,drop = FALSE]
	#schools.matched <- schools.raw[c(schools.treat[as.numeric(rownames(match.clean))], schools.ctrl[match.clean]),]
	
	id.idx <- which(colnames(schools.raw) == match.obj$school.id)
	
	school.bal <- balanceTable(schools.raw[-id.idx], schools.matched[-id.idx], treatment)
	
	return(list('students' = student.bal, 'schools' = school.bal))
}

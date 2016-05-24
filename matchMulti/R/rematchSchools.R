rematchSchools <-
function(match.out, students, school.fb = NULL, verbose = FALSE, keep.target = NULL, school.penalty = NULL, tol = 1e-3){

	#Are school fine balance constraints nested appropriately?
	if(!is.null(school.fb) && length(school.fb) > 1){
		for(i in c(1:(length(school.fb)-1))) {
			if(!all(school.fb[[i]] %in% school.fb[[i+1]])){
					stop('Each element of school.fb must contain all variables listed in previous elements')
				}	
		}	
	}
	
	if(!is.null(school.penalty) && school.penalty <= 0) stop('School penalty must be positive') 	

	if(is.null(match.out$student.matches)) stop('match.out has no first-stage information: rerun matchMulti with save.first.stage = TRUE')

	treatment <- match.out$treatment
	school.id <- match.out$school.id

	########## REMATCH ###########	

	school.match <- matchSchools(match.out$student.matches$schools.matrix, students, treatment, school.id, school.fb, school.penalty, verbose, tol = tol) 

	#if keep.target is provided, iterate until we get a match that keeps (close to) the desired number of schools
	if(!is.null(keep.target)) {
		treat.schools <- unique(students[[school.id]][students[[treatment]] == 1])
		if (keep.target <= 0 || keep.target > length(treat.schools) || keep.target %% 1 != 0) stop("keep.target must be a positive integer no greater than the total number of treated schools") 
		STARTVAL <- 1000
		MAXITER <- 1000
		SCALE_FACTOR <- 10
		STOPRULE <- 1
		ubound <- Inf
		lbound <- 0
		cur <- school.penalty
		if (is.null(cur)) {
			cur <- STARTVAL	
			school.match <- matchSchools(match.out$student.matches$schools.matrix, students, treatment, school.id, school.fb, penalty = cur, verbose, tol = tol) 
		} 
		next.match <- school.match		
		for (i in 1:MAXITER) {
			nkeep <- length(intersect(treat.schools, next.match[,1]))
			if (nkeep == keep.target) {
				if (verbose) print(paste('Iteration',i,': reached target number of schools with penalty of',cur))
				break
			}
			if (nkeep > keep.target) {
				if (verbose) print(paste('Iteration',i,': too many schools retained with penalty of',cur))			
			ubound <- cur
				cur <- lbound + (ubound - lbound)/2
			}
			if (nkeep < keep.target) {
				if (verbose) print(paste('Iteration',i,': not enough schools retained with penalty of',cur))			
				lbound <- cur
				if (is.finite(ubound)) {
					cur <- lbound + (ubound - lbound)/2
				} else {
					cur <- SCALE_FACTOR*cur
				}
			}
			if(ubound - lbound < STOPRULE) break
			next.match <- matchSchools(match.out$student.matches$schools.matrix, students, treatment, school.id, school.fb, cur, verbose,tol = tol) 
		}
		school.match <- next.match
	}
	
	########### OUTPUT ###########

	out.match <- assembleMatch(match.out$student.matches$student.matches, school.match, school.id, treatment)
	
	#record dropped schools and students
	drop.obj <- list()
	dropped.schools <- setdiff(unique(students[[school.id]]), school.match)
	treat.table <- table(students[[school.id]], students[[treatment]])
	treat.table <- treat.table[match(as.character(dropped.schools),rownames(treat.table)),]
	drop.z <- apply(treat.table, 1, function(x) which(x >0 )-1)
	drop.obj$schools.t <- dropped.schools[drop.z == 1]
	drop.obj$schools.c <- dropped.schools[drop.z == 0]
	
 	student.count.df <- data.frame('school.id' = c(students[[school.id]],out.match[[school.id]]), 'before.after' = c(rep(0,nrow(students)), rep(1, nrow(out.match))))
 	schooltab <- table(student.count.df$school.id, student.count.df$before.after)
 	schooltab[,2] <- schooltab[,1] - schooltab[,2]
 	colnames(schooltab) <- c('Student Count', 'Dropped')
	drop.obj$students.by.school <- schooltab
	
	out.obj <- match.out
	out.obj$matched <- out.match
	out.obj$school.match <- school.match
	out.obj$dropped <- drop.obj
	out.obj
}

matchMulti <-
function(data, treatment, school.id, match.students = TRUE, student.vars = NULL, school.fb = NULL, verbose = FALSE, student.penalty.qtile = 0.05, min.keep.pctg = 0.8, school.penalty = NULL, save.first.stage = TRUE, tol = 1e-3){
	
	students <- data
	##### Validate input #####
	
	#are variable names present in data frame?
	if (!(treatment %in% colnames(students))) {
		stop(paste('Treatment variable', treatment,'not found'))
	}
	if (!(school.id %in% colnames(students))) {
		stop(paste('School ID variable', school.id,'not found'))
	}	
	if (!is.null(student.vars) && !all(student.vars %in% colnames(students))) {
		bad.vars <- student.vars[!(student.vars %in% colnames(students))]
		stop(paste('Student variable', bad.vars,'not found\n'))
	}
	if (!is.null(school.fb) && !all(unlist(school.fb) %in% colnames(students))) {
		bad.vars <- unique(unlist(school.fb)[!(unlist(school.fb) %in% colnames(students))])
		stop(paste('School variable', bad.vars,'not found\n'))
	}	
	
	#do schools nest within treatment categories?
	treat.tab <- table(students[[school.id]], students[[treatment]]) 
	if(any(apply(treat.tab, 1, min) > 0)) {
		stop('Some schools contain both treated and control students')
	}
	
	#Are school fine balance constraints nested appropriately?
	if(!is.null(school.fb) && length(school.fb) > 1){
		for(i in c(1:(length(school.fb)-1))) {
			if(!all(school.fb[[i]] %in% school.fb[[i+1]])){
					stop('Each element of school.fb must contain all variables listed in previous elements')
				}	
		}	
	}
	
	if(student.penalty.qtile < 0 || student.penalty.qtile > 1) stop('Student penalty quantile must be in [0,1]') 

	if(!is.null(school.penalty) && school.penalty <= 0) stop('School penalty must be positive') 	
	
	if(match.students && is.null(student.vars)) stop('Cannot match students unless student variables are specified')
	
	
	########### MATCHING ###########
		
	#student matches in all school pairings	
	student.matches <- matchStudents(students, treatment, school.id, match.students, student.vars, verbose, student.penalty.qtile, min.keep.pctg)

	school.match <- matchSchools(student.matches$schools.matrix, students, treatment, school.id, school.fb, school.penalty, verbose, tol = tol) 


	########### OUTPUT ###########

	out.match <- assembleMatch(student.matches$student.matches, school.match, school.id, treatment)
	
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
	
	out.obj <- list('raw' = students, 'matched' = out.match, 'school.match' = school.match,'dropped' = drop.obj, 'school.id' = school.id, 'treatment' = treatment)
	
	if(save.first.stage){ 
		out.obj$student.matches <- student.matches
	}
	out.obj
}

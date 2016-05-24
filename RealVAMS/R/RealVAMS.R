RealVAMS <-
function(score.data,outcome.data,persistence = "CP", school.effects=FALSE, REML = TRUE,
score.fixed.effects=formula(~as.factor(year) + 0),
outcome.fixed.effects=formula(~1), max.iter.EM = 10,
outcome.family = binomial(link = "probit"),tol1 = 1e-7,
max.PQL.it = 30, pconv = .Machine$double.eps*1e9,
var.parm.hessian = TRUE, verbose = TRUE){

control <- list(persistence=persistence,school.effects=school.effects, REML=REML,score.fixed.effects=score.fixed.effects,outcome.fixed.effects=outcome.fixed.effects,
max.iter.EM=max.iter.EM,outcome.family=outcome.family,tol1=tol1,max.PQL.it=max.PQL.it,pconv=pconv,verbose=TRUE)
Z_mat<-score.data
if(control$school.effects) Z_mat$schoolID<-factor(Z_mat$schoolID)
B.mat<-outcome.data
control$NR<-FALSE
control$max.iter.EM<-c(control$max.iter.EM,5)
control$EM.to.NR<-50
control$hes.method<-"simple"
control$hessian<-var.parm.hessian
Z_mat <- Z_mat[order(Z_mat$year, Z_mat$teacher), ]
max.PQL.it <- control$max.PQL.it
huge.flag <- FALSE
if (class(try(na.fail(Z_mat[, !(names(Z_mat) %in% c("teacher", "y"))]), silent = TRUE)) == "try-error") {
   cat("*Error: NA values present.\n*NA values are allowed for the 'teacher; and 'y' variables, but no others.\n*Please remove these observations from your data frame.\n")
   flush.console()
   return(0)
}
     if (identical(sort(unique(outcome.data$r)),c(0,1)) | class(outcome.data$r) != "integer") {
        cat("*Error: outcome.data$r should be a numeric column of 0's and 1's, with 1 representing the positive class")
        flush.console()
        return(0)
    }
if (class("student.side") != "character" | class("persistence") != "character") {
    cat("*Error: student.side and persistence must be characters (using quotation marks)")
    flush.console()
    return(0)
}
original_year <- Z_mat$year
tyear <- factor(Z_mat$year, ordered = TRUE)
pyear <- original_year
nyear.score <- nlevels(tyear)
for (i in 1:nyear.score) {
    pyear[Z_mat$year == levels(tyear)[i]] <- i
}
Z_mat$year <- factor(pyear, ordered = TRUE)
key <- unique(cbind(as.numeric(original_year), as.numeric(pyear)))
key <- key[order(key[, 2]), ]
control$key <- key
for (j in 1:nyear.score) {
    Z_mat[Z_mat$year == levels(Z_mat$year)[j], ]$teacher <- paste(Z_mat[Z_mat$year == levels(Z_mat$year)[j], ]$teacher, "(year", levels(Z_mat$year)[j], ")", sep = "")
}

    if (persistence=="CP"|persistence=="VP") {
        res <- vp_cp(Z_mat=Z_mat, B.mat=B.mat, control=control)
    }    else {
        cat("Error in specification of student.side or persistence\n")
        return(0)
    }


    res<-c(res,list(key=key))
#    res$teach.effects$teacher_year<-key[match(res$teach.effects$teacher_year,key[,2]),1]
#    res$teach.effects$effect_year<-key[match(res$teach.effects$effect_year,key[,2]),1]
    class(res) <- "RealVAMS"
    return(res)

}

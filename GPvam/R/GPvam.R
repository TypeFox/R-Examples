GPvam <-
function (vam_data, fixed_effects = formula(~as.factor(year) + 
    0), student.side = "R", persistence="GP", max.iter.EM = 1000, 
    tol1 = 1e-07, hessian = FALSE, hes.method = "simple", verbose = TRUE) 
{
    control<-list(max.iter.EM=max.iter.EM,tol1=tol1,hessian=hessian,hes.method=hes.method,verbose=verbose,persistence=persistence)
    Z_mat <- vam_data
    if (class(try(na.fail(Z_mat[, !(names(Z_mat) %in% c("teacher", 
        "y"))]), silent = TRUE)) == "try-error") {
        cat("*Error: NA values present.\n*NA values are allowed for the 'teacher; and 'y' variables, but no others.\n*Please remove these observations from your data frame.\n")
        flush.console()
        return(0)
    }
    if (student.side == "G" & persistence!="GP") {
    cat("*Error: student.side=\"G\" may only be paired with persistence=\"GP\"")
        flush.console()
        return(0)
    }
    if (class("student.side")!="character"|class("persistence")!="character") {
    cat("*Error: student.side and persistence must be characters (using quotation marks)")
        flush.console()
        return(0)
    }
    original_year<-Z_mat$year
    tyear<- factor(Z_mat$year,ordered=TRUE)
    pyear<-original_year
    nyear <- nlevels(tyear)
    for(i in 1:nyear){
     pyear[Z_mat$year==levels(tyear)[i]]<-i
    }
    Z_mat$year<-factor(pyear,ordered=TRUE)
    key<-unique(cbind(as.numeric(original_year),as.numeric(pyear)))
    key<-key[order(key[,2]),]
    control$key<-key
    for (j in 1:nyear) {
        Z_mat[Z_mat$year == levels(Z_mat$year)[j], ]$teacher <- paste(Z_mat[Z_mat$year == 
            levels(Z_mat$year)[j], ]$teacher, "(year", levels(Z_mat$year)[j], ")", sep = "")
    }
    if (student.side == "R" & persistence=="GP") {
        res <- GP.un(Z_mat, fixed_effects, control)
    }
    else if (student.side == "G"& persistence=="GP") {
        res <- GP.csh(Z_mat, fixed_effects, control)
    } 
    else if (student.side == "R"& persistence=="rGP") {
        res <- rGP.un(Z_mat, fixed_effects, control)
    } 
    else if (student.side == "R"& (persistence=="VP"|persistence=="CP"|persistence=="ZP")) {
        res <- VP.CP.ZP.un(Z_mat, fixed_effects, control)
    } 
    else {
        cat("Error in specification of student.side or persistence\n")
        return(0)
    }
    #res$teach.effects$year
    res<-c(res,list(key=key))
    res$teach.effects$teacher_year<-key[match(res$teach.effects$teacher_year,key[,2]),1]
    res$teach.effects$effect_year<-key[match(res$teach.effects$effect_year,key[,2]),1]
    class(res) <- "GPvam"
    return(res)
}

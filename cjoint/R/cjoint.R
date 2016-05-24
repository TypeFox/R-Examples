## cjoint: An R Package for estimating Average Marginal Component-specific Effects from conjoint survey experiments
## July, 2015
## Anton Strezhnev, Elissa Berwick, Jens Hainmueller, Daniel Hopkins, Teppei Yamamoto

#############################
# Clustered standard errors
#############################

cluster_se_glm <- function(model, cluster){
    
    #  Drop unused cluster indicators, if cluster var is a factor
    if (class(cluster) == "factor") {
        cluster <- droplevels(cluster)
    }
    
    if (nrow(model.matrix(model)) != length(cluster)) {
        stop("check your data: cluster variable has different N than model - you may have observations with missing data") 
    }
    
    M <- length(unique(cluster))
    N <- length(cluster)           
    K <- model$rank
    
    ## if(M<50) {
    ##     warning("Fewer than 50 clusters, variances may be unreliable (could try block bootstrap instead).")
    ## }
   
    dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
    uj  <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum));
    rcse.cov <- dfc * sandwich(model, meat. = crossprod(uj)/N)
    return(rcse.cov)
}

###################################################################
## function for removing ALL whitespaces
## from elements of vector "vec"
## trailing and leading stripped
## white space in middle replaced with "repl"
## also removes special characters that mess up grep and formula
###################################################################

clean.names <- function(vec,repl="_") {
    #define meta characters to be removed (excl. ^, *)
    meta <- c("[","$","{","(","\\","+",")","|","?","<",">",".")
    #add escapes
    meta2 <- sapply(meta,USE.NAMES = F,function(x) paste(c("\\",x),collapse=""))
    #add some other special characters to be removed
    special <- c(meta2,"!",",",";","'","`","/","-")
    #make regexp version
    special.regexp <- paste(c("[",special,"]"),collapse = "")
    # clean names in vec
    sapply(vec,USE.NAMES = F,function(x) {
        #remove special characters
        x <- gsub(special.regexp,"",x)
        #remove leading and trailing whitespaces
        x <- gsub("^\\s+|\\s+$", "", x)
        #replace middle space with repl
        gsub("\\s+",repl,x)
    })
}

#########################
## amce function
#########################

amce <- function(formula, data, design="uniform", respondent.varying = NULL, subset=NULL, respondent.id=NULL, cluster=TRUE, na.ignore=FALSE, weights=NULL, baselines = NULL) {
    
###### Formula and Variables

    #we will split the formula into separate lists
    #unique_vars = all VARIABLES in formula
    #respondent_vars = VARIABLES varying by respondent
    #profile_vars = VARIABLES varying by profile
    #orig_effects = all EFFECTS in formula 
    #profile_effects = profile varying EFFECTS in formula
    #user inputted names for above end with "_user"
    
##### Parse formula, clean variable and input names

    formula_user <- formula
    #all variables in formula
    formula_char_user <- all.vars(formula)
    #lists with original names of variables and levels
    user_names <- list()
    user_levels <- list()
    for (char in formula_char_user) {
        user_names[[clean.names(char)]] <- char
        if (class(data[[char]]) == "factor") {
            old_names <- names(user_levels)
            user_levels <- c(user_levels,levels(data[[char]]))
            names(user_levels) <- c(old_names,clean.names(levels(data[[char]])))
        }
    }

    #make sure no duplicates after spaces and special characters removed
    formula_char <- clean.names(formula_char_user)    
    #if this makes for non-unique names, stop
    if(length(unique(formula_char)) != length(formula_char)) {
        stop("Error: Variable names must be unique when whitespace and meta-characters are removed. Please rename.")
    }
        
    #separate dependent and independent variables and clean
    y_var <- clean.names(formula_char_user[1])
    #identify ALL original effects; will add in missing base terms automatically
    orig_effects <- clean.names(attr(terms(formula_user),"term.labels"))
    #formula sorting part I: sort non-interaction terms and put them first
    orig_effects <- c(sort(orig_effects[!grepl(":",orig_effects)]), orig_effects[grepl(":",orig_effects)])
    #combine with "+"
    vars_plus <- paste(orig_effects,collapse = " + ")
    #then remake formula 
    form <- formula(paste(c(y_var,vars_plus),collapse = "~"))
    orig_effects <- attr(terms(form),"term.labels")

    #find missing base terms
    full_terms <- attr(terms(formula(paste(y_var,paste(sapply(orig_effects,function(x) gsub(":","*",x)),collapse=" + "),sep=" ~ "))),"term.labels")
    # add in any missing base terms for interactions
    missing_terms <- full_terms[!is.element(full_terms,orig_effects)]
    if (length(missing_terms > 0)) {
        orig_effects <- c(orig_effects,missing_terms)
        warning("Missing base terms for interactions added to formula")
    }

    #formula sorting redux: sort non-interaction terms and put them first
    orig_effects <- c(sort(orig_effects[!grepl(":",orig_effects)]), orig_effects[grepl(":",orig_effects)])
    #combine with "+"
    vars_plus <- paste(orig_effects,collapse = " + ")
    #then remake formula 
    form <- formula(paste(c(y_var,vars_plus),collapse = "~"))
    orig_effects <- attr(terms(form),"term.labels")

    #unique variables only (no interactions)
    unique_vars <- clean.names(rownames(attr(terms(form),"factor"))[-1])
    #respondent variables
    respondent_vars <- clean.names(respondent.varying)
    #profile variables
    profile_vars <- unique_vars[!is.element(unique_vars,respondent_vars)]

    #identify the REQUESTED profile effects and respondent effects (if any)
    if (length(respondent_vars) > 0) {
        #identify profile only effects
        profile_effects <- unlist(sapply(orig_effects,USE.NAMES = F,function(x) {
            y <- strsplit(x,":")[[1]]
            if (!any(is.element(y,respondent_vars))) x
        }))
        #terms containing a respondent var
        resp_only <- unlist(sapply(orig_effects,USE.NAMES = F, function(x) {
            y <- strsplit(x,":")[[1]]
            if(any(is.element(y,respondent_vars))) x
        }))
        #things that respondent vary is interacted with
        resp_mod <- unlist(sapply(resp_only,USE.NAMES = F,function(x) {
            y <- strsplit(x,":")[[1]]
            vars <- y[!is.element(y,respondent_vars)]
            if (length(vars) > 0) paste(vars,collapse = ":")
        }))
        resp_effects <- c(resp_mod,resp_only)  
    } else {
        profile_effects <- orig_effects
        resp_effects <- NULL
    }

    ### Extra name cleaning

    #cleaning additional inputs
    if (!is.null(respondent.id)) respondent.id <- clean.names(respondent.id)
    if (!is.null(weights)) weights <- clean.names(weights)
    if (!is.null(baselines)) {
        names(baselines) <- clean.names(names(baselines))
        baselines <- lapply(baselines,function(x) clean.names(x))
    }

    #cleaning within data 
    colnames(data) <- clean.names(colnames(data))
    data <- data.frame(data) #in case of dplyr etc.
    for (var in colnames(data)) {
        if (class(data[[var]]) == "factor") {
            clean.labels <- clean.names(levels(data[[var]]))
            if (length(unique(clean.labels)) != length(clean.labels)) {
                stop (paste("Error: levels of variable", var, "when whitespace and meta-characters are removed. Please rename."))
            }
            data[[var]] <- factor(data[[var]],levels=levels(data[[var]]),labels=clean.names(levels(data[[var]])))
        }
    }
    
#######  Sanity Checks Re: Data

    # Are variables in data?
    for(var in formula_char) {
        if(!(var %in% colnames(data))) {
            stop(paste("Error:", var, "not in 'data'"))
        }
    }

    # Make sure non-respondent varying are factors
    for (var in profile_vars) {
        if (class(data[[var]]) != "factor") {
            data[[var]] <- as.factor(data[[var]])
            warning(paste(c("Warning: ",var," changed to factor"),collapse=""))
        }
    }
    
    # Is there missing data?
    if(na.ignore == FALSE){
      for(variab in formula_char){
        if (sum(is.na(data[[variab]])) != 0 ){
          stop(paste("Error:", variab, "has missing values in 'data'"))
        }
      }
    }

    # Is the respondent varying characteristic even in the formula obj?
    if (!is.null(respondent_vars)) {
        for (var in respondent_vars) {
            found <- 0
            for (formulavars in formula_char) {
                if (var == formulavars) {
                    found <- 1
                }
            }
            if (found == 0) {
                stop(paste("Error:", var, "is specified in respondent.varying, but is not in the formula"))
            }
        }
    }

    # Check whether outcome variable is a binary 0-1 or numeric
    if (!is.numeric(data[[y_var]]) & !is.integer(data[[y_var]])) {
        stop(paste("Error:", y_var, "is not numeric or integer"))
    }

    # Are the user-supplied desired baselines in the data?
    if (!is.null(baselines)) {
        for(var in names(baselines)) {
            if(!(baselines[[var]] %in% data[[var]])) {
                stop(paste("Error: user supplied baseline",baselines[[var]],"is not a level of",var))
            }
        }      
    }

#######  Sanity Checks Re: clustering

    # If respondent.id is NULL make sure cluster is FALSE
    if(is.null(respondent.id) & cluster == TRUE) {
        warning("respondent.id is NULL - setting cluster to FALSE. Please specify a respondent.id variable if you want to estimate clustered standard errors")
        cluster <- FALSE
    }

    # If (cluster = TRUE), make sure respondent.id is specified and in data
    if (cluster == TRUE) {
        if (is.null(respondent.id)) {
            stop('Error: Must specify a respondent.id if cluster = TRUE')
        } else if (!(respondent.id %in% colnames(data))) {
            stop('Error: respondent.id not in data')
        }
    }

######## Sanity checks re: weights

    if(!is.null(weights)) {
        #is the weight var in the data matrix?
        if (!weights %in% colnames(data)) {
            stop('Error: weights not in data')
        }
        #are weights uniform? if so turn off weights
        if(length(unique(data[[weights]])) == 1) {
            weights <- NULL 
        }
    }

##### Sanity Checks re: design matrix
    
    # If design is already conjointDesign object, proceed to relevant sanity checks
    if (class(design) == "conjointDesign") {
        # Remove whitespaces etc from dimension names of design array 
        names(dimnames(design$J)) <- clean.names(names(dimnames(design$J)))
        dimnames(design$J) <- lapply(dimnames(design$J),function(x) clean.names(x))  
        #and design dependencies
        names(design$depend) <- clean.names(names(design$depend))
        design$depend <- lapply(design$depend,function(x) clean.names(x))  
        #Now check to make sure profile varying attributes are in conjointDesign
        for (eff in profile_vars) {   
            if (!(eff %in% names(dimnames(design$J)))) {
                stop(paste("Error:", var, "not in 'design' object"))
            }
        }      
        #Check to make sure conjointDesign attributes are in data and level names match
        for (eff in names(dimnames(design$J))) {
            if (!(eff %in% colnames(data))){
                stop(paste("Error: attribute", eff, "in 'design' object is not in 'data'"))
            } else {
        # Check all level names for the attribute in dataset appear in design
                for (lev in levels(as.factor(data[[eff]]))) {
                    if (!(lev %in% dimnames(design$J)[[eff]])) {
                        stop(paste("Error: factor level", lev, "of attribute", eff, "not in 'design' object"))
                    }
                }
            }
        }    
    } else if (design == "uniform") {    
        # else if design == "uniform", create J-dimensional array 
        design <- list()        
        # Determine dimensions
        # And create J matrix with uniform probabilities across all vars
        design.dim <- vector(length=length(profile_vars))
        dim_list <- list()
        for (i in 1:length(profile_vars)) {
            design.dim[i] <- length(unique(data[[profile_vars[i]]]))
            dim_list[[i]] <- levels(factor(data[[profile_vars[i]]]))
        }
        names(dim_list) <- profile_vars
        design$J <- array(1/prod(design.dim), dim=design.dim, dimnames=dim_list)
        design$depend <- compute_dependencies(design$J)       
    } else {
         #if neither uniform nor conjointDesign, error
        stop('Error: argument \'design\' must be a valid character string ("uniform") or a conjointDesign object')   
    }

####### Subsetting data    

    if (is.null(subset)) {
        data <- data 
    } else {
        if (class(subset) == "logical") {
            if (length(subset) == nrow(data)) {
                data <- subset(data, subset) 
            } else {
                warning("Warning: invalid argument to 'subset' - must be the same length as the number of rows in data")
            }
        } else {
            warning("Warning: invalid argument to 'subset' - must be a logical")
        }
    }

###### Adjust baselines if given    

    if (!is.null(baselines)) {
        for (var in names(baselines)) {
            data[[var]] <- factor(data[[var]])
            data[[var]] <- relevel(data[[var]], baselines[[var]])
        } 
    }

####### Adding relevant interaction terms to model

    #If there are any dependencies-- only for profile-varying!
    if(any(profile_vars %in% names(design$depend))) {
        #initialize full interaction set
        depend_vars <- c()     
        #loop over effects with dependencies
        for(eff in profile_vars[profile_vars %in% names(design$depend)]) {
            #initialize interaction set for given variable
            inter <- c()
            #identify higher order occurences of variable
            #make sure it's just that variable, not followed or begun by "_"
            eff_all <- grep(paste(c(":",eff,"(?!_)","|",eff,":(?<!_)"),collapse=""),
                            orig_effects,value=T,perl=T)
            #if you find some, break up, sort and replace ":" with "*"
            if (length(eff_all) > 0) {
                eff_all <- sapply(strsplit(eff_all,":"),function(x) paste(sort(x),collapse="*"))
            }
            #combine with lower order (main effect)
            eff_all <- c(eff,eff_all)
            #for each occurrence, create interaction
            inter <- sapply(eff_all,USE.NAMES = F,function(x) {
                #get conditioning set
                T_r <- design$depend[[eff]]
                #make factors
                for (t in T_r){
                    data[[t]] <- as.factor(data[[t]])
                }
                #combine name and dependency
                T_r_d <- c(x,T_r)
                #make interaction term
                paste(T_r_d,collapse="*")
            })
            #add to list
            depend_vars <- c(depend_vars,inter)
        }      
        #drop repeats
        depend_vars <- unique(depend_vars)
        #add to formula
        form_full <- formula(paste(c(form,depend_vars),collapse = " + "))
    } else {
      form_full <- form
  }

    #all variables to be run
    all_run_vars <- attr(terms(form_full),"term.labels")
   #formula sorting redux: sort non-interaction terms and put them first
    all_run_vars <- c(sort(all_run_vars[!grepl(":",all_run_vars)]), all_run_vars[grepl(":",all_run_vars)])
    #combine with "+"
    vars_plus <- paste(all_run_vars,collapse = " + ")
    #then remake formula 
    form_full<- formula(paste(c(y_var,vars_plus),collapse = "~"))
    all_run_vars <- attr(terms(form_full),"term.labels")
    
####### If there are respondent varying terms, split into two formulas
######## One contains only profile effects
######## Second is full formula

    if (length(respondent_vars) > 0) {
        ### profile only formula
        #remove those involving respondent things
        prof_only <- unlist(sapply(all_run_vars,function(x) {
            y <- strsplit(x,":")[[1]]
            if(!any(is.element(y,respondent_vars))) x
        }))
        prof_only_plus <- paste(prof_only,collapse = " + ")
        #formula with profile only
        form_prof <- paste(all.vars(form_full)[1],prof_only_plus,sep=" ~ ")
        form_prof <- formula(form_prof)        
    } else {
        #otherwise use full formula
        form_prof <- form_full
    }
    all_prof <- attr(terms(form_prof),"term.labels")
    if (any(!is.element(all_prof,all_run_vars))) {
        warning("Warning: mismatch of term names between full formula and profile formula")
    }
    
####### Running OLS

    #run model(s)-- if using weights, use tools from "survey" package
    #otherwise run usual lm function
    if (is.null(weights)) {
        lin.mod.prof <- lm(form_prof, data=data)
        if (length(respondent_vars) > 0) {
            lin.mod.full <- lm(form_full, data=data)
        } else {
            lin.mod.full <- NULL
        }
    } else {
        if (cluster) {
            out.design <- svydesign(ids = data[[respondent.id]], weights = data[[weights]], data=data)
        } else {
            out.design <- svydesign(ids = ~0, weights = data[[weights]], data=data)
        }
        lin.mod.prof <- svyglm(form_prof, data=data, design = out.design)
        if (length(respondent_vars) > 0) {
            lin.mod.full <- svyglm(form_full, data=data, design = out.design)
        } else {
            lin.mod.full <- NULL
        }
    }

    #If there's missing data - flag it
    #if (na.ignore == TRUE & !is.null(lin.mod.full$na.action)) {
     #   stop(paste("Error: Observations with missing data in 'data'"))
    #}
  
    #Get sample size
    sample_size_prof <- length(lin.mod.prof$residuals)
    if (length(respondent.varying) > 0) {
        sample_size_full <- length(lin.mod.full$residuals)
    } else {
        sample_size_full <- NULL
    }
    
    #Compute vcov of OLS
    if (is.null(weights) & cluster == TRUE) {
        #clusters but no weights
        vcov_mat_prof <- cluster_se_glm(lin.mod.prof, data[[respondent.id]])
        if (length(respondent.varying) > 0) {
            vcov_mat_full <- cluster_se_glm(lin.mod.full, data[[respondent.id]])
        } else {
            vcov_mat_full <- NULL
        }
    } else if (!is.null(weights)) {
        #weights with or without cluster
        vcov_mat_prof <- vcov(lin.mod.prof)
        if (length(respondent_vars) > 0) {
            vcov_mat_full <- vcov(lin.mod.full)
        } else {
            vcov_mat_full <- NULL
        }        
    } else {
    #Not clustered or weighted
        vcov_mat_prof <- vcovHC(lin.mod.prof,type="HC2")
        if (length(respondent.varying) > 0) {
            vcov_mat_full <- vcovHC(lin.mod.full,type="HC2")
        } else {
            vcov_mat_full <- NULL
        }
    }

######### Extract Effects from the profile-vars only linear model

# proposed nomenclature here:
# effect = attribute in question, which has "effect levels"
# depends = attributes it depends on, each of which has "depend levels"

    #Make R CMD check happy
    J_baseline <- NULL
    J_effect <- NULL
    
    #before we start, make a blank call to the design array J
    J_call <- Quote(design$J[])
    J_call <- J_call[c(1, 2, rep(3, length(dim(design$J))))]

    #warnings counter
    warn_i <- 0
    
############## loop over unique profile vars only (AMCE and ACIE); interactions etc. below

    #blank list for output
    estimates <- list()
    #re-sort profile effects
    profile_effects <- c(sort(profile_effects[!grepl(":",profile_effects)]), profile_effects[grepl(":",profile_effects)])
    #initialize list for weighted cross-terms
    covariance_list <- list()
    #blank matrix of var probs
    varprob_mat <- matrix(0,nrow(vcov_mat_prof),ncol(vcov_mat_prof))
    colnames(varprob_mat) <- rownames(varprob_mat) <- colnames(vcov_mat_prof)
        
    for(i in 1:length(profile_effects)) {

        #split into sections if it's an interaction
        substrings <- strsplit(profile_effects[i], "[:*]", perl=TRUE)[[1]]

        #administrative loop to find levels
        all_levels <- list()
        all_levels_coefs <- list()
        for(effect in substrings) {
            #get all level names and coefficient names-- sans baseline!!!
            all_levels[[effect]] <- levels(data[[effect]])[-1]
            all_levels_coefs[[effect]] <- sapply(all_levels[[effect]], function(x) {
                paste(effect, x, sep="")
            })
        }
  
        #find all combinations of level names-- add as FIRST column
        levels <- expand.grid(all_levels, stringsAsFactors = FALSE)
        #make level combos in first column
        levels <- cbind(apply(levels,1,function(x) paste(x,collapse=":")),levels)
        colnames(levels) <- c("name",substrings)

        #and all combinations of actual coefficient names
        coefs <- expand.grid(all_levels_coefs, stringsAsFactors = FALSE)
        coefs <- apply(coefs,1,function(x) paste(x,collapse=":"))
       
        # Initialize the results 
        results <- matrix(nrow=2, ncol = nrow(levels))
        if (length(substrings) > 1) {
            rownames(results) <- c("ACIE", "Std. Error")
        } else {
            rownames(results) <- c("AMCE", "Std. Error")
        }
        colnames(results) <- as.character(levels[,1])

         #### find extra times when this effect is mentioned
        all_depends <- unlist(sapply(all_prof,USE.NAMES = F,function(x) {
            y <- strsplit(x,":")[[1]]
            if (all(is.element(substrings,y))) x
        }))
        # remove the actual term
        all_depends <- all_depends[-is.element(all_depends,profile_effects[i])]

 #### loop over every combination of levels of component effects
        for(j in 1:nrow(levels)) {
                
            #figure out which level of inter we're doing
            effect_level <- as.character(levels[j,1])
            effect_level_coef <- coefs[j]

            #get its beta and var-cov matrix
            initial_beta <- coefficients(lin.mod.prof)[effect_level_coef]
            if (effect_level_coef %in% colnames(vcov_mat_prof)) {
                initial_var <- vcov_mat_prof[effect_level_coef, effect_level_coef]
            } else {
                initial_var <- NA
            }

            #if interaction,make sure there is baseline support for this level combination
            if (!is.na(initial_beta) & !is.na(initial_var) & length(substrings) > 1) {
                for (effect1 in substrings) {
                    #get effect base
                    effect_base1 <- levels(data[[effect1]])[1]
                    #subset data to that level
                    base.subset <- data[which(data[[effect1]] == effect_base1),]                       
                    #loop over other profile-varying vars in interaction to subset further
                    for(effect in substrings[!(substrings %in% effect1)]) {
                        base.subset <- base.subset[which(base.subset[[effect]] == as.character(levels[j,effect])),]
                    }                        
                    #if there's no support left, change beta and var to NA
                    if (nrow(base.subset) == 0) {
                        initial_beta <- initial_var <- NA
                        #and give a warning that you had to do it
                        warn_i <- warn_i + 1
                    }
                }
            }
               
            # If initial_beta and initial_variance are not NA (are valid level combination)
            # and there are dependent variables to incorporate
            if (!is.na(initial_beta) & !is.na(initial_var) & length(all_depends) > 0) {
                #get the slice of design array J associated with baseline and inter level
                J_effect_call <- J_base_call <-  J_call
                for(effect in substrings) {
                    #identify its baseline and modify baseline call accordingly
                    base <- levels(data[[effect]])[1]
                    effect_index <- which(names(dimnames(design$J)) == effect)
                    J_base_call[effect_index + 2] <- base                            
                    #identify level of each effect and modify inter call accordingly
                    level <- levels[j,effect]
                    J_effect_call[effect_index + 2] <- level
                }
                eval(call("<-", Quote(J_baseline), J_base_call))
                eval(call("<-", Quote(J_effect), J_effect_call))

                # Initialize some vectors to store interactions and probabilities
                interaction_probabilities <- c()
                interaction_names <- c()
                covariance_names <- c()
                covariance_probs <- c()

                #### loop over dependencies for all components of interaction
                for(k in 1:length(all_depends)) {

                    #attribute effect is dependent on
                    depend <- all_depends[[k]]
                    #figure out what levels of what variables are involved
                    substrings_d <- strsplit(depend,":")[[1]]
                    substrings_d <- substrings_d[!is.element(substrings_d,substrings)]
                    all_depend_coefs <- list()
                    for (sub in substrings_d) {
                        all_depend_coefs[[sub]] <-  sapply(levels(data[[sub]]), function(x) paste(sub, x,sep=""))
                    }
                    all_depend_levels <- expand.grid(all_depend_coefs)
                    substrings_l <- strsplit(effect_level_coef,":")[[1]]
                    for (l in length(substrings_l):1) {
                        all_depend_levels <- cbind(rep(substrings_l[l], nrow(all_depend_levels)), all_depend_levels)
                    }
                    colnames(all_depend_levels)[1:length(substrings_l)] <- substrings
                    all_depend_levels <- all_depend_levels[sort(colnames(all_depend_levels))]
                    all_depend_level_coefs <- apply(all_depend_levels, 1, function(x) paste(x,collapse=":"))
                
                    #baseline support for depend attribute level 
                    if (!(is.null(dim(J_baseline)))){
                        baseline_support <- apply(J_baseline,substrings_d,sum)
                    } else {
                        baseline_support <- J_baseline
                    }
                    baseline_support[baseline_support != 0] <- 1

                    #probs for depend attribute levels WITH baseline support
                    if (!is.null(dim(J_effect))) {
                        joint_prob <- apply(J_effect, substrings_d, sum)*baseline_support
                    } else {
                        joint_prob <- J_effect*baseline_support
                    }
                    #make it a vector
                    joint_prob <- as.vector(joint_prob)
                    names(joint_prob) <- all_depend_level_coefs

                    all_depend_level_coefs <- all_depend_level_coefs[!is.na(lin.mod.prof$coefficients[all_depend_level_coefs])]
                    varprob_mat[effect_level_coef,all_depend_level_coefs] <- as.numeric(joint_prob[all_depend_level_coefs])/as.numeric(sum(joint_prob))
                    
                    ##### if all_depend_level_coefs is 1 or longer (because R doesn't parse for 1:0 correctly)
                    if (length(all_depend_level_coefs)){
                      ##### loop over levels of depends attribute
                      #if present baselines will omit automatically because coefs are NA
                      for (z in 1:length(all_depend_level_coefs)) {
  
                          #coefficient name that goes with this effect level & depend level
                          depend_level_coef <- all_depend_level_coefs[z] 
                          #calculate probabilities for this effect and depend level 
                          var_prob <- joint_prob[depend_level_coef]
                          var_prob <- as.numeric(var_prob)/as.numeric(sum(joint_prob))
                                  
                          #now add interaction beta and variance to initial 
                          if (!is.na(lin.mod.prof$coefficients[depend_level_coef])) {
                              if (!is.na(vcov_mat_prof[depend_level_coef,depend_level_coef]) & !is.na(vcov_mat_prof[effect_level_coef, depend_level_coef])) {
                                  # add weighted beta to initial_beta
                                  initial_beta <- initial_beta + var_prob*lin.mod.prof$coefficients[depend_level_coef]
                                  # add weighted variance + covariance terms too
                                  initial_var <- initial_var + (var_prob^2)*vcov_mat_prof[depend_level_coef, depend_level_coef] +  2*(var_prob)*vcov_mat_prof[effect_level_coef, depend_level_coef]
                                          
                                   # add probabilities and names to compute covariances
                                  interaction_probabilities <- c(interaction_probabilities, var_prob)
                                  interaction_names <- c(interaction_names, depend_level_coef)
                                  #and across different variables
                                  covariance_probs <- c(covariance_probs,var_prob)
                                  covariance_names <- c(covariance_names,depend_level_coef)
  
                              }
                          } #end if that added beta & var, cov
                          
                      } #end loop over levels of dependent attribute
                    }
                } #end for loop over different dependent attributes

                # after going through all levels of the depends and all dependent attributes
                # add remaining covariance terms to the parameter variance 
                if (length(interaction_probabilities) > 1) {
                    #loop over all depend attributes 1 to N-1
                    for (x in 1:(length(interaction_probabilities) - 1)) {
                        #loop over depend attributes one ahead of previous to the end
                        for (y in (x+1):(length(interaction_probabilities))) {
                            initial_var <- initial_var + 2*interaction_probabilities[x]*interaction_probabilities[y]*vcov_mat_prof[interaction_names[x], interaction_names[y]]
                        }
                    }
                } #end if has more than 1 depend levels and/or depend attributes

                #add names of depend levels and their var probs to list
                #if probs not null that is
                if (!is.null(covariance_probs)) {
                    covariance_list[[effect_level_coef]] <- data.frame(covariance_names,covariance_probs)
                }
            
            } #end if has valid beta, var, dependencies

            # Store effect and standard error estimates
            results[1,j] <- initial_beta
            if (!is.na(initial_var)) {
                results[2,j] <- sqrt(initial_var)
            } else {
                results[2,j] <- NA
            }
            
        } #end for loop over all level combinations

         # combine estimates + SEs into single matrix - store in list
        estimates[[profile_effects[i]]] <- results
            
    } #end for loop over profile effects      

### fix var-cov matrix
    
    ## #add in single sum corrections
    vcov_prof <- varprob_mat %*% vcov_mat_prof + t(varprob_mat %*% vcov_mat_prof) + vcov_mat_prof
  
    ## ## too slow!!
    ## vcov_test <- vcov_prof
    ## non_zero_rows <- unlist(sapply(rownames(varprob_mat),function(x) if(sum(varprob_mat[x,]) != 0) x))
    ## function(x,y,vcov) {
    ##     out <- vcov[x,y]
    ##     return(out)
    ## }
    ## cov.ij <- Vectorize(cov.ij,vectorize.args = c("x","y"))
    ## function(e1,e2,varprob,vcov) {        
    ##     non_zero_e1 <- names(varprob[e1,varprob[e1, ] != 0])
    ##     non_zero_e2 <- names(varprob[e2,varprob[e2, ] != 0])
    ##     out <- sum(outer(non_zero_e1,non_zero_e2,function(x,y) varprob[e1,x]*varprob[e2,y]*cov.ij(x,y,vcov_mat_prof)))
    ##     return(out)
    ## }
    ## weights.ij <- Vectorize(weights.ij,vectorize.args = c("e1","e2"))
    ## test <- outer(non_zero_rows,non_zero_rows,function(x,y) weights.ij(x,y,varprob_mat))
    
    #final modifications for var-cov matrix (double sum corrections)
    #these only exist when both variables have depends terms
    #so only modify previously modified variables
    if (length(covariance_list) > 0) {
            #loop over each modified coefficient
        for (x in 1:length(covariance_list)) {
            var1 <- names(covariance_list)[x]
            names1 <- as.character(covariance_list[[var1]][,1])
            probs1 <- covariance_list[[var1]][,2]
                #loop over all other modified coefficients, so i != j
            for (y in 1:length(covariance_list)) {
                var2 <- names(covariance_list)[y]
                names2 <- as.character(covariance_list[[var2]][,1])
                probs2 <- covariance_list[[var2]][,2]
                #loop over each one of interaction names for var1
                for (z in 1:length(names1)) {
                    vcov_prof[var1,var2] <- vcov_prof[var1,var2] + sum(probs1[z]*probs2*vcov_mat_prof[names1[z],names2])  
                }
            } 
        }
    }

    #determine term names for profile effects (no depends) to keep
    profile_effects_plus <- paste(profile_effects,collapse=" + ")
    profile_effects_form <- formula(paste(c(y_var,profile_effects_plus),collapse = " ~ "))   
    profile_effects_terms <- colnames(model.matrix(profile_effects_form,data))   
    profile_effects_terms <- profile_effects_terms[profile_effects_terms %in% colnames(vcov_mat_prof)]
    vcov_prof <- vcov_prof[profile_effects_terms,profile_effects_terms]

######### Extract Effects from the full model (if have respondent interactions)

# proposed nomenclature here:
# effect = attribute in question, which has "effect levels"
# depends = attributes it depends on, each of which has "depend levels"
# inters = attributes in interaction terms each of which has "inter levels"

 #if there are any respondent effects
    if (length(respondent_vars) > 0) {

        #blank list for output
        conditional.estimates <- list()
        #re-sort respondent effects
        resp_effects <- c(sort(resp_effects[!grepl(":",resp_effects)]), resp_effects[grepl(":",resp_effects)])
        #initialize list for weighted cross-terms
        covariance_list <- list()
        #blank matrix for var probs
        varprob_mat <- matrix(0,nrow(vcov_mat_full),ncol(vcov_mat_full))
        colnames(varprob_mat) <- rownames(varprob_mat) <- colnames(vcov_mat_full)

        #loop over respondent-related effects
        for (i in 1:length(resp_effects)) {
            
            #split into component effects, if interaction
            substrings <- strsplit(resp_effects[i], "[:*]", perl=TRUE)[[1]]

            ## start by finding levels
            #administrative loop over components of interaction
            all_levels <- list()
            all_levels_coefs <- list()
            for(effect in substrings) {
                #if it's not a factor, only has the 1 "level" and coefficient name stays
                if (class(data[[effect]]) != "factor") {
                    all_levels[[effect]] <- effect
                    all_levels_coefs[[effect]] <- effect
                } else {
                    #if it is a factor, get all level names and coefficient names-- sans baseline!!!
                    all_levels[[effect]] <- levels(data[[effect]])[-1]
                    all_levels_coefs[[effect]] <- sapply(all_levels[[effect]],
                                                         function(x) {
                                                             paste(effect, x, sep="")
                                                         })
                }
            }
                
            #find all combinations of level names-- add as FIRST column
            levels <- expand.grid(all_levels, stringsAsFactors = FALSE)
            levels <- cbind(apply(levels,1,function(x) paste(x,collapse=":")),levels)
            colnames(levels) <- c("name",substrings)

            #and all combinations of coefficient names
            coefs <- expand.grid(all_levels_coefs, stringsAsFactors = FALSE)
            coefs <- apply(coefs,1,function(x) paste(x,collapse=":"))
       
            # Initialize the results 
            results <- matrix(nrow=2, ncol = nrow(levels))
            rownames(results) <- c("Conditional Estimate", "Std. Error")
            colnames(results) <- levels[,1]

            #### find extra times when this effect is mentioned in full formula
            # only if anything related to profile var is involved
            if (any(substrings %in% profile_vars)) {
                all_depends <- unlist(sapply(all_run_vars,USE.NAMES = F,function(x) {
                    y <- strsplit(x,":")[[1]]
                    if (all(is.element(substrings,y))) x
                }))
                # remove the actual term
                all_depends <- all_depends[!is.element(all_depends,resp_effects[i])]
                # remove any that involve any other respondent varying terms
                resp.other <- respondent_vars[!is.element(respondent_vars,substrings)]
                all_depends <- unlist(sapply(all_depends,function(x) {   
                    sub_depends <- strsplit(x,":")[[1]]
                    if (all(!is.element(sub_depends,resp.other))) x
                }))
            } else {
                #no profile vars, no depends
                all_depends <- c()
            }
                              
            #### loop over every combination of levels of component effects
            for(j in 1:nrow(levels)) {
                
                #figure out which level of inter we're doing
                effect_level <- as.character(levels[j,1])
                effect_level_coef <- coefs[j]

                #get its beta and var-cov matrix
                initial_beta <- coefficients(lin.mod.full)[effect_level_coef]
                if (effect_level_coef %in% colnames(vcov_mat_full)) {
                    initial_var <- vcov_mat_full[effect_level_coef, effect_level_coef]
                } else {
                    initial_var <- NA
                }

                #make sure there is baseline support for this level combination
                if (!is.na(initial_beta) & !is.na(initial_var)) {                   
                    for (effect1 in substrings[substrings %in% profile_vars]) {
                        #get effect base
                        effect_base1 <- levels(data[[effect1]])[1]
                        #subset data to that level
                        base.subset <- data[which(data[[effect1]] == effect_base1),]                       
                        #loop over other profile-varying vars in interaction to subset further
                        for(effect in substrings[substrings %in% profile_vars][!(substrings[substrings %in% profile_vars] %in% effect1)]) {
                            base.subset <- base.subset[which(base.subset[[effect]] == as.character(levels[j,effect])),]
                        }                        
                        #if there's no support left, change beta and var to NA
                        if (nrow(base.subset) == 0) {
                            initial_beta <- initial_var <- NA
                            #and give a warning that you had to do it
                            warn_i <- warn_i + 1
                        }
                    }
                }
               
                # If initial_beta and initial_variance are not NA and there are depends
                # proceed to add to beta and var
                if (!is.na(initial_beta) & !is.na(initial_var) & length(all_depends) > 0) {

                    #get the slice of design array J associated with baseline and inter level
                    #profile variables only!
                    J_effect_call <- J_base_call <-  J_call
                    for(effect in substrings[substrings %in% profile_vars]) {
                        #identify its baseline and modify baseline call accordingly
                        base <- levels(data[[effect]])[1]
                        effect_index <- which(names(dimnames(design$J)) == effect)
                        J_base_call[effect_index + 2] <- base                           
                        #identify level of each effect and modify inter call accordingly
                        level <- levels[j,effect]
                        J_effect_call[effect_index + 2] <- level
                    }
                    eval(call("<-", Quote(J_baseline), J_base_call))
                    eval(call("<-", Quote(J_effect), J_effect_call))

                    # Initialize some vectors to store interactions and probabilities
                    interaction_probabilities <- c()
                    interaction_names <- c()
                    covariance_names <- c()
                    covariance_probs <- c()
  
                    #### loop over dependencies for all components of effect
                    for(k in 1:length(all_depends)) {

                        #attribute effect is dependent on
                        depend <- all_depends[[k]]
                       #figure out what levels of what variables are involved
                        substrings_d <- strsplit(depend,":")[[1]]
                        substrings_d <- substrings_d[!is.element(substrings_d,substrings)]
                        all_depend_coefs <- list()
                        for (sub in substrings_d) {
                            all_depend_coefs[[sub]] <-  sapply(levels(data[[sub]]), function(x) paste(sub, x,sep=""))
                        }
                        all_depend_levels <- expand.grid(all_depend_coefs)
                        substrings_l <- strsplit(effect_level_coef,":")[[1]]
                        for (l in length(substrings_l):1) {
                            all_depend_levels <- cbind(rep(substrings_l[l], nrow(all_depend_levels)), all_depend_levels)
                        }
                        colnames(all_depend_levels)[1:length(substrings_l)] <- substrings
                        ####put terms together in proper order
                        all_depend_levels <- all_depend_levels[sort(colnames(all_depend_levels))]
                        all_depend_level_coefs <- apply(all_depend_levels, 1, function(x) paste(x,collapse=":"))
                            
                       #baseline support for depend attribute level in inter
                        if (!(is.null(dim(J_baseline)))){
                            baseline_support <- apply(J_baseline,substrings_d,sum)
                        } else {
                            baseline_support <- J_baseline
                        }
                        baseline_support[baseline_support != 0] <- 1

                        #support for depend attribute levels WITH baseline support
                        if (!is.null(dim(J_effect))) {
                            joint_prob <- apply(J_effect, substrings_d, sum)*baseline_support
                        } else {
                            joint_prob <- J_effect*baseline_support
                        }
                        #make it a vector
                        joint_prob <- as.vector(joint_prob)
                        names(joint_prob) <- all_depend_level_coefs

                        all_depend_level_coefs <- all_depend_level_coefs[!is.na(lin.mod.full$coefficients[all_depend_level_coefs])]
                        varprob_mat[effect_level_coef,all_depend_level_coefs] <- as.numeric(joint_prob[all_depend_level_coefs])/as.numeric(sum(joint_prob))
                        
                        ## If there are non-null # of depend-level-coefs (because R doesn't handle for 1:0 correctly)
                        if (length(all_depend_level_coefs) != 0){
                          ##### loop over levels of depends attribute
                          #baselines will omit automatically because coefs are NA
                          for (z in 1:length(all_depend_level_coefs)) {
  
                              #coefficient name that goes with this effect level & depend level
                              depend_level_coef <- all_depend_level_coefs[z] 
                              #calculate probabilities for this effect and depend level 
                              var_prob <- joint_prob[depend_level_coef]
                              var_prob <- as.numeric(var_prob)/as.numeric(sum(joint_prob))
                                  
                              #now add interaction beta and variance to initial 
                              if (!is.na(lin.mod.full$coefficients[depend_level_coef])) {
                                  if (!is.na(vcov_mat_full[depend_level_coef,depend_level_coef]) & !is.na(vcov_mat_full[effect_level_coef, depend_level_coef])) {
                                      # add probabilities to initial_beta
                                      initial_beta <- initial_beta +
                                          var_prob*lin.mod.full$coefficients[depend_level_coef]
                                      # add variance + covariance terms too
                                      initial_var <- initial_var + (var_prob^2)*vcov_mat_full[depend_level_coef,depend_level_coef] +  2*(var_prob)*vcov_mat_full[effect_level_coef, depend_level_coef]
                                          
                                      # add probabilities and names to compute covariances
                                      interaction_probabilities <- c(interaction_probabilities, var_prob)
                                      interaction_names <- c(interaction_names, depend_level_coef)
                                      #and across different variables
                                      covariance_probs <- c(covariance_probs,var_prob)
                                      covariance_names <- c(covariance_names,depend_level_coef)
                                          
                                  }
                              } #end if that added beta & var, cov       
                          } #end loop over levels of dependent attribute
                      } #end for loop over different dependent attributes
                    }
                    # after going through all levels of the depends and all dependent attributes
                    # add remaining covariance terms to the parameter variance 
                    if (length(interaction_probabilities) > 1) {
                            #loop over all depend attributes 1 to N-1
                        for (x in 1:(length(interaction_probabilities) - 1)) {
                                #loop over depend attributes one ahead of previous to the end
                            for (y in (x+1):(length(interaction_probabilities))) {
                                initial_var <- initial_var + 2*interaction_probabilities[x]*interaction_probabilities[y]*vcov_mat_full[interaction_names[x], interaction_names[y]]
                            }
                        }
                    } #end if has more than 1 depend levels and/or depend attributes

                    #add names of depend levels and their var probs to list
                    #if there are any that is
                    if (!is.null(covariance_probs)) {
                        covariance_list[[effect_level_coef]] <- data.frame(covariance_names,covariance_probs)
                    }
                } #end if initial beta and var are NA, has depends
                
                # Store effect and standard error estimates
                results[1,j] <- initial_beta
                if (!is.na(initial_var)) {
                    results[2,j] <- sqrt(initial_var)
                } else {
                    results[2,j] <- NA
                }
                
            } #end for loop over all level combinations

            # combine estimates + SEs into single matrix - store in list
            conditional.estimates[[resp_effects[i]]] <- results
            
        } #end for loop over respondent related effects      

        ## #add in single sum corrections
        vcov_resp <- varprob_mat %*% vcov_mat_full + t(varprob_mat %*% vcov_mat_full) + vcov_mat_full

       #final modifications for var-cov matrix
       #these only exist when both variables have depends terms
       #so only modify previously modified variables
        if (length(covariance_list) > 0) {
            #loop over each modified coefficient
            for (x in 1:length(covariance_list)) {
                var1 <- names(covariance_list)[x]
                names1 <- as.character(covariance_list[[var1]][,1])
                probs1 <- covariance_list[[var1]][,2]
                #loop over all other modified coefficients, so i != j
                for (y in 1:length(covariance_list)) {
                    var2 <- names(covariance_list)[y]
                    names2 <- as.character(covariance_list[[var2]][,1])
                    probs2 <- covariance_list[[var2]][,2]
                    #loop over each one of interaction names for var1
                    for (z in 1:length(names1)) {
                        vcov_resp[var1,var2] <- vcov_resp[var1,var2] + sum(probs1[z]*probs2*vcov_mat_full[names1[z],names2])  
                    }
                } 
            }
        }

        ###terms to keep (no depends)
        resp_effects_plus <- paste(resp_effects,collapse=" + ")
        resp_effects_form <- formula(paste(c(y_var,resp_effects_plus),collapse = " ~ "))
        resp_effects_terms <- colnames(model.matrix(resp_effects_form,data))
        resp_effects_terms <- resp_effects_terms[resp_effects_terms %in% colnames(vcov_mat_full)]
        vcov_resp <- vcov_resp[resp_effects_terms,resp_effects_terms]
        
    } #end if there are any respondent related effects
      
############  create conjoint object for output
    
    output <- list()
    class(output) <- c("amce")
    
    #saving things for unconditional estimates
    output$estimates <- estimates
    #saving profile attributes
    output$attributes <- dimnames(design$J)
    #save original profile-only vcov matrix
    output$vcov.prof <- vcov_prof
    #save sample size used for unconditional estimates
    output$samplesize_prof <- sample_size_prof
    #save style edited formula (no depends)
    output$formula <- form

    #final warning tally
    if (warn_i > 0) {
        warning(paste("Warning: ",warn_i," interaction levels lacked support at baseline, effects undefined unless alternative baseline is provided."))
    }
    
    #saving things for conditional estimates
    if (length(respondent.varying) > 0) {     
        output$cond.estimates <- conditional.estimates
        output$vcov.resp <- vcov_resp
        output$samplesize_full <- sample_size_full
        #save style edited formula (no depends), only resp-related
        output$cond.formula <- resp_effects_form
    }
    
    # Save baselines of unique (main) effects (if factor) to "baselines"
    # If continuous save summary information to "continuous"
    output$baselines <- list()
    output$continuous <- list()
    for (k in unique_vars) {
        if (class(data[[k]]) == "factor") {
            output$baselines[[k]] <- levels(data[[k]])[1]
        } else {
             output$continuous[[k]] <- quantile(model.matrix(form,data)[,k], probs=c(0.25,0.5,0.75), na.rm=T)
            #output$continuous[[k]] <- quantile(data[[k]],  probs=c(0.25,0.5,0.75),na.rm=T)
        }
    }

    #save number of respondents if ID given
    if (!is.null(respondent.id)) {
        output$numrespondents <- length(unique(data[[respondent.id]]))
    } else {
        output$numrespondents <- NULL
    }

    #save respondent variables if given
    if (!is.null(respondent.varying)) {
        output$respondent.varying <- respondent_vars
    } else {
        output$respondent.varying <- NULL
    }

    #save weights as output (if any)
    if (!is.null(weights)) {
        output$weights <- subset(data, select = weights)
    } else {
        output$weights <- NULL
    }

    #save original names
    output$user.names <- user_names
    output$user.levels <- user_levels
    #save the original data
    output$data <- data
    return(output)
}

############################################################
## summary function for results of main AMCE function
############################################################

# Function for summarizing output from main amce function
# LIST given to "covariate.values" contains VECTORS at which ...
# ... conditional effects will be calculated; default is quantiles...
# ... in the case of continuous, levels otherwise
# ... can be given manual names by naming entries
# Note that must be NAMED LIST, entry name is the respondent.varying effect

summary.amce <- function(object, covariate.values=NULL,  ...) {
    amce_obj <- object

######################### administrative section
    
    # Initialize list to store summary object
    summary_results <- list()
    # Create header of data.frame
    header <- c("Attribute", "Level", "Estimate", "Std. Err", "z value", "Pr(>|z|)", " ")

    #more than 1 resp variable warning
    if (length(amce_obj$respondent.varying) > 1) {
        warning("Warning: More than 1 respondent-varying characteristic detected. Conditional estimates will be presented varying only one characteristic at time while the other is held at 0 (if continuous) or its base (if factor).")
    }

    #checks on user-supplied covariate values
    if(!is.null(covariate.values)) {
        #clean variable names
        names(covariate.values) <- clean.names(names(covariate.values))
        for (i in 1:length(covariate.values)) {
            #make sure they appear AMCE object
            if (!names(covariate.values)[i] %in% amce_obj$respondent.varying) {
                stop(paste(c("Error: variable",names(covariate.values)[i],"is not a respondent-varying characteristic in AMCE object."),collapse=" "))
            }
            #if they do appear and are factors, clean them and make sure valid
            if (names(covariate.values)[i] %in% names(amce_obj$baselines)) {
                covariate.values[[i]] <- clean.names(covariate.values[[i]])
                if (any(!covariate.values[[i]] %in% levels(amce_obj$data[[names(covariate.values)[i]]]))) {
                    stop(paste("Error: some level(s) of variable",names(covariate.values)[i],"do not appear in data."))
                }
            }
        }
    }
 
##################### reporting unconditional results

    #all attribute estimates
    all_prof <- names(amce_obj$estimates) 
    #get AMCE only
    all_amce <- grep(":",all_prof,value=T,invert=T)
    #get ACIE only
    all_acie <- grep(":",all_prof,value=T,invert=F)

     #How many AMCE?
     namce <-  sum(sapply(all_amce,function(x) length(unlist(amce_obj$estimates[[x]]))/2))
    # Create results matrix for AMCE only
    summary_results[["amce"]] <- matrix(nrow= namce, ncol=length(header))
    colnames(summary_results[["amce"]]) <- header
    summary_results[["amce"]] <- as.data.frame(summary_results[["amce"]])
    #amce index
    amce_i <- 1
    
    #summary results matrix for ACIE only
    if (length(all_acie) > 0) {
        #number affected
        nacie <- sum(sapply(all_acie,function(x) length(unlist(amce_obj$estimates[[x]]))/2))
        #make results matrix
        summary_results[["acie"]] <- matrix(nrow= nacie, ncol=length(header))
        colnames(summary_results[["acie"]]) <- header
        summary_results[["acie"]] <- as.data.frame(summary_results[["acie"]])
        #acie index
        acie_i <- 1
    }

    #objects for storing baselines
    baselines_amce <- c()
    baselines_acie <- c()
    #objects for storing effect names
    names_amce <- c()
    names_acie <- c()

    # Loop over non-respondent varying attributes, which are all factors by assumption
    for (effect in all_prof) {
        #split terms
        variates <- strsplit(effect, ":")[[1]]
        # Figure out the baseline level(s)
        lev_list <- c()
        print_names <- c()
        for (var in variates) {
            lev_list <- c(lev_list, amce_obj$user.levels[[amce_obj$baselines[[var]]]])
            print_names <- c(print_names,amce_obj$user.names[[var]])
        }
        print_level <- paste(lev_list,sep="",collapse=":")
        print_effect <- paste(print_names,sep="",collapse=":")       
        # If ACIE
        if (grepl(":",effect)) {
            #set entry name
            entry_name <- "acie"
            #report baselines and names
            baselines_acie <- c(baselines_acie,print_level)           
            names_acie <- c(names_acie,print_effect)
            #which index?
            index <- acie_i
        } else {
            entry_name <- "amce"
            # report baselines and names
            baselines_amce <- c(baselines_amce,print_level)
            names_amce <- c(names_amce, print_effect)
            #which index?
            index <- amce_i
        }   
         # Append results to the estimates dataframe
        for (p in 1:ncol(amce_obj$estimates[[effect]])) {
            variate_levels <- strsplit(colnames(amce_obj$estimates[[effect]])[p],":")[[1]]
            lev_list <- c()
            for (lev in variate_levels) {
                lev_list <- c(lev_list, amce_obj$user.levels[[lev]])
            }
            print_level <- paste(lev_list,sep="",collapse=":")
            summary_results[[entry_name]][index,1] <- print_effect
            summary_results[[entry_name]][index,2] <- print_level
            summary_results[[entry_name]][index,3] <- amce_obj$estimates[[effect]][1,p]
            summary_results[[entry_name]][index,4] <- amce_obj$estimates[[effect]][2,p]
            zscr <- amce_obj$estimates[[effect]][1,p]/amce_obj$estimates[[effect]][2,p]
            summary_results[[entry_name]][index,5] <- zscr
            pval <- 2*pnorm(-abs(zscr))
            summary_results[[entry_name]][index,6] <- pval            
                # Stars!
            if (!is.na(pval)) {
                if (pval < .001) {
                    summary_results[[entry_name]][index,7] <- "***"
                } else if (pval < .01) {
                    summary_results[[entry_name]][index,7] <- "**"
                } else if (pval < .05) {
                    summary_results[[entry_name]][index,7] <- "*"
                } else {
                    summary_results[[entry_name]][index,7] <- ""
                }
            } else {
                summary_results[[entry_name]][index,7] <- ""
            }
            index <- index + 1
        }
        #advance appropriate index
        if (grepl(":",effect)) {
            acie_i <- index
        } else {
            amce_i <- index
        }
    }

    #save as data frame and save baselines
    summary_results[["amce"]] <- as.data.frame(summary_results[["amce"]])
    summary_results[["baselines_amce"]] <- data.frame("Attribute" = names_amce, "Level" = baselines_amce)

    #if any acie save those too
    if (grepl(":",effect)) {
        summary_results[["acie"]] <- as.data.frame(summary_results[["acie"]])
        summary_results[["baselines_acie"]] <- data.frame("Attribute" = names_acie, "Level" = baselines_acie)
    }

################ reporting conditional results

## all appearances by the modified variable (say "var1")
## have already had the beta and var, cov for any other appearances added in
## (such as "dependency:var1:respondent.var" for "var1:respondent.var")
## So the only things that need to be grabbed are the profile term
## and the interaction with the respondent varying term (var1 and var1:respondent.var)
    
    # If there are respondent varying attributes, add estimates at levels/quantiles
    if (length(amce_obj$cond.estimates) > 0) {

        # Extract vectors of betas for all coefficients
        beta_vector <- do.call(cbind,amce_obj$cond.estimates)[1,]
        #get correct coefficient names
        all_coef <- c()
        for (var in names(amce_obj$cond.estimates)) {
            prof_vars <- strsplit(var, ":")[[1]]
            for (lev in colnames(amce_obj$cond.estimates[[var]])) {
                prof_levs <- strsplit(lev, ":")[[1]]
                orig_profs <- c()
                for (v in 1:length(prof_vars)) {
                    #continuous component keeps simple name
                    if (prof_vars[v] == prof_levs[v]) {
                        orig_profs[v] <- prof_vars[v]
                    } else  {
                        #factor component has 2 parts, var and level
                        orig_profs[v] <- paste(c(prof_vars[v], prof_levs[v]), collapse="")
                    }
                }
                orig_prof <- paste(orig_profs,collapse=":")
                all_coef <- c(all_coef,orig_prof)
            }
        }
        names(beta_vector) <- all_coef
        #also get true levels for all factors in the cond.formula
        xlevs <- sapply(all.vars(amce_obj$cond.formula)[-1][all.vars(amce_obj$cond.formula)[-1] %in% names(amce_obj$baselines)],function(x) levels(amce_obj$data[[x]]), simplify = F)

        #empty vectors for table key
        tab_name_amce <- c()
        tab_var_amce  <- c()
        tab_val_amce  <- c()
        tab_name_acie <- c()
        tab_var_acie <- c()
        tab_val_acie <- c()
        
        # Loop over respondent-varying characteristics
        for (effect in amce_obj$respondent.varying) {

            #how to print the effect name
            print_effect <- amce_obj$user.names[[effect]]
            # identify all REQUESTED terms involving effect
            all_req_vars <- attr(terms(amce_obj$formula),"term.labels")
            all_resp <- unlist(sapply(all_req_vars,function(x) {
                y <- strsplit(x,":")[[1]]
                if (any(y == effect)) x                
            }))
            #figure out profile attributes these refer to
            all_mod <- unlist(sapply(all_resp,function(x) {
                subs <- strsplit(x,":")[[1]]
                subs <- subs[which(is.element(subs,all_prof))]
                if (length(subs) > 0) paste(subs,collapse=":")
            }))
            #just unique ones
            all_mod <- unique(all_mod)
            #make sure there are some
            if (length(all_mod) == 0) {
                stop(paste(c("respondent characteristic",effect,"not interacted with profile attributes, no interpretation"),collapse=" "))
            }
            
            #### split modified terms into AMCE or ACIE
            mod_amce <- grep(":",all_mod,value = T,invert = T)
            mod_acie <- grep(":",all_mod,value = T,invert = F)
            #how many AMCE are affected by this respondent characteristic?
            namce <- sum(sapply(mod_amce,function(x) length(unlist(amce_obj$estimates[[x]]))/2))
            #how many ACIE are affected by this respondent characteristic? 
            if (length(mod_acie) > 0) {                    
                nacie <- sum(sapply(mod_acie,function(x) length(unlist(amce_obj$estimates[[x]]))/2))
            }

            #get raw levels default for this effect if no values given or install given names
            raw_levels <- c()
            if (is.null(covariate.values[[effect]])) {
                if (effect %in% names(amce_obj$baselines)) {
                    # if it's a factor get levels, adding baseline
                    raw_levels <- colnames(amce_obj$cond.estimates[[effect]])
                    raw_levels <- c(amce_obj$baselines[[effect]],raw_levels)
                } else {
                    # otherwise get summary information from "continuous"
                    raw_levels <- amce_obj$continuous[[effect]]
                }
            } else {
                raw_levels <- covariate.values[[effect]]
            }

            ######## Now loop over levels/quantiles of "effect"
            for (i in 1:length(raw_levels)) {

                # get the appropriate model matrix for prediction
                # first, make a data_dummy data matrix
                data_dummy <- amce_obj$data
                # if continuous, setting value of "effect" using fake data matrix
                # also name of level; if continuous just effect name
                if (effect %in% names(amce_obj$continuous)) {
                    #name of "level" is just effect name
                    resp_lev <- effect
                    #as is original variable name
                    orig_resp <- effect
                    #how to print level name
                    if (!is.null(names(raw_levels))) {
                        print_resp_lev <- names(raw_levels)[i]
                    } else {
                        print_resp_lev <- raw_levels[i]
                    }
                } else {
                    #if resp var is a factor set level
                    resp_lev <- raw_levels[i]
                    # make original var name
                    orig_resp <- paste(c(effect,resp_lev),collapse="")
                    #printed name from user input
                    print_resp_lev <- amce_obj$user.levels[[resp_lev]]
                }
                #set value of in data
                data_dummy[[effect]] <- raw_levels[i]
                
                #### prep AMCE results
                # make results matrices
                entry_name_amce <- paste(c(effect,i,"amce"),collapse="")
                # make empty results matrix
                summary_results[[entry_name_amce]] <- matrix(NA,nrow = namce, ncol=length(header))
                colnames(summary_results[[entry_name_amce]]) <- header
                summary_results[[entry_name_amce]] <- as.data.frame(summary_results[[entry_name_amce]])
                amce_i <- 1
                #edit table key
                tab_name_amce <- c(tab_name_amce,entry_name_amce)
                tab_var_amce <- c(tab_var_amce, effect)
                tab_val_amce <- c(tab_val_amce, print_resp_lev)

                #### prep ACIE results, if any
                if (length(mod_acie) > 0) {
                    #entry name
                    entry_name_acie <- paste(c(effect,i,"acie"),collapse="")
                    #results matrix
                    summary_results[[entry_name_acie]] <- matrix(NA,nrow= nacie, ncol=length(header))
                    colnames(summary_results[[entry_name_acie]]) <- header
                    summary_results[[entry_name_acie]] <- as.data.frame(summary_results[[entry_name_acie]])
                    acie_i <- 1
                    #edit table key
                    tab_name_acie <- c(tab_name_acie,entry_name_acie)
                    tab_var_acie <- c(tab_var_acie, effect)
                    tab_val_acie <- c(tab_val_acie, print_resp_lev)
                }

                # loop over all attribute effects modified by "effect"
                for (prof_var in all_mod) {

                    #get name to print
                    prof_vars <- strsplit(prof_var, ":")[[1]]
                    print_vars <- c()
                    for (var in prof_vars) {
                        print_vars <- c(print_vars,amce_obj$user.names[[var]])  
                    }
                    print_prof_var <- paste(c(print_vars),collapse=":")
                    # Set entry to add to 
                    if (grepl(":",prof_var)) {
                        #set entry name
                        entry_name <- entry_name_acie
                        #which index?
                        index <- acie_i
                    } else {
                        entry_name <- entry_name_amce
                        #which index?
                        index <- amce_i
                    }

                    #loop over the associated betas 
                    for (p in 1:ncol(amce_obj$cond.estimates[[prof_var]])) {

                        #name of level we're modifying
                        prof_lev <- colnames(amce_obj$cond.estimates[[prof_var]])[p]
                        prof_levs <- strsplit(prof_lev, ":")[[1]]
                        #how to print it
                        print_levs <- c()
                        for (lev in prof_levs) {
                            print_levs <- c(print_levs, amce_obj$user.levels[[lev]])
                        }
                        print_prof_level <- paste(print_levs,sep="",collapse=":")
                        #put together var and level for original name
                        orig_profs <- c()
                        for (v in 1:length(prof_vars)) {
                            orig_profs[v] <- paste(c(prof_vars[v], prof_levs[v]), collapse="")
                            data_dummy[[prof_vars[v]]] <- prof_levs[v]
                        }
                        orig_prof <- paste(orig_profs,collapse=":")
                        
                        #use modified data_dummy to make model matrix, preserving old levels
                        pred_mat <- model.matrix(amce_obj$cond.formula,data_dummy,xlev = xlevs)
                        #all column names containing original profile var name
                        pred_cols <- grep(paste(c("(?=.*", orig_resp, ")(?=.*", orig_prof, ")"), collapse=""),colnames(pred_mat) ,value=T,perl=T)                            
                        #make sure it contains no unrelated vars
                        pred_cols <- unlist(sapply(pred_cols,function(x) {
                            subs <- strsplit(x,":")[[1]]
                            subs <- subs[!is.element(subs,orig_profs)]
                            subs <- subs[!grepl(orig_resp,subs)]
                            if(length(subs) < 1) x
                        }))
                        #add in base modified attribute
                        pred_cols <- c(orig_prof,pred_cols)
                        #but no duplicates
                        pred_cols <- unique(pred_cols)

                        #calculation for coefficient and SE
                        beta_inter <- pred_mat[1,pred_cols] %*% beta_vector[pred_cols]
                        if (!is.na(beta_inter)) {
                            #gather covariance terms; other occurences dealt with in main fxn
                            all_cov <- c()
                            for (a in 1:length(pred_cols)) {
                                for (b in 1:length(pred_cols)) {
                                    all_cov <- c(all_cov,pred_mat[1,pred_cols[a]]*pred_mat[1,pred_cols[b]]*amce_obj$vcov.resp[pred_cols[a], pred_cols[b]])
                                }
                            }
                            var_inter <- sum(all_cov)
                            #modified zscr and pval
                            se_inter <- sqrt(var_inter)
                            zscr <- beta_inter/se_inter
                            pval <- 2*pnorm(-abs(zscr))
                        } else {
                            se_inter <- zscr <- pval <- NA
                        }

                        #write results
                        summary_results[[entry_name]][index,1] <- print_prof_var
                        summary_results[[entry_name]][index,2] <- print_prof_level
                        summary_results[[entry_name]][index,3] <- beta_inter
                        summary_results[[entry_name]][index,4] <- se_inter
                        summary_results[[entry_name]][index,5] <- zscr
                        summary_results[[entry_name]][index,6] <- pval
                        # Stars!
                        if (!is.na(beta_inter)) {
                            if (pval < .001) {
                                summary_results[[entry_name]][index,7] <- "***"
                            } else if (pval < .01) {
                                summary_results[[entry_name]][index,7] <- "**"
                            } else if (pval < .05) {
                                summary_results[[entry_name]][index,7] <- "*"
                            } else {
                                summary_results[[entry_name]][index,7] <- ""
                            }
                        } else {
                            summary_results[[entry_name]][index,7] <- ""
                        }
                            index <- index + 1

                    } #end loop over levels of profile var
                    #advance appropriate index
                    if (grepl(":",effect)) {
                        acie_i <- index
                    } else {
                        amce_i <- index
                    }
                } #end loop over modified profile vars

                #save AMCE results as data frame
                summary_results[[entry_name_amce]] <- as.data.frame(summary_results[[entry_name_amce]])
                #and ACIE, if any
                if (length(mod_acie) > 0) {
                    summary_results[[entry_name_acie]] <- as.data.frame(summary_results[[entry_name_acie]])
                }
                    
            } #end loop over respondent var levels (lev_list)
                    
        } #end loop over respondent varying characteristics

        #save results as data frame and save baselines
         summary_results[["table_values_amce"]] <- data.frame("Table Name" = tab_name_amce, "Level Name" = tab_var_amce, "Level Value" = tab_val_amce)
        summary_results[["table_values_amce"]] <- apply(summary_results[["table_values_amce"]],c(1,2),function(x) as.character(x))
      
         summary_results[["table_values_acie"]] <- data.frame("Table Name" = tab_name_acie, "Level Name" = tab_var_acie, "Level Value" = tab_val_acie)
        summary_results[["table_values_acie"]] <- apply(summary_results[["table_values_acie"]],c(1,2),function(x) as.character(x))
 
        } else {
        summary_results[["table_values_amce"]] <- NULL
        summary_results[["table_values_acie"]] <- NULL
    }
  
    # Save sample size(s)
    summary_results[["samplesize_estimates"]] <- amce_obj$samplesize_prof
    if (!is.null(amce_obj$samplesize_full)) {
        summary_results[["samplesize_resp"]] <- amce_obj$samplesize_full
    } else {
        summary_results[["samplesize_resp"]] <- NULL
    }

    # If there's a respondent number, add that as well
    if (!is.null(amce_obj$numrespondents)) {
        summary_results[["respondents"]] <- amce_obj$numrespondents
    } else {
        summary_results[["respondents"]] <- NULL
    }

    # Set class
    class(summary_results) <- c("summary.amce")
  
    # Return
    return(summary_results)
}

############################################################
## print summary function for results of main AMCE function
############################################################

print.summary.amce <- function(x, digits=5, ...) {
    summary_result <- x

    #basic print for AMCE
    cat("------------------------------------------\n")
    cat("Average Marginal Component Effects (AMCE):\n")
    cat("------------------------------------------\n")
    print(summary_result$amce, digits=digits, row.names=F)
    cat("---\n")
    cat(paste("Number of Obs. = ", summary_result$samplesize_estimates, sep=""))
    cat("\n")
    cat("---\n")   
    if (!is.null(summary_result$respondents)) {
        cat(paste("Number of Respondents = ", summary_result$respondents, sep=""))
        cat("\n")
        cat("---\n")
    }
    
    #add extra tables for AMCE interactions with respondent varying
    if (!is.null(summary_result$table_values_amce)) {
        if (nrow(summary_result$table_values_amce) > 0) {
            for (i in 1:nrow(summary_result$table_values_amce)) {
                cat("------------------------------------------------------------\n")
                cat(paste(c("Conditional AMCE's","(",summary_result$table_values_amce[i,2],"=", summary_result$table_values_amce[i,3],"):\n"),collapse=" "))
                cat("------------------------------------------------------------\n")
                print(summary_result[[summary_result$table_values_amce[i,1]]], digits=digits, row.names=F)
                cat("---\n")
                cat(paste("Number of Obs. = ", summary_result$samplesize_resp, sep=""))
                cat("\n")            
                if (!is.null(summary_result$respondents)) {
                    cat(paste("Number of Respondents = ", summary_result$respondents, sep=""))
                    cat("\n")
                }
                cat("---\n")  
                cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05")
                cat("\n")
                cat("\n") 
            }
        }
    }

    #print AMCE baselines
    cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05")
    cat("\n")
    cat("\n")
    cat("--------------------\n")
    cat("AMCE Baseline Levels:\n")
    cat("--------------------\n")
    print(summary_result$baselines_amce, row.names=F)
    cat("\n")
    cat("\n")


    #Tables for UNCONDITIONAL interaction
    if (!is.null(summary_result$acie)) {
        cat("---------------------------------------------\n")
        cat("Average Component Interaction Effects (ACIE):\n")
        cat("---------------------------------------------\n")
        print(summary_result$acie, digits=digits, row.names=F)
        cat("---\n")
        cat(paste("Number of Obs. = ", summary_result$samplesize_estimates, sep=""))
        cat("\n")
        if (!is.null(summary_result$respondents)) {
            cat(paste("Number of Respondents = ", summary_result$respondents, sep=""))
            cat("\n")
            
        }
    }
    
    #add extra tables for ACIE interactions with respondent varying
    if (!is.null(summary_result$table_values_acie)) {
        if (nrow(summary_result$table_values_acie) > 0) {
            for (i in 1:nrow(summary_result$table_values_acie)) {
                cat("------------------------------------------------------------\n")
                cat(paste(c("Conditional ACIE's","(",summary_result$table_values_acie[i,2],"=", summary_result$table_values_acie[i,3],"):\n"),collapse=" "))
                cat("------------------------------------------------------------\n")
                print(summary_result[[summary_result$table_values_acie[i,1]]], digits=digits, row.names=F)
                cat("---\n")
                cat(paste("Number of Obs. = ", summary_result$samplesize_resp, sep=""))
                cat("\n")              
                if (!is.null(summary_result$respondents)) {
                    cat(paste("Number of Respondents = ", summary_result$respondents, sep=""))
                    cat("\n") 
                }
                cat("---\n")
                cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05")
                cat("\n")
                cat("\n")
            }             
        }
    }

    #baselines for ACIE
    if (!is.null(summary_result$acie)) {
        cat("---\n")
        cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05")
        cat("\n")
        cat("\n")
        cat("--------------------\n")
        cat("ACIE Baseline Levels:\n")
        cat("--------------------\n")
        print(summary_result$baselines_acie, row.names=F)
        cat("\n")
        cat("\n")
    }

}

########################################################
# plot amce function
#######################################################

# default will return single plot with point estimates and CI's
# to facet plots, give "facet.name" the name of the variable to facet by...
# ... by default variable facets will be ALL level combinations...
# ... or if a continuous variable is given, will use quantiles...
# ... to customize, directly give "facet.levels" a named LIST ...
# ... with desired values of each variable entered as DATA FRAME with desired names
# "display" takes one of: "all", "unconditional", "conditional"
# with no given facet and no respondent vars, this choice is irrelevant (unconditional only)
# with no given facet and respondent vars, default is all but can choose just unconditional or interaction
# with a facet given (respondent or otherwise), similarly can choose

plot.amce <- function(x, main="", xlab="Change in E[Y]", ci=.95, colors=NULL, xlim=NULL, breaks=NULL, labels=NULL, attribute.names = NULL, level.names = NULL, label.baseline = TRUE, text.size=11, text.color = "black", point.size = .5, dodge.size=0.9, plot.theme = NULL, plot.display = "all", facet.names = NULL, facet.levels = NULL, ...) {
    
    # You need ggplot2
    amce_obj <- x
    ylim <- xlim
    
    # Make R CMD check happy
    pe <- NULL
    se <- NULL
    group <- NULL
    lower <- NULL
    upper <- NULL
    
############################## basic set-up: get attributes and levels

    # Extract raw attribute names from the amce_obj$estimates object
    raw_attributes <- names(amce_obj$estimates)
    # Extract raw levels 
    raw_levels <- list()
    for(m in names(amce_obj$estimates)) {
        raw_levels[[m]] <- colnames(amce_obj$estimates[[m]])
    }
    
    # Determine baseline level for each effect estimate in raw_levels and append to beginning of each vector in raw_levels
    for (effect in names(raw_levels)) {
        effect_elements <- strsplit(effect, ":")[[1]]
        baseline_interactions <- c()
        for (elem in effect_elements) {
            #if effect element is a factor variable get baseline
            if (elem %in% names(amce_obj$baselines)) {
                baseline_interactions <- c(baseline_interactions, amce_obj$baselines[[elem]])
            } else {
                #otherwise just add name of continuous variable
                baseline_interactions <- c(baseline_interactions,elem)
            }
        }
        interaction_str <- paste(baseline_interactions,sep="",collapse=":")
        raw_levels[[effect]] <- c(interaction_str, raw_levels[[effect]]) 
    }

############################## Incorporate and adjust user-input
 
    # Convert ci to z-score
    if (ci < 1 & ci > 0) {
        zscr <- qnorm(1- ((1-ci)/2))
    } else {
        cat("Invalid confidence interval -- Defaulting to 95%")
        zscr <- qnorm(1- ((1-.95)/2))
    }
    
    # Sanity check user-provided attribute.names against AMCE objects
    if (!is.null(attribute.names)) {
        attribute.names <- unique(attribute.names)
        if (length(attribute.names) != length(raw_attributes)) {
            cat(paste("Error: The number of unique elements in attribute.names ", length(attribute.names), " does not match the attributes in amce object for which estimates were obtained: ", paste(raw_attributes,collapse=", "), "\n", sep=""))
            cat("Defaulting attribute.names to attribute names in AMCE object\n")
            attribute.names <- NULL
        }
    }

    # Sanity check user-provided level.names against AMCE object
    if (!is.null(level.names)) {
        names(level.names) <- clean.names(names(level.names))
        for (name in names(level.names)) {
            if (name %in% names(raw_levels)) {
                if (length(level.names[[name]]) != length(raw_levels[[name]])) {
                    cat(paste("Error: level.names lengths do not match levels for attribute ", name, "\n",sep=""))
                    cat(paste("Defaulting level.names for attribute ", name, " to level names in AMCE object", "\n",sep=""))
                    level.names[[name]] <- NULL
                }   
            } else {
                cat(paste("Error: level.names entry ",name," not in AMCE object. Removing level.names for attribute.","\n",sep=""))
                level.names[[name]] <- NULL
            }
        }
    }

    # valid plot.display option?
    plot.display.opts <- c("all","unconditional","interaction")
    if (!is.element(plot.display,plot.display.opts)) {
        stop(paste(c("Error-- plot.display must be once of: ",paste(plot.display.opts,collapse=", ")),collapse=" "))
    }

    #clean facet names; if levels but no names? level names are facets
    if (!is.null(facet.names)) {
        facet.names <- clean.names(facet.names)
    } else if (!is.null(facet.levels)) {
        facet.names <- clean.names(names(facet.levels))
    }

    #check that they are in AMCE object
    facet.names.check <- c()
    for (facet.name in facet.names) {
        if (!facet.name %in% names(amce_obj$estimates) & !facet.name %in% names(amce_obj$cond.estimates)) {
            stop(paste(c("Error-- cannot find facet name",facet.name,"in AMCE object output."),collapse=" "))
        } else {
            facet.names.check <- c(facet.names.check,facet.name)
        }
    }
    facet.names <- facet.names.check

    #if no facets but there are respondent varying characteristics, use those
    if ((is.null(facet.names)) & (length(amce_obj$respondent.varying) > 0) & (plot.display != "unconditional")) {
        facet.names <- amce_obj$respondent.varying
    }

    #no facet name or resp var, must be unconditional
    if (is.null(facet.names) & plot.display == "interaction") {
        warning("Warning: no facet name or respondent varying characteristic provided to calculate conditional estimates. Will display unconditional only")
        plot.display <- "unconditional"
    }

    #unconditional but facet names given? remove facet names
    if(plot.display == "unconditional" & !is.null(facet.names)) {
        warning("Warning-- plot display is set to unconditional, facet names will be ignored")
        facet.names <- NULL
        facet.levels <- NULL
    }

    #check and clean facet levels if provided
    if (!is.null(facet.levels)) {
        #clean names of facet levels
        names(facet.levels) <- clean.names(names(facet.levels))
        #clean actual levels
        for (facet.name in names(facet.levels)) {
            #if it's a factor, clean up level names
            if (facet.name %in% names(amce_obj$baselines)) {
                facet.levels[[facet.name]] <- clean.names(facet.levels[[facet.name]])
                #make sure that if it's profile-varying, there's more than base
                if (facet.name %in% names(amce_obj$estimates) && is.element(amce_obj$baselines[[facet.name]],facet.levels[[facet.name]])) {
                    stop (paste(c("Error: Facet level \"",as.character(amce_obj$baselines[[facet.name]]), "\" is the baseline level of a profile varying attribute. Please provide alternative facet level or use defaults."), collapse=""))               
                }
                #names from user input if none provided
                if (is.null(names(facet.levels[[facet.name]]))) {
                    names(facet.levels[[facet.name]]) <- sapply(facet.levels[[facet.name]], USE.NAMES = F, function(x) amce_obj$user.levels[[x]])
                }
            } else if (is.null(names(facet.levels[[facet.name]]))) {
                names(facet.levels[[facet.name]]) <- as.character(facet.levels[[facet.name]])
            }
        }
    }

################################### Compile estimates into plottable objects

    #blank data frame for plot data
    d <- data.frame(pe=c(), se=c(), upper=c(), lower=c(), var=c(), group=c(), facet=c())
  
############# Unconditional estimates

    #only display if plot.display == all or unconditional
    if (plot.display != "interaction") {
        #if plot.display == all, add unconditional facet name (not needed for unconditional only)
        if (plot.display == "all") {
            uncond.facet.name <- "Unconditional"
        } else {
            uncond.facet.name <- NA
        }
        #if plot = all and there are non-respondent varying facet names
        #remove them from raw attributes
        if (plot.display == "all" && !is.null(facet.names)) {
            attr_remove <- c()
            for (facet.name in facet.names[!is.element(facet.names, amce_obj$respondent.varying)]) {
                attr_remove1 <- raw_attributes[grepl(":",raw_attributes)]
                attr_remove1 <- attr_remove1[grepl(facet.name,attr_remove1)]
                attr_remove <- c(attr_remove,attr_remove1)
            }
            raw_attributes <- raw_attributes[!is.element(raw_attributes,attr_remove)]
        }       
        #loop over raw attribute names
        for (i in 1:length(raw_attributes)) {
            #get raw attribute name
            attr_name <- raw_attributes[i]
            #get attribute name to print
            if (!is.null(attribute.names)) {
                print_attr_name <- attribute.names[which(names(amce_obj$estimates) == raw_attributes[i])]
            } else {
                variates <- strsplit(attr_name,":")[[1]]
                var_list <- c()
                for (var in variates) {
                    var_list <- c(var_list, amce_obj$user.names[[var]])
                }
                print_attr_name <- paste(var_list,collapse=":")
            }
            #set up basic group header and add to plot
            d_head <- data.frame(pe=NA, se=NA, upper=NA, lower=NA, var=paste(print_attr_name, ":", sep=""), group="<NA>",facet=uncond.facet.name)
            d <- rbind(d,d_head)    
            #iterate over levels
            for (j in 1:length(raw_levels[[attr_name]])) {
                #raw level name
                level_name <- raw_levels[[attr_name]][j]
                #get level name to print
                if (!is.null(level.names) && (attr_name %in% names(level.names))) {
                    print_level_name <- level.names[[attr_name]][j]
                } else {
                    levels <- strsplit(level_name,":")[[1]]
                    lev_list <- c()
                    for (lev in levels) {
                        lev_list <- c(lev_list, amce_obj$user.levels[[lev]])
                    }
                    print_level_name <- paste(lev_list,collapse=":") 
                }
                #if on the first level
                if (j == 1) {
                    if (label.baseline) {
                        print_level_name <- paste("(Baseline = ",print_level_name,")",sep="")
                    }
                    #get the baseline and print a blank line
                    d_lev <- data.frame(pe=NA, se=NA, upper=NA, lower=NA, var=paste("   ", print_level_name,sep=""), group=print_attr_name, facet=uncond.facet.name) 
                } else {
                    #retrieve estimate and SE
                    val_pe <- amce_obj$estimates[[attr_name]][1,level_name]
                    val_se <- amce_obj$estimates[[attr_name]][2,level_name]       
                    #calculate bounds
                    upper_bnd <- val_pe + zscr*val_se
                    lower_bnd <- val_pe - zscr*val_se 
                    #make line to add to plot data
                    d_lev <- data.frame(pe=val_pe, se=val_se, upper=upper_bnd, lower=lower_bnd, var=paste("   ", print_level_name,sep=""), group=print_attr_name, facet=uncond.facet.name)       
                } #end if a baseline
                #add to plot
                d <- rbind(d,d_lev)
            } #end loop over levels
        } #end loop over non-facet related attribute names
    } #end if plot.display == all or plot.display == conditional
        
############# Conditional estimates

    #Only if plot.display is all or conditional and we got a facet name from somehere
    if (plot.display != "unconditional" & !is.null(facet.names)) {

        #if user didn't give any levels, make blank list
        if (is.null(facet.levels)) facet.levels <- list()

        #loop over facets
        for (facet.name in facet.names) {    

            #how to print it
            print_facet_name <- amce_obj$user.names[[facet.name]]

            #### get objects for estimates
            if (facet.name %in% amce_obj$respondent.varying) {
                #if respondent varying betas from conditional
                estimate.source <- amce_obj$cond.estimates
                #also get original levels for all factors in the cond.formula
                xlevs <- sapply(all.vars(amce_obj$cond.formula)[-1] [all.vars(amce_obj$cond.formula)[-1] %in% names(amce_obj$baselines)], function(x) levels(amce_obj$data[[x]]), simplify=F)
            } else { #if a facet var
                #draw from unconditional estimates
                estimate.source <- amce_obj$estimates
            }
            # Extract vectors of betas for all coefficients
            beta.vector <- do.call(cbind,estimate.source)[1,]
            # SE's
            se.vector <- do.call(cbind,estimate.source)[2,]
            #get correct coefficient names
            all.coef <- c()
            for (var in names(estimate.source)) {
                prof_vars <- strsplit(var, ":")[[1]]
                for (lev in colnames(estimate.source[[var]])) {
                    prof_levs <- strsplit(lev, ":")[[1]]
                    orig_profs <- c()
                    for (v in 1:length(prof_vars)) {
                        #continuous component keeps simple name
                        if (prof_vars[v] == prof_levs[v]) {
                            orig_profs[v] <- prof_vars[v]
                        } else {
                            #factor component has 2 parts, var and level
                            orig_profs[v] <- paste(c(prof_vars[v], prof_levs[v]), collapse="")
                        }
                    }
                    orig_prof <- paste(orig_profs,collapse=":")
                    all.coef <- c(all.coef,orig_prof)
                }
            }
            names(beta.vector) <- all.coef
            names(se.vector) <- all.coef

            #### get levels of facet name
            #if not provided by user, get defaults
            if (!is.element(facet.name,names(facet.levels))) {
               #if it's a factor, default facet levels are all levels
                if (facet.name %in%  names(amce_obj$baselines)) {
                    #if NOT respondent varying get levels and names from ESTIMATES
                    if (facet.name %in% names(amce_obj$estimates)) {
                        facet.levels[[facet.name]] <- colnames(amce_obj$estimates[[facet.name]])
                        #names from user input
                        names(facet.levels[[facet.name]]) <- sapply(facet.levels[[facet.name]], USE.NAMES = F,function(x) amce_obj$user.levels[[x]])
                    } else {
                        #get levels and names from COND.ESTIMATES
                        facet.levels[[facet.name]] <- colnames(amce_obj$cond.estimates[[facet.name]])
                        #add in baseline!!
                        facet.levels[[facet.name]] <- c(amce_obj$baselines[[facet.name]], facet.levels[[facet.name]])
                        #names from user input
                        names(facet.levels[[facet.name]]) <- sapply(facet.levels[[facet.name]], USE.NAMES = F,function(x) amce_obj$user.levels[[x]])
                    }
                } else if (facet.name %in% names(amce_obj$continuous)) {
                    #if it's continuous, default is quantiles
                    facet.levels[[facet.name]] <-  amce_obj$continuous[[facet.name]]
                    #names stay as values, eg 25%
                    names(facet.levels[[facet.name]]) <- names(amce_obj$continuous[[facet.name]])
                } 
            }
            
            #### identify all REQUESTED terms involving facet name
            all_req_vars <- attr(terms(amce_obj$formula),"term.labels")
            all_mod <- unlist(sapply(all_req_vars,function(x) {
                y <- strsplit(x,":")[[1]]
                if (any(y == facet.name)) x
            }))
            #figure out profile attributes these refer to
            all_mod <- unlist(sapply(all_mod,function(x) {
                subs <- strsplit(x,":")[[1]]
                subs <- subs[is.element(subs,names(amce_obj$estimates))]
                subs <- subs[subs != facet.name]
                if (length(subs) > 0) paste(subs,collapse=":")
            }))
            #make sure there are some
            if (length(all_mod) == 0) {
                stop(paste(c("Error: Facet variable",facet.name,"not interacted with profile attributes"),collapse=" "))
            }
            #just unique ones
            all_mod <- unique(all_mod)

           #for each actual facet level make new set of plot data
            for (k in 1:length(facet.levels[[facet.name]])) {

                # get the appropriate model matrix for prediction
                # first, make a data.dummy data matrix
                data.dummy <- amce_obj$data
                # set level
                if (facet.name %in% names(amce_obj$continuous)) {
                    #name of "level" is just facet name
                    resp_lev <- facet.name
                    #as is original variable name
                    orig_facet <- facet.name
                } else {
                    # if facet is a factor, get level
                    resp_lev <- facet.levels[[facet.name]][k]
                    # make original var name
                    orig_facet <- paste(c(facet.name,resp_lev),collapse="")
                }
                #set level in data
                data.dummy[[facet.name]] <- facet.levels[[facet.name]][k]
                #how to print facet level
                if (is.element(facet.name,names(amce_obj$estimates))) {
                    #if ACIE
                    print_facet_level <- paste(c("ACIE", paste(c(print_facet_name, names(facet.levels[[facet.name]])[k]), collapse = " = ")), collapse = "\n")
                } else {
                    #if conditional on respondent varying
                    print_facet_level <- paste(c("Conditional on",paste(c(print_facet_name, names(facet.levels[[facet.name]])[k]), collapse = " = ")), collapse = "\n")
                }

                #loop over variables to be modified
                for (attr_name in all_mod) {
                    #split into components
                    variates <- strsplit(attr_name,":")[[1]]
                    #get attribute name to print
                    var_list <- c()
                    for (var in variates) {
                        var_list <- c(var_list, amce_obj$user.names[[var]])
                    }
                    if (!is.null(attribute.names)) {
                        print_attr_name <- attribute.names[which(names(amce_obj$estimates) == attr_name)]
                    } else {                       
                        print_attr_name <- paste(var_list,collapse=":")
                    }
                    #set up header to reflect base (non-facet) category
                    d_head <- data.frame(pe=NA, se=NA, upper=NA, lower=NA, var=paste(print_attr_name, ":", sep=""), group="<NA>", facet=print_facet_level)
                     #add new header
                    d <- rbind(d, d_head)

                    #iterate over levels of modified variable
                    for (p in 1:length(raw_levels[[attr_name]])) {
                        #raw level name
                        level_name <- raw_levels[[attr_name]][p]
                        #split it up
                        levels <- strsplit(level_name,":")[[1]]
                        #get original name of modified var and modify data dummy
                        orig_mods <- c()
                        lev_list <- c()
                        for (lev in 1:length(levels)) {
                            orig_mods[lev] <- paste(c(variates[lev], levels[lev]), collapse="")
                            data.dummy[[variates[lev]]] <- levels[lev]
                            lev_list <- c(lev_list, amce_obj$user.levels[[levels[lev]]])
                        }
                        orig_mod <- paste(orig_mods,collapse=":")
                        #get level name to print
                        if (!is.null(level.names) && (attr_name %in% names(level.names))) {
                            print_level_name <- level.names[[attr_name]][p]
                        } else {
                            print_level_name <- paste(lev_list,collapse=":") 
                        }
 
                       #get the baseline of modified var and make a blank line
                        if (p == 1) {
                            if (label.baseline) {
                                print_level_name <- paste("(Baseline = ",print_level_name,")",sep="")
                            }
                            d_lev <- data.frame(pe=NA, se=NA, upper=NA, lower=NA, var=paste("   ",  print_level_name,sep=""), group=print_attr_name, facet= print_facet_level)        
                        } else {

                            #figure out original interaction name
                            orig_inter <- grep(paste(c("(?=.*", orig_mod, ")(?=.*", orig_facet,")"),collapse=""),names(beta.vector),value=T,perl=T)              
                            #make sure it contains no other unrelated variables
                            orig_inter <- unlist(sapply(orig_inter,function(x) {
                                subs <- strsplit(x,":")[[1]]
                                subs <- subs[!is.element(subs,orig_mods)]
                                subs <- subs[!grepl(orig_facet,subs)]
                                if(length(subs) < 1) x
                            }))

                            ### if NOT respondent varying, straight grab
                            if (!facet.name %in% amce_obj$respondent.varying) {
                                #retrieve estimate and SE
                                val_pe <- beta.vector[orig_inter]
                                val_se <- se.vector[orig_inter]
                                #calculate bounds
                                upper_bnd <- val_pe + zscr*val_se
                                lower_bnd <- val_pe - zscr*val_se 
                            } else {
                                #add in modified var to orig_inter
                                pred.cols <- c(orig_mod,orig_inter)
                                #use modified data.dummy to make model matrix, keep old levels
                                pred.mat <- model.matrix(amce_obj$cond.formula, data.dummy, xlev = xlevs)
                                ##### do actual calculation for beta and variance
                                beta.inter <- pred.mat[1,pred.cols] %*% beta.vector[pred.cols]
                                if (!is.na(beta.inter)) {
                                    # for SE, first obtain relevant covariance
                                    all.cov <- c()
                                    for (a in 1:length(pred.cols)) {
                                        for (b in 1:length(pred.cols)) {
                                            all.cov <- c(all.cov, pred.mat[1,pred.cols[a]]*pred.mat[1,pred.cols[b]]*amce_obj$vcov.resp[pred.cols[a], pred.cols[b]])
                                        }
                                    }
                                    var.inter <- sum(all.cov)
                                    val_pe <- beta.inter
                                    val_se <- sqrt(var.inter)
                                    #calculate bounds
                                    upper_bnd <- val_pe + zscr*val_se
                                    lower_bnd <- val_pe - zscr*val_se 
                                } else {
                                    val_pe <- NA
                                    val_se <- NA
                                    upper_bnd <- NA
                                    lower_bnd <- NA
                                }
                            }
                            #make line to add to plot data
                            d_lev <- data.frame(pe=val_pe, se=val_se, upper=upper_bnd, lower=lower_bnd, var=paste("   ", print_level_name,sep=""), group=print_attr_name, facet=print_facet_level)       
                        } #end if p = 1

                        #add level data to plot data
                        d  <- rbind(d, d_lev)
                      
                    } #end loop over levels of modified var
                } #end loop over all modified vars
            } #end loop over level of facetted variable
        } #end loop over facets
    } else {
        #if there are no facets or plot.display is unconditional, remove that column
        d <- d[,1:6]
    }
    
#################    format "d" dataframe
    
  # Set Y bounds
    if (is.null(ylim)) {
        max_upper <- max(d$upper, na.rm=T) + .05
        min_lower <- min(d$lower, na.rm=T) - .05
        ylim <- c(min_lower, max_upper)
        d[is.na(d)] <- max_upper + 100
    } else {
        d[is.na(d)] <- max(ylim) + 100
    }

  # Make group factors <NA> actually NA
    d$group[d$group == "<NA>"] <- NA
    #same with facet
    if(!is.null(facet.names)) d$facet[d$facet == "<NA>"] <- NA

    #in case of duplicate level names
    realvar <- c()
    for(i in 1:nrow(d)) {
        var <-  as.character(d[i,"var"])
        group <- as.character(d[i,"group"])
        group <- ifelse(group == "<NA>","",group)
        realvar <- c(realvar,paste(c(group,var),collapse = ""))
    }
    d$realvar <- realvar
    
  # Reverse factor ordering
    d$realvar <- factor(d$realvar,levels=unique(d$realvar)[length(d$realvar):1])
    #d$var <- factor(d$var,levels=unique(d$var)[length(d$var):1])
    #make facet into factor, if it exists
    if (!is.null(facet.names)) {
        d$facet <- factor(d$facet,levels=unique(d$facet))
    }
    
########## plot output
                                        
    p = ggplot(d, aes(y=pe, x=realvar, colour=group))
    p = p + coord_flip(ylim = ylim)  
    p = p + geom_hline(yintercept = 0,size=.5,colour="black",linetype="dotted") 
    p = p + geom_pointrange(aes(ymin=lower,ymax=upper),position=position_dodge(width=dodge.size),size=point.size)

    #add facetting
    if (!is.null(facet.names)) {
        p = p + facet_wrap(~ facet)
    } 
    
  # If breaks and labels Null, use default
    if (is.null(breaks) & is.null(labels)) {
        p = p + scale_y_continuous(name=xlab)
    } else if (is.null(breaks) & !is.null(labels)) {
        p = p + scale_y_continuous(name=xlab, labels=labels)
    } else if (!is.null(breaks) & is.null(labels)) {
        p = p + scale_y_continuous(name=xlab, breaks=breaks)
    } else if (!is.null(breaks) & !is.null(labels)) {
        p = p + scale_y_continuous(name=xlab, breaks=breaks, labels=labels)
    }
  
    fix.xlabs <- as.character(d$var)[!duplicated(d$realvar)]
    p = p + scale_x_discrete(name="",labels=fix.xlabs[length(fix.xlabs):1])
  # If there's a title,add it
    if (!is.null(main)) {
        if (main != "") {
            p = p + ggtitle(main)
        }
    }
  # If no colors specified, use default
    if (is.null(colors)) {
        p = p + scale_colour_discrete(" ") 
    } else if (is.vector(colors)) {
    # Make manual palette
        cPal <- rep(colors, ceiling(length(unique(d$group))/length(colors)))
    # Use manual palette
        p = p + scale_colour_manual(values=cPal)
    } else {
        cat("Error: 'colors' must be a vector. Using default colors\n")
        p = p + scale_colour_discrete(" ") 
    }

  # colour scheme
    # if no theme specified, use default
    if (is.null(plot.theme)){
        
        theme_bw1 <- function(base_size = text.size, base_family = "") {
            theme_grey(base_size = base_size, base_family = base_family) %+replace%
            theme(axis.text.x = element_text(size = base_size*.9, colour = text.color,  hjust = .5 , vjust=1),axis.text.y = element_text(size = base_size, colour = text.color, hjust = 0 , vjust=.5 ), axis.ticks = element_line(colour = "grey50"),axis.title.y =  element_text(size = base_size,angle=90,vjust=.01,hjust=.1),plot.title = element_text(face = "bold"),legend.position = "none")
        }
        
        p = p + theme_bw1()
        print(p)
        
    } else if (is.null(class(plot.theme)))  {
      
      cat("Error: 'plot.theme' is not a valid ggplot theme object. Using default theme\n")
      theme_bw1 <- function(base_size = text.size, base_family = "") {
          theme_grey(base_size = base_size, base_family = base_family) %+replace%
          theme(axis.text.x = element_text(size = base_size*.9, colour = text.color,  hjust = .5 , vjust=1),axis.text.y = element_text(size = base_size, colour = text.color, hjust = 0 , vjust=.5 ), axis.ticks = element_line(colour = "grey50"),axis.title.y =  element_text(size = base_size,angle=90,vjust=.01,hjust=.1),plot.title = element_text(face = "bold"),legend.position = "none")
      }
      
      p = p + theme_bw1()
      print(p)
      
  } else if (class(plot.theme)[1] != "theme") {
      
      cat("Error: 'plot.theme' is not a valid ggplot theme object. Using default theme\n")
      theme_bw1 <- function(base_size = text.size, base_family = "") {
          theme_grey(base_size = base_size, base_family = base_family) %+replace%
          theme(axis.text.x = element_text(size = base_size*.9, colour = text.color,  hjust = .5 , vjust=1),axis.text.y = element_text(size = base_size, colour = text.color, hjust = 0 , vjust=.5 ), axis.ticks = element_line(colour = "grey50"),axis.title.y =  element_text(size = base_size,angle=90,vjust=.01,hjust=.1),plot.title = element_text(face = "bold"),legend.position = "none")
      }
      
      p = p + theme_bw1()
      print(p)
        
    # otherwise use the user-passed theme   
    } else {
        p = p + plot.theme
        print(p)
    }
      
}

#####################
# dependencies
######################

compute_dependencies <- function(J, tol=1e-14){
  # Get attribute names
  attribute_names <- names(dimnames(J))
  # If only one attribute, no dependence
  if (length(attribute_names) == 1){
    dependency_list <- list()
    dependency_list[[attribute_names[1]]] <- c()
    return(dependency_list)
  }else{
  # Create list for each attribute_name
  dependency_list <- list()
  for (k in attribute_names){
    dependency_list[[k]] <- c()
  }
  # Loop over each pair of attributes - figure out if they're independent
  for (i in 1:(length(attribute_names)-1)){
    for(j in (i+1):length(attribute_names)){
      attr1 <- attribute_names[i]
      attr2 <- attribute_names[j]
      cross_tab <- apply(J, c(attr1,attr2), sum)
      # Standardize
      sums <- apply(cross_tab, 1, sum)
      cross_tab_std <- cross_tab/sums
      # Compute similarities
      is_equal = TRUE
      r_1 <- cross_tab_std[1,]
      if (nrow(cross_tab_std) > 1){
        for (m in 2:nrow(cross_tab_std)){
          if (any(as.vector(r_1) - as.vector(cross_tab_std[m,]) > tol)){
            is_equal <- FALSE 
          }
        }
      }
      
      # If not the same, append to dependency dictionary
      if (!is_equal){
        dependency_list[[attr1]] <- c(dependency_list[[attr1]], attr2)
        dependency_list[[attr2]] <- c(dependency_list[[attr2]], attr1)
      } 
    }
  }
  return(dependency_list)
  }
}

##############################
# Make design matrix 
#############################

makeDesign <- function(type="file", J=NULL, filename=NULL, attribute.levels = NULL, constraints=NULL, level.probs=NULL, tol=1e-14){
  ## Initialize conjointDesign object
  design.obj <- NULL
  ## If type="file", then try to load from file generated by Conjoint SDT
  if (type == 'file'){
    if (is.null(filename)){
      cat("Error: Must provide a valid filename argument for type = 'file'\n")
      stop()
    }
    # Make a connection and read the lines
    connection <- file(filename, open="r")
    file_lines <- readLines(connection)
    close(connection)
    attr_index <- which(file_lines == "Attributes")
    weight_index <- which(file_lines == "Weights")
    restriction_index <- which(file_lines == "Restrictions")
    
    attributes <- file_lines[(attr_index+1):(weight_index-1)] 
    weight <- file_lines[(weight_index+1):(restriction_index-1)]
    if (restriction_index+1 != length(file_lines)){
      constr <- file_lines[(restriction_index+1):length(file_lines)]
    }else{
      constr <- NULL
    }
    attribute.levels <- list()
    for (attrstr in attributes){
      attributename <- strsplit(attrstr, ":")[[1]][1]
      levels <- strsplit(strsplit(attrstr, ":")[[1]][2], ",")[[1]]
      attribute.levels[[attributename]] <- levels
    }
    level.probs <- list()
    for (probstr in weight){
      attributename <- strsplit(probstr, ":")[[1]][1]
      weights <- strsplit(strsplit(probstr, ":")[[1]][2], ",")[[1]]
      level.probs[[attributename]] <- as.numeric(weights)
    }
    if (is.null(constr) != TRUE){
      constraints <- list()
      for (i in 1:length(constr)){
        allconstraints <- strsplit(constr[i], ";")[[1]]
        constraints[[i]] <- list()
        for (m in allconstraints){
          attributename <- strsplit(m, ":")[[1]][1]
          constrained_levels <- strsplit(strsplit(m, ":")[[1]][2], ",")[[1]]
          constraints[[i]][[attributename]] <- constrained_levels
        }
      }
    }else{
      constraints <- NULL
    }
  }
  ## if type = "array", check whether J is a valid array, then create conjointDesign object  
  if (type == 'array'){
    if (sum(J) != 1){
      cat("Error: Profile assignment probability array invalid: Does not sum to 1\n")
    }else{
      design.obj$J <- J
      design.obj$dependence <- compute_dependencies(J)
    }
  }else if (type == 'constraints' | type == 'file'){
    ## if type = "constraints"
    if (is.null(attribute.levels) | is.list(attribute.levels) != TRUE){
     cat("Error: Must provide a valid list() object in attribute.levels argument for type = 'constraints'\n")
    }
    # Calculate number of dimensions
    dimensions <- c()
    for (attr in names(attribute.levels)){
      dimensions <- c(dimensions, length(attribute.levels[[attr]]))
    }
    # Initialize 
    J_mat <- array(NA, dim=dimensions, dimnames = attribute.levels)

    # Fill in constrained cells with 0 probability
    for (cstr in constraints){
      # Save the names of the constrained attributes
      constraint_names <- names(cstr)
      # Construct a call to select relevant rows
      select_call <- Quote(J_mat[])
      select_call <- select_call[c(1, 2, rep(3, length(dim(J_mat))))]
      for (f in 1:length(constraint_names)){
        name <- constraint_names[f]
        index <- which(names(dimnames(J_mat)) == name)
        select_call[index+2] <- cstr[f]
      }
      # Make a Call
      eval(call("<-", select_call, 0))
    }
    
    # If no randomization weights, then uniform marginal randomization
    if (is.null(level.probs)){
      cell_prob <- 1/sum(is.na(J_mat))
      J_mat[is.na(J_mat)] <- cell_prob
    }else{
      # Normalize level.probs to sum to 1
      for (attr in names(level.probs)){
        level.probs[[attr]] <- level.probs[[attr]]/sum(level.probs[[attr]])
      }
      
      # If no names in level.probs, assume they're in order
      for (attr in names(level.probs)){
        if (is.null(names(level.probs[[attr]]))){
          names(level.probs[[attr]]) <- attribute.levels[[attr]]
        }
      }
      # If randomization weights specified: more calculations
      unconstrained_probs <- J_mat
      unconstrained_probs[TRUE] <- 1
      # For each Attribute
      for (attr in names(dimnames(J_mat))){
        # For each level in each Attribute
        for (level in attribute.levels[[attr]]){
          # Get a marginal probability
          marg_prob <- level.probs[[attr]][[level]]
        
          # Extract the rows pertaining to that marginal probability
          select_call <- Quote(unconstrained_probs[])
          select_call <- select_call[c(1, 2, rep(3, length(dim(J_mat))))]
          
          index <- which(names(dimnames(J_mat)) == attr)
          select_call[index+2] <- level
          # Make a Call
          
          eval(call("<-", select_call, eval(call("*",select_call, marg_prob))))
        }
      }
      missing_prob <- sum(unconstrained_probs[is.na(J_mat)==FALSE])
      increase_prob <- unconstrained_probs*1/(1-missing_prob)

      J_mat[is.na(J_mat)] <- increase_prob[is.na(J_mat)]
    }

    
    
    design.obj$J <- J_mat
    design.obj$dependence <- compute_dependencies(J_mat, tol)
    
  }else{
    cat("Invalid type argument: Must be either 'file', 'array', or 'constraints")
  }

  # Return design.obj if valid
  ## Make it a conjointDesign object 
  class(design.obj) <- "conjointDesign"
  return(design.obj)
}

###########################
# read in qualtrics data
###########################

read.qualtrics <- function(filename, responses, covariates=NULL, respondentID=NULL){
  # If nothing in "responses"

  if (is.null(responses)){
    stop("responses cannot be NULL")
    return(NULL)
  }
  
  # Load CSV Results
  qualtrics_results <- read.csv(filename, stringsAsFactors=F)

  # Extract variable names/question names
  var_names <- as.character(qualtrics_results[1,])
  q_names <- colnames(qualtrics_results)
  qualtrics_results <- qualtrics_results[2:nrow(qualtrics_results),]
  colnames(qualtrics_results) <- var_names
  
  # Find the attribute names
  attr_name_cols <- var_names[grep("F-[0-9]+-[0-9]+(?!-)", var_names, perl=TRUE)]
  
  # If no attributes fit the description
  if (length(attr_name_cols) == 0){
    stop("Error: Cannot find any columns designating attributes and levels. Please make sure the input file originated from a Qualtrics survey designed using the Conjoint SDT")
    return(NULL)
  }
  
  # Check whether attribute columns are empty or not
  for (attr_column in attr_name_cols){
    if (is.null(unique(qualtrics_results[,attr_column]))){
      stop(paste("Error, attribute column ", attr_column, " has no attribute names - recommend deleting this column"))
    }else if (unique(qualtrics_results[,attr_column])[1] == "" & length(unique(qualtrics_results[,attr_column])) == 1){
      stop(paste("Error, attribute column ", attr_column, " has no attribute names - recommend deleting this column"))
    }
  }
  
  # Parse to matrix
  attr_name_matrix <- matrix(unlist(strsplit(attr_name_cols,"-")),nrow=3,ncol=length(attr_name_cols))
  colnames(attr_name_matrix) <- attr_name_cols
  attr_name_matrix <- attr_name_matrix[2:nrow(attr_name_matrix),]
  attr_name_matrix <- as.data.frame(t(attr_name_matrix))

  num_tasks <-unique(as.integer(attr_name_matrix[,1]))
  
  # If number of responses doesn't match num_tasks
  if (length(num_tasks) != length(responses)){
    cat("Error: Number of response columns doesn't equal number of tasks in data")
    return(NULL)
  }
  
  # Find the level names
  level_name_cols <- var_names[grep("F-[0-9]+-[0-9]+-[0-9]", var_names, perl=TRUE)]

  # Check whether level columns are empty or not
  for (lev_column in level_name_cols){
    if (is.null(unique(qualtrics_results[,lev_column]))){
      stop(paste("Error, level column ", lev_column, " has no attribute names - recommend deleting this column"))
    }else if (unique(qualtrics_results[,lev_column])[1] == "" & length(unique(qualtrics_results[,lev_column])) == 1){
      stop(paste("Error, level column ", lev_column, " has no attribute names - recommend deleting this column"))
    }
  }  
  
  # Convert to matrix
  level_name_matrix <- matrix(unlist(strsplit(level_name_cols,"-")),nrow=4,ncol=length(level_name_cols))
  colnames(level_name_matrix) <- level_name_cols
  level_name_matrix <- level_name_matrix[2:nrow(level_name_matrix),]
  level_name_matrix <- as.data.frame(t(level_name_matrix))

  # If respondentID is not null
  if (!is.null(respondentID)){
    respondent_index <- qualtrics_results[,which(q_names %in% respondentID)]

  }else{
    respondent_index <- 1:nrow(qualtrics_results)
  }
  
  # Get the response rows
  if (is.character(responses[1])){
    response_vars <- which(q_names %in% responses)
  }else{
    response_vars <- responses
  }

  # Initialize output dataframe
  out_data_dataset <- NULL
  
  # Parse each row (this is probably really slow...)
  for (r in 1:nrow(qualtrics_results)){
    # Get the attributes
    skipRow = FALSE
    # If no attributes, skip the row
    attribute_refs_1 <- rownames(attr_name_matrix[attr_name_matrix[,1] == 1,])
    attribute_vector_1 <- qualtrics_results[r,attribute_refs_1]

    if (is.na(attribute_vector_1[1])){
      skipRow = TRUE
    }else if (attribute_vector_1[1] == ""){
      skipRow = TRUE
    }
    if (skipRow != TRUE){
    # Extract a covariate vector
    if (!is.null(covariates)){
      covariate_index <- which(q_names %in% covariates)
      covnames <- q_names[covariate_index]
      unit_cov <- qualtrics_results[r,covariate_index]

    }else{
      unit_cov <- c()
    }
    # For each question
    for (k in num_tasks){
       attribute_refs <- rownames(attr_name_matrix[attr_name_matrix[,1] == k,])
       
       attribute_vector <- qualtrics_results[r,attribute_refs]
       
       num_profiles <- as.integer(unique(level_name_matrix[,2]))


       selec_num <- qualtrics_results[r,response_vars[k]]

       # For each profile
       for (j in num_profiles){
         profile_ref <- rownames(level_name_matrix[level_name_matrix[,1] == k&level_name_matrix[,2] == j,])
         profile_levels <- qualtrics_results[r,profile_ref]


         names(profile_levels) <- attribute_vector
         
         if (is.na(as.integer(selec_num))){
           selec <- NA 
         }else if (as.integer(selec_num) == as.integer(j)){
           selec <- 1
         }else{
           selec <- 0
         }
         
         if (!is.null(covariates)){
           row_vec <- data.frame(r,respondent_index[r], k, j, profile_levels, selec, unit_cov)
  
           header <- as.vector(unlist(c("respondentIndex", "respondent","task","profile",attribute_vector, "selected", covnames)))
          
           colnames(row_vec) <- header
         }else{
           row_vec <- data.frame(r,respondent_index[r], k, j, profile_levels, selec)
           
           header <- as.vector(unlist(c("respondentIndex", "respondent","task","profile",attribute_vector, "selected")))
           
           colnames(row_vec) <- header
         }
         if (is.null(out_data_dataset)){
           out_data_dataset <- row_vec
         }else{
           out_data_dataset <- rbind(out_data_dataset, row_vec)
         }

       }
    }
    }
  }
  # Do some post-processing
  for (m in attribute_vector){
    out_data_dataset[[m]] <- as.factor(out_data_dataset[[m]])
  }
  out_data_dataset$respondentIndex <- as.integer(out_data_dataset$respondentIndex)
  out_data_dataset$selected <- as.integer(out_data_dataset$selected)
  out_data_dataset$task <- as.integer(out_data_dataset$task)
  out_data_dataset$profile <- as.integer(out_data_dataset$profile)
  
  # Return dataset
  return(out_data_dataset)
}





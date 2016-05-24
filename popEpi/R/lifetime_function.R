# #####################################################################################################################
# #Seppa K and Hakulinen T (2009) Mean and median survival times of cancer patients should be corrected for 
# #informative censoring. Journal of Clinical Epidemiology 62: 1095-1102.
# #
# #This R script calculates the crude and the bias-reduced estimates of the mean and the median survival times. 
# #The script has been tested using R version 2.15.2 for Windows and the Epi package version 1.1.40.
# #####################################################################################################################
# 
# lifetime <- function(patients, age_classes,interval_length,cut_point,end,r, popmort) {   # function lifetime
#   
#   ## INTENTION: only really intended to be tested against since Karri put a lot
#   ## of effort into this.
#   ## USAGE: E.g.
#   # 
#   # patients <- read.table("some_file.csv")
#   # attach(patients)
#   # # Variables:
#   # # diag_time = time of diagnosis in calender years
#   # # diag_age = age at diagnosis in years
#   # # d = censoring indicator, d=0 for a censored and d=1 for a deceased individual 
#   # # t = observed follow-up time in months
#   # 
#   # # Annual survival probabilities of Finnish women stratified by age (1-year age classes: 0-99 years (rows)) and 
#   # # calender time (1-year caleder periods: 1951-2005 (columns))
#   # popmort <- scan("some_popmort_file")
#   # popmort <- matrix(popmort,nrow=100)
#   # 
#   # library(Epi)  # Lexis() and splitLexis() functions (included in the Epi package, version 1.1.40) 
#   # 
#   # 
#   # lt <- lifetime(
#   #   patients = patients,
#   #   age_classes=c(0,25,45,65,75,100),	# Age classes: 0-24, 25-44, 45-64, 65-74, 75-99.
#   #   interval_length=6,			# Length of the interval in months: 6 months.
#   #   cut_point=c(9,9,9,9,9,9),		# Extrapolation is started after 9 years of survival both in the crude analysis (the first coordinate)
#   #   # and also in the analysis of every age class (the rest of the coordinates).
#   #   end=100,					# Extrapolation is continued for 100 years.
#   #   r=1,						# Constant interval specific relative survival ratio assumed in the extrapolation 
#   #   popmort = popmort
#   #   # (r<1, if patients have a persistent excess risk of death due to cancer)
#   # )
#   # 
#   # rm(popmort)
#   # detach(patients)
#   
#   diag_age <- evalq(diag_age, envir = patients)
#   diag_time <- evalq(diag_time, envir = patients)
#   diag_age <- evalq(diag_age, envir = patients)
#   
#   weights <- NULL
#   Results <- matrix(rep(NA,length(age_classes)*4),ncol=length(age_classes))
#   S <- matrix(rep(NA,(min(cut_point)+end)*12/interval_length*length(age_classes)),ncol=length(age_classes))
#   S_star <- matrix(rep(NA,(min(cut_point)+end)*12/interval_length*length(age_classes)),ncol=length(age_classes))
#   var_S <- matrix(rep(NA,(min(cut_point)+end)*12/interval_length*length(age_classes)),ncol=length(age_classes))
#   var_E <- NULL
#   
#   ################################
#   #Crude and age specific results
#   ################################
#   age_c <- NULL
#   for(age_c in 1:length(age_classes)) {
#     cat("Starting computation using age group ", age_c, 
#         "/",length(age_classes), "... \n", sep = "")
#     if (age_c == 1) {	#crude results
#       group <- subset(patients, diag_age >= age_classes[age_c] & diag_age < age_classes[length(age_classes)])
#       x <- cut_point[age_c]
#     }
#     if (age_c > 1) {	#age specific results
#       group <- subset(patients, diag_age >= age_classes[age_c-1] & diag_age < age_classes[age_c])
#       x <- cut_point[age_c]
#     }
#     weights[age_c] <- nrow(group)
#     
#     #####################
#     #Observed life table 
#     #####################
#     attach(group,warn.conflicts=F)
#     
#     lexis <- NULL
#     lex <- NULL
#     lexis <- Lexis(entry=rep(0,weights[age_c]), duration=t, exit.status = d)
#     cat("Splitting... \n")
#     lex <- splitLexisDT(lexis, breaks = seq(0,240,interval_length),
#                         timeScale = timeScales(lexis)[1], drop = FALSE)
#     # lex <- splitLexis(lexis, breaks = seq(0,240,interval_length))
#     setDF(lex)
#     cat("done splitting. \n")
#     
#     intervals<-ceiling(max(t/interval_length))
#     l <- NULL
#     w <- NULL
#     l_eff <- NULL
#     deaths <- NULL
#     p <- NULL
#     
#     cat("looping over intervals... \n")
#     setDT(lex)
#     l <- lex[.(entry = (1:(intervals+1)-1)*interval_length), 
#              .(obs = .N), on = "entry", by = .EACHI]$obs
#     deaths <- lex[.(lex.Xst = 1, entry = (1:(intervals)-1)*interval_length), 
#                   .(obs = .N), on = c("lex.Xst","entry"), by = .EACHI]$obs
#     # for(j in 1:(intervals+1)) l[j] <- sum(lex$entry == (j-1)*interval_length)
#     # for(j in 1:(intervals))	deaths[j] <- sum(lex$entry == (j-1)*interval_length & lex$lex.Xst == 1 )
#     
#     w <- -diff(l)-deaths
#     l <- l[-(intervals+1)]
#     l_eff <- l - w/2
#     p <- 1 - deaths/l_eff
#     
#     cumul <- cumprod(p)
#     
#     ###################################################
#     #Interval specific expected survival probabilities
#     ###################################################
#     cat("Interval-specific exp survs... \n")
#     p_star <- matrix(rep(NA,(end+x)*12/interval_length*nrow(group)),ncol=(end+x)*12/interval_length)
#     for(i in 1:((end+x)*12/interval_length)) {
#       rows <- ceiling(diag_age+(interval_length/12)*(i-1))
#       columns <- floor(diag_time+(interval_length/12)*(i-1)-1950)
#       rows[rows>nrow(popmort)] <- nrow(popmort)
#       columns[columns>ncol(popmort)] <- ncol(popmort)
#       p_star[,i] <- popmort[cbind(rows,columns)]^(interval_length/12)
#     }
#     
#     # Cumulative survival proportions of a comparable group in the general population
#     s_e <- c(1,apply(apply(p_star,1,cumprod),1,mean))
#     
#     ###########################################################
#     #Mean and median estimates, if extrapolation is NOT needed
#     ###########################################################
#     if (l[intervals]==deaths[intervals] & x*12/interval_length>intervals) {
#       cat("mean and median estimates (no extrapolation)... \n")
#       
#       s_p <- c(1,cumul,rep(0,x*12/interval_length-intervals))
#       
#       # Mean and median estimates
#       E <- (0.5+sum(s_p[-1]))*interval_length/12
#       E_star <- (0.5+sum(s_e[-1],na.rm = T))*interval_length/12
#       
#       time <- seq(0,length(s_p)-1,1)
#       md <- (max(time[s_p >= 0.5]) + (min(s_p[s_p >= 0.5]) - 0.5)/(min(s_p[s_p >= 0.5]) - max(s_p[s_p < 0.5])))*interval_length/12
#       md_star <- (max(time[s_e >= 0.5]) + (min(s_e[s_e >= 0.5]) - 0.5)/(min(s_e[s_e >= 0.5]) - max(s_e[s_e < 0.5])))*interval_length/12
#       
#       S[,age_c] <- c(s_p,rep(0,nrow(S)-length(s_p)))
#       S_star[,age_c] <- s_e[1:nrow(S_star)]
#       
#       Results[,age_c] <- c(E,E_star,md,md_star)
#       
#       # Variance estimates
#       f <- NULL
#       for(k in 1:(intervals-2)) {
#         f[k] <- c(1,cumul)[k]*( 1 + sum(cumprod(p[(k+1):(intervals-1)])) )
#       }
#       f[intervals-1] <- c(1,cumul)[intervals-1]
#       f[intervals] <- c(1,cumul)[intervals]
#       
#       var_p <- p*(1-p)/(l-w/2)
#       var_E[age_c] <- sum(f^2*var_p)
#       
#       var_S[1,age_c] <- 0
#       var_S[2:intervals,age_c] <- cumul[1:(intervals-1)]^2 * cumsum( deaths/(l_eff*(l_eff-deaths)) )[1:(intervals-1)]
#       var_S[(intervals+1):nrow(var_S),age_c] <- 0
#       
#     } #if extrapolation is NOT needed
#     
#     #######################################################
#     #Mean and median estimates, if extrapolation is needed
#     #######################################################
#     if (x*12/interval_length<=intervals) { 
#       
#       cat("mean and median estimates (using extrapolation)... \n")
#       
#       # Cumulative survival probabilities for a patient alive at the beginning of (x*12/interval_length+1)th interval
#       s_x_star <- r^(seq(1,end*12/interval_length)*interval_length/12) * 
#         apply(matrix(p_star[lex$lex.id[lex$entry==x*12],(x*12/interval_length+1):((end+x)*12/interval_length)],ncol=end*12/interval_length),1,cumprod) 
#       
#       # Expected survival proportions for the whole patient group alive at the beginning of the (x*12/interval_length+1)th interval
#       s_x_star_group <- apply(s_x_star,1,mean)  
#       
#       # Extrapolated cumulative survival proportions for a patient group (extrapolated from x to x + end years)
#       s_p <- c(1,cumul[1:(x*12/interval_length)],cumul[(x*12/interval_length)]*s_x_star_group)
#       
#       # Mean and median estimates
#       E <- (0.5+sum(s_p[-1]))*interval_length/12
#       E_star <- (0.5+sum(s_e[-1],na.rm = T))*interval_length/12
#       
#       time <- seq(0,nrow(S)-1,1)	
#       md <- (max(time[s_p >= 0.5]) + (min(s_p[s_p >= 0.5]) - 0.5)/(min(s_p[s_p >= 0.5]) - max(s_p[s_p < 0.5])))*interval_length/12
#       md_star <- (max(time[s_e >= 0.5]) + (min(s_e[s_e >= 0.5]) - 0.5)/(min(s_e[s_e >= 0.5]) - max(s_e[s_e < 0.5])))*interval_length/12
#       
#       S[,age_c] <- s_p[1:nrow(S)]
#       S_star[,age_c] <- s_e[1:nrow(S_star)]
#       
#       Results[,age_c] <- c(E,E_star,md,md_star)
#       
#       # Variance estimates
#       f <- NULL
#       for(k in 1:(x*12/interval_length-2)) {
#         f[k] <- c(1,cumul)[k]*( 1 + sum(cumprod(p[(k+1):(x*12/interval_length-1)])) + prod(p[(k+1):(x*12/interval_length)])*(1+sum(s_x_star_group)) )
#       }
#       f[x*12/interval_length-1] <- c(1,cumul)[x*12/interval_length-1]*( 1 + p[x*12/interval_length]*(1+sum(s_x_star_group)) )
#       f[x*12/interval_length] <- c(1,cumul)[x*12/interval_length]*(1+sum(s_x_star_group))
#       
#       var_p <- p*(1-p)/(l-w/2)
#       var_E[age_c] <- sum(f^2*var_p[1:(x*12/interval_length)])
#       
#       var_S[1,age_c] <- 0
#       var_S[2:(x*12/interval_length+1),age_c] <- (cumul^2 * cumsum( deaths/(l_eff*(l_eff-deaths)) ) )[1:(x*12/interval_length)]
#       var_S[(x*12/interval_length+2):nrow(var_S),age_c] <- s_x_star_group[-length(s_x_star_group)]^2 * (cumul^2 * cumsum( deaths/(l_eff*(l_eff-deaths)) ) )[x*12/interval_length]
#       
#     } #if extrapolation is needed
#     
#   } #for age classes
#   
#   ########################
#   #Bias-reduced estimates
#   ########################
#   
#   cat("overall bias-reduced estimates... \n")
#   
#   # Bias-reduced mean
#   E_weighted <- NULL
#   for(i in 1:2) {
#     E_weighted <- c( E_weighted,weighted.mean(Results[i,2:length(age_classes)],weights[2:length(age_classes)]) )
#   }
#   
#   # Bias-reduced cumulative survival proportions
#   S_weighted <- NULL
#   for(i in 1:nrow(S)) {
#     S_weighted <- c( S_weighted,weighted.mean(S[i,2:length(age_classes)],weights[2:length(age_classes)]) )
#   }
#   
#   # Bias-reduced median
#   time <- seq(0,nrow(S)-1,1)
#   md_bias_reduced <- ( max(time[S_weighted >= 0.5]) + (min(S_weighted[S_weighted >= 0.5]) - 0.5)/
#                          (min(S_weighted[S_weighted >= 0.5]) - max(S_weighted[S_weighted < 0.5])) )*interval_length/12
#   
#   # Standard errors of the mean
#   SE_E <- sqrt( (1/weights[1]^2)*sum( (weights[2:length(age_classes)]^2)*var_E[2:length(age_classes)] ) )*interval_length/12
#   SE_E_crude <- sqrt(var_E[1])*interval_length/12
#   
#   # Standard errors and 95% CI's for the cumulative survival proportions
#   SE_s_p_ages <- matrix(rep(NA,(min(cut_point)+end)*12/interval_length*(length(age_classes)-1)),ncol=(length(age_classes)-1))
#   for(i in 2:length(age_classes)) {
#     SE_s_p_ages[,i-1] <- weights[i]^2*var_S[,i]
#   }
#   SE_s_p <- sqrt( (1/weights[1]^2)*apply( SE_s_p_ages, 1, sum ) )
#   SE_s_p_crude <- sqrt( var_S[,1] )
#   
#   S_l <- S_weighted - qnorm(0.975)*SE_s_p
#   S_u <- S_weighted + qnorm(0.975)*SE_s_p
#   S_l_crude <- S[,1] - qnorm(0.975)*SE_s_p_crude
#   S_u_crude <- S[,1] + qnorm(0.975)*SE_s_p_crude
#   
#   # 95% CI's for the crude and the bias-reduced mean
#   CI_E_crude <-  c(Results[1,1]-qnorm(0.975)*SE_E_crude,Results[1,1]+qnorm(0.975)*SE_E_crude)
#   CI_E_weighted <- c(E_weighted[1]-qnorm(0.975)*SE_E,E_weighted[1]+qnorm(0.975)*SE_E)
#   
#   # The first and the last point where the 95% CI of the cumulative survival proportion covers the value of 50% 
#   CI_md_crude <- c(min(time[S_l_crude < 0.5])*interval_length/12,max(time[S_u_crude > 0.5])*interval_length/12)
#   CI_md_bias_reduced <- c(min(time[S_l < 0.5])*interval_length/12,max(time[S_u > 0.5])*interval_length/12)
#   
#   ########################
#   #Printing the estimates
#   ########################
#   
#   table <- rbind(t(Results[,2:length(age_classes)]),Results[,1],c(E_weighted,md_bias_reduced,Results[4,1]))
#   
#   stratification <- NULL
#   for(i in 1:(length(age_classes)-1)) {
#     stratification <- c(stratification,paste(age_classes[i],"-",age_classes[i+1]-1,sep=""))
#   }
#   
#   list(
#     matrix(cbind(round(cbind(table,table[,2]-table[,1],table[,4]-table[,3]),1),round(cbind((table[,2]-table[,1])/table[,2],(table[,4]-table[,3])/table[,4]),2)),ncol=8,
#            dimnames=list(c(stratification,"Crude","Bias-reduced"),c("E","E*","Md","Md*","E*-E","Md*-Md","(E*-E)/E*","(Md*-Md)/Md*")) ),
#     matrix(round(rbind(CI_E_crude,CI_E_weighted,CI_md_crude,CI_md_bias_reduced),1),ncol=2,
#            dimnames=list(c("Crude mean","Bias-reduced mean","Crude median","Bias-reduced median"),c("95% CI","")) )
#   )
#   
# } #function lifetime
# 

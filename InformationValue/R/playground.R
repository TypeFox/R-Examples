# # Functions for binning continuous variable
# # 
# # Create funcs to fine classing and coarse classing of factor variables. 
# # 
# # Add Confusion Matrix
# # 
# # Add KS statistic
# # 
# # Func for WOE for all vars in df
# # 
# # Func to get IV summary for all vars in a df
# # 
# # plot IVs
# 
# 
# # library(InformationValue)
# inputData <- read.csv("http://rstatistics.net/wp-content/uploads/2015/09/adult.csv")
# inp <- read.csv("http://rstatistics.net/wp-content/uploads/2015/09/adult.csv")
# 
# head(inputData)
# 
# factor_vars <- c ("WORKCLASS", "EDUCATION", "MARITALSTATUS", "OCCUPATION", "RELATIONSHIP", "RACE", "SEX", "NATIVECOUNTRY")
# 
# factor_var <- factor_vars[1]
# for(factor_var in factor_vars){
#   inputData[[factor_var]] <- WOE(X=inputData[, factor_var], Y=inputData$ABOVE50K)
# }
# 
# #   AGE  WORKCLASS FNLWGT  EDUCATION EDUCATIONNUM MARITALSTATUS OCCUPATION
# # 1  39  0.1608547  77516  0.7974104           13    -1.8846680  -0.713645
# # 2  50  0.2254209  83311  0.7974104           13     0.9348331   1.084280
# # 3  38 -0.1278453 215646 -0.5201257            9    -1.0030638  -1.555142
# # 4  53 -0.1278453 234721 -1.7805021            7     0.9348331  -1.555142
# # 5  28 -0.1278453 338409  0.7974104           13     0.9348331   0.943671
# # 6  37 -0.1278453 284582  1.3690863           14     0.9348331   1.084280
# 
# #   RELATIONSHIP        RACE        SEX CAPITALGAIN CAPITALLOSS HOURSPERWEEK
# # 1    -1.015318  0.08064715  0.3281187        2174           0           40
# # 2     0.941801  0.08064715  0.3281187           0           0           13
# # 3    -1.015318  0.08064715  0.3281187           0           0           40
# # 4     0.941801 -0.80794676  0.3281187           0           0           40
# # 5     1.048674 -0.80794676 -0.9480165           0           0           40
# # 6     1.048674  0.08064715 -0.9480165           0           0           40
# 
# #   NATIVECOUNTRY ABOVE50K WORKCLASS_WOE EDUCATION_WOE MARITALSTATUS_WOE
# # 1    0.02538318        0        0.2561        0.7772           -1.8079
# # 2    0.02538318        0        0.2050        0.7772            0.9162
# # 3    0.02538318        0       -0.1407       -0.4454           -1.0337
# # 4    0.02538318        0       -0.1407       -1.8525            0.9162
# # 5    0.11671564        0       -0.1407        0.7772            0.9162
# # 6    0.02538318        0       -0.1407        1.4590            0.9162
# 
# #   OCCUPATION_WOE RELATIONSHIP_WOE RACE_WOE SEX_WOE NATIVECOUNTRY_WOE AGE_WOE
# # 1        -0.7563          -0.9691   0.0823  0.3207            0.0311  0.4016
# # 2         1.1321           0.9246   0.0823  0.3207            0.0311  0.4016
# # 3        -1.3650          -0.9691   0.0823  0.3207            0.0311  0.4016
# # 4        -1.3650           0.9246  -0.7781  0.3207            0.0311  0.4016
# # 5         0.9655           0.9700  -0.7781 -0.9266           -0.4714 -0.8000
# # 6         1.1321           0.9700   0.0823 -0.9266            0.0311  0.4016
# 
# #   FNLWGT_WOE EDUCATIONNUM_WOE HOURSPERWEEK_WOE CAPITALGAIN_WOE
# # 1    -0.0569           1.7019           0.7471         -0.2037
# # 2    -0.0569           1.7019          -1.5126         -0.2037
# # 3    -0.1646          -0.4267           0.7471         -0.2037
# # 4    -0.1646          -1.5335           0.7471         -0.2037
# # 5    -0.1646           1.7019           0.7471         -0.2037
# # 6    -0.1646           1.7019           0.7471         -0.2037
# 
# all_factor_vars <- c(factor_vars, paste0(factor_vars, "_WOE"))
# 
# all_iv <- data.frame(VARS=factor_vars, IV=numeric(length(factor_vars)), STRENGTH=character(length(factor_vars)), stringsAsFactors = F)
# factor_var <- all_factor_vars[2]
# for (factor_var in factor_vars){
#   all_iv[all_iv$VARS == factor_var, "IV"] <- InformationValue::IV(X=inputData[, factor_var], Y=inputData$ABOVE50K)
#   all_iv[all_iv$VARS == factor_var, "STRENGTH"] <- attr(InformationValue::IV(X=inputData[, factor_var], Y=inputData$ABOVE50K), "howgood")
# }
# 
# all_iv <- all_iv[order(-all_iv$IV), ]
# 
# library(ggplot2)
# ggplot(all_iv, aes(x=reorder(VARS, IV), y=IV, fill=STRENGTH)) + geom_bar(stat = "identity") + coord_flip() + theme(legend.position="none") + labs(x="", y="Information Value", title="Information Value")
# 
# #>           VARS         IV            STRENGTH
# #>   RELATIONSHIP 1.53560810   Highly Predictive
# #>  MARITALSTATUS 1.33882907   Highly Predictive
# #>     OCCUPATION 0.77622839   Highly Predictive
# #>      EDUCATION 0.74105372   Highly Predictive
# #>            SEX 0.30328938   Highly Predictive
# #>      WORKCLASS 0.16338802   Highly Predictive
# #>  NATIVECOUNTRY 0.07939344 Somewhat Predictive
# #>           RACE 0.06929987 Somewhat Predictive
# 
# # a <- iv.mult(inputData[, c(factor_vars, "ABOVE50K")], y="ABOVE50K", summary = T)
# 
# #>   AGE  WORKCLASS FNLWGT  EDUCATION EDUCATIONNUM MARITALSTATUS OCCUPATION
# #> 1  39  0.1608547  77516  0.7974104           13    -1.8846680  -0.713645
# #> 2  50  0.2254209  83311  0.7974104           13     0.9348331   1.084280
# #> 3  38 -0.1278453 215646 -0.5201257            9    -1.0030638  -1.555142
# #> 4  53 -0.1278453 234721 -1.7805021            7     0.9348331  -1.555142
# #> 5  28 -0.1278453 338409  0.7974104           13     0.9348331   0.943671
# #> 6  37 -0.1278453 284582  1.3690863           14     0.9348331   1.084280
# 
# #>   RELATIONSHIP        RACE        SEX CAPITALGAIN CAPITALLOSS HOURSPERWEEK
# #> 1    -1.015318  0.08064715  0.3281187        2174           0           40
# #> 2     0.941801  0.08064715  0.3281187           0           0           13
# #> 3    -1.015318  0.08064715  0.3281187           0           0           40
# #> 4     0.941801 -0.80794676  0.3281187           0           0           40
# #> 5     1.048674 -0.80794676 -0.9480165           0           0           40
# #> 6     1.048674  0.08064715 -0.9480165           0           0           40
# 
# #>   NATIVECOUNTRY ABOVE50K
# #> 1    0.02538318        0
# #> 2    0.02538318        0
# #> 3    0.02538318        0
# #> 4    0.02538318        0
# #> 5    0.11671564        0
# #> 6    0.02538318        0
# 
# # KS Statistic
# data("ActualsAndScores")
# 
# ks_table <- function(actuals, predictedScores){
#   # sort the actuals and predicred scores and create 10 groups.
#   dat <- data.frame(actuals, predictedScores)
#   dat <- dat[order(-dat$predictedScores), ]
#   rows_in_each_grp <- round(nrow(dat)/10)
#   first_9_grps <- rep(1:9, each=rows_in_each_grp) 
#   last_grp <- rep(10, nrow(dat)-length(first_9_grps))
#   grp_index <- c(first_9_grps, last_grp)
#   dat <- cbind(grp_index, dat)
#   
#   # init the ks_table and make the columns.
#   ks_tab <- data.frame(rank=1:10, total_pop=as.numeric(table(dat$grp_index)))
#   ks_tab[c("non_responders", "responders")] <- as.data.frame.matrix(table(dat$grp_index, dat$actuals))
#   perc_responders_tot <- sum(ks_tab$responders)/sum(ks_tab$total_pop)  # percentage of total responders.
#   ks_tab$expected_responders_by_random <- ks_tab$total_pop * perc_responders_tot  # expected responders if there was no model.
#   ks_tab$perc_responders <- ks_tab$responders/sum(ks_tab$responders)
#   ks_tab$perc_non_responders <- ks_tab$non_responders/sum(ks_tab$non_responders)
#   ks_tab$cum_perc_responders <- cumsum(ks_tab$perc_responders)
#   ks_tab$cum_perc_non_responders <- cumsum(ks_tab$perc_non_responders)
#   ks_tab$difference <- ks_tab$cum_perc_responders - ks_tab$cum_perc_non_responders
#   return(ks_tab)
# }
# 
# # ks_table(a, p)
# 
# ks_stat <- function(actuals, predictedScores){
#   # the max of ks_table$difference
#   return(round(max(ks_table(actuals=actuals, predictedScores = predictedScores)$difference), 4))
# }
# 
# # ks_stat(a, p)
# 
# ks_plot <- function(actuals, predictedScores){
#   rank <- 0:10
#   model <- c(0, ks_table(actuals = actuals, predictedScores = predictedScores)$cum_perc_responders)*100
#   random <- seq(0, 100, 10)
#   df <- data.frame(rank, random, model)
#   df_stack <- stack(df, c(random, model))
#   df_stack$rank <- rep(rank, 2)
#   df_stack$delta <- df_stack$values[12:22]-df_stack$values[1:11]
#   
#   print(ggplot2::ggplot(df_stack, aes(x=rank, y=values, colour=ind, label=paste0(round(values, 2), "%"))) + geom_line(size=1.25) + labs(x="rank", y="Percentage Responders Captured", title="KS Plot") +
#           theme(plot.title = element_text(size=20, face="bold")) + geom_text(aes(y=values+4)))
# }
# 
# # ks_plot(a, p)



# R CMD build InformationValue
# R CMD check InformationValue_1.1.2.tar.gz --as-cran
# R CMD rd2pdf InformationValue
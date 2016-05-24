#' get_statistics_from_dataFrame
#' 
#' @import pwr plyr MASS ggplot2
#' @description 
#'    This function computes Cohen's f, f2 and w, adjusted p-value from GLM quasi-Poisson, negative binomial and Normal distribution.
#' @param
#' df_contrast
#'    A data frame that consists of 'ID' column and expression profile (columns after 'ID' column).
#'    'ID' column should be unique. Column names after 'ID' column should be unique.
#'    Only positive numbers are allowed in expression data. Here is an example.
#' \tabular{rrrrr}{
#'   ID \tab Y500U100_001 \tab Y500U100_002 \tab Y500U200_001 \tab Y500U200_002 \cr
#'   YKL060C \tab 151 \tab 195 \tab 188 \tab 184 \cr
#'   YDR155C \tab 154 \tab 244 \tab 237 \tab 232 \cr
#'   YOL086C \tab  64 \tab 89 \tab 128 \tab 109 \cr
#'   YJR104C \tab 161 \tab 155 \tab 158 \tab 172 \cr
#'   YGR192C \tab 157 \tab 161 \tab 173 \tab 175 \cr
#'   YLR150W \tab 96 \tab 109 \tab 113 \tab 115 \cr
#'   YPL037C \tab 23 \tab 28 \tab 27 \tab 27 \cr
#'   YNL007C \tab 53 \tab 58 \tab 64 \tab 63 \cr
#'   YBR072W \tab 52 \tab 53 \tab 54 \tab 44 \cr
#'   YDR418W_1 \tab 76 \tab 53 \tab 62 \tab 74 \cr
#'   }
#' @param
#' df_group
#'    A data frame that consists of 'Col_Name' and 'Group' columns
#'    This parameter is to match experiment groups to expression profiles of df_contrast.
#'    'Col_Name' should be corresponding to column names of expression profile of df_contrast.
#'    'Group' columns have experiment informaion of columns in expression profile of df_contrast.  Here is an example. See the example of df_contrast together.
#'    \tabular{rr}{
#'    Col_Name \tab Group \cr
#'    Y500U100_001 \tab U100 \cr
#'    Y500U100_002 \tab U100 \cr
#'    Y500U200_001 \tab U200 \cr
#'    Y500U200_002 \tab U200 \cr
#'    }
#' @param
#' padj
#'    Choose one of these c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). 
#'    "fdr" is default option. The option is same to \code{\link[stats]{p.adjust}}.
#' @return
#' A list that consists of the following items:
#'    \tabular{ll}{
#'    $data_table \tab A data frame that have statistics for each IDs\cr
#'    $min_rep \tab Common number of replicates in your group information. \cr
#'    $max_rep \tab Maximum number of replicates in your group information. \cr
#'    $nt \tab The number of total experiments in your expression profile. \cr
#'    $ng \tab The number of groups in your group information. \cr
#'    $method_pvalue_adjustment \tab The selected method for p-value adjustment \cr
#'    }
#'    \tabular{ll}{
#'    data_table's elements \tab  \cr
#'    Cohens_W \tab Cohen's w \cr
#'    Cohens_F \tab Cohen's f\cr
#'    Cohens_F2 \tab Cohen's f2\cr
#'    Max_FC \tab Maximum fold change among all the possible group pairs\cr
#'    QP_Pval_adjusted \tab Adjusted p-value from GLM quasi-Poisson \cr
#'    NB_Pval_adjusted \tab Adjusted p-value from GLM negative binomial \cr
#'    Normal_Pval_adjusted \tab Adjusted p-value from Normal ANOVA \cr
#'    }
#' @examples
#' library(selfea)
#' 
#' ## Test selfea for single protein expression
#' values <- c(6,8,10,29,26,22)
#' groups <- c("U200","U200","U200","U600","U600","U600")
#' experiments <- c("exp1","exp2","exp3","exp4","exp5","exp6")
#' 
#' df_expr <- data.frame(ID="Protein_1",exp1=6,exp2=8,exp3=10,exp4=29,exp5=26,exp6=22)
#' df_group <- data.frame(Col_Name=experiments,Group=groups)
#' list_result <- get_statistics_from_dataFrame(df_expr,df_group)
#' top_table(list_result)
#' 
#' ## For this example we will import Gregori data
#' ## Josep Gregori, Laura Villareal, Alex Sanchez, Jose Baselga, Josep Villanueva (2013).
#' ## An Effect Size Filter Improves the Reproducibility 
#' ## in Spectral Counting-based Comparative Proteomics.
#' ## Journal of Proteomics, DOI http://dx.doi.org/10.1016/j.jprot.2013.05.030')
#' 
#' ## Description:
#' ## Each sample consists in 500ng of standard yeast lisate spiked with 
#' ## 100, 200, 400 and 600fm of a mix of 48 equimolar human proteins (UPS1, Sigma-Aldrich).
#' ## The dataset contains a different number of technical replimessagees of each sample
#' 
#' ## import Gregori data
#' data(example_data1)
#' df_contrast <- example_data
#' df_group <- example_group
#' 
#' ## Get statistics through 'get_statistics_from_dataFrame' function
#' list_result <- get_statistics_from_dataFrame(df_contrast,df_group)
#' 
#' ## Get significant features (alpha >= 0.05 and power >= 0.90)
#' significant_qpf <- top_table(list_result,pvalue=0.05,power_desired=0.90,method='QPF')
#' 
#' @export
#'
get_statistics_from_dataFrame <- function(df_contrast, df_group, padj = 'fdr')
{
  # Check if df_contrast is data frame
  if (!is.data.frame(df_contrast)) {
    stop('df_contrast should be data frame')
  }
  
  # Check if df_group is data frame
  if (!is.data.frame(df_group)) {
    stop('df_group should be data frame')
  }

  # Check if the first column of contrast file is "ID" column
  if (!colnames(df_contrast)[1] == "ID"){
    message("\n'ID' column is not found, so first column is considered as ID column")
    dataset.ID <- df_contrast[,1]
  } else {
    dataset.ID <- df_contrast$ID
  }

  message("Expression data dimension:\t",nrow(df_contrast)," IDs x ",ncol(df_contrast)-1," Columns")
  
  adj_options <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if (padj %in% adj_options) {
    message("You chose \"",padj,"\" for p-value adjustment.\n")
  } else {
    stop('Wrong options for p-value adjustment. It should be one of holm, hochberg, hommel, bonferroni, BH, BY, fdr and none')
  }

  # Change factor vector into character vector
  if(is.factor(dataset.ID))
    dataset.ID <- as.character(dataset.ID)
  
  # Check group data frame has two columns
  if (length(colnames(df_group)) != 2) {
    stop("The number of columns of 'Group' from your input should have two columns. Please check the manual.\n")
  } 
  
  # Check if the first column of group file is "Col_Name"
  if (!colnames(df_group)[1] == "Col_Name"){
    message("'Col_Name' column is not found in your Group information, so first column is considered as column names of contrast file")
    colName_group <- df_group[,1]
  } else {
    colName_group <- df_group$Col_Name
  }
  
  # Change factor vector into character vector
  if(is.factor(colName_group))
    colName_group <- as.character(colName_group)
  
  group <- df_group$Group
  # Check if the second column of group file is "Group"
  if (!colnames(df_group)[2] == "Group"){
    message("'Group' column is not found in your Group information, so second column is considered as group names")
    group <- df_group[,2]
  } else {
    group <- df_group$Group
  }
  
  # Change factor vector into character vector
  if(is.factor(group))
    group <- as.character(group)
  
  dataset.expr <- df_contrast[,seq(2,length(df_contrast))]
  colName_expr <- colnames(dataset.expr)
  
  # Check if column names in contrast file and group file are same
  if (length(colName_group) != length(colName_expr)) {
    stop("The number of experiments between group and expression files are not same.  Also, the order of column names should be same.")
  }
  num_col <- length(colName_group)
  for(i in 1:num_col) {
    if (colName_group[i] != colName_expr[i]) {
      stop("There is different column names between contrast file and group file. Check your group file and try again. The order of column names should be same.")
    }
  }
  
  # Count replicates of each groups
  message("Group dimension:")
  exp_groups <- unique(group)
  num_groups <- length(exp_groups)
  num_reps <- c()
  for (j in 1:num_groups){
    group_name <- exp_groups[j]
    num_replicates <- length(group[group==group_name])
    num_reps <- c(num_reps,num_replicates)
    message("\t\t\t\tGroup ",group_name," has ",num_replicates," replicates")
  }

  list_return <- list(data_table = glm_anova(dataset.expr, dataset.ID, group, padj), 
                      min_rep = min(num_reps),
                      max_rep = max(num_reps),
                      nt = num_col,
                      ng = num_groups,
                      method_pvalue_adjustment=padj)
  return (list_return)
}


#' get_statistics_from_file
#' 
#' @import pwr plyr MASS ggplot2
#' @description 
#'    This function computes Cohen's f, f2 and w, adjusted p-value from GLM quasi-Poisson, negative binomial and Normal distribution.
#' @param 
#' file_expr
#'    a CSV type file, comma (,) seperated file format, that has unique "ID" at the first column and expression data for the corresponding ID.  Here is an short example.
#'    \tabular{l}{
#'    ID,Y500U100_001,Y500U100_002,Y500U200_001,Y500U200_002 \cr
#'    YKL060C,151,195,221,201 \cr
#'    YDR155C,154,244,190,187 \cr
#'    YOL086C,64,89,116,119 \cr
#'    }
#' 
#' @param 
#' file_group
#'    a CSV type file, comma (,) seperated file format, that consists of "Col_Name", column names of "file_expr" parameter, and "Group" information of the corresponding column name.
#'    The order of "Col_Name" column have to be same to order of columns in "file_expr".  Here is an example.  See also the example above.
#'    \tabular{l}{
#'    Col_Name,Group \cr
#'    Y500U100_001,U100 \cr
#'    Y500U100_002,U100 \cr
#'    Y500U200_001,U200 \cr
#'    Y500U200_002,U200 \cr
#'    }
#'    
#' @param
#' padj
#'    Choose one of these c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). 
#'    "fdr" is default option. The option is same to \code{\link[stats]{p.adjust}}.
#' @return 
#' A list that consists of the following items:
#'    \tabular{ll}{
#'    $data_table \tab A data frame that have statistics for each IDs\cr
#'    $min_rep \tab Common number of replicates in your group information. \cr
#'    $max_rep \tab Maximum number of replicates in your group information. \cr
#'    $nt \tab The number of total experiments in your expression profile. \cr
#'    $ng \tab The number of groups in your group information. \cr
#'    $method_pvalue_adjustment \tab The selected method for p-value adjustment \cr
#'    }
#'    \tabular{ll}{
#'    data_table's elements \tab  \cr
#'    Cohens_W \tab Cohen's w \cr
#'    Cohens_F \tab Cohen's f\cr
#'    Cohens_F2 \tab Cohen's f2\cr
#'    Max_FC \tab Maximum fold change among all the possible group pairs\cr
#'    QP_Pval_adjusted \tab Adjusted p-value from GLM quasi-Poisson \cr
#'    NB_Pval_adjusted \tab Adjusted p-value from GLM negative binomial \cr
#'    Normal_Pval_adjusted \tab Adjusted p-value from Normal ANOVA \cr
#'    }
#' @examples
#' library(selfea)
#' 
#' ## For this example we will import Gregori data
#' ## Josep Gregori, Laura Villareal, Alex Sanchez, Jose Baselga, Josep Villanueva (2013).
#' ## An Effect Size Filter Improves the Reproducibility 
#' ## in Spectral Counting-based Comparative Proteomics.
#' ## Journal of Proteomics, DOI http://dx.doi.org/10.1016/j.jprot.2013.05.030')
#' 
#' ## Description:
#' ## Each sample consists in 500ng of standard yeast lisate spiked with 
#' ## 100, 200, 400 and 600fm of a mix of 48 equimolar human proteins (UPS1, Sigma-Aldrich).
#' ## The dataset contains a different number of technical replimessagees of each sample
#' 
#' ## Import Gregori data
#' data(example_data1)
#' df_contrast <- example_data
#' df_group <- example_group
#' 
#' ## Write Gregori data to use 'get_statistics_from_file' function
#' write.csv(df_contrast,"expression.csv",row.names=FALSE)
#' write.csv(df_group,"group.csv",row.names=FALSE)
#' 
#' ## Get statistics
#' list_result <- get_statistics_from_file("expression.csv","group.csv","fdr")
#' 
#' ## Get significant features (alpha >= 0.05 and power >= 0.90)
#' significant_qpf <- top_table(list_result,pvalue=0.05,power_desired=0.90,method='QPF')
#' 
#' @export
get_statistics_from_file <- function(file_expr = '',file_group = '', padj = 'fdr')
{
  if (file_expr == '') {
    stop('file_expr parameter is empty. Please put valid file path\n')

  } else if (!file.exists(file_expr)){
    stop(file_expr,' does not exist.  Please check the file path\n')
    
  } else if (file_group == '') {
    stop('file_group parameter is empty. Please put valid file path\n')
    
  } else if (!file.exists(file_group)){
    stop(file_group,' does not exist.  Please check the file path\n')
    
  } else {
    message('Expression dataset:\t\t',file_expr)
    message('Group dataset:\t\t\t',file_group)
    
    df_contrast <- read.csv(file_expr, check.names=F)
    df_group <- read.csv(file_group, check.names=F)
  }
  
  ### Use get_statistics_from_dataFrame
  return (get_statistics_from_dataFrame(df_contrast, df_group, padj))
}

#' calculate_cohen_f2
#' 
#' @import pwr plyr MASS ggplot2
#' @description 
#'    Calculate Cohen's f2. Followed formulars at wikipages (https://en.wikipedia.org/wiki/Effect_size , https://en.wikipedia.org/wiki/Coefficient_of_determination)
#' @param 
#' model_glm
#'    GLM model generated by 'glm' function
#' @param
#' df_aov
#'    A data frame containing groups in 'Run' column and values in 'SC' column
#' @return
#' Cohen's f2 (an effect size for linear models)
calculate_cohen_f2 <- function(model_glm, df_aov)
{
  SStot <- sum((df_aov$SC - mean(df_aov$SC))^2)
  SSres <- sum((df_aov$SC - model_glm$fitted.values)^2)
  R2 <- 1 - (SSres/SStot)
  return(R2/(1-R2))
}


#' glm_anova
#' 
#' @import pwr plyr MASS ggplot2
#' @description 
#'    Calculate P-values from ANOVA using Normal, Quasi-Poisson and Negative Binomial distribution and Cohen's effect sizes
#' @param 
#' dataset.expr
#'    A data frame that has column names for distinguishing experiments and numerical values for expression levels
#' @param
#' dataset.ID
#'    A vector of the obtained expression profile's ID column
#' @param
#' group
#'    A data frame that consists of 'Col_Name' and 'Group' obtained from the user file through get_statistics_from_file.
#' @param
#' padj
#'    Choose one of these c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"). 
#'    "fdr" is default option.
#' @return
#' A data frame containing ID, Cohen's W, Cohen's F, Max fold change,
#'    GLM Negative Binomial P-value, GLM Quasi-Poisson P-value and ANOVA with Normal P-value.
glm_anova <- function(dataset.expr, dataset.ID, group, padj ='fdr')
{
  pvals_nor <- c()
  pvals_qp <- c()
  pvals_nb <- c()
  pvals_logN <- c()

  cohen_f2_qp <- c()
  cohen_f2_logN <- c()
  
  cohen_f <- c()
  cohen_f_qp <- c()
  
  cohen_w <- c()
  max_fcs <- c()
  ID_return <- c()
  SC <- NULL
  NR <- nrow(dataset.expr)
  for(i in 1:NR) {
    sc <- as.numeric(dataset.expr[i,])
    if (sum(sc) <= 0) {
      message ("In your data ", dataset.ID[i], " has only zero values. Selfea can't use it, so skip this data.")
    } else {
      ID_return <- c(ID_return, dataset.ID[i])
      df_aov <- data.frame(SC=sc,Run=group)

      ### Gaussian
      anv <- aov(SC ~ Run, df_aov)
      sum_anv <- summary(anv)
      pvals_nor <- c(pvals_nor,sum_anv[[1]][["Pr(>F)"]][[1]])
      
      ### ANOVA calculation
      ### http://www.itl.nist.gov/div898/handbook/prc/section4/prc434.htm
      # CM <- sum(sc)^2/length(sc)
      # SST <- (sum(sc[1:4])^2)/4 + (sum(sc[5:10])^2)/6 + (sum(sc[11:13])^2)/3 + (sum(sc[14:19])^2)/6 - CM
      # SStot <- sum(sc^2) - CM
      # SSE <- SStot - SST
      
      SST <- sum_anv[[1]][["Sum Sq"]][[1]]
      SSE <- sum_anv[[1]][["Sum Sq"]][[2]]
      MST <- SST/sum_anv[[1]][["Df"]][[1]]
      MSE <- SSE/sum_anv[[1]][["Df"]][[2]]
      anv_f <- MST/MSE
      
      ## Log Normal
      sc4logN <- replace(sc, sc==0, 0.0001)
      df_aov_4logN <- data.frame(SC=sc4logN,Run=group)
      model_glm <- glm(SC~Run,data=df_aov_4logN,family=gaussian(link = "log"))
      model_glm0 <- glm(SC~1,data=df_aov_4logN,family=gaussian(link = "log"))
      p.value <- anova(model_glm,model_glm0,test="F")$"Pr(>F)"[2]
      pvals_logN <- c(pvals_logN,p.value)
      cohen_f2_logN <- c(cohen_f2_logN,calculate_cohen_f2(model_glm, df_aov_4logN))
      
      ### QuasiPoisson
      model_glm <- glm(SC~Run,data=df_aov,family=quasipoisson(link = "log"))
      model_glm0 <- glm(SC~1,data=df_aov,family=quasipoisson(link = "log"))
      p.value <- anova(model_glm,model_glm0,test="F")$"Pr(>F)"[2]
      pvals_qp <- c(pvals_qp,p.value)
      cohen_f2_qp <- c(cohen_f2_qp,calculate_cohen_f2(model_glm, df_aov))
      
      # https://en.wikipedia.org/wiki/Mean_squared_error
      MSE_qp <- sum((df_aov$SC - model_glm$fitted.values)^2)/model_glm$df.residual
      
      ### GLM using Negative Binomial Dist.
      model_glm <- MASS::glm.nb(SC~Run,data=df_aov)
      model_glm0 <- MASS::glm.nb(SC~1,data=df_aov)
      p.value <- anova(model_glm,model_glm0)$"Pr(Chi)"[2]
      pvals_nb <- c(pvals_nb,p.value)
      
      ### Cohen's f
      numerator <- 0
      N <- nrow(df_aov)
      # http://www.r-tutor.com/elementary-statistics/numerical-measures/variance
      # https://onlinecourses.science.psu.edu/stat501/node/254
      grand_mean <- mean(df_aov$SC)
      df_run2mean <- plyr::ddply(df_aov, "Run", plyr::summarize, mean=mean(SC))
      df_run2sd <- plyr::ddply(df_aov, "Run", plyr::summarize, sd=sd(SC))
      df_run2nrow <- plyr::ddply(df_aov, "Run", nrow)
      for (j in 1:nrow(df_run2mean)) {
        p = df_run2nrow[j,2]/N
        numerator <- numerator + p*(grand_mean-df_run2mean[j,2])^2
      }
      cohen_f     <- c(cohen_f,sqrt(numerator/MSE))
      cohen_f_qp  <- c(cohen_f_qp,sqrt(numerator/MSE_qp))
      
      ### Cohen's w
      num_runs <- length(sc)
      p0 <- c(rep(1/num_runs,num_runs))
      p1 <- df_aov$SC/sum(df_aov$SC)
      cohen_w <- c(cohen_w,pwr::ES.w1(p0,p1))
      
      ### Max fold change
      max_fc <- 0
      df_run2sum <- plyr::ddply(df_aov, "Run", plyr::summarize, sum=sum(SC))
      for (k in 1:nrow(df_run2sum)) {
        for (l in 1:nrow(df_run2sum)) {
          group_numer <- as.character(df_run2sum$Run[k])
          group_denom <- as.character(df_run2sum$Run[l])
          
          denom <- df_run2sum$sum[df_run2sum$Run == group_denom]
          numer <- df_run2sum$sum[df_run2sum$Run == group_numer]
          fc <- numer/denom
          
          # Code for debugging
          #cat(i)
          #cat('\t')
          #cat(df_run2sum$sum[df_run2sum$Run == group_numer])
          #cat('\t')
          #cat(df_run2sum$sum[df_run2sum$Run == group_denom])
          #cat('\t')
          #cat(fc)
          #cat('\n')
          
          if (fc == 'NaN') {
            max_fc <- max_fc
          } else if (fc > max_fc) {
            max_fc <- fc
          }
          
        }
      }
      max_fcs <- c(max_fcs, max_fc)
    }
  }
  
  adj.pval_qp <- p.adjust(pvals_qp,method=padj,n=length(pvals_qp))
  adj.pvals_nb <- p.adjust(pvals_nb,method=padj,n=length(pvals_nb))
  adj.pval_nor <- p.adjust(pvals_nor,method=padj,n=length(pvals_nor))
  adj.pval_logN <- p.adjust(pvals_logN,method=padj,n=length(pvals_logN))
  
  df_out <- data.frame(Protein=ID_return, 
                       Cohens_W=cohen_w, 
                       Cohens_F=cohen_f, 
                       Cohens_F2=cohen_f2_qp, 
                       Max_FC=max_fcs,
                       QP_Pval_adjusted=adj.pval_qp, 
                       NB_Pval_adjusted=adj.pvals_nb, 
                       Normal_Pval_adjusted=adj.pval_nor)
  return (df_out)
}

#' top_table
#' 
#' @import pwr plyr MASS ggplot2
#' @description 
#'    Get IDs that pass two filters, p-value and effect-size.  This top_table will make a significant list that is less than p-value and greater than effect-size.  Effect-size are calculated by obtained power level.
#'    This function requires four parameters. ex) top_table(input_data,pvalue=0.05,power_desired=0.90,method='QPF')
#' @param
#' input_list
#'    The list should be produced by 'get_statistics_from_file' or 'get_statistics_from_dataFrame' function.
#'    See \code{\link{get_statistics_from_file}} and \code{\link{get_statistics_from_dataFrame}} for more information.
#'    It consists of the following items:
#'    \tabular{ll}{
#'    $data_table \tab A data frame that have statistics for each IDs\cr
#'    $min_rep \tab Common number of replicates in your group information. \cr
#'    $max_rep \tab Maximum number of replicates in your group information. \cr
#'    $nt \tab The number of total experiments in your expression profile. \cr
#'    $ng \tab The number of groups in your group information. \cr
#'    }
#' @param 
#' pvalue
#'    p-value should be ranged between 0 to 1.
#'    default is 0.05.
#' @param 
#' power_desired
#'    Give the statistical power you desired for output significant list
#' @param 
#' method
#'    Choose statistics method you want to use for making significant list
#'    \tabular{ll}{
#'    "QPF" \tab combination of Quasi-Poisson and Cohen's f. Default. \cr
#'    "QPF2" \tab combination of Quasi-Poisson and Cohen's f2. \cr
#'    "QPFC" \tab combination of Quasi-Poisson and Fold change. \cr
#'    "NBW" \tab combination of Negative Binomial and Cohen's w. \cr
#'    "NBF2" \tab combination of Negative Binomial and Cohen's f2. \cr
#'    "NBFC" \tab combination of Negative Binomial and Fold change. \cr
#'    "NORF" \tab combination of ANOVA with normal distribution and Cohen's f. \cr
#'    "NORFC" \tab combination of ANOVA with normal distribution and Fold change. \cr
#'    }
#' @param
#' FC_threshold
#'    Fold change you want to use. Default is 2.
#' @return 
#' A list containing the follow items and a scatter plot that x-axis is effect size and y-axis is probability.  
#' Vertical line the plot is minimum effect size and horizontal line is maximum probability threshold.
#' Red dots means insignificant, while blue dots are significant.
#'    \tabular{ll}{
#'    top_table \tab a data frame that have calculated statistics for top table IDs \cr
#'    minimum_effect_size \tab Minimum effect size threshold \cr
#'    selected_effect_size_filter \tab The selected effect size filter \cr
#'    minimum_power \tab Minimum statistical power in the top_table \cr
#'    selected_model \tab The selected probability model for calculating p-value \cr
#'    alpha \tab Maximum adjusted p-value \cr
#'    method_pvalue_adjustment \tab The selected method for p-value adjustment \cr
#'    num_group \tab The number of groups used for generating the top_table \cr
#'    common_replicates \tab The number of common replicates. \cr
#'    num_columns \tab The number of columns (samples or experiments) \cr
#'    }
#'    \tabular{ll}{
#'    top_table's elements \tab  \cr
#'    Cohens_W \tab Cohen's w \cr
#'    Cohens_F \tab Cohen's f\cr
#'    Cohens_F2 \tab Cohen's f2\cr
#'    Max_FC \tab Maximum fold change among all the possible group pairs\cr
#'    QP_Pval_adjusted \tab Adjusted p-value from GLM quasi-Poisson \cr
#'    NB_Pval_adjusted \tab Adjusted p-value from GLM negative binomial \cr
#'    Normal_Pval_adjusted \tab Adjusted p-value from Normal ANOVA \cr\cr
#'    }
#' @examples
#' library(selfea)
#' 
#' ## Test selfea for single protein expression
#' values <- c(6,8,10,29,26,22)
#' groups <- c("U200","U200","U200","U600","U600","U600")
#' experiments <- c("exp1","exp2","exp3","exp4","exp5","exp6")
#' 
#' df_expr <- data.frame(ID="Protein_1",exp1=6,exp2=8,exp3=10,exp4=29,exp5=26,exp6=22)
#' df_group <- data.frame(Col_Name=experiments,Group=groups)
#' list_result <- get_statistics_from_dataFrame(df_expr,df_group)
#' top_table(list_result)
#' 
#' ## For this example we will import Gregori data
#' ## Josep Gregori, Laura Villareal, Alex Sanchez, Jose Baselga, Josep Villanueva (2013).
#' ## An Effect Size Filter Improves the Reproducibility 
#' ## in Spectral Counting-based Comparative Proteomics.
#' ## Journal of Proteomics, DOI http://dx.doi.org/10.1016/j.jprot.2013.05.030')
#' 
#' ## Description:
#' ## Each sample consists in 500ng of standard yeast lisate spiked with 
#' ## 100, 200, 400 and 600fm of a mix of 48 equimolar human proteins (UPS1, Sigma-Aldrich).
#' ## The dataset contains a different number of technical replimessagees of each sample
#' 
#' ## import Gregori data
#' data(example_data1)
#' df_contrast <- example_data
#' df_group <- example_group
#' 
#' ## Get statistics through 'get_statistics_from_dataFrame' function
#' list_result <- get_statistics_from_dataFrame(df_contrast,df_group)
#' 
#' ## Get significant features (alpha >= 0.05 and power >= 0.90)
#' significant_qpf <- top_table(list_result,pvalue=0.05,power_desired=0.90,method='QPF')
#' 
#' @export
#' 
top_table <- function(input_list,pvalue=0.05,power_desired=0.90,method='QPF',FC_threshold = 2)
{
  input_data <- input_list$data_table
  min_rep <- input_list$min_rep
  max_rep <- input_list$max_rep
  nt <- input_list$nt
  ng <- input_list$ng
  padj <- input_list$method_pvalue_adjustment
  
  # check input_data
  if (!is.data.frame(input_data)) {
    stop('top_table function needs data frame produced by get_statistics_from_file function.')
  }
  if (length(grep("Protein",colnames(input_data))) == 0) {
    stop('Protein column is missing.')
  }
  if (length(grep("NB_Pval_adjusted",colnames(input_data))) == 0) {
    stop('NB_Pval_adjusted column is missing.')
  }
  if (length(grep("Cohens_W",colnames(input_data))) == 0) {
    stop('Cohens_W column is missing.')
  }
  if (length(grep("Cohens_F",colnames(input_data))) == 0) {
    stop('Cohens_F column is missing.')
  }
  if (length(grep("Max_FC",colnames(input_data))) == 0) {
    stop('Max_FC column is missing.')
  }
  if (length(grep("QP_Pval_adjusted",colnames(input_data))) == 0) {
    stop('QP_Pval_adjusted column is missing.')
  }
  if (length(grep("Normal_Pval_adjusted",colnames(input_data))) == 0) {
    stop('Normal_Pval_adjusted column is missing.')
  }
  
  # check pvalue
  if (pvalue > 1) {
    stop('P-value threshould cannot be bigger than 1.')
  } else if (pvalue < 0) {
    stop('P-value threshould should be positive number.')
  }
  
  # check power_desired
  if (power_desired > 1) {
    stop('Statistical power cannot be bigger than 1.')
  } else if (power_desired < 0) {
    stop('Statistical power should be positive number.')
  }
  
  
  df_return <- data.frame()
  min_es <- 0
  method_PVALUE <- ""
  method_ES <- ""
  
  if (method == "QPF") {
    #########################################
    ### Quasi-Poisson model and Cohen's f ###
    #########################################
    tmp <- pwr::pwr.anova.test(k = ng, n = min_rep, sig.level = pvalue, power = power_desired)
    min_es <- as.double(sprintf("%.4f",tmp$f))

    df_analysis <- data.frame(ID=as.character(input_data$Protein), 
                       PVALUE=input_data$QP_Pval_adjusted, ES=input_data$Cohens_F)
    method_PVALUE <- "Adjusted P-value of Quasi-Poisson Model"
    method_ES <- "Cohen's f"
  
    
    }else if (method == "QPF2") {
      ##########################################
      ### Quasi-Poisson model and Cohen's f2 ###
      ##########################################
      # https://people.richland.edu/james/lecture/m170/ch13-1wy.html
      # u = numerator df (# of group-1) and v = denumerator df (# of exp - # of group)
      tmp <- pwr::pwr.f2.test(u = ng-1, v = nt-ng, sig.level = pvalue, power = power_desired)
      min_es <- as.double(sprintf("%.4f",tmp$f2))
      
      df_analysis <- data.frame(ID=as.character(input_data$Protein), 
                                PVALUE=input_data$QP_Pval_adjusted, ES=input_data$Cohens_F2)
      method_PVALUE <- "Adjusted P-value of Quasi-Poisson Model"
      method_ES <- "Cohen's f2"
      
      
    } else if (method == "NBW") {
    ############################################
    ### Negative Binomial modeland Cohen's w ###
    ############################################
    tmp <- pwr::pwr.chisq.test(N = nt, df = ng - 1 , sig.level = pvalue, power = power_desired)
    min_es <- as.double(sprintf("%.4f",tmp$w))
    
    df_analysis <- data.frame(ID=as.character(input_data$Protein), 
                              PVALUE=input_data$NB_Pval_adjusted, ES=input_data$Cohens_W)
    method_PVALUE <- "Adjusted P-value of Negative Binomial Model"
    method_ES <- "Cohen's w"

    
    }else if (method == "NBF2") {
      #################################################
      ### Negative Binomial modeland and Cohen's f2 ###
      #################################################
      # https://people.richland.edu/james/lecture/m170/ch13-1wy.html
      # u = numerator df (# of group-1) and v = denumerator df (# of exp - # of group)
      tmp <- pwr::pwr.f2.test(u = ng-1, v = nt-ng, sig.level = pvalue, power = power_desired)
      min_es <- as.double(sprintf("%.4f",tmp$f2))
      
      df_analysis <- data.frame(ID=as.character(input_data$Protein), 
                                PVALUE=input_data$NB_Pval_adjusted, ES=input_data$Cohens_F2)
      method_PVALUE <- "Adjusted P-value of Negative Binomial Model"
      method_ES <- "Cohen's f2"
      
      
    }else if (method == "LogNF2") {
      stop ("LogNF2 is not provided anymore.")
      #######################################
      ### Log Normal model and Cohen's f2 ###
      #######################################
      # https://people.richland.edu/james/lecture/m170/ch13-1wy.html
      # u = numerator df (# of group-1) and v = denumerator df (# of exp - # of group)
      tmp <- pwr::pwr.f2.test(u = ng-1, v = nt-ng, sig.level = pvalue, power = power_desired)
      min_es <- as.double(sprintf("%.4f",tmp$f2))
      
      df_analysis <- data.frame(ID=as.character(input_data$Protein), 
                                PVALUE=input_data$LogNormal_Pval_adjusted, ES=input_data$Cohens_F2)
      method_PVALUE <- "Adjusted P-value of Log-Normal Model"
      method_ES <- "Cohen's f2"
      
      
    } else if (method == "NORF") {
    ########################################
    ### Normal ANOVA model and Cohen's f ###
    ########################################
    tmp <- pwr::pwr.anova.test(k = ng, n = min_rep, sig.level = pvalue, power = power_desired)
    min_es <- as.double(sprintf("%.4f",tmp$f))
    
    df_analysis <- data.frame(ID=as.character(input_data$Protein), 
                              PVALUE=input_data$Normal_Pval_adjusted, ES=input_data$Cohens_F)
    method_PVALUE <- "Adjusted P-value of Normal ANOVA Model"
    method_ES <- "Cohen's f"


  } else if (method == "QPFC") {
    ###########################################
    ### Quasi-Poisson model and Fold change ###
    ###########################################

    ### Min. power and Cohen's f calculation
    df_tmp <- subset(input_data,input_data$QP_Pval_adjusted < pvalue)
    df_tmp <- subset(df_tmp, df_tmp$Max_FC >= FC_threshold)
    min_f <- min(df_tmp$Cohens_F)
    tmp <- pwr::pwr.anova.test(k = ng, n = min_rep, sig.level = pvalue, f = min_f)
    
    power_desired <- as.double(sprintf("%.2f",tmp$power))
    min_es <- FC_threshold
    
    df_analysis <- data.frame(ID=as.character(input_data$Protein), 
                              PVALUE=input_data$QP_Pval_adjusted, ES=input_data$Max_FC)
    method_PVALUE <- "Adjusted P-value of Quasi-Poisson Model"
    method_ES <- "Fold change"
    
    message ("Since you chose Fold Change, 'power_desired' is ignored.")


  } else if (method == "NBFC") {
    ###############################################
    ### Negative binomial model and Fold change ###
    ###############################################
    
    ### Min. power and Cohen's w calculation
    df_tmp <- subset(input_data,input_data$NB_Pval_adjusted < pvalue)
    df_tmp <- subset(df_tmp, df_tmp$Max_FC >= FC_threshold)
    min_w <- min(df_tmp$Cohens_W)
    tmp <- pwr::pwr.chisq.test(N = nt, df = ng - 1 , sig.level = pvalue, w = min_w)
    
    power_desired <- as.double(sprintf("%.2f",tmp$power))
    min_es <- FC_threshold
    
    df_analysis <- data.frame(ID=as.character(input_data$Protein), 
                              PVALUE=input_data$NB_Pval_adjusted, ES=input_data$Max_FC)
    method_PVALUE <- "Adjusted P-value of negative binomial Model"
    method_ES <- "Fold change"

    message ("Since you chose Fold Change, 'power_desired' is ignored.")

  
  } else if (method == "NORFC") {
    ##########################################
    ### Normal ANOVA model and Fold change ###
    ##########################################
    
    ### Min. power and Cohen's f calculation
    df_tmp <- subset(input_data,input_data$Normal_Pval_adjusted < pvalue)
    df_tmp <- subset(df_tmp, df_tmp$Max_FC >= FC_threshold)
    min_f <- min(df_tmp$Cohens_F)
    tmp <- pwr::pwr.anova.test(k = ng, n = min_rep, sig.level = pvalue, f = min_f)
    
    power_desired <- as.double(sprintf("%.2f",tmp$power))
    min_es <- FC_threshold
    
    df_analysis <- data.frame(ID=as.character(input_data$Protein), 
                              PVALUE=input_data$Normal_Pval_adjusted, ES=input_data$Max_FC)
    method_PVALUE <- "Adjusted P-value of Normal ANOVA Model"
    method_ES <- "Fold change"
    
    message ("Since you chose Fold Change, 'power_desired' is ignored.")
        
  } else {
    stop ("Please choose one between 'QPF', 'QPFC', 'NBW', 'NBFC', 'NORF', 'NORFC'.")
  }
  
  # find significant IDs by filtering by p-value and effect size
  df_significant <- subset(df_analysis,(df_analysis$PVALUE < pvalue & df_analysis$ES >= min_es))

  # find insignificant IDs
  df_insig <- subset(df_analysis,(df_analysis$PVALUE >= pvalue | df_analysis$ES < min_es))
  
  # Draw a scatter plot
  input_data_plot <- rbind(df_significant,df_insig)
  cat <- c(rep('Significant',nrow(df_significant)), rep('Insignificant',nrow(df_insig)))
  df_es <- c(rep(min_es,nrow(df_significant)), rep(min_es,nrow(df_insig)))
  df_pv <- c(rep(pvalue,nrow(df_significant)), rep(pvalue,nrow(df_insig)))
  input_data_plot <- data.frame(x=input_data_plot$ES,
                                y=input_data_plot$PVALUE,
                                Category=cat,
                                min_ES=df_es,
                                max_pvalue=df_pv)
  draw_scatter_plots(input_data_plot,pvalue,min_es,power_desired,method_ES,method_PVALUE)

  # Make a list to return
  colnames(df_significant) <-c("ID",method_PVALUE,method_ES)
  list_return <- list(top_table = df_significant, 
                      minimum_effect_size = min_es,
                      selected_effect_size_filter = method_ES,
                      minimum_power = power_desired,
                      selected_model = method_PVALUE,
                      alpha = pvalue,
                      method_pvalue_adjustment = padj,
                      num_group = ng,
                      common_replicates = min_rep,
                      num_columns = nt)
  return (list_return)
  
}


#' draw_scatter_plots
#' 
#' @import pwr plyr MASS ggplot2
#' @description 
#' Draw a scatterplot to show how significant IDs are distinguished from the total
#' @param
#' input_data_frame
#'    A data frame that consists of 'x' (P-value), 'y' (Effect size),'cat' (significant or not).
#' @param 
#' max_pvalue
#'    P-value threshold
#' @param 
#' power_desired
#'    Give the statistical power you desired for output significant list
#' @param
#' min_ES
#'    Effect size filter threshold
#' @param
#' x_label
#'    Label of X axis
#' @param
#' y_label
#'    Label of y axis
#' @return 
#' A scatter plot
#'
draw_scatter_plots <- function(input_data_frame,max_pvalue,min_ES,power_desired,x_label,y_label)
{
  x <- NULL
  y <- NULL
  Category <- NULL
  p <- ggplot2::ggplot(input_data_frame, aes(x,y,colour=Category))
  p<- p  + theme_bw() + scale_colour_brewer(palette="Set1") + 
    geom_point() + 
    geom_vline(aes(xintercept=min_ES)) + 
    geom_hline(aes(yintercept=max_pvalue)) +
    #annotate("text", x = min_ES, y = max(input_data_frame$y), label = min_ES) +
    #annotate("text", x = max(input_data_frame$x), y = max_pvalue, label = max_pvalue) +
    labs(x = x_label, y = y_label) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          legend.position="none")
  print(p)
  #ggsave(filename="scatter_plot.pdf")
  #message("The scatter plot was generated and saved at ",getwd(),".")
  #message("In the plot, vertical line is effect size threshold and horizontal line is p-value threshold.")
}



#' ttest_cohens_d
#' 
#' @import pwr plyr MASS ggplot2
#' @description 
#' Fulfill Welch Two Sample t-test (\code{\link{t.test}}) and calculate Cohen's d as well as determine significance by p-value and effect size threshold.
#' @param
#' values
#'    A scalar vector.  Length of both of two vectors, values and groups, should be same.
#' @param 
#' groups
#'    Experiment groups for the vector 'values'. Length of both of two vectors, values and groups, should be same.  
#'    The number of groups is not limited to two, such as group <- c('A','A','A','B','B').
#' @param 
#' power
#'    Give the statistical power you desired for output significant list
#' @param
#' alpha
#'    P-value threshold
#' @param
#' alternative
#'    Choose one of these c("two.sided", "less", "greater"). Default is "two.sided".
#' @param
#' paired
#'    if two groups are paired, set it to TRUE.  Default is FALSE.
#' @param
#' var.equal
#'    if two groups are assumed to have same variance, set it to TRUE.  Default is FALSE.
#' @return 
#' A list containing the followings:
#'    \tabular{ll}{
#'    observed_pvalue \tab Calculated P-value from T-test \cr
#'    observed_cohens_d \tab Calculated Cohen's f\cr
#'    threshold_cohens_d \tab Cohen's d threshold at the desired power\cr
#'    threshold_pvalue \tab Desired p-value threshold\cr
#'    flag_pvalue \tab TRUE=passed the pvalue threshold, FALSE=not\cr
#'    flag_cohens_d \tab TRUE=passed the Cohen's d threshold, FALSE=not \cr
#'    power_desired \tab Statistical power in you input parameters \cr
#'    method \tab 'Welch Two Sample t-test' \cr
#'    alternative \tab alternative option in you input parameters \cr
#'    paired \tab paired option in you input parameters \cr
#'    var.equal \tab var.equal option in you input parameters \cr
#'    }
#' @examples
#' library(selfea)
#' 
#' values <- c(8,10,8,8,11,29,26,22,27,26)
#' groups <- c("U200","U200","U200","U200","U200","U600","U600","U600","U600","U600")
#' list_result <- ttest_cohens_d (values, groups, 0.05, 0.90)
#' 
#' @export
#' 
ttest_cohens_d <- function(values, groups, alpha = 0.05, power = 0.90, alternative = 'two.sided', paired = FALSE, var.equal = FALSE) {
  if (length(values) != length(groups)){
    stop('Length of both of two vectors, values and groups, should be same.')
  }
  if (alpha > 1){
    stop('P-value threshold (alpha) should be between 0 and 1')
  }
  if (alpha < 0){
    stop('P-value threshold (alpha) should be between 0 and 1')
  }
  if (power < 0){
    stop('Statistical power (power) should be between 0 and 1')
  }
  if (power > 1){
    stop('Statistical power (power) should be between 0 and 1')
  }
  
  SC <- NULL
  data <- data.frame(SC=values,Run=groups)
  t_test <- t.test(SC~Run, data=data, alternative=alternative, paired=paired, var.equal = var.equal)
  pval <- t_test$p.value
  df_run2mean <- plyr::ddply(data, "Run", plyr::summarize, mean=mean(SC))
  df_run2sd <- plyr::ddply(data, "Run", plyr::summarize, sd=sd(SC))
  df_run2nrow <- plyr::ddply(data, "Run", nrow)
  
  n1 <- as.numeric(df_run2nrow[1,2])
  n2 <- as.numeric(df_run2nrow[2,2])
  sd1 <- as.numeric(df_run2sd[1,2])
  sd2 <- as.numeric(df_run2sd[2,2])
  sd_pooled <- sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
  cohen_d <-abs(df_run2mean[1,2] - df_run2mean[2,2])/sd_pooled
  
  tmp <- pwr.t2n.test(n1 = n1, n2= n2, sig.level = alpha, power = power)
  power_threshold <- tmp$d
  flag_power_threshold <- TRUE
  if (power_threshold > cohen_d) {
    flag_power_threshold <- FALSE
  }
  flag_alpha <- TRUE
  if (pval > alpha) {
    flag_alpha <- FALSE
  }
  list_return <- list(observed_pvalue = pval,
                      observed_cohens_d = cohen_d,
                      threshold_cohens_d=power_threshold,
                      threshold_pvalue=alpha,
                      flag_pvalue=flag_alpha,
                      flag_cohens_d=flag_power_threshold,
                      power_desired=power,
                      method='Welch Two Sample t-test',
                      alternative=alternative,
                      paired=paired,
                      var.equal=var.equal)
}

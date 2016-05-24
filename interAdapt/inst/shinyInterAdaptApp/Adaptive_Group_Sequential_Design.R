## This file computes power, expected sample size, and expected trial duration for the adaptive design and standard designs, for the interAdapt software tool.

## R Library for computing multivariate normal distribution function (for use in setting efficacy boundaries)
library(mvtnorm)
## R Library for constructing tables
library(xtable)

## User Controlled Parameters that are set with Sliders (with initial values and allowed ranges of each variable)

## Subpopulation 1 proportion (Range: 0 to 1)
p1_user_defined <- 0.33

## Note: throughout, we denote the treatment arm by A=1 and control arm by A=0. 
## We represent p_{1c} (probability of successful outcome Y=1 under assignment to control, for subpopulation s) by ps0, for each s = 1, 2.
## Similarly, we represent p_{1t} (probability of successful outcome Y=1 under assignment to treatment, for subpopulation s) by ps1, for each s = 1, 2.

## Probability outcome = 1 under control:
## for Subpopulation 1 (Range: 0 to 1)
p10_user_defined <- 0.25
## for Subpopulation 2 (Range: 0 to 1)
p20_user_defined <- 0.20

## Prob. outcome = 1 under treatment, at alternative:
## for Subpopulation 1 (Range: 0 to 1)
p11_user_defined<- 0.25 + 0.125 
## The user does not input the corresponding value for Subpopulation 2, since this quantity is varied (it is the quantity on the horizontal axis in the plots in the Performance section of the user interface)

## Alpha allocation
# Desired familywise type I error rate (one-sided) (Range: 0 to 1)
alpha_FWER_user_defined <- 0.025
# Proportion of alpha_FWER_user_defined to test of H0C (Range: 0 to 1)
alpha_H0C_proportion_user_defined <- 0.09


## Adaptive Design Per-stage Sample Sizes
per_stage_sample_size_combined_adaptive_design_user_defined <- 150 #(Range: 0 to 1000)
per_stage_sample_size_when_only_subpop_1_enrolled_adaptive_design_user_defined  <- 311 #(Range: 0 to 1000)

## End of list of User Controlled Parameters that are set with Sliders

## Additional Parameters (expected to be less frequently changed than parameters set by sliders) to be input using textboxes under PARAMETERS Tab

## Parameters used in all designs (adaptive and standard)
# Group Sequential Boundary Parameter (exponent) (Range: -0.5 to 0.5)
Delta <- (-0.5)
# Number simulated trials used to compute power, expected sample size, and expected trial duration for plots and tables
iter <- 7000  # Range: 1 to 500,000
# Time limit for primary function 
time_limit<-45 #(seconds, range from 5 to 60)
# Number stages in trial design
total_number_stages <- 5 # Range 1:20
# Enrollment rate for combined population (patients per year)
enrollment_rate_combined_population <- 420
# Delay from enrollment to primary outcome observed in years
delay_from_enrollment_to_primary_outcome <- 1/2 
# Horizontal axis range in plots, in terms of mean treatment effect on risk difference scale:
lower_bound_treatment_effect_subpop_2 <- (-0.2)  # range (-1,1)
upper_bound_treatment_effect_subpop_2 <- (0.2)  # range (-1,1)
## Parameters used only by adaptive design
last_stage_subpop_2_enrolled_adaptive_design <- 3  #(Range: 1 to total_number_stages)
# Stopping boundary proportionality constant for subpopulation 2 in adaptive design (z-statistic scale)
subpopulation_2_stopping_boundary_proportionality_constant_adaptive_design <- 0 #(Range: -10 to 10)
# Futility boundary proportionality constant for H01 for adaptive design (z-statistic scale)
H01_futility_boundary_proportionality_constant_adaptive_design <- 0 #(Range: -10 to 10)
## Parameters used only by standard designs
# Futility boundary proportionality constant for standard design enrolling combined population (z-statistic scale)
H0C_futility_boundary_proportionality_constant_standard_design <- -0.1 #(Range: -10 to 10)
# Per stage sample size for standard design enrolling combined population
per_stage_sample_size_combined_standard_design_H0C <- 90 #(Range: 0 to 1000)
# Futiltiy boundary proportionality constant for standard design enrolling subpopulation 1 only
H01_futility_boundary_proportionality_constant_standard_design <- -0.1 ##(Range: -10 to 10)
# Per stage sample size for standard design enrolling on subpopulation 1 
per_stage_sample_size_combined_standard_design_H01 <- 106 #(Range: 0 to 1000)

## List of global variables not directly accessible by user; these are functions of above variables, but are updated every time table_constructor function is called
H0C_efficacy_boundary_proportionality_constant_standard_design <- 2.04
H01_efficacy_boundary_proportionality_constant_standard_design <- 2.04
H0C_efficacy_boundary_proportionality_constant_adaptive_design <- 2.54
H01_efficacy_boundary_proportionality_constant_adaptive_design <- 2.12
subpop_1_efficacy_boundaries_adaptive_design <- H01_efficacy_boundary_proportionality_constant_adaptive_design*((1:total_number_stages)/total_number_stages)^Delta
# List of treatment effect values (on risk difference scale) at which power, expected sample size, and expected duration will be evaluated 
risk_difference_list <- seq(lower_bound_treatment_effect_subpop_2,upper_bound_treatment_effect_subpop_2,length=7)

# Construct stopping boundaries for adaptive design based on user-input proportionality constants:
subpopulation_2_stopping_boundaries_adaptive_design <- c(subpopulation_2_stopping_boundary_proportionality_constant_adaptive_design*(((1:(last_stage_subpop_2_enrolled_adaptive_design-1))/(last_stage_subpop_2_enrolled_adaptive_design-1))^Delta),rep(Inf,total_number_stages-last_stage_subpop_2_enrolled_adaptive_design+1))
# Compute subpopulation 1 cumulative sample size vector for adaptive design
if(total_number_stages>last_stage_subpop_2_enrolled_adaptive_design){
	subpop_1_sample_size_vector <- c((1:last_stage_subpop_2_enrolled_adaptive_design)*per_stage_sample_size_combined_adaptive_design_user_defined*p1_user_defined,(last_stage_subpop_2_enrolled_adaptive_design*per_stage_sample_size_combined_adaptive_design_user_defined*p1_user_defined) + (1:(total_number_stages-last_stage_subpop_2_enrolled_adaptive_design))*per_stage_sample_size_when_only_subpop_1_enrolled_adaptive_design_user_defined)
} else {
	subpop_1_sample_size_vector <- c((1:last_stage_subpop_2_enrolled_adaptive_design)*per_stage_sample_size_combined_adaptive_design_user_defined*p1_user_defined)
}
subpop_1_futility_boundaries_adaptive_design <- c(H01_futility_boundary_proportionality_constant_adaptive_design*(subpop_1_sample_size_vector[1:total_number_stages-1]/subpop_1_sample_size_vector[total_number_stages-1])^Delta,H01_efficacy_boundary_proportionality_constant_adaptive_design)

# Construct stopping boundaries for standard designs based on user-input proportionality constants:
combined_pop_futility_boundaries_standard_design_H0C <-c(H0C_futility_boundary_proportionality_constant_standard_design*(((1:(total_number_stages-1))/(total_number_stages-1))^Delta),Inf)
subpop_1_futility_boundaries_standard_design_H0C <-c(rep(Inf,total_number_stages-1),Inf)

combined_pop_futility_boundaries_standard_design_H01 <-c(H01_futility_boundary_proportionality_constant_standard_design*(((1:(total_number_stages-1))/(total_number_stages-1))^Delta),Inf)
subpop_1_futility_boundaries_standard_design_H01 <-c(rep(Inf,total_number_stages-1),Inf)

# Function to compute power for given design and data generating distribution.
# The data generating distribution is defined by the signal to noise ratio (SNR) for subpopulation 1 patients: SNR_subpop_1
# and the signal to noise ratio (SNR) for subpopulation 2 patients: SNR_subpop_2. 
# SNR for a given population s (s=1 or 2) is defined as the mean difference (p_{st}-p_{sc}) divided by the corresponding standard deviation. Formulas for these quantities are given below.
# The inputs to the function are:
#  design_type: may be "adaptive" or "standard"
#  p1: proportion of population in subpopulation 1 (denoted pi_1 in the documentation)
#  total_number_stages: total number of stages in the trial (denoted K in the documentation)
#  standard_design_futility_boundaries: only used for standard designs
#  subpop_1_futility_boundaries: vector (l_{1,1},...l_{1,K}), where each l_{1,k} is the subpopulation 1 futility boundary defined in the documentation
#  subpop_2_futility_boundaries: vector (l_{2,1},...l_{2,K}), where each l_{2,k} is the subpopulation 2 futility boundary defined in the documentation
#  n1: per-stage sample size for combined population during each stage at or before last_stage_subpop_2_enrolled_adaptive_design (i.e., for each stage k <= k*)
#  n2: per-stage sample size for subpopulation 1 during each stage from last_stage_subpop_2_enrolled_adaptive_design through last stage (i.e., for each stage k > k*)
#  e_AD_C: the proportionality constant for efficacy boundaries for the combined population in the adaptive design, denoted e_{AD,C} in documentation
#  e_AD_1: the proportionality constant for efficacy boundaries for subpopulation 1 in the adaptive design, denoted e_{AD,1} in documentation
#  In the following, r1 and r2 are the probabilities of being randomized to treatment, for subpopulation 1 and subpopulation 2, respectively; these are set throughout to equal 1/2, so as to achieve 1:1 randomization to treatment vs. control. 
#  SNR_subpop_1: (p11-p10)/sqrt(p11*(1-p11)/r1+p10*(1-p10)/(1-r1)) 
#  SNR_subpop_2:	(p21-p20)/sqrt(p21*(1-p21)/r2+p20*(1-p20)/(1-r2)) 
#  outcome_variance_subpop_1: p11*(1-p11)/r1+p10*(1-p10)/(1-r1)
#  outcome_variance_subpop_2: p21*(1-p21)/r2+p20*(1-p20)/(1-r2)

get_power <- function(design_type="adaptive",p1,total_number_stages=5,k,standard_design_futility_boundaries=rep(-Inf,total_number_stages),subpop_1_futility_boundaries=rep(-Inf,total_number_stages),subpop_2_futility_boundaries=rep(-Inf,total_number_stages),n1,n2,e_AD_C,e_AD_1,SNR_subpop_1,SNR_subpop_2,outcome_variance_subpop_1,outcome_variance_subpop_2){
	p2 <- (1-p1) # proportion of population in subpopulation 2
	# Enrollment rate subpop. 1 (patients per year)
	enrollment_rate_subpop_1 <- p1_user_defined*enrollment_rate_combined_population
	# Enrollment rate subpop. 2 (patients per year)
	enrollment_rate_subpop_2 <- (1-p1_user_defined)*enrollment_rate_combined_population 

	# These are sample sizes of number who have final outcomes available at each interim analysis
	
	combined_population_stagewise_sample_sizes <- c(rep(n1,k),rep(n2,total_number_stages-k))
	subpop_1_stagewise_sample_sizes <- c(rep(p1*n1,k),rep(n2,total_number_stages-k))
	subpop_2_stagewise_sample_sizes <- combined_population_stagewise_sample_sizes - subpop_1_stagewise_sample_sizes

	# Build up vectors of cumulative sample sizes for combined population and each subpopulation, respectively, based on stagewise sample sizes 
	cum_sample_sizes_combined_population <- rep(0,total_number_stages)
	cum_sample_sizes_subpop_1 <- rep(0,total_number_stages)
	cum_sample_sizes_subpop_2 <- rep(0,total_number_stages)
	for(stage in 1:total_number_stages)
	{
		if(stage==1){cum_sample_sizes_subpop_1[stage] <- subpop_1_stagewise_sample_sizes[1]} else {cum_sample_sizes_subpop_1[stage] <- cum_sample_sizes_subpop_1[stage-1] + subpop_1_stagewise_sample_sizes[stage]}
		if(stage==1){cum_sample_sizes_subpop_2[stage] <- subpop_2_stagewise_sample_sizes[1]} else {cum_sample_sizes_subpop_2[stage] <- cum_sample_sizes_subpop_2[stage-1] + subpop_2_stagewise_sample_sizes[stage]}
	}
	cum_sample_sizes_combined_population <- cum_sample_sizes_subpop_2 + cum_sample_sizes_subpop_1
 
	combined_population_efficacy_boundaries <- e_AD_C*c((cum_sample_sizes_combined_population[1:k]/cum_sample_sizes_combined_population[k])^Delta,rep(Inf,total_number_stages-k))
	subpop_1_efficacy_boundaries <- e_AD_1*(cum_sample_sizes_subpop_1/cum_sample_sizes_subpop_1[total_number_stages])^Delta
        
## Get list of sample sizes corresponding to each interim analysis
	all_relevant_subpop_1_sample_sizes <- sort(unique(c(cum_sample_sizes_subpop_1)))
	all_relevant_subpop_2_sample_sizes <- sort(unique(c(cum_sample_sizes_subpop_2)))

## generate z-statistic increments
	Z_subpop_1_increment <- array(0,c(length(all_relevant_subpop_1_sample_sizes),iter)) 
	Z_subpop_1_increment[1,] <- rnorm(iter)+SNR_subpop_1*sqrt(all_relevant_subpop_1_sample_sizes[1])
	if(length(all_relevant_subpop_1_sample_sizes)>1)
	{	for(i in 2:length(all_relevant_subpop_1_sample_sizes))
		{
			Z_subpop_1_increment[i,] <- rnorm(iter)+SNR_subpop_1*sqrt(all_relevant_subpop_1_sample_sizes[i]-all_relevant_subpop_1_sample_sizes[i-1])
		}
	}
	Z_subpop_2_increment <- array(0,c(length(all_relevant_subpop_2_sample_sizes),iter)) 
	Z_subpop_2_increment[1,] <- rnorm(iter)+SNR_subpop_2*sqrt(all_relevant_subpop_2_sample_sizes[1])
	if(length(all_relevant_subpop_2_sample_sizes)>1)
	{
		for(i in 2:length(all_relevant_subpop_2_sample_sizes))
		{
			Z_subpop_2_increment[i,] <- rnorm(iter)+SNR_subpop_2*sqrt(all_relevant_subpop_2_sample_sizes[i]-all_relevant_subpop_2_sample_sizes[i-1])
		}
	}
	
## generate partial sums of increments
## Construct cumulative z-statistics:
	# First for subpop_1 
	Z_subpop_1_partial_weighted_sum_of_increments <- Z_subpop_1_increment
	if(length(all_relevant_subpop_1_sample_sizes)>1)
	{
		for(i in 2:length(all_relevant_subpop_1_sample_sizes))
		{
			Z_subpop_1_partial_weighted_sum_of_increments[i,] <- 
		((sqrt(all_relevant_subpop_1_sample_sizes[i-1]/all_relevant_subpop_1_sample_sizes[i])*Z_subpop_1_partial_weighted_sum_of_increments[i-1,])		
			+ (sqrt((all_relevant_subpop_1_sample_sizes[i]-all_relevant_subpop_1_sample_sizes[i-1])/all_relevant_subpop_1_sample_sizes[i])*Z_subpop_1_increment[i,]))
		}
	}
	Z_subpop_1_cumulative <- array(0,c(total_number_stages,iter))
	for(i in 1:total_number_stages){
		index <- which(all_relevant_subpop_1_sample_sizes==cum_sample_sizes_subpop_1[i])
		Z_subpop_1_cumulative[i,] <- Z_subpop_1_partial_weighted_sum_of_increments[index,]
	}
	# For subpopulation 2
	Z_subpop_2_partial_weighted_sum_of_increments <- Z_subpop_2_increment
	if(length(all_relevant_subpop_2_sample_sizes)>1)
	{
		for(i in 2:length(all_relevant_subpop_2_sample_sizes))
		{
			Z_subpop_2_partial_weighted_sum_of_increments[i,] <- 
		((sqrt(all_relevant_subpop_2_sample_sizes[i-1]/all_relevant_subpop_2_sample_sizes[i])*Z_subpop_2_partial_weighted_sum_of_increments[i-1,])		
			+ (sqrt((all_relevant_subpop_2_sample_sizes[i]-all_relevant_subpop_2_sample_sizes[i-1])/all_relevant_subpop_2_sample_sizes[i])*Z_subpop_2_increment[i,]))
		}
	}
	Z_subpop_2_cumulative <- array(0,c(total_number_stages,iter))
	for(i in 1:total_number_stages){
		index <- which(all_relevant_subpop_2_sample_sizes==cum_sample_sizes_subpop_2[i])
		Z_subpop_2_cumulative[i,] <- Z_subpop_2_partial_weighted_sum_of_increments[index,]
	}
	# Define combined_population population z-statistics
	variance_component1 <- (p1^2)*outcome_variance_subpop_1/cum_sample_sizes_subpop_1
	if(p2!=0){variance_component2 <- (p2^2)*outcome_variance_subpop_2/cum_sample_sizes_subpop_2}else{variance_component2 <- 0*variance_component1}
	correlation_Z_subpop_1_with_Z_combined_population <- sqrt(variance_component1/(variance_component1+variance_component2)) 
	correlation_Z_subpop_2_with_Z_combined_population <- sqrt(variance_component2/(variance_component1+variance_component2))
	Z_combined_population_cumulative <- (correlation_Z_subpop_1_with_Z_combined_population*Z_subpop_1_cumulative + correlation_Z_subpop_2_with_Z_combined_population*Z_subpop_2_cumulative)
	
## Determine outcomes of each simulated trial
if(design_type=="adaptive"){
    # record if efficacy boundary ever crossed, for each of H0C and H01:
 	ever_cross_H0C_efficacy_boundary <- rep(0,iter)
	ever_cross_H01_efficacy_boundary <- rep(0,iter)
	# indicator of stopping all enrollment, and of stopping only subpopulation 2, respectively:
	all_stopped <- rep(0,iter)
    subpop_2_stopped <- rep(0,iter)
    # indicators of rejecting null hypotheses:
	reject_H0C <- rep(0,iter)
    reject_H01 <- rep(0,iter)
    # record stage (just) after which trial stops
	final_stage_subpop_1_enrolled_up_through <- rep(total_number_stages,iter)
    final_stage_subpop_2_enrolled_up_through <- rep(total_number_stages,iter)
	for(stage in 1:total_number_stages)
	{
		  #below, k represents k^* from paper
          if(stage <= k){ever_cross_H0C_efficacy_boundary <- ifelse(Z_combined_population_cumulative[stage,]>combined_population_efficacy_boundaries[stage],1,ever_cross_H0C_efficacy_boundary)} # since always stop H0C testing after stage k
          ever_cross_H01_efficacy_boundary <- ifelse(Z_subpop_1_cumulative[stage,]>subpop_1_efficacy_boundaries[stage],1,ever_cross_H01_efficacy_boundary)
		# Step 1 of algorithm: Determine if any new events where H0C rejected for efficacy:
         if(stage <= k){reject_H0C <- ifelse((!all_stopped) & (!subpop_2_stopped) & Z_combined_population_cumulative[stage,]>combined_population_efficacy_boundaries[stage],1,reject_H0C)}

         reject_H01 <- ifelse((!all_stopped) & Z_subpop_1_cumulative[stage,]>subpop_1_efficacy_boundaries[stage],1,reject_H01)
          
         all_stopped <- ifelse(reject_H0C | reject_H01 | (Z_subpop_1_cumulative[stage,]< subpop_1_futility_boundaries[stage]),1,all_stopped)          

         subpop_2_stopped <- ifelse(all_stopped | (Z_subpop_2_cumulative[stage,] < subpop_2_futility_boundaries[stage]),1,subpop_2_stopped)

        # force subpop 2 stop at stage k if not yet stopped already
        if(stage==k){subpop_2_stopped <- rep(1,iter)}
		
		# record at what stage each subpop. stopped
        final_stage_subpop_1_enrolled_up_through <- ifelse((final_stage_subpop_1_enrolled_up_through==total_number_stages) & (all_stopped==1),stage,final_stage_subpop_1_enrolled_up_through)
        final_stage_subpop_2_enrolled_up_through <- ifelse((final_stage_subpop_2_enrolled_up_through==total_number_stages) & (subpop_2_stopped==1),stage,final_stage_subpop_2_enrolled_up_through)
	}
return(c(
mean(cum_sample_sizes_subpop_1[final_stage_subpop_1_enrolled_up_through]+cum_sample_sizes_subpop_2[final_stage_subpop_2_enrolled_up_through]), # expected sample size
mean(cum_sample_sizes_subpop_1[final_stage_subpop_1_enrolled_up_through])/(p1*enrollment_rate_combined_population), # expected duration
mean(reject_H0C), # power to reject H0C
mean(reject_H01), # power to reject H01
mean(reject_H0C | reject_H01) # power to reject H01 or H0C
))
} else if(design_type=="standard"){  
    # record if efficacy boundary ever crossed:
 	ever_cross_H0C_efficacy_boundary <- rep(0,iter)
	# indicator of stopping all enrollment, and of stopping only subpopulation 2, respectively:
	all_stopped <- rep(0,iter)
    # indicator of rejecting null hypotheses:
	reject_H0C <- rep(0,iter)
    # record stage (just) after which trial stops
	final_stage_enrolled_up_through <- rep(total_number_stages,iter)
	for(stage in 1:total_number_stages)
	{
		  #below, k represents k^* from paper
          ever_cross_H0C_efficacy_boundary <- ifelse(Z_combined_population_cumulative[stage,]>combined_population_efficacy_boundaries[stage],1,ever_cross_H0C_efficacy_boundary)
		# Step 1 of algorithm: Determine if any new events where H0C rejected for efficacy:
         reject_H0C <- ifelse((!all_stopped) & Z_combined_population_cumulative[stage,]>combined_population_efficacy_boundaries[stage],1,reject_H0C)          
         all_stopped <- ifelse(reject_H0C | (Z_combined_population_cumulative[stage,]< standard_design_futility_boundaries[stage]),1,all_stopped)          		
		# record at what stage each subpop. stopped
        final_stage_enrolled_up_through <- ifelse((final_stage_enrolled_up_through==total_number_stages) & (all_stopped==1),stage,final_stage_enrolled_up_through)
	}
return(c(
mean(cum_sample_sizes_subpop_1[final_stage_enrolled_up_through]+cum_sample_sizes_subpop_2[final_stage_enrolled_up_through]), # expected sample size
mean(reject_H0C) # power to reject H0C
))
}}

#Binary search to find smallest proportionality constants e_AD_C and e_AD_1 such that worst-case familywise Type I error is at most alpha (set by user) for adaptive design
get_adaptive_efficacy_boundaries <- function(alpha_FWER,alpha_H0C,outcome_variance_subpop_1,
outcome_variance_subpop_2)
{
	p1 <- p1_user_defined
	p2 <- 1-p1
	k_star <- last_stage_subpop_2_enrolled_adaptive_design
ss<- rep(per_stage_sample_size_combined_adaptive_design_user_defined,k_star)
if(k_star>1)
{
	for(i in 2:k_star){
	ss[i] <- ss[i-1]+per_stage_sample_size_combined_adaptive_design_user_defined
	}
}
	cov_matrix <- diag(k_star)
for(i in 1:k_star){for(j in 1:k_star) cov_matrix[i,j] <- sqrt(min(ss[i],ss[j])/max(ss[i],ss[j]))}

boundary_vector_with_unit_proportionality_constant <- ((1:k_star)/k_star)^Delta

OF_prop_constant_upper_bnd <- 10
OF_prop_constant_lower_bnd <- 0
while(OF_prop_constant_upper_bnd-OF_prop_constant_lower_bnd > 0.000001)
{
	OF_prop_constant_midpt <- mean(c(OF_prop_constant_lower_bnd,OF_prop_constant_upper_bnd))
	type_I_error <- 1-(pmvnorm(lower=rep(-Inf,k_star),upper=OF_prop_constant_midpt*boundary_vector_with_unit_proportionality_constant,mean=rep(0,k_star),sigma=cov_matrix))
	if(type_I_error < alpha_H0C) OF_prop_constant_upper_bnd <- OF_prop_constant_midpt else OF_prop_constant_lower_bnd <- OF_prop_constant_midpt
}
H0C_prop_const <- OF_prop_constant_midpt
OF_final_boundaries_H0C <- OF_prop_constant_midpt*boundary_vector_with_unit_proportionality_constant

# Generate candidate efficacy boundaries for subpopulation 1
# First construct cumulative sample sizes for subpopulation 1
cum_sample_sizes_subpop_1 <- c(p1_user_defined*ss,rep(0,total_number_stages-k_star))
if(total_number_stages-k_star >0)
{
	for(i in (k_star+1):total_number_stages){
	cum_sample_sizes_subpop_1[i] <-  cum_sample_sizes_subpop_1[i-1]+ per_stage_sample_size_when_only_subpop_1_enrolled_adaptive_design_user_defined
	}
}
# Construct covariance matrix of Z_{C,1},dots,Z_{C,k_star},Z_{S,1},dots,Z_{S,total_number_stages}
cov_matrix_full <- cov_matrix_full_check <- diag(k_star+total_number_stages)
# upper left set to cov_matrix for Z_{C,j}:
cov_matrix_full[1:k_star,1:k_star] <- cov_matrix
# build lower right for Z_{S,j} 
for(i in 1:total_number_stages){for(j in 1:total_number_stages) cov_matrix_full[i+k_star,j+k_star] <- sqrt(min(cum_sample_sizes_subpop_1[i],cum_sample_sizes_subpop_1[j])/max(cum_sample_sizes_subpop_1[i],cum_sample_sizes_subpop_1[j]))}

cum_sample_sizes_subpop_2 <- c(ss*(1-p1_user_defined),rep(ss[k_star]*(1-p1_user_defined),total_number_stages-k_star))

# build upper right and lower left parts of covariance matrix Z_{C,i}Z_{S,j}
for(i in 1:k_star){
	for(j in 1:total_number_stages){
	cov_matrix_full[j+k_star,i] <-cov_matrix_full[i,j+k_star] <- sqrt((min(cum_sample_sizes_subpop_1[i],cum_sample_sizes_subpop_1[j])/max(cum_sample_sizes_subpop_1[i],cum_sample_sizes_subpop_1[j]))*(p1*outcome_variance_subpop_1/(p1*outcome_variance_subpop_1+p2*outcome_variance_subpop_2)))}}
# Do binary search to set efficacy threshold for subpopulation 1.
boundary_vector_with_unit_proportionality_constant_subpop_1 <- (cum_sample_sizes_subpop_1/cum_sample_sizes_subpop_1[total_number_stages])^Delta
#((1:total_number_stages)/total_number_stages)^Delta

OF_prop_constant_upper_bnd <- 10
OF_prop_constant_lower_bnd <- 0
while(OF_prop_constant_upper_bnd-OF_prop_constant_lower_bnd > 0.000001)
{
	OF_prop_constant_midpt <- mean(c(OF_prop_constant_lower_bnd,OF_prop_constant_upper_bnd))
	type_I_error <- 1-(pmvnorm(lower=rep(-Inf,k_star+total_number_stages),upper=c(OF_final_boundaries_H0C,OF_prop_constant_midpt*boundary_vector_with_unit_proportionality_constant_subpop_1),mean=rep(0,k_star+total_number_stages),sigma=cov_matrix_full))
	
	if(type_I_error < alpha_FWER) OF_prop_constant_upper_bnd <- OF_prop_constant_midpt else OF_prop_constant_lower_bnd <- OF_prop_constant_midpt
}
H01_prop_const <- OF_prop_constant_midpt
OF_final_boundaries_H01 <- OF_prop_constant_midpt*boundary_vector_with_unit_proportionality_constant_subpop_1
return(c(H0C_prop_const,H01_prop_const,boundary_vector_with_unit_proportionality_constant_subpop_1))
}

##
## Construct table of values to be displayed in user interface
##
table_constructor <- function(){

setTimeLimit(time_limit) # stops computation if taking greater than time_limit

error_counter <- 0 # checks if errors encountered
k <- last_stage_subpop_2_enrolled_adaptive_design
p1 <- p1_user_defined
p2 <- (1-p1)
p11 <- p11_user_defined
p10 <- p10_user_defined
p20 <- p20_user_defined

# Probability randomized to control Arm
# for Subpop. 1 (Range: 0 to 1)
r1 <- 1/2 
# for Subpop. 2 (Range: 0 to 1)
r2 <- 1/2

futility_boundaries_standard_design_H0C <<-c(H0C_futility_boundary_proportionality_constant_standard_design*(((1:(total_number_stages-1))/(total_number_stages-1))^Delta),Inf)

futility_boundaries_standard_design_H01 <<-c(H01_futility_boundary_proportionality_constant_standard_design*(((1:(total_number_stages-1))/(total_number_stages-1))^Delta),Inf)

#Placeholders: these are not used in algorithm
subpop_1_futility_boundaries_standard_design_H0C <- rep(-Inf,total_number_stages)
subpop_1_futility_boundaries_standard_design_H01 <- rep(-Inf,total_number_stages)
subpop_2_futility_boundaries_standard_design_H0C <- rep(-Inf,total_number_stages)
subpop_2_futility_boundaries_standard_design_H01 <- rep(-Inf,total_number_stages)
combined_population_stopping_boundaries_adaptive_design <- rep(-Inf,total_number_stages)

outcome_variance_subpop_1_under_null <- p10*(1-p10)/r1+p10*(1-p10)/(1-r1)
outcome_variance_subpop_2_under_null <- p20*(1-p20)/r2+p20*(1-p20)/(1-r2)
prop_consts <- get_adaptive_efficacy_boundaries(alpha_FWER_user_defined,alpha_H0C_proportion_user_defined*alpha_FWER_user_defined,outcome_variance_subpop_1_under_null,outcome_variance_subpop_2_under_null)
H0C_efficacy_boundary_proportionality_constant_adaptive_design <<- prop_consts[1] 
H01_efficacy_boundary_proportionality_constant_adaptive_design <<- prop_consts[2] 
adaptive_design_boundary_vector_with_unit_proportionality_constant_subpop_1 <- prop_consts[3:(3+total_number_stages-1)]

## Compute efficacy boundary for standard design enrolling combined population and testing only H0C, such that worst-case Type I error is at most alpha
ss <- 1:total_number_stages
cov_matrix <- diag(total_number_stages)
for(i in 1:total_number_stages){for(j in 1:total_number_stages) cov_matrix[i,j] <- sqrt(min(ss[i],ss[j])/max(ss[i],ss[j]))}
boundary_vector_with_unit_proportionality_constant <- ((1:total_number_stages)/total_number_stages)^Delta
OF_prop_constant_upper_bnd <- 10
OF_prop_constant_lower_bnd <- 0
while(OF_prop_constant_upper_bnd-OF_prop_constant_lower_bnd > 0.000001)
{
	OF_prop_constant_midpt <- mean(c(OF_prop_constant_lower_bnd,OF_prop_constant_upper_bnd))
	type_I_error <- 1-(pmvnorm(lower=rep(-Inf,total_number_stages),upper=c(OF_prop_constant_midpt*boundary_vector_with_unit_proportionality_constant),mean=rep(0,total_number_stages),sigma=cov_matrix))

	if(type_I_error < alpha_FWER_user_defined) OF_prop_constant_upper_bnd <- OF_prop_constant_midpt else OF_prop_constant_lower_bnd <- OF_prop_constant_midpt
}
H0C_efficacy_boundary_proportionality_constant_standard_design <<- OF_prop_constant_midpt
H01_efficacy_boundary_proportionality_constant_standard_design <<- OF_prop_constant_midpt

H0C_efficacy_boundaries <- H0C_efficacy_boundary_proportionality_constant_adaptive_design*c(((1:k)/k)^Delta,rep(Inf,total_number_stages-k))

subpop_1_efficacy_boundaries_adaptive_design <<- H01_efficacy_boundary_proportionality_constant_adaptive_design*adaptive_design_boundary_vector_with_unit_proportionality_constant_subpop_1

subpopulation_2_stopping_boundaries_adaptive_design <<- c(subpopulation_2_stopping_boundary_proportionality_constant_adaptive_design*(((1:(last_stage_subpop_2_enrolled_adaptive_design-1))/(last_stage_subpop_2_enrolled_adaptive_design-1))^Delta),rep(Inf,total_number_stages-last_stage_subpop_2_enrolled_adaptive_design+1))
# Compute subpop. 1 cumulative sample size vector
if(total_number_stages>last_stage_subpop_2_enrolled_adaptive_design){
	subpop_1_sample_size_vector <- c((1:last_stage_subpop_2_enrolled_adaptive_design)*per_stage_sample_size_combined_adaptive_design_user_defined*p1_user_defined,(last_stage_subpop_2_enrolled_adaptive_design*per_stage_sample_size_combined_adaptive_design_user_defined*p1_user_defined) + (1:(total_number_stages-last_stage_subpop_2_enrolled_adaptive_design))*per_stage_sample_size_when_only_subpop_1_enrolled_adaptive_design_user_defined)
} else {
	subpop_1_sample_size_vector <- c((1:last_stage_subpop_2_enrolled_adaptive_design)*per_stage_sample_size_combined_adaptive_design_user_defined*p1_user_defined)
}
subpop_1_futility_boundaries_adaptive_design <<- c(H01_futility_boundary_proportionality_constant_adaptive_design*(subpop_1_sample_size_vector[1:total_number_stages-1]/subpop_1_sample_size_vector[total_number_stages-1])^Delta,H01_efficacy_boundary_proportionality_constant_adaptive_design)

print_ss_and_boundaries_flag <- 1
SNR_subpop_1 <- (p11-p10)/sqrt(p11*(1-p11)/r1+p10*(1-p10)/(1-r1))
outcome_variance_subpop_1 <- p11*(1-p11)/r1+p10*(1-p10)/(1-r1)

risk_difference_list <<- sort(unique(c(seq(max(c(min(c(lower_bound_treatment_effect_subpop_2,upper_bound_treatment_effect_subpop_2,0)),-p20)),min(c(max(c(lower_bound_treatment_effect_subpop_2,upper_bound_treatment_effect_subpop_2,0)),1-p20)),length=10))))

standard_combined_population_df <- array(0,c(length(risk_difference_list),4))
standard_subpop_1_only_df <- array(0,c(length(risk_difference_list),3))
adaptive_df <- array(0,c(length(risk_difference_list),5))
overrun_df <- array(0,c(length(risk_difference_list),3))
counter_combined_population <- 1
counter_subpop_1 <- 1
counter_adaptive <- 1
 
for(percent_benefit_subpop_2 in rev(risk_difference_list))
{
	p21 <- p20 + percent_benefit_subpop_2 
	#next two lines ensure p21 is in the range  (.001 , .999)
	p21<-max(p21,.001)
	p21<-min(p21,.999)
	SNR_subpop_2 <- (p21-p20)/sqrt(p21*(1-p21)/r2+p20*(1-p20)/(1-r2)) 
    outcome_variance_subpop_2 <- p21*(1-p21)/r2+p20*(1-p20)/(1-r2)

power_vec <- get_power(design_type="adaptive",p1=p1_user_defined,total_number_stages,k=last_stage_subpop_2_enrolled_adaptive_design,combined_population_stopping_boundaries_adaptive_design,subpop_1_futility_boundaries_adaptive_design,subpopulation_2_stopping_boundaries_adaptive_design,n1=per_stage_sample_size_combined_adaptive_design_user_defined,n2=per_stage_sample_size_when_only_subpop_1_enrolled_adaptive_design_user_defined,e_AD_C=H0C_efficacy_boundary_proportionality_constant_adaptive_design,e_AD_1=H01_efficacy_boundary_proportionality_constant_adaptive_design,SNR_subpop_1=SNR_subpop_1,SNR_subpop_2=SNR_subpop_2,outcome_variance_subpop_1,outcome_variance_subpop_2)

adaptive_df[counter_adaptive,] <-  power_vec # Expected sample size, Expected Duration, Power to Reject H0C, Power to reject H01, Power to reject at least one of H0C or H01
overrun_df[counter_adaptive,1] <- 0
counter_adaptive <- counter_adaptive + 1

power_vec <- get_power(design_type="standard",p1=p1_user_defined,total_number_stages=total_number_stages,k=total_number_stages,futility_boundaries_standard_design_H0C,subpop_1_futility_boundaries=subpop_1_futility_boundaries_standard_design_H0C,subpop_2_futility_boundaries=subpop_2_futility_boundaries_standard_design_H0C,n1=per_stage_sample_size_combined_standard_design_H0C,n2=0,e_AD_C=H0C_efficacy_boundary_proportionality_constant_standard_design,e_AD_1=Inf,SNR_subpop_1=SNR_subpop_1,SNR_subpop_2=SNR_subpop_2,outcome_variance_subpop_1,outcome_variance_subpop_2)
standard_combined_population_df[counter_combined_population,] <- c(power_vec[1],power_vec[1]/enrollment_rate_combined_population,power_vec[2],0)
overrun_df[counter_combined_population,2] <- 0
counter_combined_population <- counter_combined_population +1

power_vec <- get_power(design_type="standard",p1=1,total_number_stages=total_number_stages,k=total_number_stages,futility_boundaries_standard_design_H01,subpop_1_futility_boundaries=subpop_1_futility_boundaries_standard_design_H01,subpop_2_futility_boundaries=subpop_2_futility_boundaries_standard_design_H01,n1=per_stage_sample_size_combined_standard_design_H01,n2=0,e_AD_C=H01_efficacy_boundary_proportionality_constant_standard_design,e_AD_1=Inf,SNR_subpop_1=SNR_subpop_1,SNR_subpop_2=SNR_subpop_2,outcome_variance_subpop_1,outcome_variance_subpop_2)
standard_subpop_1_only_df[counter_subpop_1,] <- c(power_vec[1],power_vec[1]/(p1*enrollment_rate_combined_population),power_vec[2])
overrun_df[counter_subpop_1,3] <- 0

counter_subpop_1 <- counter_subpop_1 +1

}
return(data.frame(cbind(rev(risk_difference_list),adaptive_df,standard_combined_population_df,standard_subpop_1_only_df,overrun_df)))
}

## Construct Power Curve Plot for display in user interface
power_curve_plot <- function()
{
plot(0,type="n",xlim=c(min(risk_difference_list),max(risk_difference_list)),ylim=c(0,1),main="Power versus Average Treatment Effect in Subpopulation 2",xlab="Avg. Treatment Effect on Risk Difference Scale in Subpopulation 2",ylab="Power")

lines(x=rev(risk_difference_list),y=table1[,6],lty=1,col=1,lwd=3)
# H0C adaptive
lines(x=rev(risk_difference_list),y=table1[,4],lty=2,col=1,lwd=3)
# H01 adaptive
lines(x=rev(risk_difference_list),y=table1[,5],lty=3,col=1,lwd=3)

# H0C standard
lines(x=rev(risk_difference_list),y=table1[,9],lty=4,col=3,lwd=3)

# H01 standard
lines(x=rev(risk_difference_list),y=table1[,13],lty=5,col=4,lwd=3)
ltext<-rep(NA,5)
ltext[1]<-expression(paste("Adaptive, Power to Reject at Least one of H"[0][C]," or H"[0][1]))
ltext[2]<-expression(paste("Adaptive, Power to Reject H"[0][C]))
ltext[3]<-expression(paste("Adaptive, Power to Reject H"[0][1]))
ltext[4]<-expression(paste("Standard Design Total Pop., Power to Reject H"[0][C],""))
ltext[5]<-expression(paste("Standard Design Subpop. 1 Only, Power to Reject H"[0][1],""))
legend("bottomright",legend=ltext,lty=c(1,2,3,4,5),col=c(1,1,1,3,4),lwd=c(3,3,3,3,3))
}

## Expected Sample Size Plot
expected_sample_size_plot <- function()
{
min_ess <- 0
max_ess <- max(c(table1[,2],table1[,7],table1[,11]))
plot(0,type="n",xlim=c(min(risk_difference_list),max(risk_difference_list)),ylim=c(min_ess,max_ess),main="Expected Sample Size versus Average Treatment Effect in Subpopulation 2",xlab="Avg. Treatment Effect on Risk Difference Scale in Subpopulation 2",ylab="Expected Sample Size")

# adaptive
lines(x=rev(risk_difference_list),y=table1[,2],lty=1,col=1,lwd=3)

# H0C standard
lines(x=rev(risk_difference_list),y=table1[,7],lty=2,col=3,lwd=3)

# H01 standard
lines(x=rev(risk_difference_list),y=table1[,11],lty=3,col=4,lwd=3)
legend("bottomright",legend=c("Adaptive Design","Standard Design Total Pop.","Standard Design Subpop. 1 Only"),lty=c(1,2,3),col=c(1,3,4),lwd=c(3,3,3))

}

## Expected Duration Plot
expected_duration_plot <- function()
{
min_dur <- 0
max_dur <- max(c(table1[,3],table1[,8],table1[,12]))
plot(0,type="n",xlim=c(min(risk_difference_list),max(risk_difference_list)),ylim=c(min_dur,max_dur),main="Expected Duration versus Average Treatment Effect in Subpopulation 2",xlab="Avg. Treatment Effect on Risk Difference Scale in Subpopulation 2",ylab="Expected Duration in Years")

# adaptive
lines(x=rev(risk_difference_list),y=table1[,3],lty=1,col=1,lwd=3)

# H0C standard
lines(x=rev(risk_difference_list),y=table1[,8],lty=2,col=3,lwd=3)

# H01 standard
lines(x=rev(risk_difference_list),y=table1[,12],lty=3,col=4,lwd=3)
legend("bottomright",legend=c("Adaptive Design","Standard Design Total Pop.","Standard Design Subpop. 1 Only"),lty=c(1,2,3),col=c(1,3,4),lwd=c(3,3,3))

}

## Construct diagram of efficacy and futility boundaries for adaptive design
boundary_adapt_plot <-function()
{
adapt_boundary_mat<- t(adaptive_design_sample_sizes_and_boundaries_table()[[1]][c("H0C Efficacy Boundaries u(C,k) for z-statistics Z(C,k)", "Boundaries l(2,k) for Z(2,k) to Stop Subpop. 2 Enrollment", "H01 Efficacy Boundaries u(1,k) for Z(1,k)", "Boundaries l(1,k) for Z(1,k) to Stop All Enrollment"), ])
adapt_boundary_mat[adapt_boundary_mat==Inf]<-NA
fancyTitle<-expression(atop('Decision Boundaries for Sequential Test of Combined Population','Null Hypothesis ( H'[0][C]~') and Subpopulation 1 Null Hypothesis ( H'[0][1]~')'))
matplot(adapt_boundary_mat,type='o',main=fancyTitle,pch=c(0,1,2,3),col=c('red','red','blue','blue'),lty=2,cex=1.5, xlab='Stage',ylab='Boundaries on Z-score scale')
ltext<-rep(NA,4)
ltext[1]<-expression(paste("H"[0][C]," Efficacy Boundaries"))
ltext[2]<-expression(paste("Boundaries l(2,k) for Z(2,k) to Stop Subpop. 2"))
ltext[3]<-expression(paste('H'[0][1],' Efficacy Boundaries'))
ltext[4]<-expression(paste("Boundaries l(1,k) for Z(1,k) to Stop All Enrollment"))
legend('topright',ltext,pch=c(0,1,2,3),col=c('red','red','blue','blue'),lty=2,cex=1.25) 
}

## Construct diagram of efficacy and futility boundaries for design enrolling combined population
boundary_standard_H0C_plot <-function()
{
H0C_boundary_mat<- t(standard_H0C_design_sample_sizes_and_boundaries_table()[[1]][c("H0C Efficacy Boundaries for z-statistics Z(C,k)", "H0C Futility Boundaries for z-statistics Z(C,k)"), ])
fancyTitle<-expression(atop('Decision Boundaries for Sequential Test of', 'Combined Population Null Hypothesis ( H'[0][C]~')'))
matplot(H0C_boundary_mat,type='o',main=fancyTitle,lty=2,pch=c(0,1),col='red', xlab='Stage',ylab='Boundaries on Z-score scale',cex=1.5)
ltext<-rep(NA,2)
ltext[1]<-expression(paste('H'[0][C],' Efficacy Boundaries'))
ltext[2]<-expression(paste('H'[0][C],' Futility Boundaries'))
legend('topright',ltext,lty=2,pch=c(0,1),col='red',cex=1.25)
}

## Construct diagram of efficacy and futility boundaries for design enrolling only subpopulation 1
boundary_standard_H01_plot <-function()
{
H01_boundary_mat<- t(standard_H01_design_sample_sizes_and_boundaries_table()[[1]][c("H01 Efficacy Boundaries for z-statistics Z(1,k)", "H01 Futility Boundaries for z-statistics Z(1,k)"), ])
fancyTitle<-expression(atop('Decision Boundaries for Sequential Test of', 'Combined Population Null Hypothesis ( H'[0][1]~')'))
matplot(H01_boundary_mat,type='o',main=fancyTitle,lty=2,pch=c(0,1),col='blue', xlab='Stage',ylab='Boundaries on Z-score scale',cex=1.5)
ltext<-rep(NA,2)
ltext[1]<-expression(paste('H'[0][1],' Efficacy Boundaries'))
ltext[2]<-expression(paste('H'[0][1],' Futility Boundaries'))
legend('topright',ltext,lty=2,pch=c(0,1),col='blue',cex=1.25)
}

## Construct table displaying performance of each design, including the following: power, expected sample size, and expected trial duration
performance_table <- function()
{
output_df <- cbind(as.matrix(table1))
output_df_formatted <- cbind(output_df[,1],output_df[,2],output_df[,3],100*output_df[,4],100*output_df[,5],100*output_df[,6],output_df[,7],output_df[,8],100*output_df[,9],
#100*output_df[,10],
output_df[,11],output_df[,12],100*output_df[,13])
colnames(output_df_formatted) <- c("Subpop.2 Tx. Effect","AD:Sample Size","AD:DUR","AD:Power H0C","AD:Power H01","AD:Power H0C or H01","SC:Sample Size","SC:DUR","SC:Power H0C","SS:Sample Size","SS:DUR","SS:Power H01")
return(list(output_df_formatted,digits=c(0,2,0,1,0,0,0,0,1,0,0,1,0),caption=paste0("Comparison of avg sample size, avg duration (DUR), and power (as a percent), for the following designs: the Adaptive Design (AD), the Standard Design Enrolling Combined Population (SC), and the Standard Design Enrolling Subpop. 1 Only (SS). All designs strongly control the familywise Type I error rate at ",alpha_FWER_user_defined,".")))
}

## Format the performance table
transpose_performance_table<-function(ptable){
	#make a matrix of digit values for each cell
	#take off the 1st entry (=0) from the digits vector (for the head col names)
	#and add a zero column for the new row names
	dpt<-dim(ptable[[1]])
	digitMatPre<-matrix(ptable$digits[-1],nrow=dpt[2],ncol=dpt[1],byrow=FALSE)
	digitMat<-cbind(0,digitMatPre) #add a column of zeros for the row names
	orderedTab<-t(ptable[[1]])[,(dpt[1]:1)] #reverse the col order
	outTab<-data.frame(orderedTab)
	#We need to say include.colnames=FALSE in our custom renderTable function
	return(list( outTab, digits=digitMat, caption=ptable$caption ) )
}

## Construct table displaying efficacy and futility boundaries for adaptive design
adaptive_design_sample_sizes_and_boundaries_table <- function()
{
k <- last_stage_subpop_2_enrolled_adaptive_design
p1 <- p1_user_defined
p2 <- 1-p1

H0C_efficacy_boundaries <- H0C_efficacy_boundary_proportionality_constant_adaptive_design*c(((1:k)/k)^Delta,rep(NA,total_number_stages-k))

subpopulation_2_stopping_boundaries_adaptive_design <<- c(subpopulation_2_stopping_boundary_proportionality_constant_adaptive_design*(((1:(last_stage_subpop_2_enrolled_adaptive_design-1))/(last_stage_subpop_2_enrolled_adaptive_design-1))^Delta),Inf,rep(NA,total_number_stages-last_stage_subpop_2_enrolled_adaptive_design))
# Compute subpop. 1 cumulative sample size vector
if(total_number_stages>last_stage_subpop_2_enrolled_adaptive_design){
	subpop_1_sample_size_vector <- c((1:last_stage_subpop_2_enrolled_adaptive_design)*per_stage_sample_size_combined_adaptive_design_user_defined*p1_user_defined,(last_stage_subpop_2_enrolled_adaptive_design*per_stage_sample_size_combined_adaptive_design_user_defined*p1_user_defined) + (1:(total_number_stages-last_stage_subpop_2_enrolled_adaptive_design))*per_stage_sample_size_when_only_subpop_1_enrolled_adaptive_design_user_defined)
} else {
	subpop_1_sample_size_vector <- c((1:last_stage_subpop_2_enrolled_adaptive_design)*per_stage_sample_size_combined_adaptive_design_user_defined*p1_user_defined)
}
subpop_1_futility_boundaries_adaptive_design <<- c(H01_futility_boundary_proportionality_constant_adaptive_design*(subpop_1_sample_size_vector[1:total_number_stages-1]/subpop_1_sample_size_vector[total_number_stages-1])^Delta,subpop_1_efficacy_boundaries_adaptive_design[total_number_stages])

if(k<total_number_stages){
row1 <- c(p1*per_stage_sample_size_combined_adaptive_design_user_defined*(1:k),p1*per_stage_sample_size_combined_adaptive_design_user_defined*k+per_stage_sample_size_when_only_subpop_1_enrolled_adaptive_design_user_defined*(1:(total_number_stages-k)))

row2 <- c(p2*per_stage_sample_size_combined_adaptive_design_user_defined*(1:k),rep(p2*per_stage_sample_size_combined_adaptive_design_user_defined*k,total_number_stages-k))

row3 <- c(per_stage_sample_size_combined_adaptive_design_user_defined*(1:k),per_stage_sample_size_combined_adaptive_design_user_defined*k+per_stage_sample_size_when_only_subpop_1_enrolled_adaptive_design_user_defined*(1:(total_number_stages-k)))
}else{
row1 <- c(p1*per_stage_sample_size_combined_adaptive_design_user_defined*(1:k))

row2 <- c(p2*per_stage_sample_size_combined_adaptive_design_user_defined*(1:k))

row3 <- c(per_stage_sample_size_combined_adaptive_design_user_defined*(1:k))
}

H0C_efficacy <-  H0C_efficacy_boundaries

subpopulation_2_stopping_boundaries_adaptive_design_copy <- subpopulation_2_stopping_boundaries_adaptive_design
H01_efficacy <- subpop_1_efficacy_boundaries_adaptive_design
H01_futility <- subpop_1_futility_boundaries_adaptive_design

output_df <- rbind(row1,row2,row3,H0C_efficacy,subpopulation_2_stopping_boundaries_adaptive_design_copy,H01_efficacy,H01_futility)
row.names(output_df) <- c("Cumulative Sample Size Subpop. 1","Cumulative Sample Size Subpop. 2","Cumulative Sample Size Combined Pop.","H0C Efficacy Boundaries u(C,k) for z-statistics Z(C,k)","Boundaries l(2,k) for Z(2,k) to Stop Subpop. 2 Enrollment","H01 Efficacy Boundaries u(1,k) for Z(1,k)","Boundaries l(1,k) for Z(1,k) to Stop All Enrollment")
colnames(output_df) <- 1:total_number_stages
dig_array <- array(0,c(7,(total_number_stages+1)))

dig_array[4:7,2:(total_number_stages+1)] <- ifelse(is.na(output_df[4:7,]),0,ifelse(output_df[4:7,]==0,0,2))
return(list(output_df,digits=dig_array,caption="Cumulative Sample Sizes and Decision Boundaries for Adaptive Design. Each column corresponds to a stage. All thresholds are given on the z-statistic scale."))

}

## Construct table displaying efficacy and futility boundaries for standard design enrolling combined population
standard_H0C_design_sample_sizes_and_boundaries_table <- function()
{
k<-total_number_stages
p1 <- p1_user_defined
p2 <- 1-p1

H0C_efficacy_boundaries <- H0C_efficacy_boundary_proportionality_constant_standard_design*c(((1:k)/k)^Delta)

futility_boundaries_standard_design_H0C <<-c(H0C_futility_boundary_proportionality_constant_standard_design*(((1:(total_number_stages-1))/(total_number_stages-1))^Delta),H0C_efficacy_boundaries[total_number_stages])

row1 <- c(p1*per_stage_sample_size_combined_standard_design_H0C*(1:k))

row2 <- c(p2*per_stage_sample_size_combined_standard_design_H0C*(1:k))

row3 <- c(per_stage_sample_size_combined_standard_design_H0C*(1:k))

H0C_efficacy <-  H0C_efficacy_boundaries

H0C_futility <- futility_boundaries_standard_design_H0C


output_df <- rbind(row1,row2,row3,H0C_efficacy,H0C_futility)
row.names(output_df) <- c("Cumulative Sample Size Subpop. 1","Cumulative Sample Size Subpop. 2","Cumulative Sample Size Combined Pop.","H0C Efficacy Boundaries for z-statistics Z(C,k)","H0C Futility Boundaries for z-statistics Z(C,k)")
colnames(output_df) <- 1:total_number_stages
dig_array <- array(0,c(5,(total_number_stages+1)))

dig_array[4:5,2:(total_number_stages+1)] <- ifelse(is.na(output_df[4:5,]),0,ifelse(output_df[4:5,]==0,0,2))
return(list(output_df,digits=dig_array,caption="Cumulative Sample Sizes and Decision Boundaries for Standard Design Enrolling Combined Population. Each column corresponds to a stage. All thresholds are given on the z-statistic scale."))


}

## Construct table displaying efficacy and futility boundaries for standard design enrolling subpopulation 1 only
standard_H01_design_sample_sizes_and_boundaries_table <- function()
{
k<-total_number_stages
p1 <- p1_user_defined
p2 <- 1-p1

H01_efficacy_boundaries <- H01_efficacy_boundary_proportionality_constant_standard_design*c(((1:k)/k)^Delta)

futility_boundaries_standard_design_H01 <<-c(H01_futility_boundary_proportionality_constant_standard_design*(((1:(total_number_stages-1))/(total_number_stages-1))^Delta),H01_efficacy_boundaries[total_number_stages])


row1 <- c(per_stage_sample_size_combined_standard_design_H01*(1:k))

H01_efficacy <-  H01_efficacy_boundaries

H01_futility <- futility_boundaries_standard_design_H01

output_df <- rbind(row1,H01_efficacy,H01_futility)
row.names(output_df) <- c("Cumulative Sample Size","H01 Efficacy Boundaries for z-statistics Z(1,k)","H01 Futility Boundaries for z-statistics Z(1,k)")
colnames(output_df) <- 1:total_number_stages
dig_array <- array(0,c(3,(total_number_stages+1)))
#dig_array[4:7,] <- array(2,c(4,(total_number_stages+1)))
dig_array[2:3,2:(total_number_stages+1)] <- ifelse(is.na(output_df[2:3,]),0,ifelse(output_df[2:3,]==0,0,2))
return(list(output_df,digits=dig_array,caption="Cumulative Sample Sizes and Decision Boundaries for Standard Design enrolling only Subpopulation 1. Each column corresponds to a stage. All thresholds are given on the z-statistic scale."))
}

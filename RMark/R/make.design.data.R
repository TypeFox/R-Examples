#' Create design dataframes for MARK model specification
#' 
#' For each type of parameter in the analysis model (e.g, p, Phi, r), this
#' function creates design data based on the parameter index matrix (PIM) in
#' MARK terminology. Design data relate parameters to the sampling and data
#' structures; whereas \code{data} relate to the object(animal) being sampled.
#' Design data relate parameters to time, age, cohort and group structure.  Any
#' variables in the design data can be used in formulas to specify the model in
#' \code{\link{make.mark.model}}.
#' 
#' After processing the data, the next step is to create the design data for
#' building the models which is done with this function.  The design data are
#' different than the capture history data that relates to animals.  The types
#' of parameters and design data are specific to the type of analysis.  For
#' example, consider a CJS analysis that has parameters Phi and p.  If there
#' are 4 occasions, there are 3 cohorts and potentially 6 different Phi and 6
#' different p parameters for a single group.  The format for each parameter
#' information matrix (PIM) in MARK is triangular.  RMark uses the all
#' different formulation for PIMS by default, so the PIMs would be
#' \preformatted{ Phi p 1 2 3 7 8 9 4 5 10 11 6 12 } If you chose
#' pim.type="time" for each parameter in "CJS", then the PIMS are structured as
#' \preformatted{ Phi p 1 2 3 4 5 6 2 3 5 6 3 6 } That structure is only useful
#' if there is only a single release cohort represented by the PIM.  If you
#' choose this option and there is more than one cohort represented by the PIM
#' then it will restrict the possible set of models that can be represented.
#' 
#' Each of these parameters relates to different times, different cohorts (time
#' of initial release) and different ages (at least in terms of time since
#' first capture).  Thus we can think of a data frame for each parameter that
#' might look as follows for Phi for the all different structure:
#' \preformatted{ Index time cohort age 1 1 1 0 2 2 1 1 3 3 1 2 4 2 2 0 5 3 2 1
#' 6 3 3 0 } With this design data, one can envision models that describe Phi
#' in terms of the variables time, cohort and age.  For example a time model
#' would have a design matrix like: \preformatted{ Int T2 T3 1 1 0 0 2 1 1 0 3
#' 1 0 1 4 1 1 0 5 1 0 1 6 1 0 1 } Or a time + cohort model might look like
#' \preformatted{ Int T2 T3 C2 C3 1 1 0 0 0 0 2 1 1 0 0 0 3 1 0 1 0 0 4 1 1 0 1
#' 0 5 1 0 1 1 0 6 1 0 1 0 1 } While you could certainly develop these designs
#' manually within MARK, the power of the R code rests with the internal R
#' function \code{\link{model.matrix}} which can take the design data and
#' create the design matrix from a formula specification such as \code{~time}
#' or \code{~time+cohort} alleviating the need to create the design matrix
#' manually.  While many analyses may only need age, time or cohort, it is
#' quite possible to extend the kind of design data, to include different
#' functions of these variables or add additional variables such as effort.
#' One could consider design data for p as follows: \preformatted{ Index time
#' cohort age effort juvenile 7 1 1 1 10 1 8 2 1 2 5 0 9 3 1 3 15 0 10 2 2 1 5
#' 1 11 3 2 2 15 0 12 3 3 1 15 1 } The added columns represent a time dependent
#' covariate (effort) and an age variable of juvenile/adult. With these design
#' data, it is easy to specify different models such as \code{~time},
#' \code{~effort}, \code{~effort+age} or \code{~effort+juvenile}.
#' 
#' With the simplest call:
#' 
#' \code{ddl=make.design.data(proc.example.data)}
#' 
#' the function creates default design data for each type of parameter in the
#' model as defined by \code{proc.example.data$model}.  If
#' \code{proc.example.data} was created with the call in the first example of
#' \code{\link{process.data}}, the model is "CJS" (the default model) so the
#' return value is a list with 2 data frames: one for Phi and one for p. They
#' can be accessed as \code{ddl$Phi} (or \code{ddl[["Phi"]]}) and \code{ddl$p}
#' (or \code{ddl[["p"]]}) or as \code{ddl[[1]]} and \code{ddl[[2]]}
#' respectively.  Using the former notation is easier and safer because it is
#' independent of the ordering of the parameters in the list. For this example,
#' there are 16 groups and each group has 21 Phi parameters and 21 p
#' parameters. Thus, there are 336 rows (parameters) in the design data frame
#' for both Phi and p and thus a total of 772 parameters.
#' 
#' The default fields in each dataframe are \code{group}, \code{cohort},
#' \code{age}, \code{time}, \code{Cohort}, \code{Age}, and \code{Time}. The
#' first 4 fields are factor variables, whereas \code{Cohort}, \code{Age} and
#' \code{Time} are numeric covariate versions of \code{cohort}, \code{age}, and
#' \code{time} shifted so the first value is always zero. In addition, for
#' closed capture heterogeneity models a factor variable \code{mixture} is
#' included. If \code{groups} were created in the call to
#' \code{\link{process.data}}, then the factor variables used to create the
#' \code{groups} are also included in the design data for each type of
#' parameter.  If one of the grouping variables is an age variable it is named
#' \code{initial.age.class} to recognize explicitly that it represents a static
#' initial age and to avoid naming conflicts with \code{age} and \code{Age}
#' variables which represent dynamic age variables of the age of the animal
#' through time.  Non-age related grouping variables are added to the design
#' data using their names in \code{data}.  For example if
#' \code{proc.example.data} is defined using the first example in
#' \code{process.data}, then the fields \code{sex}, \code{initial.age.class}
#' (in place of \code{age} in this case), and \code{region} are added to the
#' design data in addition to the \code{group} variable that has 16 levels. The
#' levels of the \code{group} variable are created by pasting (concatenating)
#' the values of the grouping factor in order.  For example, M11 is sex=M, age
#' class=1 and region=1.
#' 
#' By default, the factor variables \code{age}, \code{time}, and \code{cohort}
#' are created such that there is a factor level for each unique value.  By
#' specfying values for the argument \code{parameters}, the values of
#' \code{age}, \code{time}, and \code{cohort} can be binned (put into
#' intervals) to reduce the number of levels of each factor variable. The
#' argument \code{parameters} is specified as a list of lists. The first level
#' of the list specifies the parameter type and the second level specifies the
#' variables (\code{age}, \code{time}, or \code{cohort}) that will be binned
#' and the cutpoints (endpoints) for the intervals.  For example, if you
#' expected that survival may change substantially to age 3 (i.e. first 3 years
#' of life) but remain relatively constant beyond then, you could bin the ages
#' for survival as 0,1,2,3-8.  Likewise, as well, you could decide to bin time
#' into 2 intervals for capture probability in which effort and expected
#' capture probability might be constant within each interval.  This could be
#' done with the following call using the argument \code{parameters}:
#' 
#' \preformatted{ddl=make.design.data(proc.example.data,
#' parameters=list(Phi=list(age.bins=c(0,0.5,1,2,8)),
#' p=list(time.bins=c(1980,1983,1986))))}
#' 
#' In the above example, \code{age} is binned for Phi but not for p; likewise
#' \code{time} is binned for p but not for Phi. The bins for age were defined
#' as 0,0.5,1,2,8 because the intervals are closed ("]" - inclusive) on the
#' right by default and open ("(" non-inclusive) on the left, except for the
#' first interval which is closed on the left.  Had we used 0,1,2,8, 0 and 1
#' would have been in a single interval.  Any value less than 1 and greater
#' than 0 could be used in place of 0.5. Alternatively, the same bins could be
#' specified as:
#' 
#' \preformatted{ddl=make.design.data(proc.example.data,
#' parameters=list(Phi=list(age.bins=c(0,1,2,3,8)),
#' p=list(time.bins=c(1980,1984,1986))),right=FALSE)}
#' 
#' To create the design data and maintain flexibility, I recommend creating the
#' default design data and then adding other binned variables with the function
#' \code{\link{add.design.data}}.  The 2 binned variables defined above can be
#' added to the default design data as follows:
#' 
#' \preformatted{ ddl=make.design.data(proc.example.data)
#' ddl=add.design.data(proc.example.data,ddl,parameter="Phi",type="age",
#' bins=c(0,.5,1,2,8),name="agebin")
#' ddl=add.design.data(proc.example.data,ddl,parameter="p",type="time",
#' bins=c(1980,1983,1986),name="timebin") }
#' 
#' Adding the other binned variables to the default design data, allows models
#' based on either time, timebin, or Time for p and age, agebin or Age for Phi.
#' Any number of additional binned variables can be defined as long as they are
#' given unique names. Note that R is case-specific so \code{~Time} specifies a
#' model which is a linear trend over time ((e.g. Phi(T) or p(T) in MARK)
#' whereas \code{~time} creates a model with a different factor level for each
#' \code{time} in the data (e.g. Phi(t) or p(t) in MARK) and \code{~timebin}
#' creates a model with 2 factor levels 1980-1983 and 1984-1986.
#' 
#' 
#' Some circumstances may require direct manipulation of the design data to
#' create the needed variable when simple binning is not sufficient or when the
#' design data is a variable other than one related to \code{time}, \code{age},
#' \code{cohort} or \code{group} (e.g., effort index).  This can be done with
#' any of the vast array of R commands.  For example, consider a situation in
#' which 1983 and 1985 were drought years and you wanted to develop a model in
#' which survival was different in drought and non-drought years.  This could
#' be done with the following commands:
#' 
#' \code{ddl$Phi$drought=0}
#' 
#' \code{ddl$Phi$drought[ddl$phi$time==1983 | ddl$Phi$time==1985]= 1}
#' 
#' The first command creates a variable named drought in the Phi design data
#' and initializes it with 0.  The second command changes the drought variable
#' to 1 for the years 1983 and 1985. The single brackets [] index a data frame,
#' matrix or vector.  In this case \code{ddl$Phi$drought} is a vector and
#' \code{ddl$Phi$time==1983 | ddl$Phi$time==1985} selects the values in which
#' time is equal (==) to 1983 or ("|") 1985.  A simpler example might occur if
#' we want to create a function of one of the continuous variables.  If we
#' wanted to define a model for p that was a function of age and age squared,
#' we could add the age squared variable as:
#' 
#' \code{ddl$p$Agesq=ddl$p$Age^2}
#' 
#' Any of the fields in the design data can be used in formulae for the
#' parameters.  However, it is important to recognize that additional variables
#' you define and add to the design data are specific to a particular type of
#' parameter (e.g., Phi, p, etc). Thus, in the above example, you could not use
#' Agesq in a model for Phi without also adding it to the Phi design data.  As
#' described in \code{\link{make.mark.model}}, there is actually a simpler way
#' to add simple functions of variables to a formula without defining them in
#' the design data.
#' 
#' 
#' The above manipulations are sufficient if there is only one or two variables
#' that need to be added to the design data.  If there are many covariates that
#' are time(occasion)-specific then it may be easier to setup a dataframe with
#' the covariate data and use \code{\link{merge_design.covariates}}.
#' 
#' 
#' The fields that are automatically created in the design data depends on the
#' model.  For example, with models such as "POPAN" or any of the "Pradel"
#' models, the PIM structure is called square which really means that it is a
#' single row and all the rows are the same length for each group.  Thus,
#' rectangular or row may have been a better descriptor.  Regardless, in this
#' case there is no concept of a cohort within the PIM which is equivalent to a
#' row within a triangular PIM for "CJS" type models. Thus, for parameters with
#' "Square" PIMS the cohort (and Cohort) field is not generated.  The cohort
#' field is also not created if \code{pim.type="time"} for "Triangular" PIMS,
#' because that structure has the same structure for each row (cohort) and
#' adding cohort effects would be inappropriate.
#' 
#' 
#' For models with "Square" PIMS or \code{pim.type="time"} for "Triangular"
#' PIMS, it is possible to create a cohort variable by defining the cohort
#' variable as a covariate in the capture history data and using it as a
#' variable for creating groups.  As with all grouping variables, it is added
#' to the design data.  Now the one caution with "Square" PIMS is that they are
#' all the same length.  Imagine representing a triangular PIM with a set of
#' square PIMS with each row being a cohort.  The resulting set of PIMS is now
#' rectangular but the lower portion of the matrix is superfluous because the
#' parameters represent times prior to the entry of the cohort, assuming that
#' the use of cohort is to represent a birth cohort.  This is only problematic
#' for these kinds of models when the structure accomodates age and the concept
#' of a birth cohort.  The solution to the problem is to delete the design data
#' for the superfluous parameters after is is created (see warning below).
#' For example, let us presume that you used cohort with 3 levels 
#' as a grouping variable for a model with "Square" PIMS which has 3 occasions.  
#' Then, the PIM structure would look as follows for Phi: 
#' \preformatted{ Phi 1 2 3 4 5 6 7 8 9 }.  If
#' each row represented a cohort that entered at occasions 1,2,3 then
#' parameters 4,7,8 are superfluous or could be thought of as representing
#' cells that are structural zeros in the model because no observations can be
#' made of those cohorts at those times.
#' 
#' After creating the design data, the unneeded rows can be deleted with R
#' commands or you can use the argument \code{remove.unused=TRUE}.  As an
#' example, a code segment might look as follows if \code{chdata} was defined
#' properly: \preformatted{
#' mydata=process.data(chdata,model="POPAN",groups="cohort")
#' ddl=make.design.data(mydata) ddl$Phi=ddl$Phi[-c(4,7,8),] } If cohort and
#' time were suitably defined an easier solution especially for a larger
#' problem would be \preformatted{
#' ddl$Phi=ddl$Phi[as.numeric(ddl$Phi$time)>=as.numeric(ddl$Phi$cohort),] }
#' Which would only keep parameters in which the time is the same or greater
#' than the cohort.  Note that time and cohort would be factor variables and <
#' and > do not make sense which is the reason for the \code{as.numeric} which
#' translates the factor to a numeric ordering of factors (1,2,...) but not the
#' numeric value of the factor level (e.g., 1987,1998).  Thus, the above
#' assumes that both time and cohort have the same factor levels. The design
#' data is specific to each parameter, so the unneeded parameters need to be
#' deleted from design data of each parameter.
#' 
#' However, all of this is done automatically by setting the argument
#' \code{remove.unused=TRUE}.  It functions differently depending on the type
#' of PIM structure.  For models with "Triangular" PIMS, unused design data are
#' determined based on the lack of a release cohort.  For example, if there
#' were no capture history data that started with 0 and had a 1 in the second
#' position ("01.....") that would mean that there were no releases on occasion
#' 2 and row 2 in the PIM would not be needed so it would be removed from the
#' design data.  If \code{remove.unused=TRUE} the design data are removed for
#' any missig cohorts within each group. For models with "Square" PIMS, cohort
#' structure is defined by a grouping variable.  If there is a field named
#' "cohort" within the design data, then unused design data are defined to
#' occur when time < cohort.  This is particularly useful for age structured
#' models which define birth cohorts.  In that case there will be sampling
#' times prior to the birth of the cohort which are not relevant and should be
#' treated as "structural zeros" and not as a zero based on stochastic events.
#' 
#' If the design data are removed, when the model is constructed with
#' \code{\link{make.mark.model}}, the argument \code{default.fixed} controls
#' what happens with the real parameters defined by the missing design data.
#' If \code{default.fixed=TRUE}, then the real parameters are fixed at values
#' that explain with certainty the observed data (e.g., p=0).  That is
#' necessary for models with "Square" PIMS (eg, POPAN and Pradel models) that
#' include each capture-history position in the probability calculation.  For
#' "Triangular" PIMS with "CJS" type models, the capture(encounter) history
#' probability is only computed for occasions past the first "1", the release.
#' Thus, when a cohort is missing there are no entries and the missing design
#' data are truly superfluous and \code{default.fixed=FALSE} will assign the
#' missing design data to a row in the design matrix which has all 0s.  That
#' row will show as a real parameter of (0.5 for a logit link) but it is not
#' included in any parameter count and does not affect any calculation. The
#' advantage in using this approach is that the R code recognizes these and
#' displays blanks for these missing parameters, so it makes for more readable
#' output when say every other cohort is missing. See
#' \code{\link{make.mark.model}} for more information about deleted design data
#' and what this means to development of model formula.
#' 
#' For design data of "Multistrata" models, additional fields are added to
#' represent strata.  A separate PIM is created for each stratum for each
#' parameter and this is denoted in the design data with the addition of the
#' factor variable \code{stratum} which has a level for each stratum.  In
#' addition, for each stratum a dummy variable is created with the name of the
#' stratum (\code{strata.label})and it has value 1 when the parameter is for
#' that stratum and 0 otherwise.  Using these variables with the interaction
#' operator ":" in formula allows more flexibility in creating model structure
#' for some strata and not others.  All "Multistrata" models contain "Psi"
#' parameters which represent the transitions from a stratum to all other
#' strata.  Thus if there are 3 strata, there are 6 PIMS for the "Psi"
#' parameters to represent transition from A to B, A to C, B to A, B to C, C to
#' A and C to B.  The "Psi" parameters are represented by multimonial logit
#' links and the probability of remaining in the stratum is determined by
#' substraction.  To represent these differences, a factor variable
#' \code{tostratum} is created in addition to \code{stratum}. Likewise, dummy
#' variables are created for each stratum with names created by pasting "to"
#' and the strata label (e.g., toA, toB etc). Some examples of using these
#' variables to create models for "Psi" are given in
#' \code{\link{make.mark.model}}.
#' 
#' \preformatted{ 
#' 
#' ######WARNING######## 
#' Deleting design data for mlogit parameters like Psi in the multistate
#' model can fail if you do things like delete certain transitions.  It is better
#' to add the field fix. It should be assigned the value NA for parameters that
#' are estimated and a fixed real value for those that are fixed. Here is an example
#' with the mstrata data example:
#' 
#' data(mstrata)
#' # deleting design data approach to fix Psi A to B to 0 (DON'T use this approach) 
#' dp=process.data(mstrata,model="Multistrata")
#' ddl=make.design.data(dp)
#' ddl$Psi=ddl$Psi[!(ddl$Psi$stratum=="A" & ddl$Psi$tostratum=="B"),]
#' ddl$Psi
#' summary(mark(dp,ddl,output=FALSE),show.fixed=TRUE)
#' #new approach using fix to set Phi=1 for time 2 (USE this approach)
#' ddl=make.design.data(dp)
#' ddl$Psi$fix=NA
#' ddl$Psi$fix[ddl$Psi$stratum=="A" & ddl$Psi$tostratum=="B"]=0
#' ddl$Psi
#' summary(mark(dp,ddl,output=FALSE),show.fixed=TRUE)
#' }
#' @param data Processed data list; resulting value from process.data
#' @param parameters Optional list containing a list for each type of parameter
#' (list of lists); each parameter list is named with the parameter name (eg
#' Phi); each parameter list can contain vectors named age.bins,time.bins and
#' cohort.bins \tabular{ll}{ \code{subtract.stratum} \tab a vector of strata
#' letters (one for each strata) \cr \tab that specifies the tostratum that is
#' computed by subtraction \cr \tab for mlogit parameters like Psi\cr
#' \code{age.bins} \tab bins for binning ages\cr \code{time.bins} \tab bins for
#' binning times\cr \code{cohort.bins} \tab bins for binning cohorts\cr
#' \code{pim.type} \tab either "all" for all different, "time" for column time
#' structure, or \cr \tab "constant" for all values the same within the PIM\cr}
#' @param remove.unused If TRUE, unused design data are deleted; see details
#' below
#' @param right If TRUE, bin intervals are closed on the right
#' @param common.zero if TRUE, uses a common begin.time to set origin (0) for
#' Time variable defaults to FALSE for legacy reasons but should be set to TRUE
#' for models that share formula like p and c with the Time model
#' @return The function value is a list of data frames. The list contains a
#' data frame for each type of parameter in the model (e.g., Phi and p for
#' CJS).  The names of the list elements are the parameter names (e.g., Phi).
#' The structure of the dataframe depends on the calling arguments and the
#' model & data structure as described in the details above.
#' @author Jeff Laake
#' @export
#' @seealso \code{\link{process.data}},\code{\link{merge_design.covariates}},
#' \code{\link{add.design.data}}, \code{\link{make.mark.model}},
#' \code{\link{run.mark.model}}
#' @keywords utility
#' @examples
#' 
#' 
#' data(example.data)
#' proc.example.data=process.data(example.data)
#' ddl=make.design.data(proc.example.data)
#' ddl=add.design.data(proc.example.data,ddl,parameter="Phi",type="age",
#'   bins=c(0,.5,1,2,8),name="agebin")
#' ddl=add.design.data(proc.example.data,ddl,parameter="p",type="time",
#'   bins=c(1980,1983,1986),name="timebin")
#' 
#' 
#' 
make.design.data <-                                                      
function(data,parameters=list(),remove.unused=FALSE,right=TRUE,common.zero=FALSE)
{
#------------------------------------------------------------------------------------------------------
# make.design.data -  creates a design dataframe that is used to construct the design matrix for mark
#                     in make.mark.model
#
# Arguments:
#
#    data             - data list after using process.data
#    parameters       - list with an element for each parameter
#                       each element is a list with age.bins, time.bins and cohort.bins
#                          age.bins         - bins for grouping ages
#                          time.bins        - bins for grouping times
#                          cohort.bins      - bins for grouping cohorts
#                          pim.type         - type of pim structure "all","time","constant"
#                          subtract.stratum - for each stratum, the one to compute by subtraction (for Psi only)
#                          time.varying     - vector of field names that are time varying for this parameter
#                          fields           - vector of field names to be included in design data that are not time varying
#    remove.unused    - if TRUE, unused design data are removed; for triangular
#                       pims, unused design data are determined based on lack of
#                       ch for a particular row (cohort) of a group;  for square
#                       pims. if there is a cohort field in the design data, then
#                       it excludes any design data in which cohort < time.
#
#    common.zero      - if TRUE, uses a common begin.time to set origin (0) for Time variable
#                      defaults to FALSE for legacy reasons but should be set to TRUE
#                      for models that share formula like p and c with the Time model
#
# Value:
#    full.design.data - list of design data frames for each type of parameter in the model
#
#
# Functions used: setup.parameters, compute.design.data, valid.parameters, setup.model
#
#----------------------------------------------------------------------------------------------------
remove.unused.occasions=function(data,ddl)
{
#
# Check validity of parameter list; stop if not valid
#
		parameter="p"
		if(!valid.parameters(data$model,parameter)) stop()
		parameters=setup.parameters(data$model,parameters=NULL,data$nocc,check=FALSE,
				number.of.groups=dim(data$freq)[2])
		if(!parameters[[parameter]]$type%in%c("Triang","STriang"))stop("\nDoes not work for parameters with non-triangular PIM\n")
		ch=data$data$ch
		if(data$model=="Multistrata")
			ch=gsub("[1-9 a-z A-Z]","1",ch)
#
#    Loop over groups
#
		number.of.groups=dim(data$group.covariates)[1]
		if(is.null(number.of.groups))
		{
			number.of.groups=1
			chsplit=ch
		}
		else
			ch.split=split(ch,data$data$group)
		for(j in 1:number.of.groups)
		{
			chmat=matrix(as.numeric(unlist(strsplit(ch.split[[j]],split=vector(length=0)))),ncol=nchar(ch[1]),byrow=TRUE)
			exclude.occ=(1:dim(chmat)[2])[colSums(chmat)==0]
			if(number.of.groups==1)
				ddl[[parameter]]=ddl[[parameter]][!ddl[[parameter]]$time%in%levels(ddl[[parameter]]$time)[exclude.occ-1],]
			else
			{
				group=levels(ddl[[parameter]]$group)[j]
				ddl[[parameter]]=ddl[[parameter]][(ddl[[parameter]]$group!=group) |
								(ddl[[parameter]]$group==group & !ddl[[parameter]]$time%in%levels(ddl[[parameter]]$time)[exclude.occ-1]),]
			}
		}
		return(ddl)
}
#### start of make.design.data
#
#
#
if(!is.list(data))
	stop("data argument is not a processed data list")
else
	if(!"data"%in%names(data)|!"model"%in%names(data))
		stop("data argument is not a processed data list")
#
# Check validity of parameter list; stop if not valid
#
  if(!valid.parameters(data$model,parameters)) stop()
#
#  Setup model and parameters
#
  par.list=setup.parameters(data$model,check=TRUE)
  parameters=setup.parameters(data$model,parameters,data$nocc,check=FALSE,
          number.of.groups=dim(data$freq)[2])
  parameters=parameters[par.list]
  model.list=setup.model(data$model,data$nocc,data$mixtures)
# If reverse, set remove.unused=TRUE  
  if(data$reverse)
  {
	  remove.unused=TRUE
	  temp=parameters[["Psi"]]$subtract.stratum
	  if(!is.null(temp) && any(temp!=data$strata.labels)) stop("Cannot set subtract.stratum and use reverse\n")
  }
#
# Create a data matrix for the each parameter in the model with age, year and cohort for each index
# This data matrix (design.data) is used below to create the design matrix from the formulas
# If age,cohort or year bins are given, use those.  Otherwise each is treated as a factor 
# wihtout binning.
#
# 10 Jan 06 ; added pim.type argument in call to compute.design.data
#
   full.design.data=vector("list",length=length(parameters))
   pimtypes=vector("list",length=length(parameters))
   anyTriang=FALSE
   anySquare=FALSE
   for(i in 1:length(parameters))
   {
   #
   # For mixtures, multistrata and robust designs set up values for input to
   # compute.design.data
   #
	 limits=NULL
     if(data$mixtures==1)
     {
        if(!is.null(parameters[[i]]$mix)) 
		{
			parameters[[i]]$mix=FALSE
            parameters[[i]]$rows=1
		}
     }
     if(!is.null(parameters[[i]]$bystratum) && parameters[[i]]$bystratum)
     {
        strata.labels=data$strata.labels
        nstrata=data$nstrata
		if(!is.null(parameters[[i]]$subset) && nchar(parameters[[i]]$subset)>0)
		{
			limits=strsplit(parameters[[i]]$subset,"")[[1]]
			if(limits[1]=="0")
			{
				strata.labels=c("0",strata.labels)
				nstrata=nstrata+1
			} 
		#	else
		#	{
		#		strata.labels=strata.labels[as.numeric(limits[1]):nstrata]
		#		nstrata=length(strata.labels)
		#	}
		}
        if(!is.null(parameters[[i]]$tostrata) && parameters[[i]]$tostrata)
        {
           if(!is.null(parameters[[i]]$subtract.stratum))
              subtract.stratum=parameters[[i]]$subtract.stratum
           else
              subtract.stratum=strata.labels
           tostrata=TRUE
        }
        else
        {
			if(data$model%in%c("RDMSOpenMCSeas","RDMSOpenMisClass","RDMSMisClass") & names(parameters)[i]%in%c("pi","Omega"))
			{
			   subtract.stratum=data$strata.labels[nstrata]
		    } else
			{
			   subtract.stratum=NULL
	        }
            tostrata=FALSE
        }
     } 
     else
     {
        subtract.stratum=NULL
        strata.labels=NULL
        nstrata=1
        tostrata=FALSE
     }
     if(!model.list$robust) parameters[[i]]$secondary=FALSE
#
#    Compute design data for this parameter if conditions are valid
#    mod 27 June 2011 -- if data structure (too few occasions) is such that no parameters can be estimated it does not create the design data
     if(is.na(parameters[[i]]$num)||(parameters[[i]]$num+data$nocc)>0)
	 {
		 sub.stratum=0
		 if(!is.null(parameters[[i]]$sub.stratum))sub.stratum=parameters[[i]]$sub.stratum
		 # Special code for parameter Phi0 in RDMSOccRepro model
		 #
		 if(data$model=="RDMSOccRepro" & names(parameters)[i]=="Phi0")
		 {
			 design.data=expand.grid(stratum=data$strata.labels,group=1:ncol(data$freq))      
			 if(ncol(data$freq)>1)
			 {
				 ix=grep("age",names(data$group.covariates))
				 cnames=names(data$group.covariates)
				 if(length(ix)!=0)
					 if(names(data$group.covariates)[ix]=="age")
					 {
						 cnames[ix]="initial.age.class"
						 names(data$group.covariates)=cnames
					 }
				 gc=data.frame(data$group.covariates[design.data$group,]) 
				 names(gc)=cnames 
				 row.names(gc)=NULL
				 design.data=cbind(design.data,gc)	  
			 }
		 } else
			 
             design.data=compute.design.data(data,parameters[[i]]$begin,parameters[[i]]$num,
                         parameters[[i]]$type,parameters[[i]]$mix,parameters[[i]]$rows,
                         parameters[[i]]$pim.type,parameters[[i]]$secondary, nstrata,
                         tostrata,strata.labels,subtract.stratum,common.zero=common.zero,
						 sub.stratum=sub.stratum,limits=limits)
         if(!is.null(parameters[[i]]$mix) && parameters[[i]]$mix)design.data$mixture=as.factor(design.data$mixture)
         if(parameters[[i]]$secondary)
		 {
			 session.labels=data$begin.time+cumsum(c(0,data$time.intervals[data$time.intervals>0]))
			 design.data$session=factor(session.labels[design.data$session])
		 }
         design.data$group=as.factor(design.data$group)
         if(!is.null(data$group.covariates))
            levels(design.data$group)=apply(data$group.covariates,1,paste,collapse="")
         if(!is.null(design.data$cohort))
            if(is.null(parameters[[i]]$cohort.bins))
               design.data$cohort=factor(design.data$cohort,levels=unique(levels(factor(design.data$cohort))))
            else
               design.data$cohort=cut(design.data$cohort,parameters[[i]]$cohort.bins,include.lowest=TRUE,right=right)
         if(!is.null(design.data$age))
         if(is.null(parameters[[i]]$age.bins))
           design.data$age=factor(design.data$age,levels=unique(levels(factor(design.data$age))))
         else
           design.data$age=cut(design.data$age,parameters[[i]]$age.bins,include.lowest=TRUE,right=right)
         if(!is.null(design.data$time))
         # mod 30 Sept 09 to remove unused time factor levels
         if(is.null(parameters[[i]]$time.bins))
            design.data$time=factor(design.data$time,levels=unique(levels(factor(design.data$time))))
         else
            design.data$time=cut(design.data$time,parameters[[i]]$time.bins,include.lowest=TRUE,right=right)
         if(model.list$closed | model.list$robust )
         {
            if(names(parameters)[i]=="p" )
            {
			   if(!is.null(parameters[[i]]$share)) design.data$c=0
               design.data$age=NULL
               design.data$Age=NULL
            }
            if(names(parameters)[i]=="c")
            {
               design.data$c=1
               design.data$age=NULL
               design.data$Age=NULL
            }
            if(names(parameters)[i]=="N" | (names(parameters)[i]=="pi" & !is.null(parameters[[i]]$mix)))
            {
               design.data$age=NULL
               design.data$Age=NULL
               design.data$time=NULL
               design.data$Time=NULL
            }
         }
		 if(data$model=="RDMSOccRepro")
		 {
			 if(names(parameters)[i]=="R")
				 design.data=design.data[order(design.data$tostratum,design.data$stratum),]
			 else
				 if(names(parameters)[i]=="Delta")
				     design.data=design.data[order(design.data$stratum,design.data$tostratum,design.data$session),]
		 }
         full.design.data[[i]]=cbind(par.index=1:nrow(design.data),model.index=1:nrow(design.data),design.data)
		 row.names(full.design.data[[i]])=1:nrow(full.design.data[[i]])
         pimtypes[[i]]=list(pim.type=parameters[[i]]$pim.type)
	     if(!is.null(subtract.stratum))pimtypes[[i]]$subtract.stratum=subtract.stratum
         if(parameters[[i]]$type%in%c("Triang","STriang")&&parameters[[i]]$pim.type=="all")anyTriang=TRUE
         if(parameters[[i]]$type =="Square")anySquare=TRUE
	  }
   }
   names(full.design.data)=names(parameters)
   null.design.data=sapply(full.design.data,is.null)
   parameters=parameters[!null.design.data]
   full.design.data=full.design.data[!null.design.data]
#  add model indices
   prev=0
   for(i in 1:length(full.design.data))
   {
	   full.design.data[[i]]$model.index=full.design.data[[i]]$par.index+prev
	   prev=max(full.design.data[[i]]$model.index)
   }	   
#
#  Remove unused design data
#
   if(remove.unused)
   {
      ch=data$data$ch
      if(data$model=="Multistrata")
        ch=gsub("[A-Z a-z 1-9]","1",ch)
      if(anyTriang)
      {
#
#    Loop over groups
#
         number.of.groups=dim(data$group.covariates)[1]
         if(is.null(number.of.groups))number.of.groups=1
         for(j in 1:number.of.groups)
         {
           remove.cohort=NULL
           for(k in 1:data$nocc)
           {
              if(k>1)
                 first.0=paste(rep("0",k-1),collapse="")
              else
                 first.0=""
              if(number.of.groups==1)
              {
                 if(!any(substr(ch,1,k)==paste(first.0,"1",sep="")))
                    remove.cohort=c(remove.cohort,k)
              }
              else
                 if(!any(substr(ch[data$data$group==j],1,k)==paste(first.0,"1",sep="")))
                    remove.cohort=c(remove.cohort,k)
           }
           for(i in 1:length(parameters))
           {
              if(parameters[[i]]$type %in%c("Triang","STriang")&&parameters[[i]]$pim.type=="all")
              {
                 if(number.of.groups==1)
					 full.design.data[[i]]=full.design.data[[i]][!(full.design.data[[i]]$occ.cohort%in%remove.cohort),]
#				 full.design.data[[i]]=full.design.data[[i]][!(as.numeric(full.design.data[[i]]$cohort)%in%remove.cohort),]
                 else
                 {
#          modified 7 Apr 08 to handle different begin.times between groups
				full.design.data[[i]]=full.design.data[[i]][!(as.numeric(full.design.data[[i]]$group)==j &
									full.design.data[[i]]$occ.cohort%in%remove.cohort),]
#				full.design.data[[i]]=full.design.data[[i]][!(as.numeric(full.design.data[[i]]$group)==j &
#                             as.numeric(factor(full.design.data[[i]]$cohort,levels=unique(full.design.data[[i]]$cohort[as.numeric(full.design.data[[i]]$group)==j ])))%in%remove.cohort),]
#          modified 10 Aug to remove unused levels created in removing cohorts                    
                    full.design.data[[i]]$cohort=factor(full.design.data[[i]]$cohort)
                    full.design.data[[i]]$age=factor(full.design.data[[i]]$age)
                    full.design.data[[i]]$time=factor(full.design.data[[i]]$time)
                 }
              }
           }                             
        }
     }
#    if reverse Multistrata model, remove design data for S,Psi and p for added occasions/intervals
     if(data$reverse)
	 {
		 full.design.data[["S"]]=full.design.data[["S"]][!full.design.data[["S"]]$occ%in%seq(1,data$nocc-1,2),]
		 full.design.data[["p"]]=full.design.data[["p"]][!full.design.data[["p"]]$occ%in%seq(1,data$nocc-1,2),]
		 full.design.data[["Psi"]]=full.design.data[["Psi"]][!full.design.data[["Psi"]]$occ%in%seq(2,data$nocc-1,2),]
	 }
     if(anySquare)
     {
        for(i in 1:length(parameters))
        {
           if(parameters[[i]]$type =="Square"&is.null(parameters[[i]]$leave.unused))
           {
              time=full.design.data[[i]]$time
              cohort=full.design.data[[i]]$cohort
              full.design.data[[i]]=full.design.data[[i]][as.numeric(levels(time)[as.numeric(time)])
                      >= as.numeric(levels(cohort)[as.numeric(cohort)]),]
           }
        }
     }
#    drop any unused factor levels after removing design data
     for(i in 1:length(parameters))full.design.data[[i]]=droplevels(full.design.data[[i]])
   }
#  Delete occ.cohort which is only used to remove unused cohorts if any
   if(data$reverse)
      for(i in 1:length(parameters))
	      full.design.data[[i]]$occ.cohort=NULL
#  make pim type assignments and return results
   pimtypes=pimtypes[!null.design.data]
   names(pimtypes)=names(parameters)
   full.design.data$pimtypes=pimtypes
   return(full.design.data)
}

        
         

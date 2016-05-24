
#' Fibtele
#' 
#' 
#' 
#' @docType data
#' @name fibtele
#' @usage fibtele
#' @format A data frame with 147 observations on the following 35
#' variables. The first ten variables are segmentation variables. 
#' The rest of the variables refer to five latent concepts: 1) \code{Image}=Image, 2) 
#'   \code{Qual.spec}=Specific Quality, 3) \code{Qual.gen}=Generic Quality, 4) 
#' \code{Value}=Value, 5) \code{Satis}=Satisfaction. 
#'  Variables description
#' \itemize{
#' \item{\code{Image}}: Generic students perception of ICT schools: (internationally recognized, 
#'      ranges of courses, leader in research).
#' \item{\code{Qual.spec}}: Perception about the achieved quality on the specific skills in the school.
#' \item{\code{Qual.gen}}: Perception about achieved quality on the  generic skills in 
#' the school (abilities in solving problem, communication skills).
#' \item{\code{Value}}: The advantage or profit that the alumni may draw from the school 
#' degree (well paid job, motivated job, prospectives in improvement and promotion).
#' \item{\code{Satis}}: Degree of alumni satisfaction about the formation in school respect to 
#' their actual work conditions.
#' }
#'
#' Manifest variables description
#'
#'\itemize{
#'\item{\code{ima1}}{MV:It is the best college to study IE}
#'\item{\code{ima2}}{MV:It is internationally recognized}
#'\item{\code{ima3}}{MV:It has a wide range of courses}
#'\item{\code{ima4}}{MV:The Professors are good}
#'\item{\code{ima5}}{MV:Facilities and equipment are good}
#'\item{\code{ima6}}{MV:It is leader in research}
#'\item{\code{ima7}}{MV:It is well regarded by the companies}
#'\item{\code{ima8}}{MV:It is oriented to new needs and technologies}
#'\item{\code{quaf1}}{MV:Basic skills}
#'\item{\code{quaf2}}{MV:Specific Technic skills}
#'\item{\code{quaf3}}{MV:Applied skills}
#'\item{\code{qutr1}}{MV:Achieved abilities in solving problem}
#'\item{\code{qutr2}}{MV:Training in business management}
#'\item{\code{qutr3}}{MV:The written and oral communication skills}
#'\item{\code{qutr4}}{MV:Planning and time management acquired}
#'\item{\code{qutr5}}{MV:Team-work skills}
#'\item{\code{val1}}{MV:It has allowed me to find a well paid job}
#'\item{\code{val2}}{MV:I have good prospectives in improvement and promotion}
#'\item{\code{val3}}{MV:It has allowed me to find a job that motivates me}
#'\item{\code{val4}}{MV:The training received is the basis on which I will develope my career}
#'\item{\code{sat1}}{MV:I am satisfied with the training received}
#'\item{\code{sat2}}{MV:I am satisfied with my current situation}
#'\item{\code{sat3}}{MV:I think I will have a good career}
#'\item{\code{sat4}}{MV:What do you think is the prestige of your work}
#'}
#' 
#' Segmentation Variables description 
#'\itemize{
#'\item{\code{Career}}{a factor with levels \code{EI} \code{ETS} \code{TEL}}
#'\item{\code{Gender}}{a factor with levels \code{female} \code{male}}      
#'\item{\code{Age}}{a factor with levels \code{25-26years} \code{27-28years} \code{29-30years} \code{31years+}}   
#'\item{\code{Studying}}{a factor with levels \code{no.stud} \code{yes.stud}}    
#'\item{\code{Contract}} {a factor with levels \code{fix.cont} \code{other.cont} \code{temp.cont}} 
#'\item{\code{Salary}}{a factor with levels \code{18k} \code{>45k} \code{25k} \code{35k} \code{45k}}  
#'\item{\code{Firmtype}}{a factor with levels \code{priva} \code{publi}}       
#'\item{\code{Accgrade}}{a factor with levels \code{7-8accnote} \code{ accnote<7} \code{accnote>8}}      
#'\item{\code{Grade}}{a factor with levels \code{<6.5note} \code{>7.5note} \code{6.5-7note} \code{7-7.5note}}     
#'\item{\code{Startwork}}{a factor with levels \code{after.grad} \code{befor.grad}}
#'}



#'@references Lamberti, G. (2014) \emph{Modeling with Heterogeneity.} PhD Dissertation.
#' @source Laboratory of Information Analysis and Modeling (LIAM). 
#'    Facultat de Informatica de Barcelona, Universitat Politecnica de Catalunya.
#' @keywords datasets
NULL


#' Fibtelereg
#' 
#' Fibtelereg dataset 
#' 
#' @docType data
#' @name fibtelereg
#' @usage fibtelereg
#' @format A data frame with 147 observations on the following 18 variables. The first ten variables 
#' are segmentation variables. The rest of the variables refer to five variables 1) 
#' \code{Image} = Image, 2) \code{Exp.spec} = Specific Expectation, 3) \code{Exp.gen} = Generic Expectation,
#' 4)\code{Qual.spec} = Specific Quality, 5) \code{Qual.gen} = Generic Quality, 6) \code{Value} = Value, 7) 
#' \code{Satis} = Satisfaction. 
#'  Variables description
#' \itemize{
#' \item{\code{Image}}: Generic students perception of ICT schools: (internationally recognized, 
#'      ranges of courses, leader in research).
#' \item{\code{Exp.spec}}: Specific Expectation on specific skills (technic or applied skills).
#' \item{\code{Exp.gen}}: Generic Expectation on generic skills (abilities in problem solving, 
#'  communication skills).
#' \item{\code{Qual.spec}}: Perception about the achieved quality on the specific skills in the school.
#' \item{\code{Qual.gen}}: Perception about achieved quality on the  generic skills in 
#' the school (abilities in solving problem, communication skills).
#' \item{\code{Value}}: The advantage or profit that the alumni may draw from the school 
#' degree (well paid job, motivated job, prospectives in improvement and promotion).
#' \item{\code{Satis}}: Degree of alumni satisfaction about the formation in school respect to 
#' their actual work conditions.
#' }
#' 
#' Segmentation Variables description 
#'\itemize{
#'\item{\code{Career}}{a factor with levels \code{EI} \code{ETS} \code{TEL}}
#'\item{\code{Gender}}{a factor with levels \code{female} \code{male}}      
#'\item{\code{Age}}{a factor with levels \code{25-26years} \code{27-28years} \code{29-30years} \code{31years+}}   
#'\item{\code{Studying}}{a factor with levels \code{no.stud} \code{yes.stud}}    
#'\item{\code{Contract}} {a factor with levels \code{fix.cont} \code{other.cont} \code{temp.cont}} 
#'\item{\code{Salary}}{a factor with levels \code{18k} \code{>45k} \code{25k} \code{35k} \code{45k}}  
#'\item{\code{Firmtype}}{a factor with levels \code{priva} \code{publi}}       
#'\item{\code{Accgrade}}{a factor with levels \code{7-8accnote} \code{ accnote<7} \code{accnote>8}}      
#'\item{\code{Grade}}{a factor with levels \code{<6.5note} \code{>7.5note} \code{6.5-7note} \code{7-7.5note}}     
#'\item{\code{Startwork}}{a factor with levels \code{after.grad} \code{befor.grad}}
#'}
#' @references Lamberti, G. (2014) \emph{Modeling with Heterogeneity.} PhD Dissertation.
#' @source Laboratory of Information Analysis and Modeling (LIAM). 
#'    Facultat de Informatica de Barcelona, Universitat Politecnica de Catalunya.
#' @keywords datasets
NULL

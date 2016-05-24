#' @title Relational goods and Leisure time dataset
#' @description This dataset consists of the results of a survey aimed at measuring the evaluation
#' of relational goods and leisure time collected in December 2014. Every participant was asked 
#' to measure on a 10 point ordinal scale his/her personal score for several relational goods 
#' (for instance, time dedicated to friends and family) and to leisure time. 
#' In addition, the survey asked respondents to self-evaluate their level of happiness by marking 
#' a sign along a horizontal line of 110 millimeters according to their feeling, with the left-most 
#' extremity standing for "extremely unhappy", and the right-most extremity corresponding to
#' the status "extremely happy". 
#' @aliases relgoods
#' @usage data(relgoods)
#' @format The description of subjects' covariates is the following:
#' \describe{
#' \item{\code{ID}}{An identification number}
#' \item{\code{Gender}}{A factor with levels: 0 = man, 1 = woman}
#' \item{\code{BirthMonth}}{A variable indicating the month of  birth of the respondent}
#' \item{\code{BirthYear}}{A variable indicating the year of  birth of the respondent}
#' \item{\code{Family}}{A variable indicating the number of members of the family}
#' \item{\code{Year.12}}{A factor with levels: 1 = if there is any child aged less than 12  in the family,
#'  0 = otherwise}
#' \item{\code{EducationDegree}}{A factor with levels: 1 = compulsory school, 2 = high school diploma,
#'  3 = Graduated-Bachelor degree, 4 = Graduated-Master degree, 5 = Post graduated}
#' \item{\code{MaritalStatus}}{A factor with levels: 1 = Unmarried, 2 = Married/Cohabitant, 
#' 3 = Separated/Divorced, 4 = Widower}
#' \item{\code{Residence}}{A factor with levels: 1 = City of Naples, 2 = District of Naples,
#'  3 = Others Campania, 4 = Others Italia, 5= Foreign countries}
#' \item{\code{Glasses}}{A factor with levels: 1 = wearing glasses or contact lenses, 0 = otherwise}
#' \item{\code{RightHand}}{A factor with levels: 1 = right-handed, 0 = left-handed}
#' \item{\code{Smoking}}{A factor with levels: 1 = smoker, 0 = not smoker}
#' \item{\code{WalkAlone}}{A factor with levels: 1 = usually walking alone, 0= usually walking in company}
#' \item{\code{Job}}{A factor with levels: 1 = Not working, 2 = Retired, 3 = occasionally, 
#' 4 = fixed-term job, 5 = permanent job}
#' \item{\code{PlaySport}}{A factor with levels: 1 = Not playing any sport, 2= Yes, individual sport, 
#' 3 = Yes, team sport}
#' \item{\code{Pets}}{A factor with levels: 1 = owning a pet, 0 = not owning any pet}
#'} 
#' 1) Respondents were asked to evaluate the following items on a 10 point Likert scale, 
#' ranging from 1 = "never, at all" to 10 = "always, a lot":
#'  \describe{
#' \item{\code{WalkOut}}{How often the respondent goes out for a walk}
#' \item{\code{Parents}}{How often respondent talks at least to one of his/her parents}
#' \item{\code{MeetRelatives}}{How often respondent meets his/her relatives}
#' \item{\code{Association}}{Frequency of involvement in volunteering or different kinds of
#'  associations/parties, etc}
#' \item{\code{RelFriends}}{Quality of respondent's relationships with friends}
#' \item{\code{RelNeighbours}}{Quality of the relationships with neighbors}
#' \item{\code{NeedHelp}}{Easiness in asking help whenever in need}
#' \item{\code{Environment}}{Level of comfort with the surrounding environment}
#' \item{\code{Safety}}{Level of safety in the streets}
#' \item{\code{EndofMonth}}{Family making ends meet}
#' \item{\code{MeetFriend}}{Number of times the respondent met his/her friends during the month 
#' preceding the interview}
#' \item{\code{Physician}}{Importance of the kindness/simpathy in the selection of respondent's physician}
#'     }
#' \describe{
#'  \item{\code{Happiness}}{Each respondent was asked to mark a sign on a 110mm horizontal line
#'   according to his/her feeling of happiness (left endpoint corresponding to completely unhappy,
#'    right-most endpoint corresponding to extremely happy}
#' }   
#' 2) Respondents were asked to score the activities for leisure time listed below, according
#'  to their involvement/degree of amusement, on a 10 point Likert scale 
#'  ranging from 1 = "At all, nothing, never" to 10 = "Totally, extremely important, always":
#' \describe{
#' \item{\code{Videogames}}{}
#' \item{\code{Reading}}{}
#' \item{\code{Cinema}}{}
#' \item{\code{Drawing}}{}
#' \item{\code{Shopping}}{}
#' \item{\code{Writing}}{}
#' \item{\code{Bicycle}}{}
#' \item{\code{Tv}}{}
#' \item{\code{StayWFriend}}{Spending time with friends}{}
#' \item{\code{Groups}}{Taking part to associations, meetings, etc.}{}
#' \item{\code{Walking}}{}
#' \item{\code{HandWork}}{Hobby, gardening, sewing, etc. }
#' \item{\code{Internet}}{}
#' \item{\code{Sport}}{}
#' \item{\code{SocialNetwork}}{}
#' \item{\code{Gym}}{}
#' \item{\code{Quiz}}{Crosswords, sudoku, etc.}
#' \item{\code{MusicInstr}}{Playing a musical instrument}
#' \item{\code{GoAroundCar}}{ Hanging out by car}
#' \item{\code{Dog}}{Walking out the dog}
#' \item{\code{GoOutEat}}{Go to restaurants/pubs}}
#' @keywords datasets
#' @details 
#' Period of data collection: December 2014 \cr
#' Mode of collection: questionnaire \cr
#' Number of observations: 2459 \cr
#' Number of subjects' covariates: 16 \cr
#' Number of analyzed items: 34 \cr
#' Warning: with missing values
#' @source \url{http://www.labstat.it/home/wp-content/uploads/2015/09/relgoods.txt}


"relgoods"


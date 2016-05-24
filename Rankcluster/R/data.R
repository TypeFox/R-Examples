# DOCUMENTATION FOR DATA SET


#' @title rank data : APA
#' @docType data
#' @aliases APA
#' @name APA
#' @format A list containing :
#' \describe{
#'   \item{data}{ matrix of size 5738x5 containing the 5738 observed full ranks in ranking representation.
#'                The ranking representation r=(r_1,...,r_m) contains the  ranks assigned to the objects, and means that the ith object is in r_ith position.
#'                
#'               For example, if the ranking representation of a rank is (4,3,1,2,5), it means that judge ranks the first object in 4th position, second object in 3rd position, ...}
#'   
#'   \item{frequency}{matrix of size 120x6. Each row corresponds to one of the different observed rank.
#'                   The first fifth columns contains the observed ranks (ordering representation) and the sixth column
#'                    contains the frequency of observation.}
#'   
#'   \item{m}{ vector with the size of the ranks (5 here).}
#' }
#' 
#' @description  This dataset contains the 5738 full rankings resulting from the American Psychological Association (APA) presidential election of 1980.
#' For this election, members of APA had to rank five candidates in order of preference. 
#' 
#' For information, a total of 15449 votes have been registred for this election, but only the 5738 full rankings are reported in the APA dataset. Candidates A and C were research psychologists, candidates D and E were clinical psychologists and candidate B was a community psychologist.
#' 
#' 
#' @source "Group representations in probability and statistics", P. Diaconis, 1988.
#' 
#' @examples 
#'  data(APA)
#'  
#' @keywords datasets
NULL


#' @title rank data : big4
#' @docType data
#' @aliases big4
#' @name big4
#' @format A list containing :
#' \describe{
#'   \item{data}{A matrix of size 21*8 containing the 21 Premier League seasons. Each row corresponding to one ranking (ranking representation).
#'               
#'              The ranking representation r=(r_1,...,r_m) contains the  ranks assigned to the objects, and means that the ith object is in r_ith position.
#'              
#'               For example, if the ranking representation of a rank is (4,3,1,2,5), it means that judge ranks the first object in 4th position, second object in 3rd position, ...}
#'   \item{frequency}{matrix of size 21*9. Each row corresponds to one of the 21 different observed rankings, and the last column contains the observation frequency.}
#'   \item{m}{the size of the rankings (m=c(4,4) ).}
#' }
#' @description  This dataset is composed of the rankings (in ranking notation) of the "Big Four" English football teams (A: Manchester, B: Liverpool, C: Arsenal, D: Chelsea) to the English Championship (Premier League) and according to the UEFA coefficients (statistics used in Europe for ranking and seeding teams in international competitions), from 1993 to 2013.
#' 
#' In 2000-2001, Arsenal and Chelsea had the same UEFA coefficient and then are tied. UEFA ranking is (1, 4, 2, 2) for 2000-2001, what means that Manchester United is the first, Liverpool is the last, and the two intermediate positions are for Arsenal and Chelsea in an unknown order.
#' 
#' In 2009-2010, Liverpool and Arsenal have also the same UEFA coefficient, the ranking is  (1, 2, 2, 4).
#' 
#' @source   \url{http://en.wikipedia.org/wiki/Premier_League#.22Big_Four.22_dominance_.282000s.29}
#' 
#' \url{http://www.uefa.com/memberassociations/uefarankings/club/index.html}
#' 
#' @examples 
#'  data(big4)
#'  
#' @keywords datasets
NULL


#' @title Multidimensionnal partial rank data : eurovision
#' @docType data
#' @aliases eurovision
#' @name eurovision
#' @format A list containing:
#' \describe{
#'   \item{data}{ A matrix of size 34*48. Each row corresponds to the ranking representation of a multidimensionnal ranking.
#'                Columns 1 to 8 correspond to the 2007 contest, columns 9 to 18 to the 2008 contest, etc...
#'                
#'                The ranking representation r=(r_1,...,r_m) contains the  ranks assigned to the objects, and means that the ith object is in r_ith position.
#'                
#'                For example, if the ranking representation of a rank is (4,3,1,2,5), it means that judge ranks the first object in 4th position, second object in 3rd position, ...
#'                
#'   }
#'   
#'   \item{frequency}{A matrix of size 34*49 containing the differents multidimensionnal rankings. The 48 first columns are the same as in data, and the last column contains the frequency (1 for all ranks).}
#'   
#'   \item{m}{ a vector with the sizes of ranks for each dimension.}
#' }
#' @description    This dataset contains the ranking of the 8 common finalists of the Eurovision song contest from 2007 to 2012:
#' 
#' A: France, B:Germany, C:Greece, D:Romania, E:Russia, F:Spain, G:Ukrain, H:United Kingdom.
#' 
#' The number of rankings is 33, corresponding to the 33 European countries having participated to this six editions of the contest.
#' 
#' All the rankings are partial since none country has ranked this 8 countries in its 10 preferences. Missing ranking elements are zeros.
#' 
#' @source  \url{http://www.eurovision.tv} 
#' 
#' @examples 
#'  data(eurovision)
#'  
#' @keywords datasets
NULL

#' @title Multidimensionnal rank data : quiz
#' @docType data
#' @aliases quiz
#' @name quiz
#' @format A list containing:
#' \describe{
#'   \item{data}{a matrix of size 70*16. The student's answers are in row and the 16 columns correspond to the 4 rankings (for the 4 quizzes) of size 4 (ranking representation).
#' 
#'   The ranking representation r=(r_1,...,r_m) contains the  ranks assigned to the objects, and means that the ith object is in r_ith position.
#' 
#' For example, if the ranking representation of a rank is (4,3,1,2,5), it means that judge ranks the first object in 4th position, second object in 3rd position, ...}
#' \item{frequency}{a matrix of size 63*17. Each row corresponds to one of the 63 differents observed
#' rankings (ranking representation). Each row contains 4 ranks of size 4 and a last column for the frequency.}
#' \item{m}{a vector with the sizes of the ranks for each dimension.}
#' 
#' }
#' @description    This dataset contains the answers of 70 students (40 of third year and 30 of fourth year) from Polytech'Lille (statistics engineering school, France) to the four following quizzes:
#' 
#' \describe{
#'   
#'   \item{Literature Quiz}{
#'     This quiz consists of ranking four french writers according to chronological order:
#'       A=Victor Hugo, B=Moliere, C=Albert Camus, D=Jean-Jacques Rousseau.}
#'   
#'   \item{Football Quiz}{
#'     This quiz consists of ranking four national football teams according to increasing number of wins in the football World Cup: A=France, B=Germany, C=Brazil, D=Italy.}
#'   
#'   \item{Mathematics Quiz}{
#'     This quiz consists of ranking four numbers according to increasing order: 
#'       A=pi/3, B=log(1), C=exp(2), D=(1+sqrt(5))/2.}
#'   
#'   \item{Cinema Quiz}{
#'     This quiz consists of ranking four Tarentino's movies according to chronological order: 
#' A=Inglourious Basterds, B=Pulp Fiction, C=Reservoir Dogs, D=Jackie Brown.}
#' 
#' } 
#' 
#' @source   Julien Jacques
#' 
#' @examples 
#'  data(quiz)
#'  
#' @keywords datasets
NULL


#' @title rank data : sports
#' @docType data
#' @aliases sports
#' @name sports
#' @format A list containing :
#' \describe{
#'   \item{data}{a matrix containing 130 ranks of size 7 in ranking representation.
#'               
#'               The ranking representation r=(r_1,...,r_m) contains the  ranks assigned to the objects, and means that the ith object is in r_ith position.
#'               
#'               For example, if the ranking representation of a rank is (4,3,1,2,5), it means that judge ranks the first object in 4th position, second object in 3rd position, ...}
#'   
#'   \item{frequency}{a matrix with 123 differents ranks of size 7. In each row the first 7 columns correspond to one observed ranking and the last column contains the observation frequency.}
#'   \item{m}{ the size of the rankings (m=7).}
#' }
#' @description      This data set is due to Louis Roussos who asked 130 students at the
#' University of Illinois to rank seven sports according to their preference in participating:
#'  A = Baseball, B = Football, C = Basketball, D = Tennis, E = Cycling, F =
#'  Swimming, G = Jogging.
#' 
#' @source   J.I. Marden. "Analyzing and modeling rank data, volume 64 of Monographs on Statistics and Applied Probability". Chapman & Hall, London, 1995.
#' 
#' @examples 
#'  data(sports)
#'  
#' @keywords datasets
NULL

#' @title rank data : words
#' @docType data
#' @aliases words
#' @name words
#' @format A list containing :
#' \describe{
#'   \item{data}{A matrix of size 98*5 containing the 98 answers. Each row corresponding to one ranking (ranking representation).
#'               
#'               The ranking representation r=(r_1,...,r_m) contains the  ranks assigned to the objects, and means that the ith object is in r_ith position.
#'               
#'               For example, if the ranking representation of a rank is (4,3,1,2,5), it means that judge ranks the first object in 4th position, second object in 3rd position, ...}
#'  \item{frequency}{matrix of size 15*6. Each row corresponds to one of the 15 different observed rankings, and the last column contains the observation frequency.}
#'   \item{m}{the size of the rankings (m=5).}
#' }
#' @description   The data was collected under the auspices of the Graduate Record
#' Examination Board. A sample of 98 college students were asked to rank five words according to strength of association (least to most associated) with the target word "Idea":
#'   A = Thought, B = Play, C = Theory, D = Dream and E = Attention.
#' 
#' @source M.A. Fligner and J.S. Verducci. "Distance based ranking models". J. Roy. Statist. Soc. Ser. B, 48(3):359-369, 1986.
#' 
#' @examples 
#'  data(sports)
#'  
#' @keywords datasets
NULL

##' @name student.survey
##' @title Student Survey data
##' @description This data file consists of responses of graduate students in the social sciences enrolled in STA 6126 in a recent term at the University of Florida.
##' @docType data
##'
##' @format \Sexpr[stage=build,results=rd]{data(student.survey); smss:::describe_df(student.survey)}
##' \describe{
##' \item{\code{GE}}{gender}
##' \item{\code{AG}}{age in years}
##' \item{\code{HI}}{high school GPA (on a four-point scale)}
##' \item{\code{CO}}{college GPA}
##' \item{\code{DH}}{distance (in miles) of the campus from your home town}
##' \item{\code{DR}}{distance (in miles) of the classroom from your current residence}
##' \item{\code{TV}}{average number of hours per week that you watch TV}
##' \item{\code{SP}}{average number of hours per week that you participate in sports or have other physical exercise}
##'\item{\code{NE}}{number of times a week you read a newspaper}
##' \item{\code{AH}}{number of people you know who have died from AIDS or who are HIV+}
##' \item{\code{VE}}{whether you are a vegetarian}
##' \item{\code{PA}}{political affiliation (d = Democrat, r = Republican, i = independent)}
##' \item{\code{PI}}{political ideology}
##' \item{\code{RE}}{how often you attend religious services}
##' \item{\code{AB}}{opinion about whether abortion should be legal in the first three months of pregnancy}
##' \item{\code{AA}}{support affirmative action}
##' \item{\code{LD}}{belief in life after death}
##' }
##' @examples
##' data(student.survey)
##' summary(student.survey)
##' @source \url{http://www.stat.ufl.edu/~aa/social/data.html}
NULL

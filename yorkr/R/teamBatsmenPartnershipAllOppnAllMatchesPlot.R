##########################################################################################
# Designed and developed by Tinniam V Ganesh
# Date : 24 Mar 2016
# Function: teamBatsmenPartnershipAllOppnAllMatchesPlot
# This function computes the batting partnerships of a team against all oppositions and
# also the partenerships of th eopposition against this team
#
#
###########################################################################################
#' @title
#' Plots team batting partnership  all matches all oppositions
#'
#' @description
#' This function plots  the batting partnership of a team againt all oppositions in all matches
#' This function also returns a  dataframe  with the batting partnerships
#'
#' @usage
#' teamBatsmenPartnershipAllOppnAllMatchesPlot(matches,theTeam,main,plot=TRUE)
#'
#' @param matches
#' All the matches of the team against all oppositions
#'
#' @param theTeam
#' The team for which the the batting partnerships are sought
#'
#' @param main
#' The main team for which the the batting partnerships are sought

#'
#'@param plot
#' Whether the partnerships have top be rendered as a plot. If plot=FALSE the data frame is returned
#'
#' @return None or partnerships
#'
#' @references
#' \url{http://cricsheet.org/}\cr
#' \url{https://gigadom.wordpress.com/}\cr
#' \url{https://github.com/tvganesh/yorkrData}
#'
#' @author
#' Tinniam V Ganesh
#' @note
#' Maintainer: Tinniam V Ganesh \email{tvganesh.85@gmail.com}
#'
#' @examples
#' \dontrun{
#' # Get all matches for team India against all oppositions
#' d <- teamBatsmanVsBowlersAllOppnAllMatchesRept(matches,"India",rank=1,dispRows=50)
#'  #Plot the partnerships
#' teamBatsmenVsBowlersAllOppnAllMatchesPlot(d)
#'
#' #Do not plot but get the dataframe
#' e <- teamBatsmenVsBowlersAllOppnAllMatchesPlot(d,plot=FALSE)
#' }
#'
#' @seealso
#' \code{\link{teamBatsmenPartnershipAllOppnAllMatches}}\cr
#' \code{\link{teamBatsmenPartnershipOppnAllMatchesChart}}\cr
#' \code{\link{teamBatsmenVsBowlersAllOppnAllMatchesPlot}}\cr
#' \code{\link{teamBatsmenVsBowlersOppnAllMatches}}\cr
#'
#' @export
#'
teamBatsmenPartnershipAllOppnAllMatchesPlot <- function(matches,theTeam,main,plot=TRUE){
    team=batsman=nonStriker=runs=partnershipRuns=totalRuns=NULL
    a <- NULL
    a <-filter(matches,team==theTeam)

    #Get partnerships
    df <- data.frame(summarise(group_by(a,batsman,nonStriker),sum(runs)))
    names(df) <- c("batsman","nonStriker","runs")

    # Filter all rows where runs is 0. Problem when t2="Sri Lanka"
    # Sehwag and Ganguly show up as partnerships with runs=0
    #*****Check******** when the line below is removed
    df <- filter(df,runs!=0)
    df <- arrange(df,desc(runs))

    if(plot==TRUE){
        if(theTeam==main){
             plot.title <- paste(theTeam," batting partnerships")
        }else if(theTeam != main){
            plot.title <- paste(theTeam," batting partnerships against ", main)
        }
        ggplot(data=df,aes(x=batsman,y=runs,fill=nonStriker))+
            geom_bar(data=df,stat="identity") +
            xlab("Batsman") + ylab("Partnership runs") +
            ggtitle(bquote(atop(.(plot.title),
                                atop(italic("Data source:http://cricsheet.org/"),"")))) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
    } else{
        df
    }


}

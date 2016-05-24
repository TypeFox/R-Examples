#' @demoTitle baseball-aggregate
#' 
#' Demo compute aggregates
#'
#' To install and use baseball demo dataset in Aster:
#'
#' 1. download baseball.zip from
#'   https://bitbucket.org/grigory/toaster/downloads/baseball.zip
#' 2. run script to create data set in Aster
#'   sh load_baseball_data.sh -d mydbname -U username -w mypassword 
#' 3. create Aster ODBC DSN on your desktop
#'   see https://bitbucket.org/grigory/toaster/wiki/Home#markdown-header-odbc-driver-and-dns

library(toaster)
library(plyr)
library(ggplot2)

## utility input function
readlineDef <- function(prompt, default) {
  if (!is.null(prompt))
    prompt = paste0(prompt, "[", default, "]: ")
  else 
    prompt = paste0(prompt, ": ")
  
  result = readline(prompt)
  if (result == "") 
    return (default)
  else
    return (result)
}

## utility connection function
connectWithDSNToAster <- function(dsn=NULL) {
  dsn = readlineDef("Enter Aster ODBC DSN: ", dsn)
  
  tryCatch(close(conn), error=function(err) {NULL})
  
  conn = tryCatch({
    conn = odbcConnect(dsn)
    odbcGetInfo(conn)
    return (conn)
  }, error=function(err) {
    stop(paste("Can't connect to Aster - check DSN '", dsn, "'"))
  })
}

## connect to Aster first
conn = connectWithDSNToAster()

## must be connected to baseball dataset
if(!all(isTable(conn, c('pitching_enh', 'teams_enh')))) {
  stop("Must connect to baseball dataset and tables must exist.")
}

# compute team attendancy and average rank by decades
data = computeAggregates(channel = conn, tableName = "teams_enh",
                         by = c("franchid", "decadeid"),
                         aggregates = c("avg(rank) rank", "(avg(attendance)/1000)::int attendance"),
                         where="yearid >= 1980",     # attendance data became available in 1980
                         stringsAsFactors=TRUE)  
# display scatterplot of average rank vs. attendance by decade
ggplot(data) +
  geom_point(aes(x=rank, y=attendance, colour=franchid), size=4) +
  geom_text(aes(x=rank, y=attendance, label=franchid, colour=franchid), vjust=-1) +
  facet_wrap(~decadeid) + theme(legend.position="none") +
  ggtitle("Attendace vs. Rank by Decade (average)")


# compute total strikouts per team by league and decade
data = computeAggregates(channel = conn, "pitching_enh", 
               aggregates = c("avg(so) so"),
               by=c("lgid", "teamid", "decadeid"), # percent=c("lgid", "decadeid") 
               where="yearid between 1980 and 2009", stringsAsFactors=TRUE)
ggplot(data) +
  geom_histogram(aes(x=teamid, y=so, fill=teamid), stat="identity") +
  facet_grid(decadeid~lgid, scale="free_x") +
  scale_fill_manual(values = (getDiscretePaletteFactory("Set1"))(length(unique(data$teamid))),
                    guide=FALSE)


# compute total strikouts and percent using window function
data = computeAggregates(channel = conn, "pitching_enh",
               by = c("teamid", "yearid"), 
               aggregates = c("sum(r) r", 
                              "sum(r)/(sum(sum(r)) over (partition by yearid)) percent"),
               where = "yearid >= 2000 and teamid IN ('NYA','NYN','TEX','OAK','LAN','ATL','BAL','BOS','DET')", 
               stringsAsFactors=TRUE)

data = transform(data, yearid=factor(yearid)) # make yearid factor
data = data[with(data, order(teamid)), ] # reorder by teamid
data = ddply(data, .(yearid), transform, pos = cumsum(percent) - 0.5 * percent)
ggplot(data, aes(x=yearid, y=percent, fill=teamid)) +
    geom_bar(stat="identity", position="stack", width=1) +
    geom_text(aes(label = paste(teamid), y = pos), size = 3) +
  theme(legend.position = "bottom",
        axis.ticks = element_blank(),
        plot.background = element_blank()
        ,panel.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()) +  # theme_bw() +
  scale_fill_manual(values = (getDiscretePaletteFactory("Set1"))(length(unique(data$teamid))),
                    guide=guide_legend(title=NULL, nrow=2)) 
  
#' @title Create date type vectors 
#' @description This function allows you to create date type vectors.
#' @aliases Date
#' @author Javier Celigueta Munoz
#' @export Date
#' @param length Length of the date type vector
#' @return A date type vector of the indicated length
#' @examples 
#' Date(3)
#' Date(length=5)

Date=function( length = 0 ) {
  newDate = numeric( length )
  class(newDate) = "Date"
  return(newDate)
}

#' @title GanttChart
#' @description This function allows you to display a GANTT chart. It shows the duration, the sequence, the type of precedence constrain, the early and late start and finish dates and floats of all the activities.
#' @aliases GanttChart
#' @author Javier Celigueta Munoz
#' @export GanttChart
#' @import reshape2
#' @import ggplot2
#' @import grid
#' @param tasks A data.frame object with the next fields: task, start and end. All of them must be character vectors. task includes the names of the activities, start and end are the start and end dates for each activity.
#' @param information A data.frame object with the next fields: from, to, type and delay. from, to and delay must be numeric vectors while type is a character vector. from and to indicate all the sequences (the number of each activity is its position in the task vector), type indicates the type of precedence constrain for each link (SS, SF, FF or FS) and finally delay is the float for each one. 
#' @examples
#' project1=data.frame(
#'     task=c("Market Research","Concept Development","Viability Test",
#' "Preliminary Design","Process Design","Prototyping","Market Testing","Final Design",
#' "Launching"),
#'     start=c("2015-07-05","2015-07-05","2015-08-05","2015-10-05","2015-10-05","2016-02-18",
#' "2016-03-18","2016-05-18","2016-07-18"),
#'     end=c("2015-08-05","2015-08-05","2015-10-05","2016-01-05","2016-02-18","2016-03-18",
#' "2016-05-18","2016-07-18","2016-09-18"))
#' project2=data.frame(
#'     from=c(1,2,3,4,5,6,7,8),
#'     to=c(2,3,4,5,6,7,8,9),
#'     type=c("SS","FS","FS","SS","FS","FS","FS","FS"),
#'     delay=c(7,7,7,8,10,10,10,10))
#' GanttChart(project1,project2)
#' 
#' info=data.frame(
#'     task=c("Estimate market and make more exact marketing message",
#' "Design and order final package","Create press releases",
#' "Create product specification materials","Create marketing presentations",
#' "Transmit product launch details to international organization",
#' "Create sales, local, and product support groups training",
#' "Update product forecasts based on market feedback and analysis",
#' "Update launch plan based on forecast"),
#'     start=c("2015-08-20","2015-08-23","2015-08-23","2015-08-23","2015-08-23","2015-09-04",
#'     "2015-09-05","2015-08-23","2015-08-24"),
#'     end=c("2015-08-22","2015-08-29","2015-08-29","2015-09-03","2015-08-29","2015-09-05",
#'     "2015-09-17","2015-08-24","2015-08-28"))
#' details=data.frame(
#'     from=c(1,1,1,1,1,1,2,3,4,5,6,8,9),
#'     to=c(2,3,4,5,6,8,6,6,6,6,7,9,7),
#'     type=c("FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS"),
#'     delay=c(0,0,0,0,0,0,0,0,0,0,0,0,0))
#' GanttChart(info,details)

GanttChart=function(tasks,information) {
  value=NULL
  name=NULL
  timeorigin=as.Date("1970-01-01")
  term=as.Date(tasks$end)-as.Date(tasks$start)
  delay=vector(mode="numeric",length=length(tasks$task))
  for (i in seq(1,length(information$delay))) {
    delay[information$to[i]]=information$delay[i]
  }
  start.delay=as.Date(tasks$start)+delay
  end.delay=as.Date(tasks$end)+delay
  actv=rev(tasks$task)
  dfg=data.frame(
    name=factor(actv,levels=actv),
    start.date=as.Date(rev(start.delay)),
    end.date=as.Date(rev(end.delay))
  )
  mdfg=melt(dfg,id="name")
  previoustask=vector(mode="numeric",length(length(tasks$task)*length(tasks$task)))
  for (i in seq(1:length(information$from))) {
    previoustask[length(tasks$task)*information$to[i]-(length(tasks$task)-information$from[i])]=1
  }
  if (length(previoustask)!=(length(tasks$task)*length(tasks$task))) {
    for (i in seq(length(previoustask)+1,length(tasks$task)*length(tasks$task))) {
      previoustask[i]=0
    }
  }
  previous.task=rev(previoustask)
  previousmatrix=matrix(previous.task,nrow=length(actv),ncol=length(actv),byrow=TRUE)
  indexmatrix=which(previousmatrix==1,arr.ind=TRUE)
  ganttchart=ggplot(mdfg,aes(value,name)) + geom_line(size=6) + xlab(NULL) + ylab(NULL)
  for (i in seq(along=1:dim(indexmatrix)[1])) {
    if (rev(information$type)[i]=="FS") {
      ganttchart=ganttchart + geom_segment(x=as.numeric(difftime(mdfg$value[indexmatrix[i,2]+length(tasks$task)],timeorigin,units="days")),y=indexmatrix[i,2],xend=as.numeric(difftime(mdfg$value[indexmatrix[i,1]],timeorigin,units="days")),yend=indexmatrix[i,2],color="red") + geom_segment(x=as.numeric(difftime(mdfg$value[indexmatrix[i,1]],timeorigin,units="days")),y=indexmatrix[i,2],xend=as.numeric(difftime(mdfg$value[indexmatrix[i,1]],timeorigin,units="days")),yend=indexmatrix[i,1],arrow=arrow(angle=20,length=unit(0.11,"inches"),type="closed"),color="red")
    }
    if (rev(information$type)[i]=="FF") {
      ganttchart=ganttchart + geom_segment(x=as.numeric(difftime(mdfg$value[indexmatrix[i,2]+length(tasks$task)],timeorigin,units="days")),y=indexmatrix[i,2],xend=as.numeric(difftime(mdfg$value[indexmatrix[i,1]+length(tasks$task)]+1,timeorigin,units="days")),yend=indexmatrix[i,2],color="red") + geom_segment(x=as.numeric(difftime(mdfg$value[indexmatrix[i,1]+length(tasks$task)]+1,timeorigin,units="days")),y=indexmatrix[i,2],xend=as.numeric(difftime(mdfg$value[indexmatrix[i,1]+length(tasks$task)]+1,timeorigin,units="days")),yend=indexmatrix[i,1],color="red") + geom_segment(x=as.numeric(difftime(mdfg$value[indexmatrix[i,1]+length(tasks$task)]+1,timeorigin,units="days")),y=indexmatrix[i,1],xend=as.numeric(difftime(mdfg$value[indexmatrix[i,1]+length(tasks$task)],timeorigin,units="days")),yend=indexmatrix[i,1],arrow=arrow(angle=20,length=unit(0.11,"inches"),type="closed"),color="red")
    }
    if (rev(information$type)[i]=="SF") {
      ganttchart=ganttchart + geom_segment(x=as.numeric(difftime(mdfg$value[indexmatrix[i,2]],timeorigin,units="days")),y=indexmatrix[i,2],xend=as.numeric(difftime(mdfg$value[indexmatrix[i,2]],timeorigin,units="days")),yend=indexmatrix[i,2]-0.3,color="red") + geom_segment(x=as.numeric(difftime(mdfg$value[indexmatrix[i,2]],timeorigin,units="days")),y=indexmatrix[i,2]-0.3,xend=as.numeric(difftime(mdfg$value[indexmatrix[i,1]+length(tasks$task)],timeorigin,units="days")),yend=indexmatrix[i,2]-0.3,color="red") + geom_segment(x=as.numeric(difftime(mdfg$value[indexmatrix[i,1]+length(tasks$task)],timeorigin,units="days")),y=indexmatrix[i,2]-0.3,xend=as.numeric(difftime(mdfg$value[indexmatrix[i,1]+length(tasks$task)],timeorigin,units="days")),yend=indexmatrix[i,1],arrow=arrow(angle=20,length=unit(0.11,"inches"),type="closed"),color="red")
    }
    if (rev(information$type)[i]=="SS") {
      ganttchart=ganttchart + geom_segment(x=as.numeric(difftime(mdfg$value[indexmatrix[i,2]],timeorigin,units="days")),y=indexmatrix[i,2],xend=as.numeric(difftime(mdfg$value[indexmatrix[i,2]]-1,timeorigin,units="days")),yend=indexmatrix[i,2],color="red") + geom_segment(x=as.numeric(difftime(mdfg$value[indexmatrix[i,2]]-1,timeorigin,units="days")),y=indexmatrix[i,2],xend=as.numeric(difftime(mdfg$value[indexmatrix[i,2]]-1,timeorigin,units="days")),yend=indexmatrix[i,1],color="red") + geom_segment(x=as.numeric(difftime(mdfg$value[indexmatrix[i,2]]-1,timeorigin,units="days")),y=indexmatrix[i,1],xend=as.numeric(difftime(mdfg$value[indexmatrix[i,1]],timeorigin,units="days")),yend=indexmatrix[i,1],arrow=arrow(angle=20,length=unit(0.11,"inches"),type="closed"),color="red")     
    }
  }
  earlybegin=Date(length(information$to))
  earlyfinish=Date(length(information$to))
  term=as.Date(tasks$end)-as.Date(tasks$start)
  cont=0
  for (i in seq(1,length(tasks$task))) {
    if (length(which(information$to==i))!=0) {
      cont=cont+1
    }
  }
  aux1=vector(mode="numeric",length=cont)
  for (i in seq(1,length(tasks$task))) {
    if (length(which(information$to==i))!=0) {
      aux1[i]=information$to[which(information$to==i)[1]]
    }
    else {
      aux1[i]=0
    }
  }
  aux1=aux1[-which(aux1==0)]
  earlybegin=Date(max(aux1))
  earlyfinish=Date(max(aux1)) 
  for (i in aux1) {
    aux=Date(length(which(information$to==i)))
    if (length(aux)<=1) {       
      if (information$type[which(information$to==i)]=="SS") {
        earlybegin[i]=as.Date(tasks$start[information$from[which(information$to==i)]])+delay[information$from[which(information$to==i)]]
        earlyfinish[i]=earlybegin[i]+term[i]
      }
      if (information$type[which(information$to==i)]=="FS") {
        earlybegin[i]=as.Date(tasks$end[information$from[which(information$to==i)]])+delay[information$from[which(information$to==i)]]
        earlyfinish[i]=earlybegin[i]+term[i]
      }
      if (information$type[which(information$to==i)]=="FF") {
        earlyfinish[i]=as.Date(tasks$end[information$from[which(information$to==i)]])+delay[information$from[which(information$to==i)]]
        earlybegin[i]=earlyfinish[i]-term[i]
      }
      if (information$type[which(information$to==i)]=="SF") {
        earlyfinish[i]=as.Date(tasks$start[information$from[which(information$to==i)]])+delay[information$from[which(information$to==i)]]
        earlybegin[i]=earlyfinish[i]-term[i]
      }
    }
    if (length(aux)>1) {
      for (j in seq(1,length(aux))) {
        if (information$type[which(information$to==i)[j]]=="SS") {
          aux[j]=as.Date(tasks$start[information$from[which(information$to==i)[j]]])+delay[information$from[which(information$to==i)[j]]]
        }
        if (information$type[which(information$to==i)[j]]=="FS") {
          aux[j]=as.Date(tasks$end[information$from[which(information$to==i)[j]]])+delay[information$from[which(information$to==i)[j]]]
        }
        if (information$type[which(information$to==i)[j]]=="FF") {
          aux[j]=as.Date(tasks$end[information$from[which(information$to==i)[j]]])+delay[information$from[which(information$to==i)[j]]]-term[i]
        }
        if (information$type[which(information$to==i)[j]]=="SF") {
          aux[j]=as.Date(tasks$start[information$from[which(information$to==i)[j]]])+delay[information$from[which(information$to==i)[j]]]-term[i]
        }
      }
      earlybegin[i]=max(aux)
      earlyfinish[i]=earlybegin[i]+term[i]
    }
  }
  if (length(earlybegin)>cont) {
    earlybegin=earlybegin[-which(earlybegin=="1970-01-01")]
  }
  if (length(earlyfinish)>cont) {
    earlyfinish=earlyfinish[-which(earlyfinish=="1970-01-01")]
  }
  cont=0
  for (i in seq(1,length(tasks$task))) {
    if (length(which(information$from==i))!=0) {
      cont=cont+1
    }
  }
  aux2=vector(mode="numeric",length=cont)
  for (i in seq(1,length(tasks$task))) {
    if (length(which(information$from==i))!=0) {
      aux2[i]=information$from[which(information$from==i)[1]]
    }
    else {
      aux2[i]=0
    }
  }
  aux2=aux2[-which(aux2==0)]
  latebegin=Date(max(aux2))
  latefinish=Date(max(aux2))
  for (i in rev(aux2)) {
    aux=Date(length(which(information$from==i)))
    if (length(aux)<=1) {       
      if (information$type[which(information$from==i)]=="SS") {
        latebegin[i]=as.Date(tasks$start[information$to[which(information$from==i)]])+information$delay[which(information$from==i)]
        latefinish[i]=latebegin[i]+term[i]
      }
      if (information$type[which(information$from==i)]=="FS") {
        latefinish[i]=as.Date(tasks$start[information$to[which(information$from==i)]])+information$delay[which(information$from==i)]
        latebegin[i]=latefinish[i]-term[i]
      }
      if (information$type[which(information$from==i)]=="FF") {
        latefinish[i]=as.Date(tasks$end[information$to[which(information$from==i)]])+information$delay[which(information$from==i)]
        latebegin[i]=latefinish[i]-term[i]
      }
      if (information$type[which(information$from==i)]=="SF") {
        latebegin[i]=as.Date(tasks$end[information$to[which(information$from==i)]])+information$delay[which(information$from==i)]
        latefinish[i]=latebegin[i]+term[i]
      }
    }
    if (length(aux)>1) {
      for (j in seq(1,length(aux))) {
        if (information$type[which(information$from==i)[j]]=="SS") {
          aux[j]=as.Date(tasks$start[information$to[which(information$from==i)[j]]])+information$delay[which(information$from==i)[j]]+term[i]
        }
        if (information$type[which(information$from==i)[j]]=="FS") {
          aux[j]=as.Date(tasks$start[information$to[which(information$from==i)[j]]])+information$delay[which(information$from==i)[j]]
        }
        if (information$type[which(information$from==i)[j]]=="FF") {
          aux[j]=as.Date(tasks$end[information$to[which(information$from==i)[j]]])+information$delay[which(information$from==i)[j]]
        }
        if (information$type[which(information$from==i)[j]]=="SF") {
          aux[j]=as.Date(tasks$end[information$to[which(information$from==i)[j]]])+information$delay[which(information$from==i)[j]]+term[i]
        }
      }
      latefinish[i]=min(aux)
      latebegin[i]=latefinish[i]-term[i]
    }
  }
  if (length(latebegin)>cont) {
    latebegin=latebegin[-which(latebegin=="1970-01-01")]
  }
  if (length(latefinish)>cont) {
    latefinish=latefinish[-which(latefinish=="1970-01-01")]
  }
  for (i in seq(1,length(aux1))) {
    ganttchart=ganttchart + geom_point(x=as.numeric(difftime(earlybegin[i],timeorigin,units="days")),y=length(tasks$task)-aux1[i]+1+0.2,size=3,pch=25,colour="blue")
    ganttchart=ganttchart + geom_point(x=as.numeric(difftime(earlyfinish[i],timeorigin,units="days")),y=length(tasks$task)-aux1[i]+1+0.2,size=3,pch=25,colour="blue")
  }
  for (i in seq(1,length(aux2))) {
    ganttchart=ganttchart + geom_point(x=as.numeric(difftime(latebegin[i],timeorigin,units="days")),y=length(tasks$task)-aux2[i]+1-0.2,size=3,pch=24,colour="blue")
    ganttchart=ganttchart + geom_point(x=as.numeric(difftime(latefinish[i],timeorigin,units="days")),y=length(tasks$task)-aux2[i]+1-0.2,size=3,pch=24,colour="blue")
  }
  return(ganttchart)
}
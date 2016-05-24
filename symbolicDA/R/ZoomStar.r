#***********************************************************************************************************************************************
#*  
#*  (C) 2010     Andrzej Dudek    Uniwersytet Ekonomiczny we Wrocławiu
#*  
#*  Wykres typ Zoom Star dla danych symbolicznych
#*  Skrypt do książki:
#*  "Podejście symboliczne w analizie danych ekonomicznych" powstałej w ramach pro-jektu badawczego habilitacyjnego N N111 105 234
#*  
#*  Kod poniższy może być modyfikowany, kopiowany i rozprowadzany na warunkach li-cencji GPL 2 (http://gnu.org.pl/text/licencja-gnu.html), 
#*  a w szczególności pod warunkiem umieszczenia w zmodyfikowanym pliku widocznej informacji o dokonanych zmianach, wraz z datą ich dokonania. 
#*  
#***********************************************************************************************************************************************



.inn<-function(x,y){
  ret<-NULL
  for(xx in x){
    eq<-F
      for(yy in y){
      if(xx==yy){
      eq<-T
      }
    }
    ret<-c(ret,eq)
  }
  ret
}

.tickCoord<-function(r,degree,tick,labelOffset){
  c(r*tick*cos(degree)-labelOffset*cos(pi/2+degree)*sign(cos(degree)),r*tick*sin(degree)-labelOffset*sin(pi/2+degree)*sign(cos(degree)))
}

.zoomStar<-function(table.Symbolic,j,variableSelection=NULL,offset=0.2,firstTick=0.2,labelCex=.8,labelOffset=.7,tickLength=.3,histWidth=0.04,histHeight=2,rotateLabels=TRUE,variableCex=NULL){
  size<-10
  r<-size
  offset<-as.integer(offset*r)
  if(is.null(variableSelection)){
    variableSelection<-1:nrow(table.Symbolic$variables)
  }
  degreeStep<-2*pi/length(variableSelection)
  #print(variableSelection)
  plot((-size-offset):(size+offset),(-size-offset):(size+offset),type="n",xlab="Zoom star",ylab=paste("Symbolic Object",table.Symbolic$individuals[j,"name"]),axes=FALSE,pty="s")
  for(i in 1:length(variableSelection)){
    lo<-labelOffset
    degree<-degreeStep*(i-1)
    if(rotateLabels){
      srt<-degree*180/pi
    }
    else{
      srt<-0
    }
    #print(i)
    v<-variableSelection[i]
    #print(v)
    lines(c(0,r*cos(degree)),c(0,r*sin(degree)))
    if(!is.null(variableCex)){
    text(r*cos(degree)+labelOffset* sign(cos(degree)),r*sin(degree)+labelOffset*sign(sin(degree)),table.Symbolic$variables[v,"label"],cex=variableCex)
    }
    else{
    text(r*cos(degree)+labelOffset* sign(cos(degree)),r*sin(degree)+labelOffset*sign(sin(degree)),table.Symbolic$variables[v,"label"])
    }
    if(table.Symbolic$variables[v,"type"]=="IC"){
      vmin<-min(table.Symbolic$indivIC[,v,1])
      vmax<-max(table.Symbolic$indivIC[,v,2])
      v1<-table.Symbolic$indivIC[j,v,1]
      v2<-table.Symbolic$indivIC[j,v,2]
      tick<-firstTick
      lines(c(r*tick*cos(degree)-tickLength*sin(degree),r*tick*cos(degree)+tickLength*sin(degree)),c(r*tick*sin(degree)+tickLength*cos(degree),
      r*tick*sin(degree)-tickLength*cos(degree)))
      text(.tickCoord(r,degree,tick,lo)[1],.tickCoord(r,degree,tick,lo)[2],format(vmin,digits=2),srt=srt,cex=labelCex)      
      tick<-1
      lines(c(r*tick*cos(degree)-tickLength*sin(degree),r*tick*cos(degree)+tickLength*sin(degree)),c(r*tick*sin(degree)+tickLength*cos(degree),
      r*tick*sin(degree)-tickLength*cos(degree)))
      text(.tickCoord(r,degree,tick,lo)[1],.tickCoord(r,degree,tick,lo)[2],format(vmax,digits=2),srt=srt,cex=labelCex)      
      tick1<-firstTick+(1-firstTick)*(v1-vmin)/(vmax-vmin)
      tick2<-firstTick+(1-firstTick)*(v2-vmin)/(vmax-vmin)
      polygon(c(r*tick1*cos(degree)+tickLength*cos(pi/2+degree),
      r*tick1*cos(degree)-tickLength*cos(pi/2+degree),r*tick2*cos(degree)-tickLength*cos(pi/2+degree),r*tick2*cos(degree)+tickLength*cos(pi/2+degree))
      ,c(r*tick1*sin(degree)+tickLength*sin(pi/2+degree),
      r*tick1*sin(degree)-tickLength*sin(pi/2+degree),r*tick2*sin(degree)-tickLength*sin(pi/2+degree),r*tick2*sin(degree)+tickLength*sin(pi/2+degree)),
      col="lightgreen")
      text(.tickCoord(r,degree,tick1,lo)[1],.tickCoord(r,degree,tick1,lo)[2],format(v1,digits=2),srt=srt,cex=labelCex)      
      text(.tickCoord(r,degree,tick2,lo)[1],.tickCoord(r,degree,tick2,lo)[2],format(v2,digits=2),srt=srt,cex=labelCex)      
    } 
    if(table.Symbolic$variables[v,"type"]=="MN" || table.Symbolic$variables[v,"type"]=="N" ){
      varList<-table.Symbolic$detailsListNom[table.Symbolic$detailsListNom[,"details_no"]==as.numeric(as.matrix(table.Symbolic$variables)[v,"details"]),]
      if(is.null(dim(varList))){
        dim(varList)<-c(1,length(varList))
      }
      for(i in 1:nrow(varList)){
        tick<-firstTick+(1-firstTick)*(i-1)/nrow(varList)
        indivN1<-table.Symbolic$indivN[table.Symbolic$indivN[,"variable"]==v,]
        indivN2<-indivN1[indivN1[,"indiv"]==j,]
        #print(indivN2)
        #print((varList))
        #print(as.numeric(varList[i,"num"]))
        if(any(indivN2[,"value"]==as.numeric(as.matrix(varList[i,"num"])))){
          points(r*tick*cos(degree),r*tick*sin(degree),pch=21,bg="black",cex=1.3)      
          tt<-varList[i,"name"]
          #print(paste("tt",tt))
          text(.tickCoord(r,degree,tick,lo)[1],.tickCoord(r,degree,tick,lo)[2],
           substitute(underline(x), list(x=as.character(tt))),srt=srt,cex=labelCex*1.2)      

        }
        else{
          points(r*tick*cos(degree),r*tick*sin(degree),pch=21,bg="green",col="green")      
          text(.tickCoord(r,degree,tick,lo)[1],.tickCoord(r,degree,tick,lo)[2],varList[i,"name"],srt=srt,cex=labelCex)      
        }
      }
    }
    if(table.Symbolic$variables[v,"type"]=="NM"){
      varList<-table.Symbolic$detailsListNomModif[table.Symbolic$detailsListNomModif[,"details_no"]==table.Symbolic$variables[v,"details"],]
      if(is.null(dim(varList))){
        dim(varList)<-c(1,length(varList))
      }
      cl=c("lightpink3","springgreen","thistle1","paleturquoise4","olivedrab2")[sample(5)]
      for(ii in 1:nrow(varList)){
        tick<-firstTick+(1-firstTick)*(ii-1)/nrow(varList)
        indivN1<-table.Symbolic$indivNM[table.Symbolic$indivNM[,"variable"]==v,]
        indivN2<-indivN1[.inn(indivN1[,"indiv"],j),]
        if(any(indivN2[,"value"]==as.vector(varList[ii,"num"]))){
          #print(indivN2[indivN2[,"value"]==as.numeric(varList[i,"num"]),])
          #print(paste("frequency",indivN2[1,"frequency"]))
          height<-histHeight*mean(indivN2[indivN2[,"value"]==as.vector(varList[ii,"num"]),][,"frequency"])
          if(i==3 && ii==2)print(paste(i,ii,height))
          # starting 3-d box from begin or end of item position
          if(degree>pi/4 && degree<pi) {
            fromStart<-1
          }
          else{
            #print(degree)
            fromStart<-0
            #print(fromStart)
          }
          degreeOrth<-degree
          leftRight<-1
          if(degree>pi/4 && degree<pi/2){
            leftRight=-1
          }
          if(degree>pi && degree<3*pi/2){
            leftRight=-1
          }
          if(sin(degreeOrth)<.5)degreeOrth<-pi/4
          A1<-c(r*tick*cos(degree),r*tick*sin(degree))
          A2<-c(r*(tick-histWidth)*cos(degree),r*(tick-histWidth)*sin(degree))
          A3<-c(r*(tick-histWidth)*cos(degree),r*(tick-histWidth)*sin(degree)+height)
          A4<-c(r*tick*cos(degree),r*tick*sin(degree)+height)
          A5<-c(r*(tick-fromStart*histWidth)*cos(degree)+leftRight*histWidth/2*r*abs(sin(degreeOrth)),r*(tick-fromStart*histWidth)*sin(degree)+r*histWidth*abs(cos(degreeOrth)))
          A6<-c(r*(tick-fromStart*histWidth)*cos(degree)+leftRight*histWidth/2*r*abs(sin(degreeOrth)),r*(tick-fromStart*histWidth)*sin(degree)+histWidth*r*abs(cos(degreeOrth))+height  )
          A7<-c(r*(tick-(1-fromStart)*histWidth)*cos(degree)+leftRight*histWidth/2*r*abs(sin(degreeOrth)),r*(tick-(1-fromStart)*histWidth)*sin(degree)+height+histWidth*r*abs(cos(degreeOrth)))
          A1X<-A1[1]
          A1Y<-A1[2]
          A2X<-A2[1]
          A2Y<-A2[2]
          A3X<-A3[1]
          A3Y<-A3[2]
          A4X<-A4[1]
          A4Y<-A4[2]
          A5X<-A5[1]
          A5Y<-A5[2]
          A6X<-A6[1]
          A6Y<-A6[2]
          A7X<-A7[1]
          A7Y<-A7[2]

          polygon(c(A1X,A2X,A3X,A4X),
          c(A1Y,A2Y,A3Y,A4Y)
          ,col=cl)
          if(fromStart==1){
          polygon(c(A2X,A3X,A6X,A5X),
          c(A2Y,A3Y,A6Y,A5Y)
          ,col=cl)
          polygon(c(A3X,A4X,A7X,A6X),
          c(A3Y,A4Y,A7Y,A6Y)
          ,col=cl)
          }
          else{
          polygon(c(A1X,A4X,A6X,A5X),
          c(A1Y,A4Y,A6Y,A5Y)
          ,col=cl)
          polygon(c(A3X,A4X,A6X,A7X),
          c(A3Y,A4Y,A6Y,A7Y)
          ,col=cl)
          }
        }
        else{
          points(r*tick*cos(degree),r*tick*sin(degree),pch=21,bg="green",col="green")      
        }
          text(.tickCoord(r,degree,tick,lo)[1],.tickCoord(r,degree,tick,lo)[2],paste("k",ii,sep=""),srt=srt,cex=labelCex)      
      }
    }
  }
}

zoomStar<-function(table.Symbolic,j,variableSelection=NULL,offset=0.2,firstTick=0.2,labelCex=.8,labelOffset=.7,tickLength=.3,histWidth=0.04,histHeight=2,rotateLabels=TRUE,variableCex=NULL){
  .zoomStar(table.Symbolic,j,variableSelection,offset,firstTick,labelCex,labelOffset,tickLength,histWidth,histHeight,rotateLabels,variableCex)
}

zoomStar<-function(table.Symbolic,j,variableSelection=NULL,offset=0.2,firstTick=0.2,labelCex=.8,labelOffset=.7,tickLength=.3,histWidth=0.04,histHeight=2,rotateLabels=TRUE,variableCex=NULL){
  .zoomStar(table.Symbolic,j,variableSelection,offset,firstTick,labelCex,labelOffset,tickLength,histWidth,histHeight,rotateLabels,variableCex)
}

if (getRversion() >= "2.15.1") utils::globalVariables(c('positions',
                                                        'deltas',
                                                        'gammas',
                                                        'vegas',
                                                        'neutrality',
                                                        'option1',
                                                        'option2',
                                                        'tkrp'))

greekneutrality <-
function(){
  
  my.draw<-function(panel){
    
    positions<- as.numeric(strsplit(panel$positions,",")[[1]])
    deltas<- as.numeric(strsplit(panel$deltas,",")[[1]])
    gammas<- as.numeric(strsplit(panel$gammas,",")[[1]])
    vegas<-  as.numeric(strsplit(panel$vegas,",")[[1]])
    
    option1_delta<- as.numeric(strsplit(panel$option1,",")[[1]])[1]
    option1_gamma<- as.numeric(strsplit(panel$option1,",")[[1]])[2]
    option1_vega<-  as.numeric(strsplit(panel$option1,",")[[1]])[3]
    
    option2_delta<- as.numeric(strsplit(panel$option2,",")[[1]])[1]
    option2_gamma<- as.numeric(strsplit(panel$option2,",")[[1]])[2]
    option2_vega<-  as.numeric(strsplit(panel$option2,",")[[1]])[3]
    
    port_delta=sum(positions*deltas)
    port_gamma=sum(positions*gammas)
    port_vega=sum(positions*vegas)
    
    if(panel$neutrality=="Delta"){
      
      underly_posn= -port_delta
      posn1= ""
      posn2= ""
      underly_posn<-paste("Position in the underlying: $ ", underly_posn, sep="")
    }
    
    else if(panel$neutrality=="Delta and Gamma"){
      posn1<- -port_gamma/option1_gamma
      new_delta=port_delta + option1_delta*posn1
      posn1<- paste("Position in option 1: ", posn1, " contracts", sep="")
      posn2<-""
      underly_posn= -new_delta
      underly_posn<-paste("Position in the underlying: $ ", underly_posn, sep="")
    }
    
    else if(panel$neutrality=="Delta and Vega"){
      posn1<- -port_vega/option1_vega
      posn2<- ""
      new_delta=port_delta + option1_delta*posn1
      underly_posn <- -new_delta
      posn1<- paste("Position in option 1: ", posn1, " contracts", sep="")
      underly_posn<-paste("Position in the underlying: $ ", underly_posn, sep="")
    }
    
    else {
      # solve a set of eqns- Ainverse * B
      gamma_row <- c(option1_gamma,option2_gamma)
      vega_row <-  c(option1_vega,option2_vega)
      greek_matrix <- matrix(c(gamma_row,vega_row),nrow=2,byrow=T)
      B_matrix<-matrix(c(port_gamma,port_vega),byrow=T)
      posns <- solve(greek_matrix,B_matrix)
      posn1<- -posns[1]
      posn2<- -posns[2]
      new_delta=port_delta + option1_delta*posn1 + option2_delta*posn2
      posn1<- paste("Position in option 1: ", posn1, " contracts", sep="")
      posn2<- paste("Position in option 2: ", posn2, " contracts", sep="")
      underly_posn= -new_delta
      underly_posn<-paste("Position in the underlying: $ ", underly_posn, sep="")
    }
    
    
    plot(1:50, 1:50, type="n", xlab="", ylab="",
         axes=FALSE, frame = TRUE)
    text(25, 35, underly_posn,cex=1.2)
    text(25, 25, posn1,cex=1.2)
    text(25, 15, posn2,cex=1.2)
    
    
    panel
  }
  
  my.redraw <- function(panel)
  {
    rp.tkrreplot(panel, tkrp)
    panel                                                                       
  }
  
  
  
  
  my.panel <- rp.control("Hedging with Greeks")
  rp.textentry(my.panel,variable=positions,title="Positions ",initval="-1000, -500, -2000, -500",action=my.redraw)
  rp.textentry(my.panel,variable=deltas,title="Deltas      ",initval="0.5, 0.8, -0.4, 0.7",action=my.redraw)
  rp.textentry(my.panel,variable=gammas,title="Gammas ",initval="2.2, 0.6, 1.3, 1.8",action=my.redraw)
  rp.textentry(my.panel,variable=vegas,title="Vegas       ",initval="1.8, 0.2, 0.7, 1.4",action=my.redraw)
  #rp.doublebutton(my.panel,variable=discrate,step=0.25,title="Discount Rate (%  p.a.)",initval=10,range=c(1,15),showvalue=TRUE,action=my.redraw)
  rp.radiogroup(panel=my.panel,variable=neutrality, title="Type of Neutrality desired", 
                vals=c("Delta", "Delta and Gamma","Delta and Vega", "Delta, Gamma, and Vega"),action=my.redraw)
  rp.textentry(my.panel,variable=option1,title="Delta, Gamma, Vega of traded option 1 ",initval="0.6, 1.5, 0.8",action=my.redraw)
  rp.textentry(my.panel,variable=option2,title="Delta, Gamma, Vega of traded option 2 ",initval="0.1, 0.5, 0.6",action=my.redraw)
  rp.tkrplot(panel=my.panel , name=tkrp, plotfun=my.draw,hscale=2,vscale=1.5)
  
  #rp.do(my.panel, my.draw)
}

###########################################################
# This is a little GUI to help select team names
###########################################################
team.name.select <- function(newname=NULL, team.resolver, team.data, scores, type="alt.name")
{ 
  require(tcltk)
  return.stuff=FALSE #You have to hit Save or Finish to return stuff
  
  gogetem <- function() {
    tclvalue(text.var)=""
    ok=TRUE
    s <- as.integer(tkcurselection(box))
    if (!length(s)){
      ok=FALSE
      tclvalue(text.var)="\nSelect team from list\nor click Enter New Team"
    }else{ disp.name=l[s+1] }
    if(ok){
      if(type=="alt.name") team.res=rbind(team.res,data.frame(name=disp.name,alt.name=current.newname))
      if(type=="disp.name"){
        team.res$name[team.res$name==current.newname]=disp.name
      }
      #assign in the calling function
      assign("team.res",team.res,envir=sys.frame(fun.frm))
      return.stuff=TRUE
      assign("return.stuff",TRUE,envir=sys.frame(fun.frm))
      
      if(which.newname<n.newnames){ 
        which.newname=which.newname+1
        assign("which.newname",which.newname,envir=sys.frame(fun.frm))
        if(type=="alt.name") tclvalue(team.var)=paste("Bad alt.name:\n",newname[which.newname])
        if(type=="disp.name") tclvalue(team.var)=paste("Bad display name:\n",newname[which.newname])
        current.newname=newname[which.newname]
        if("venue" %in% names(scores)){
          league.names=unique(scores$venue[scores$home.team==current.newname | scores$away.team==current.newname])
        }
        tclvalue(team.leag.var)=paste(c("Leagues:",league.names), collapse="\n")
        assign("which.newname",which.newname,envir=sys.frame(fun.frm))
        assign("current.newname",current.newname,envir=sys.frame(fun.frm))
      }else{
        which.newname=n.newnames+1
        assign("which.newname",which.newname,envir=sys.frame(fun.frm))
        tclvalue(text.var)="\nNo teams left to enter\n Click save or cancel to finish.\n"
        current.newname=""
        tclvalue(team.leag.var)=""
        assign("current.newname",current.newname,envir=sys.frame(fun.frm))        
      }
    }
  }
  
  clearteam <- function() {
    tclvalue(text.var)=""
    tkselection.clear(box, 0, length(l))
  }
  #show info on selected team in panel on right
  showteam = function(){
    tclvalue(text.var)=""
    ok=TRUE
    s <- as.integer(tkcurselection(box))
    if (!length(s)){
      ok=FALSE
      tclvalue(text.var)="\nNo team selected"
    }else{ disp.name=l[s+1] }
    if(ok){
      if(tclvalue(what.to.show)=="2") disp.name=team.res$name[team.res$alt.name==disp.name]
      this.team=team.dat[team.dat$name==disp.name,]
      for(i in names(team.dat)){
        tclObj.name = as.character(get(i))
        tclvalue(tclObj.name) = this.team[[i]]
      }
      tclvalue(alt.name.text)=paste(c("Alt names used:",team.res$alt.name[team.res$name==this.team$name]), collapse="\n")
      if("venue" %in% names(scores)){
        alt.names=unique(team.res$alt.name[team.res$name==disp.name])
        league.names=unique(scores$venue[  scores$home.team==disp.name | scores$away.team==disp.name |
                                             scores$home.team %in% alt.names | scores$away.team %in% alt.names])
      }
      tclvalue(league.text)=paste(c("Leagues:",league.names), collapse="\n")
    }
  }
  enterteam=function(){
    tclvalue(text.var)=""
    bad.val=c("[\\]","[?]","@","%","[.]","[:^:]","/","&","[*]","[{]","[}]")
    
    if(tclvalue(name)=="" | any(str_detect(tclvalue(name),bad.val))){
      if(tclvalue(name)=="")  tclvalue(text.var)="\nYou need to enter\na team name."
      if(any(str_detect(tclvalue(name),bad.val)))  tclvalue(text.var)="\nDo not use punctuation\nin your team name."
    }else{
      if(tclvalue(name) %in% team.dat$name){
        for(i in names(team.dat)){
          val=tclvalue(get(i))
          if(val==i) val=""
          team.dat[team.dat$name==tclvalue(name),i]=val
        }
        if(type=="alt.name") team.res=rbind(team.res,data.frame(name=tclvalue(name),alt.name=current.newname))
        if(type=="disp.name")team.res$name[team.res$name==current.newname]=tclvalue(name)
      }else{
        new.team.data=data.frame(.tmp="")
        for(i in names(team.dat)){
          val=tclvalue(get(i))
          if(val==i) val=""
          new.team.data[i]=val
        }
        new.team.data=new.team.data[,-1]
        team.dat=rbind(team.dat,new.team.data)
        if(type=="alt.name") team.res=rbind(team.res,data.frame(name=new.team.data$name,alt.name=current.newname))
        if(type=="disp.name")team.res$name[team.res$name==current.newname]=new.team.data$name
      }
      #assign in the calling function
      assign("team.dat",team.dat,envir=sys.frame(fun.frm))
      assign("team.res",team.res,envir=sys.frame(fun.frm))
      return.stuff=TRUE
      assign("return.stuff",TRUE,envir=sys.frame(fun.frm))
      if(which.newname>=n.newnames){ 
        which.newname=n.newnames+1
        tclvalue(text.var)="\nNo teams left to enter\n Click save or cancel to finish.\n"
        assign("which.newname",which.newname,envir=sys.frame(fun.frm))
      }else{
        which.newname=which.newname+1
        if(type=="alt.name") tclvalue(team.var)=paste("Bad alt.name:\n",newname[which.newname])
        if(type=="disp.name") tclvalue(team.var)=paste("Bad display name:\n",newname[which.newname])
        current.newname=newname[which.newname]
        if("venue" %in% names(scores)){
          league.names=unique(scores$venue[scores$home.team==current.newname | scores$away.team==current.newname])
        }
        tclvalue(team.leag.var)=paste(c("Leagues:",league.names), collapse="\n")
        assign("which.newname",which.newname,envir=sys.frame(fun.frm))
        assign("current.newname",current.newname,envir=sys.frame(fun.frm))
      }
    }
  }
  
  updateteam=function(){
    tclvalue(text.var)=""
    bad.val=c("[\\]","[?]","@","%","[.]","[:^:]","/","&","[*]","[{]","[}]")
    s <- as.integer(tkcurselection(box))
    ok=TRUE
    if (!length(s)){
      ok=FALSE
      tclvalue(text.var)="\nNo team selected"
    }else{ disp.name=l[s+1] }
    if(tclvalue(what.to.show)=="2") disp.name=team.res$name[team.res$alt.name==disp.name]
    if(ok){
      if(tclvalue(name)!=disp.name){
        tclvalue(text.var)="\nYou have changed\nthe display name.\nThe corresponding\ndisplay name in the\nteam resolver is\nalso updated."
        team.res$name[team.res$name==disp.name]=tclvalue(name)
        assign("team.res",team.res,envir=sys.frame(fun.frm))
      }
      for(i in names(team.dat)){
        val=tclvalue(get(i))
        team.dat[team.dat$name==disp.name,i]=val
        #assign in the calling function
        assign("team.dat",team.dat,envir=sys.frame(fun.frm))
      }
      return.stuff=TRUE
      assign("return.stuff",TRUE,envir=sys.frame(fun.frm))
    }
  }
  
  #What to show, when radio buttons above left side are pressed
  showlist=function(){
    tclvalue(text.var)=""
    val = unique(team.dat$name)
    if(tclvalue(what.to.show)=="2") val = unique(team.res$alt.name)
    tclvalue(lvar) = paste("{",val,"}",sep="",collapse=" ")
    assign("l",val,envir=sys.frame(fun.frm))
  }
  searchteam=function(){
    tclvalue(text.var)=""
    search.val=tclvalue(search.text)
    val = unique(team.dat$name)
    if(tclvalue(what.to.show)=="2") val = unique(team.res$alt.name)
    val = val[str_detect(val,search.val)]
    tclvalue(lvar) = paste("{",val,"}",sep="",collapse=" ")
    #so that the other functions know what is listed
    assign("l",val,envir=sys.frame(fun.frm))
  }
  skipteam=function(){
    if(which.newname>n.newnames){ 
      tclvalue(team.var)="\nNo teams left.\n Use save or cancel to finish.\n"
    }else{
    #add skipped name to end of list
    newname=c(newname,newname[which.newname])
    #remove skipped name from inside list
    newname=newname[-which.newname]
    tclvalue(team.var)=newname[which.newname]
    current.newname=newname[which.newname]
    if("venue" %in% names(scores)){
      league.names=unique(scores$venue[scores$home.team==current.newname | scores$away.team==current.newname])
    }
    tclvalue(team.leag.var)=paste(c("Leagues:",league.names), collapse="\n")
    assign("newname",newname,envir=sys.frame(fun.frm))
    assign("current.newname",current.newname,envir=sys.frame(fun.frm))
    tclvalue(text.var)="\nSkipped team is added to\nend of list.  If you want\nto exit save without\nentering all teams,\nhit SAVE below.\n"
  }
  }
  
  fun.frm=sys.nframe()
  
  base <- tktoplevel()
  tkwm.title(base, "Select Team for unknown name")
  
  spec.frm <- tkframe(base,borderwidth=4)
  left.frm <- tkframe(spec.frm)
  middle.frm <- tkframe(spec.frm, borderwidth=10)
  right.frm <- tkframe(spec.frm)
  farright.frm <- tkframe(spec.frm)
  
  #LEFT COLUMN is where you select teams
  ## Two left frames:
  frame1 <- tkframe(left.frm, relief="groove", borderwidth=2)
  team.res=team.resolver #team.res is the updated one
  team.dat=team.data #team.res is the updated one
  l = unique(team.dat$name)
  lvar <- tclVar()
  tclObj(lvar) <- l
  box <- tklistbox(frame1, listvariable = lvar, selectmode = "single", width="40")
  what.to.show  <- tclVar(1)
  tkpack(tkradiobutton(frame1, command=showlist, text="team.name",
                       value=1, variable=what.to.show), anchor="w")
  tkpack(tkradiobutton(frame1, command=showlist, text="alt.name",
                       value=2, variable=what.to.show), anchor="w")
  tkpack(box)
  tkpack(tklabel(frame1, text="Use mouse or arrows\n to scroll team names"))
  #Frame 2 has the buttons
  frame2 <- tkframe(left.frm, relief="groove", borderwidth=2)  
  
  #Set up the search frame
  search.frm <- tkframe(frame2)
  left.search.frm <- tkframe(search.frm)
  right.search.frm <- tkframe(search.frm)
  search.text <- tclVar("text to search")
  entry.Name <-tkentry(left.search.frm,width="20",textvariable=search.text)
  tkpack(entry.Name)
  tkpack(tkbutton(right.search.frm, text = gettext("Search"), command = searchteam))
  
  tkpack(tkbutton(frame2, text = gettext("Clear Selection"), command = clearteam))
  
  tkpack(frame1, frame2, search.frm, fill="x")
  tkpack(search.frm, fill="x")
  tkpack(left.search.frm, right.search.frm, side="left", anchor="n")
  
  #  tkpack(frame3, frame4, fill="x")
  tkpack(left.frm, middle.frm, right.frm, farright.frm, side="left", anchor="n")
  
  
  ## RIGHT COLUMN IS THE TEAM INFO
  enterteam.frm <- tkframe(right.frm)
  tkpack(enterteam.frm, search.frm, fill="x")
  name=tclVar("") #otherwise R CMD check complains that name is not bound
  for(i in names(team.dat)){
    assign(i,tclVar(""))
    tmp=tkentry(enterteam.frm,width="40",textvariable=get(i))
    tkgrid(tklabel(enterteam.frm,text=i), tmp)
  }
  
  ## MIDDLE COLUMN IS THE BUTTONS
  tkpack(tkbutton(middle.frm, text = gettext("Show data (on right)\non selected team (on left)"), command = showteam, border=4))
  tkpack(tkbutton(middle.frm, text = gettext("<-- Use Selected Team\n(updates team resolver)"), command = gogetem))
  tkpack(tkbutton(middle.frm, text = gettext("Enter New Team -->\n(updates team data and\nteam resolver)"), command = enterteam))
  tkpack(tkbutton(middle.frm, text = gettext("Skip Team\n(return to it later)"), command = skipteam))
  tkpack(tkbutton(middle.frm, text = gettext("Update team data -->\n(updates team file\nand team resolver\nif you change name)"), command = updateteam, border=4))
  
  #Unknown name
  n.newnames=length(newname)
  team.var=tclVar()
  tclObj(team.var) <- "aalkjda;lksj"
  tclvalue(team.var) <- ""
  if(type=="alt.name" & n.newnames!=0) tclvalue(team.var)=paste("Bad alt.name:\n",newname[1])
  if(type=="disp.name" & n.newnames!=0) tclvalue(team.var)=paste("Bad display name:\n",newname[1])
  which.newname=1  #where you are at
  current.newname=newname[which.newname] #which badname is being entered
  if(is.null(current.newname)) current.newname=""
  tkpack(tklabel(middle.frm,textvariable=team.var))
  #Unknown name appears in the following venues
  team.leag.var=tclVar()
  tclObj(team.leag.var) <- "aalkjda;lksj"
  tclvalue(team.leag.var) <- ""
  if("venue" %in% names(scores)){
    league.names=unique(scores$venue[scores$home.team==current.newname | scores$away.team==current.newname])
  }
  tclvalue(team.leag.var)=paste(c("Leagues:",league.names), collapse="\n")
  tkpack(tklabel(middle.frm,textvariable=team.leag.var))
  
  #Alert text
  text.var=tclVar()
  tclObj(text.var) <- "aalkjda;lksj"
  tclvalue(text.var) <- ""
  tkpack(tklabel(middle.frm,textvariable=text.var,fg="red"))
  
  
  ## FAR RIGHT COLUMN IS More TEAM INFO
  altname.frm = tkframe(farright.frm)
  tkpack(altname.frm, fill="x")
  alt.name.text=tclVar()
  tclObj(alt.name.text) <- "aalkjda;lksj"
  tclvalue(alt.name.text) <- ""
  tkpack(tklabel(altname.frm,textvariable=alt.name.text),side="left")
  league.frm = tkframe(farright.frm)
  tkpack(league.frm, fill="x")
  league.text=tclVar()
  tclObj(league.text) <- "aalkjda;lksj"
  tclvalue(league.text) <- ""
  tkpack(tklabel(league.frm,textvariable=league.text),side="left")
  
  
  ## `Bottom frame' (on base):
  save.but <- tkbutton(base,text="SAVE (Exits the GUI now and save team resolver)", command=function(){ assign("return.stuff",TRUE,envir=sys.frame(fun.frm)); tkdestroy(base) })
  q.but <- tkbutton(base,text="Cancel (Exits the GUI without changing team resolver)", command=function(){ assign("return.stuff",FALSE,envir=sys.frame(fun.frm)); tkdestroy(base) })
  
  tkpack(spec.frm, save.but, q.but)
  
  tkwait.window(base)
  
  if(return.stuff==TRUE){
    if(which.newname>n.newnames){ skipped.teams=c()
    }else{ skipped.teams=newname[which.newname:n.newnames] }
    team.res=team.res[sort(team.res$name, index.return=TRUE)$ix,,drop=FALSE]
    team.dat=team.dat[sort(team.dat$name, index.return=TRUE)$ix,,drop=FALSE]
    return(list(team.resolver=team.res, team.data=team.dat, skipped.teams=skipped.teams, updated=TRUE)) #updated team resolver
  }
  if(return.stuff==FALSE) 
    return(list(team.resolver=team.resolver, team.data=team.data, skipped.teams=newname, updated=FALSE)) #did not update team resolver
  
}

edit_team_data = function(team.resolver, team.data, scores){
  rtn.obj=team.name.select(newname=NULL, team.resolver, team.data, scores, type="alt.name")
  return(rtn.obj)
}
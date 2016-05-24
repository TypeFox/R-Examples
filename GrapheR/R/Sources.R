#-------------------------------------------------
# Langage
#-------------------------------------------------

language<-function(lang=NULL) {
  if (!is.null(lang)) {
    Env$lang<-lang
  } else {
    Env$lang<-read.table(file.path(path.package("GrapheR"),"lang","Language.txt",fsep=.Platform$file.sep))[1,1]
  }
}


#-------------------------------------------------
# Chargement 1
#-------------------------------------------------

run.GrapheR<-function(lang=NULL,path.to.save=NULL,figurej=FALSE) {
  language(lang=lang)
  Env$path.to.save<-path.to.save
  Env$figurej <- figurej
  load.GrapheR()
}


#-------------------------------------------------
# Chargement 2 et définition des variables
#-------------------------------------------------

load.GrapheR<-function() {
  Env$img<-if (Env$lang=="en") {
    read.csv(file.path(path.package("GrapheR"),"lang","Images_en.csv",fsep=.Platform$file.sep),header=FALSE,sep=";")
    } else if (Env$lang=="fr") {
    read.csv(file.path(path.package("GrapheR"),"lang","Images_fr.csv",fsep=.Platform$file.sep),header=FALSE,sep=";")
    } else if (Env$lang=="es") {
    read.csv(file.path(path.package("GrapheR"),"lang","Images_es.csv",fsep=.Platform$file.sep),header=FALSE,sep=";")
    } else if (Env$lang=="de") {
    read.csv(file.path(path.package("GrapheR"),"lang","Images_de.csv",fsep=.Platform$file.sep),header=FALSE,sep=";")
    }
  Env$voc<-if (Env$lang=="en") {
    read.csv(file.path(path.package("GrapheR"),"lang","Language_en.csv",fsep=.Platform$file.sep),header=FALSE,as.is=1,sep=";")
    } else if (Env$lang=="fr") {
    read.csv(file.path(path.package("GrapheR"),"lang","Language_fr.csv",fsep=.Platform$file.sep),header=FALSE,as.is=1,sep=";")
    } else if (Env$lang=="es") {
    read.csv(file.path(path.package("GrapheR"),"lang","Language_es.csv",fsep=.Platform$file.sep),header=FALSE,as.is=1,sep=";")
    } else if (Env$lang=="de") {
    read.csv(file.path(path.package("GrapheR"),"lang","Language_de.csv",fsep=.Platform$file.sep),header=FALSE,as.is=1,sep=";")
    }
  Env$police<-tkfont.create(family="Arial",size=8)
  Env$police2<-tkfont.create(family="Arial",size=4)
  Env$police3<-tkfont.create(family="Arial",size=10,weight="bold")
  Env$police4<-tkfont.create(family="Courier",size=8)
  Env$police5<-tkfont.create(family="Arial",size=10)
  Env$police6<-tkfont.create(family="Arial",size=10,slant="italic")
  Env$dataset<-NULL
  Env$l.frames<-list()
  Env$l.frames$Fr1.status<-1
  Env$l.frames$Fr2.status<-0
  Env$l.frames$Fr3.status<-0
  Env$l.frames$Fr4.status<-0
  Env$l.frames$Fr5.status<-0
  Env$l.frames$Fr6.status<-0
  Env$l.fr1<-list()
  Env$l.fr2<-list()
  Env$l.fr3<-list()
  Env$l.fr4<-list()
  Env$l.fr5<-list()
  Env$l.fr6<-list()
  Env$l.fr7<-list()
  Env$l.lab<-list()
  Env$l.wdg<-list()
  Env$l.var<-list()
  Env$l.var$message<-tclVar("")
  Env$l.var$ecran<-"D"
  Env$l.var$extension<-tclVar("txt")
  Env$l.var$sepcol<-tclVar(Env$voc[5,1])
  Env$l.var$sepdec<-tclVar(Env$voc[9,1])
  Env$l.var$na<-tclVar("NA")
  Env$l.var$header<-tclVar(1)
  Env$l.var$regroup1<-tclVar(0)
  Env$l.var$regroup2<-tclVar("long")
  Env$l.var$regroup3<-tclVar("2")
  Env$l.var$var.num<-""
  Env$l.var$var.fact<-""
  Env$l.var$variable<-tclVar("")
  Env$l.var$facteur1<-tclVar("")
  Env$l.var$niveau<-tclVar("")
  Env$l.var$hist.type<-tclVar("")
  Env$l.var$encadre<-tclVar(0)
  Env$l.var$titre<-tclVar("")
  Env$l.var$titre.col<-tclVar("black")
  Env$l.var$titre.taille<-tclVar("1.5")
  Env$l.var$soustitre<-tclVar("")
  Env$l.var$graduations.col<-tclVar("black")
  Env$l.var$graduations.taille<-tclVar("1")
  Env$l.var$graduations.orient<-tclVar(Env$voc[246,1])
  Env$l.var$legendes.col<-tclVar("black")
  Env$l.var$legendes.taille<-tclVar("1")
  Env$l.var$titre.axehor<-tclVar("")
  Env$l.var$titre.axever<-tclVar("")
  Env$l.var$liminf.axehor<-tclVar("Auto")
  Env$l.var$limsup.axehor<-tclVar("Auto")
  Env$l.var$liminf.axever<-tclVar("Auto")
  Env$l.var$limsup.axever<-tclVar("Auto")
  Env$l.var$log.axehor<-tclVar(0)
  Env$l.var$log.axever<-tclVar(0)
  Env$l.var$hist.barres<-tclVar("Auto")
  Env$l.var$couleur1A<-tclVar("grey")
  Env$l.var$col.borduresA<-tclVar("black")
  Env$l.var$hist.dens<-tclVar(0)
  Env$l.var$couleur2A<-tclVar("black")
  Env$l.var$trait1<-tclVar("")
  Env$l.var$epaisseur1<-tclVar("1")
  Env$l.var$facteur.interaction<-""
  Env$l.var$box.orient<-tclVar(Env$voc[68,1])
  Env$l.var$titre.axenoms<-tclVar("")
  Env$l.var$titre.axevaleurs<-tclVar("")
  Env$l.var$liminf.axevaleurs<-tclVar("Auto")
  Env$l.var$limsup.axevaleurs<-tclVar("Auto")
  Env$l.var$log.axevaleurs<-tclVar(0)
  Env$l.var$boxmoy<-tclVar(0)
  Env$l.var$ICmediane<-tclVar(0)
  Env$l.var$varwidth<-tclVar(0)
  Env$l.var$lg.moustaches<-tclVar("1.5")
  Env$l.var$noms1<-""
  Env$l.var$moyprop<-tclVar("moy")
  Env$l.var$facteur2<-tclVar("")
  Env$l.var$proportions<-tclVar("")
  Env$l.var$prop.niveaux<-tclVar("")
  Env$l.var$nomsprop<-""
  Env$l.var$moyprop<-tclVar("moy")
  Env$l.var$facteurprop<-tclVar("")
  Env$l.var$couleur1B<-"white"
  Env$l.var$col.borduresB<-"black"
  Env$l.var$hachuresA<-tclVar("1")
  Env$l.var$hachuresB<-1
  Env$l.var$stack<-tclVar(0)
  Env$l.var$erreur<-tclVar("")
  Env$l.var$legende<-tclVar(0)
  Env$l.var$legende.pos<-tclVar("")
  Env$l.var$legende.titre<-tclVar("")
  Env$l.var$noms2<-""
  Env$l.var$plusieurs<-tclVar(0)
  Env$l.var$nomsprop.fac<-""
  Env$l.var$cam.orient<-tclVar("")
  Env$l.var$cam.start<-tclVar("0")
  Env$l.var$cam.lien<-tclVar(0)
  Env$l.var$parts.niveaux<-tclVar("")
  Env$l.var$nomsparts<-""
  Env$l.var$varX<-tclVar("")
  Env$l.var$varY<-tclVar("")
  Env$l.var$varX.prop<-tclVar("")
  Env$l.var$couleur2B<-"grey"
  Env$l.var$taille.ptsA<-tclVar("1")
  Env$l.var$taille.ptsB<-1
  Env$l.var$type.courbeA<-tclVar("")
  Env$l.var$type.courbeB<-""
  Env$l.var$trait2<-""
  Env$l.var$epaisseur2<-1
  Env$l.var$symboleA<-tclVar("1")
  Env$l.var$symboleB<-1
  Env$l.var$select<-1
  Env$l.var$droiteA<-tclVar("")
  Env$l.var$droiteB<-""
  Env$l.var$intervalA<-tclVar("")
  Env$l.var$intervalB<-""
  Env$l.var$ptlab<-tclVar(0)
  Env$l.var$sysinfo<-tclVar(0)
  Env$l.var$levels.temp <- NULL
  Env$l.var$nobar<-tclVar(0)
  Env$l.var$nw.col<-tclVar("white")
  Env$l.var$nw.lignes<-tclVar("1")
  Env$l.var$nw.colonnes<-tclVar("1")
  Env$l.var$add.abscisses<-NULL
  Env$l.var$add.hauteurs<-NULL
  Env$l.var$add.param1<-tclVar("")
  Env$l.var$add.param2<-tclVar("")
  Env$l.var$add.param3<-tclVar("")
  Env$l.var$add.trait<-tclVar("")
  Env$l.var$add.epaisseur1<-tclVar("1")
  Env$l.var$add.col1<-tclVar("black")
  Env$l.var$add.distrib<-tclVar("norm")
  Env$l.var$add.seq<-NULL
  Env$l.var$add.seq2<-NULL
  Env$l.var$add.txt<-tclVar("")
  Env$l.var$add.epaisseur2<-tclVar("1")
  Env$l.var$add.col2<-tclVar("black")
  Env$l.var$add.matrice<-NULL
  Env$l.var$fen.num<-tclVar("")
  Env$l.var$fen.larg<-tclVar("1250")
  Env$l.var$fen.type<-tclVar("jpg")
  Env$l.var$fen.res<-tclVar("150")
  Env$l.code<-list()
  if (!is.null(Env$path.to.save)) {
    Env$l.code$ask<-TRUE
    Env$l.code$save<-TRUE
    Env$l.code$folder<-Env$path.to.save
  } else {
    Env$l.code$ask<-FALSE
    Env$l.code$save<-FALSE
    Env$l.code$folder<-NULL
  }
  Env$l.code$graphsnb<-0
  Env$l.code$x.inf<-NULL
  Env$l.code$x.sup<-NULL
  Env$l.code$y.inf<-NULL
  Env$l.code$y.sup<-NULL
  ouvrir.GrapheR()
}


#-------------------------------------------------
# Affichage d'un msg dans la barre des
# msgs
#-------------------------------------------------

msg<-function(text,type) {
  if (type=="error") {
    tkconfigure(Env$l.wdg$message.wdg,foreground="red")
    tclvalue(Env$l.var$message)<-paste(Env$voc[18,1],text)
  }
  if (type=="warning") {
    tkconfigure(Env$l.wdg$message.wdg,foreground="darkgreen")
    tclvalue(Env$l.var$message)<-paste(Env$voc[160,1],text)
  }
  if (type=="info") {
    tkconfigure(Env$l.wdg$message.wdg,foreground="darkblue")
    tclvalue(Env$l.var$message)<-text
  }
  return(invisible())
}


#-------------------------------------------------
# Séparation des variables numériques et des
# facteurs
#-------------------------------------------------

variables.class<-function() {
  Env$l.var$var.num<-NULL
  for (i in 1:ncol(Env$dataset)) {
    if (is.numeric(Env$dataset[,i])) {Env$l.var$var.num<-c(Env$l.var$var.num,names(Env$dataset)[i])}
  }
  if (is.null(Env$l.var$var.num)) {Env$l.var$var.num<-""}
  Env$l.var$var.fact<-NULL
  for (i in 1:ncol(Env$dataset)) {
    if (is.factor(Env$dataset[,i])) {Env$l.var$var.fact<-c(Env$l.var$var.fact,names(Env$dataset)[i])}
  }
  if (is.null(Env$l.var$var.fact)) {Env$l.var$var.fact<-""}
}


#-------------------------------------------------
# Chargement du jeu de données - fichier externe
#-------------------------------------------------

data.load1<-function() {
  file<-tclvalue(tkgetOpenFile(filetypes=paste("{{.",tclvalue(Env$l.var$extension),"} {.",tclvalue(Env$l.var$extension),"}}",sep="")))
  if (!nchar(file)) {
    msg(text=Env$voc[19,1],type="error")
  } else {
    Env$dataset<-read.table(file,dec=ifelse(tclvalue(Env$l.var$sepdec)==Env$voc[9,1],".",","),header=ifelse(tclvalue(Env$l.var$header)==1,TRUE,FALSE),
	na.strings=tclvalue(Env$l.var$na),sep=if (tclvalue(Env$l.var$sepcol)==Env$voc[5,1]) {""} else if (tclvalue(Env$l.var$sepcol)==Env$voc[6,1]) {","} else {";"})
    Env$loading<-paste("read.table(\"",file,"\",\n  dec=",ifelse(tclvalue(Env$l.var$sepdec)==Env$voc[9,1],"\".\"","\",\""),", header=",ifelse(tclvalue(Env$l.var$header)==1,"TRUE","FALSE"),
	", na.strings=\"",tclvalue(Env$l.var$na),"\", sep=",if (tclvalue(Env$l.var$sepcol)==Env$voc[5,1]) {"\"\""} else if (tclvalue(Env$l.var$sepcol)==Env$voc[6,1]) {"\",\""} else {"\";\""},")",sep="")
    if ("var.list"%in%names(Env$l.fr2)) {
	tkdelete(Env$l.fr2$var.list,0,"end")
	for (i in 1:ncol(Env$dataset)) {tkinsert(Env$l.fr2$var.list,"end",colnames(Env$dataset)[i])}
	tkconfigure(Env$l.fr2$type.wdg,text="")
	tkconfigure(Env$l.fr2$resume.wdg,state="normal")
	tkdelete(Env$l.fr2$resume.wdg,"0.0","end")
	tkconfigure(Env$l.fr2$resume.wdg,state="disabled")
    }
    if ("var.list"%in%names(Env$l.fr3)) {
	tkdelete(Env$l.fr3$var.list,0,"end")
	for (i in 1:ncol(Env$dataset)) {tkinsert(Env$l.fr3$var.list,"end",colnames(Env$dataset)[i])}
    }
    if ("var.list"%in%names(Env$l.fr4)) {
	tkdelete(Env$l.fr4$var.list,0,"end")
	for (i in 1:ncol(Env$dataset)) {tkinsert(Env$l.fr4$var.list,"end",colnames(Env$dataset)[i])}
    }
    variables.class()
    if ("fact.wdg"%in%names(Env$l.fr5)) {
	tkconfigure(Env$l.fr5$fact.wdg,values=Env$l.var$var.fact)
	tclvalue(Env$l.var$facteur1) <- ""
	Env$l.var$levels.temp <- NULL
	tkdelete(Env$l.fr5$liste.actual,0,"end")
	tkdelete(Env$l.fr5$liste.new,0,"end")
    }
    msg(text=Env$voc[21,1],type="info")
  }
}


#-------------------------------------------------
# Chargement du jeu de données - objet R
#-------------------------------------------------

data.load2<-function() {
  if (nchar(tclvalue(tkcurselection(Env$l.fr1$obj.list)))>0) {
    tables<-NULL
    for (i in 1:length(ls(.GlobalEnv))) {
	if (is.data.frame(get(ls(.GlobalEnv)[i]))) {tables<-c(tables,ls(.GlobalEnv)[i])}
    }
    Env$dataset<-get(tables[as.numeric(tclvalue(tkcurselection(Env$l.fr1$obj.list)))+1],pos=.GlobalEnv)
    Env$loading<-tables[as.numeric(tclvalue(tkcurselection(Env$l.fr1$obj.list)))+1]
    if ("var.list"%in%names(Env$l.fr2)) {
	tkdelete(Env$l.fr2$var.list,0,"end")
	for (i in 1:ncol(Env$dataset)) {tkinsert(Env$l.fr2$var.list,"end",colnames(Env$dataset)[i])}
	tkconfigure(Env$l.fr2$type.wdg,text="")
	tkconfigure(Env$l.fr2$resume.wdg,state="normal")
	tkdelete(Env$l.fr2$resume.wdg,"0.0","end")
	tkconfigure(Env$l.fr2$resume.wdg,state="disabled")
    }
    if ("var.list"%in%names(Env$l.fr3)) {
	tkdelete(Env$l.fr3$var.list,0,"end")
	for (i in 1:ncol(Env$dataset)) {tkinsert(Env$l.fr3$var.list,"end",colnames(Env$dataset)[i])}
    }
    if ("var.list"%in%names(Env$l.fr4)) {
	tkdelete(Env$l.fr4$var.list,0,"end")
	for (i in 1:ncol(Env$dataset)) {tkinsert(Env$l.fr4$var.list,"end",colnames(Env$dataset)[i])}
    }
    variables.class()
    if ("fact.wdg"%in%names(Env$l.fr5)) {
	tkconfigure(Env$l.fr5$fact.wdg,values=Env$l.var$var.fact)
	tclvalue(Env$l.var$facteur1) <- ""
	Env$l.var$levels.temp <- NULL
	tkdelete(Env$l.fr5$liste.actual,0,"end")
	tkdelete(Env$l.fr5$liste.new,0,"end")
    }
    msg(text=Env$voc[21,1],type="info")
  } else {
    msg(text=Env$voc[20,1],type="error")
  }
}


#-------------------------------------------------
# Renommer une variable
#-------------------------------------------------

rename.variable<-function() {
  if(nchar(tclvalue(tkget(Env$l.fr3$nom.wdg)))>0) {
    if (nchar(tclvalue(tkcurselection(Env$l.fr3$var.list)))>0) {
	names(Env$dataset)[as.numeric(tclvalue(tkcurselection(Env$l.fr3$var.list)))+1]<-tclvalue(tkget(Env$l.fr3$nom.wdg))
	tkdelete(Env$l.fr3$var.list,0,"end")
	for (i in 1:ncol(Env$dataset)) {tkinsert(Env$l.fr3$var.list,"end",colnames(Env$dataset)[i])}
	tkdelete(Env$l.fr3$nom.wdg,0,"end")
	if ("var.list"%in%names(Env$l.fr2)) {
	  tkdelete(Env$l.fr2$var.list,0,"end")
	  for (i in 1:ncol(Env$dataset)) {tkinsert(Env$l.fr2$var.list,"end",colnames(Env$dataset)[i])}
	}
	if ("var.list"%in%names(Env$l.fr4)) {
	  tkdelete(Env$l.fr4$var.list,0,"end")
	  for (i in 1:ncol(Env$dataset)) {tkinsert(Env$l.fr4$var.list,"end",colnames(Env$dataset)[i])}
	}
    } else {
	msg(text=Env$voc[25,1],type="error")
    }
    variables.class()
    if ("fact.wdg"%in%names(Env$l.fr5)) {
	tkconfigure(Env$l.fr5$fact.wdg,values=Env$l.var$var.fact)
	tclvalue(Env$l.var$facteur1) <- ""
	Env$l.var$levels.temp <- NULL
	tkdelete(Env$l.fr5$liste.actual,0,"end")
	tkdelete(Env$l.fr5$liste.new,0,"end")
    }
  } else {
    msg(text=Env$voc[24,1],type="error")
  }
}


#-------------------------------------------------
# Transformer une variable en facteur
#-------------------------------------------------

convert.variable<-function(type) {
  num<-as.numeric(tclvalue(tkcurselection(Env$l.fr4$var.list)))+1
  if (type=="character") {
    Env$dataset[,num]<-factor(Env$dataset[,num])
  } else {
    if (tclvalue(Env$l.var$regroup1)==0) {
	Env$dataset[,num]<-factor(Env$dataset[,num])
    } else {
	if (tclvalue(Env$l.var$regroup1)=="long") {
	  Env$dataset[,num]<-cut(Env$dataset[,num],breaks=as.numeric(tclvalue(Env$l.var$regroup3)),labels=as.character(1:as.numeric(tclvalue(Env$l.var$regroup3))))
	} else {
	  Env$dataset[,num]<-cut(Env$dataset[,num],breaks=quantile(Env$dataset[,num],probs=seq(0,1,1/as.numeric(tclvalue(Env$l.var$regroup3))),na.rm=TRUE),labels=as.character(1:as.numeric(tclvalue(Env$l.var$regroup3))))
	}
    }
  }
  tclvalue(Env$l.var$regroup1)<-0
  tclvalue(Env$l.var$regroup2)<-"long"
  tkconfigure(Env$l.fr4$rb.noregroup,state="disabled")
  tkconfigure(Env$l.fr4$rb.regroup1,state="disabled")
  tkconfigure(Env$l.fr4$rb.regroup2,state="disabled")
  tkconfigure(Env$l.fr4$rb.regroup3,state="disabled")
  tkconfigure(Env$l.fr4$curs.wdg,state="disabled",foreground="grey")
  tkconfigure(Env$l.fr4$curs.lab,foreground="grey")
  tkconfigure(Env$l.fr4$but,state="disabled")
  variables.class()
  if ("fact.wdg"%in%names(Env$l.fr5)) {
    tkconfigure(Env$l.fr5$fact.wdg,values=Env$l.var$var.fact)
    tclvalue(Env$l.var$facteur1) <- ""
    Env$l.var$levels.temp <- NULL
    tkdelete(Env$l.fr5$liste.actual,0,"end")
    tkdelete(Env$l.fr5$liste.new,0,"end")
  }
  msg(text=Env$voc[35,1],type="info")
}


#-------------------------------------------------
# Remise à zéro des variables graphiques
#-------------------------------------------------

reinit.variables<-function() {
  tclvalue(Env$l.var$variable)<-""
  tclvalue(Env$l.var$facteur1)<-""
  tclvalue(Env$l.var$niveau)<-""
  tclvalue(Env$l.var$hist.type)<-""
  tclvalue(Env$l.var$encadre)<-0
  tclvalue(Env$l.var$titre)<-""
  tclvalue(Env$l.var$titre.col)<-"black"
  tclvalue(Env$l.var$titre.taille)<-"1.5"
  tclvalue(Env$l.var$soustitre)<-""
  tclvalue(Env$l.var$graduations.col)<-"black"
  tclvalue(Env$l.var$graduations.taille)<-"1"
  tclvalue(Env$l.var$legendes.col)<-"black"
  tclvalue(Env$l.var$graduations.orient)<-Env$voc[246,1]
  tclvalue(Env$l.var$legendes.taille)<-"1"
  tclvalue(Env$l.var$titre.axehor)<-""
  tclvalue(Env$l.var$titre.axever)<-""
  tclvalue(Env$l.var$liminf.axehor)<-"Auto"
  tclvalue(Env$l.var$limsup.axehor)<-"Auto"
  tclvalue(Env$l.var$liminf.axever)<-"Auto"
  tclvalue(Env$l.var$limsup.axever)<-"Auto"
  tclvalue(Env$l.var$hist.barres)<-"Auto"
  tclvalue(Env$l.var$couleur1A)<-"grey"
  tclvalue(Env$l.var$col.borduresA)<-"black"
  tclvalue(Env$l.var$hist.dens)<-0
  tclvalue(Env$l.var$couleur2A)<-"black"
  tclvalue(Env$l.var$trait1)<-""
  tclvalue(Env$l.var$epaisseur1)<-"1"
  Env$l.var$facteur.interaction<-""
  tclvalue(Env$l.var$box.orient)<-Env$voc[68,1]
  tclvalue(Env$l.var$titre.axenoms)<-""
  tclvalue(Env$l.var$titre.axevaleurs)<-""
  tclvalue(Env$l.var$liminf.axevaleurs)<-"Auto"
  tclvalue(Env$l.var$limsup.axevaleurs)<-"Auto"
  tclvalue(Env$l.var$log.axevaleurs)<-0
  tclvalue(Env$l.var$boxmoy)<-0
  tclvalue(Env$l.var$ICmediane)<-0
  tclvalue(Env$l.var$varwidth)<-0
  tclvalue(Env$l.var$lg.moustaches)<-"1.5"
  Env$l.var$noms1<-""
  tclvalue(Env$l.var$moyprop)<-"moy"
  tclvalue(Env$l.var$facteur2)<-""
  tclvalue(Env$l.var$proportions)<-""
  tclvalue(Env$l.var$prop.niveaux)<-""
  Env$l.var$nomsprop<-""
  tclvalue(Env$l.var$moyprop)<-"moy"
  tclvalue(Env$l.var$facteurprop)<-""
  Env$l.var$couleur1B<-"white"
  Env$l.var$col.borduresB<-"black"
  tclvalue(Env$l.var$hachuresA)<-"1"
  Env$l.var$hachuresB<-1
  tclvalue(Env$l.var$stack)<-0
  tclvalue(Env$l.var$erreur)<-""
  tclvalue(Env$l.var$legende)<-0
  tclvalue(Env$l.var$legende.pos)<-""
  tclvalue(Env$l.var$legende.titre)<-""
  Env$l.var$noms2<-""
  tclvalue(Env$l.var$plusieurs)<-0
  Env$l.var$nomsprop.fac<-""
  tclvalue(Env$l.var$cam.orient)<-""
  tclvalue(Env$l.var$cam.start)<-"0"
  tclvalue(Env$l.var$cam.lien)<-0
  tclvalue(Env$l.var$parts.niveaux)<-""
  Env$l.var$nomsparts<-""
  tclvalue(Env$l.var$varX)<-""
  tclvalue(Env$l.var$varY)<-""
  tclvalue(Env$l.var$varX.prop)<-""
  Env$l.var$couleur2B<-"grey"
  tclvalue(Env$l.var$taille.ptsA)<-"1"
  Env$l.var$taille.ptsB<-1
  tclvalue(Env$l.var$type.courbeA)<-""
  Env$l.var$type.courbeB<-""
  Env$l.var$trait2<-""
  Env$l.var$epaisseur2<-1
  tclvalue(Env$l.var$symboleA)<-"1"
  Env$l.var$symboleB<-1
  Env$l.var$select<-1
  tclvalue(Env$l.var$droiteA)<-""
  Env$l.var$droiteB<-""
  tclvalue(Env$l.var$intervalA)<-""
  Env$l.var$intervalB<-""
  tclvalue(Env$l.var$ptlab)<-0
  tclvalue(Env$l.var$sysinfo)<-0
  Env$l.var$levels.temp <- NULL
  tclvalue(Env$l.var$nobar)<-0
  Env$l.var$add.abscisses<-NULL
  Env$l.var$add.hauteurs<-NULL
  Env$l.var$add.seq<-NULL
  Env$l.var$add.seq2<-NULL
  Env$l.var$add.matrice<-NULL
}


#-------------------------------------------------
# Renommer une classe (boîte, barre, part)
#-------------------------------------------------

rename.noms1<-function(value.list,value.nom) {
  if(nchar(value.nom)>0) {
    if (nchar(value.list)>0) {
	Env$l.var$noms1[as.numeric(value.list)+1]<-value.nom
	if ("noms.list"%in%names(Env$l.fr3)) {
	  tkdelete(Env$l.fr3$noms.list,0,"end")
	  for (i in 1:length(Env$l.var$noms1)) {tkinsert(Env$l.fr3$noms.list,"end",Env$l.var$noms1[i])}
	  tkdelete(Env$l.fr3$noms.wdg,0,"end")
	}
    } else {
	msg(text=Env$voc[25,1],type="error")
    }
  } else {
    msg(text=Env$voc[24,1],type="error")
  }
}

rename.nomsprop.fac<-function(value.list,value.nom) {
  if(nchar(value.nom)>0) {
    if (nchar(value.list)>0) {
	Env$l.var$nomsprop.fac[as.numeric(value.list)+1]<-value.nom
	tkdelete(Env$l.fr3$noms.list,0,"end")
	for (i in 1:length(Env$l.var$nomsprop.fac)) {tkinsert(Env$l.fr3$noms.list,"end",Env$l.var$nomsprop.fac[i])}
	tkdelete(Env$l.fr3$noms.wdg,0,"end")
    } else {
	msg(text=Env$voc[25,1],type="error")
    }
  } else {
    msg(text=Env$voc[24,1],type="error")
  }
}

rename.nomsparts<-function() {
  if(nchar(tclvalue(tkget(Env$l.fr3$noms.wdg)))>0) {
    if (nchar(tclvalue(tkcurselection(Env$l.fr3$noms.list)))>0) {
	Env$l.var$nomsparts[as.numeric(tclvalue(tkcurselection(Env$l.fr3$noms.list)))+1]<-tclvalue(tkget(Env$l.fr3$noms.wdg))
	tkdelete(Env$l.fr3$noms.list,0,"end")
	for (i in 1:length(Env$l.var$nomsparts)) {tkinsert(Env$l.fr3$noms.list,"end",Env$l.var$nomsparts[i])}
	tkdelete(Env$l.fr3$noms.wdg,0,"end")
	if ("noms.list"%in%names(Env$l.fr4)) {
	  tkdelete(Env$l.fr4$noms.list,0,"end")
	  for (i in 1:length(Env$l.var$nomsparts)) {tkinsert(Env$l.fr4$noms.list,"end",Env$l.var$nomsparts[i])}
	}
    } else {
	msg(text=Env$voc[25,1],type="error")
    }
  } else {
    msg(text=Env$voc[24,1],type="error")
  }
}


#-------------------------------------------------
# Boîtes - couleur des boîtes
#-------------------------------------------------

col.boites<-function() {
  if (tclvalue(Env$l.var$plusieurs)==0) {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$couleur1A),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$couleur1A)<-temp
	tkconfigure(Env$l.fr4$colboites.wdg,bg=tclvalue(Env$l.var$couleur1A))
    }
  } else {
    if (nchar(tclvalue(tkcurselection(Env$l.fr4$noms.list)))>0) {
	temp<-tclvalue(tcl("tk_chooseColor",initialcolor=Env$l.var$couleur1B[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1],title=Env$voc[64,1]))
	if (nchar(temp)>0) {
	  Env$l.var$couleur1B[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1]<-temp
	  tkconfigure(Env$l.fr4$colboites.wdg,bg=Env$l.var$couleur1B[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1])
	}
    } else {
	msg(text=Env$voc[25,1],type="error")
    }
  }
}


#-------------------------------------------------
# Barplot (et boîtes) - couleur des barres et
# bordures (barres et boîtes) hachures
#-------------------------------------------------

col.barres<-function() {
  if (tclvalue(Env$l.var$plusieurs)==0) {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$couleur1A),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$couleur1A)<-temp
	tkconfigure(Env$l.fr4$colbarres.wdg,bg=tclvalue(Env$l.var$couleur1A))
    }
  } else {
    if (nchar(tclvalue(tkcurselection(Env$l.fr4$noms.list)))>0) {
	temp<-tclvalue(tcl("tk_chooseColor",initialcolor=Env$l.var$couleur1B[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1],title=Env$voc[64,1]))
	if (nchar(temp)>0) {
	  Env$l.var$couleur1B[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1]<-temp
	  tkconfigure(Env$l.fr4$colbarres.wdg,bg=Env$l.var$couleur1B[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1])
	}
    } else {
	msg(text=Env$voc[25,1],type="error")
    }
  }
}

col.bordures<-function() {
  if (tclvalue(Env$l.var$plusieurs)==0) {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$col.borduresA),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$col.borduresA)<-temp
	tkconfigure(Env$l.fr4$colbordures.wdg,bg=tclvalue(Env$l.var$col.borduresA))
    }
  } else {
    if (nchar(tclvalue(tkcurselection(Env$l.fr4$noms.list)))>0) {
	temp<-tclvalue(tcl("tk_chooseColor",initialcolor=Env$l.var$col.borduresB[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1],title=Env$voc[64,1]))
	if (nchar(temp)>0) {
	  Env$l.var$col.borduresB[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1]<-temp
	  tkconfigure(Env$l.fr4$colbordures.wdg,bg=Env$l.var$col.borduresB[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1])
	}
    } else {
	msg(text=Env$voc[25,1],type="error")
    }
  }
}

hachures<-function(num) {
  if (tclvalue(Env$l.var$plusieurs)==0) {
    tclvalue(Env$l.var$hachuresA)<-as.character(num)
    for (i in 1:9) {tkconfigure(Env$l.fr4$l.hachures[[i]],borderwidth=0)}
    tkconfigure(Env$l.fr4$l.hachures[[num]],borderwidth=2)
  } else {
    if (nchar(tclvalue(tkcurselection(Env$l.fr4$noms.list)))>0) {
	Env$l.var$hachuresB[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1]<-num
	for (i in 1:9) {tkconfigure(Env$l.fr4$l.hachures[[i]],borderwidth=0)}
	tkconfigure(Env$l.fr4$l.hachures[[num]],borderwidth=2)
    } else {
	msg(text=Env$voc[25,1],type="error")
    }
  }
}


#-------------------------------------------------
# Activer / désactiver la frame Barres d'erreur
#-------------------------------------------------

active.erreur<-function() {
  if ("type.lab"%in%names(Env$l.fr5)) {
    if (tclvalue(Env$l.var$stack)==1) {
	tkconfigure(Env$l.fr5$type.lab,foreground="grey")
	tkconfigure(Env$l.fr5$type.wdg,state="disabled")
	tkconfigure(Env$l.fr5$col.lab,foreground="grey")
	tkconfigure(Env$l.fr5$col.wdg,bg="grey")
    } else {
	tkconfigure(Env$l.fr5$type.lab,foreground="black")
	tkconfigure(Env$l.fr5$type.wdg,state="readonly")
	tkconfigure(Env$l.fr5$col.lab,foreground="black")
	tkconfigure(Env$l.fr5$col.wdg,bg=tclvalue(Env$l.var$couleur2A))
    }
  }
}


#-------------------------------------------------
# Activer / désactiver la frame Légende
#-------------------------------------------------

active.legende<-function() {
  if (tclvalue(Env$l.var$plusieurs)==0) {
    if ("legende.lab"%in%names(Env$l.fr5)) {
	tkconfigure(Env$l.fr5$legende.lab,foreground="grey")
	tkconfigure(Env$l.fr5$legende.wdg,state="disabled")
	tkconfigure(Env$l.fr5$titre.lab,foreground="grey")
	tkconfigure(Env$l.fr5$titre.wdg,state="disabled")
	tkconfigure(Env$l.fr5$position.lab,foreground="grey")
	tkconfigure(Env$l.fr5$position.wdg,state="disabled")
	tkconfigure(Env$l.fr5$noms.lab1,foreground="grey")
	tkdelete(Env$l.fr5$noms.list,0,"end")
	tkconfigure(Env$l.fr5$noms.list,state="disabled")
	tkconfigure(Env$l.fr5$noms.lab2,foreground="grey")
	tkdelete(Env$l.fr5$noms.wdg,0,"end")
	tkconfigure(Env$l.fr5$noms.wdg,state="disabled")
    }
    if ("legende.lab"%in%names(Env$l.fr6)) {
	tkconfigure(Env$l.fr6$legende.lab,foreground="grey")
	tkconfigure(Env$l.fr6$legende.wdg,state="disabled")
	tkconfigure(Env$l.fr6$titre.lab,foreground="grey")
	tkconfigure(Env$l.fr6$titre.wdg,state="disabled")
	tkconfigure(Env$l.fr6$position.lab,foreground="grey")
	tkconfigure(Env$l.fr6$position.wdg,state="disabled")
	tkconfigure(Env$l.fr6$noms.lab1,foreground="grey")
	tkdelete(Env$l.fr6$noms.list,0,"end")
	tkconfigure(Env$l.fr6$noms.list,state="disabled")
	tkconfigure(Env$l.fr6$noms.lab2,foreground="grey")
	tkdelete(Env$l.fr6$noms.wdg,0,"end")
	tkconfigure(Env$l.fr6$noms.wdg,state="disabled")
    }
  } else {
    if ("legende.lab"%in%names(Env$l.fr5)) {
	tkconfigure(Env$l.fr5$legende.lab,foreground="black")
	tkconfigure(Env$l.fr5$legende.wdg,state="normal")
	tkconfigure(Env$l.fr5$titre.lab,foreground="black")
	tkconfigure(Env$l.fr5$titre.wdg,state="normal")
	tkconfigure(Env$l.fr5$position.lab,foreground="black")
	tkconfigure(Env$l.fr5$position.wdg,state="readonly")
	tkconfigure(Env$l.fr5$noms.lab1,foreground="black")
	tkconfigure(Env$l.fr5$noms.list,state="normal")
	tkconfigure(Env$l.fr5$noms.lab2,foreground="black")
	tkconfigure(Env$l.fr5$noms.wdg,state="normal")
    }
    if ("legende.lab"%in%names(Env$l.fr6)) {
	tkconfigure(Env$l.fr6$legende.lab,foreground="black")
	tkconfigure(Env$l.fr6$legende.wdg,state="normal")
	tkconfigure(Env$l.fr6$titre.lab,foreground="black")
	tkconfigure(Env$l.fr6$titre.wdg,state="normal")
	tkconfigure(Env$l.fr6$position.lab,foreground="black")
	tkconfigure(Env$l.fr6$position.wdg,state="readonly")
	tkconfigure(Env$l.fr6$noms.lab1,foreground="black")
	tkconfigure(Env$l.fr6$noms.list,state="normal")
	tkconfigure(Env$l.fr6$noms.lab2,foreground="black")
	tkconfigure(Env$l.fr6$noms.wdg,state="normal")
    }
  }
}


#-------------------------------------------------
# Boîtes - renommer un item de légende
#-------------------------------------------------

rename.legende3<-function(value.list,value.nom) {
  if(nchar(value.nom)>0) {
    if (nchar(value.list)>0) {
	Env$l.var$noms2[as.numeric(value.list)+1]<-value.nom
	if ("noms.list"%in%names(Env$l.fr4)) {
	  tkdelete(Env$l.fr4$noms.list,0,"end")
	  for (i in 1:length(Env$l.var$noms1)) {tkinsert(Env$l.fr4$noms.list,"end",Env$l.var$noms2[i])}
	}
	tkdelete(Env$l.fr6$noms.list,0,"end")
	for (i in 1:length(Env$l.var$noms1)) {tkinsert(Env$l.fr6$noms.list,"end",Env$l.var$noms2[i])}
	tkdelete(Env$l.fr6$noms.wdg,0,"end")
    } else {
	msg(text=Env$voc[25,1],type="error")
    }
  } else {
    msg(text=Env$voc[24,1],type="error")
  }
}


#-------------------------------------------------
# Barres - renommer un item de légende
#-------------------------------------------------

rename.legende<-function(value.list,value.nom) {
  if(nchar(value.nom)>0) {
    if (nchar(value.list)>0) {
	if (tclvalue(Env$l.var$moyprop)=="moy") {
	  Env$l.var$noms2[as.numeric(value.list)+1]<-value.nom
	  if ("noms.list"%in%names(Env$l.fr6)) {
	    tkdelete(Env$l.fr6$noms.list,0,"end")
	    for (i in 1:length(Env$l.var$noms2)) {tkinsert(Env$l.fr6$noms.list,"end",Env$l.var$noms2[i])}
	    tkdelete(Env$l.fr6$noms.wdg,0,"end")
	  }
	  if ("noms.list"%in%names(Env$l.fr4)) {
	    tkdelete(Env$l.fr4$noms.list,0,"end")
	    for (i in 1:length(Env$l.var$noms2)) {tkinsert(Env$l.fr4$noms.list,"end",Env$l.var$noms2[i])}
	  }
	} else {
	  Env$l.var$nomsprop[as.numeric(value.list)+1]<-value.nom
	  if ("noms.list"%in%names(Env$l.fr4)) {
	    tkdelete(Env$l.fr4$noms.list,0,"end")
	    for (i in 1:length(Env$l.var$nomsprop)) {tkinsert(Env$l.fr4$noms.list,"end",Env$l.var$nomsprop[i])}
	  }
	  if ("noms.list"%in%names(Env$l.fr6)) {
	    tkdelete(Env$l.fr6$noms.list,0,"end")
	    for (i in 1:length(Env$l.var$nomsprop)) {tkinsert(Env$l.fr6$noms.list,"end",Env$l.var$nomsprop[i])}
	    tkdelete(Env$l.fr6$noms.wdg,0,"end")
	  }
	}
    } else {
	msg(text=Env$voc[25,1],type="error")
    }
  } else {
    msg(text=Env$voc[24,1],type="error")
  }
}


#-------------------------------------------------
# Courbe et nuage - renommer un item de légende
#-------------------------------------------------

rename.legende2<-function(value.list,value.nom) {
  if(nchar(value.nom)>0) {
    if (nchar(value.list)>0) {
	Env$l.var$noms1[as.numeric(value.list)+1]<-value.nom
	if ("noms.list"%in%names(Env$l.fr4)) {
	  tkdelete(Env$l.fr4$noms.list,0,"end")
	  for (i in 1:length(Env$l.var$noms1)) {tkinsert(Env$l.fr4$noms.list,"end",Env$l.var$noms1[i])}
	}
	if ("noms.list"%in%names(Env$l.fr5)) {
	  tkdelete(Env$l.fr5$noms.list,0,"end")
	  for (i in 1:length(Env$l.var$noms1)) {tkinsert(Env$l.fr5$noms.list,"end",Env$l.var$noms1[i])}
	}
	tkdelete(Env$l.fr6$noms.list,0,"end")
	for (i in 1:length(Env$l.var$noms1)) {tkinsert(Env$l.fr6$noms.list,"end",Env$l.var$noms1[i])}
	tkdelete(Env$l.fr6$noms.wdg,0,"end")
    } else {
	msg(text=Env$voc[25,1],type="error")
    }
  } else {
    msg(text=Env$voc[24,1],type="error")
  }
}


#-------------------------------------------------
# Barres - type moyennes ou proportions
#-------------------------------------------------

barres.moy<-function() {
  tkconfigure(Env$l.fr1$titre1,foreground="black")
  tkconfigure(Env$l.fr1$moyvar.lab,foreground="black")
  tkconfigure(Env$l.fr1$moyvar.wdg,state="readonly")
  tkconfigure(Env$l.fr1$moyfac1.lab,foreground="black")
  tkconfigure(Env$l.fr1$moyfac1.wdg,state="readonly")
  tkconfigure(Env$l.fr1$moyfac2.lab,foreground="black")
  tkconfigure(Env$l.fr1$moyfac2.wdg,state="readonly")
  tkconfigure(Env$l.fr1$titre2,foreground="grey")
  tkconfigure(Env$l.fr1$propvar.lab,foreground="grey")
  tkconfigure(Env$l.fr1$propvar.wdg,state="disabled")
  tkconfigure(Env$l.fr1$propnivx.lab,foreground="grey")
  tkconfigure(Env$l.fr1$propnivx.list,state="disabled")
  tkconfigure(Env$l.fr1$propfac.lab,foreground="grey")
  tkconfigure(Env$l.fr1$propfac.wdg,state="disabled")
  tclvalue(Env$l.var$stack)<-0
  if ("noms.list"%in%names(Env$l.fr3)) {
    tkdelete(Env$l.fr3$noms.list,0,"end")
    for (i in 1:length(Env$l.var$noms1)) {tkinsert(Env$l.fr3$noms.list,"end",Env$l.var$noms1[i])}
    tkdelete(Env$l.fr3$noms.wdg,0,"end")
  }
  if ("type.wdg"%in%names(Env$l.fr5)) {
    tclvalue(Env$l.var$erreur)<-""
    tkconfigure(Env$l.fr5$type.wdg,values=Env$voc[c(95:98),1])
  }
  if (nchar(tclvalue(tkget(Env$l.fr1$moyfac2.wdg)))>0 & tclvalue(tkget(Env$l.fr1$moyfac2.wdg))!=Env$voc[82,1]) {
    tclvalue(Env$l.var$plusieurs)<-1
    Env$l.var$couleur1B<-grey.colors(nlevels(Env$dataset[,tclvalue(Env$l.var$facteur2)]))
    Env$l.var$col.borduresB<-rep("black",nlevels(Env$dataset[,tclvalue(Env$l.var$facteur2)]))
    Env$l.var$hachuresB<-rep(1,nlevels(Env$dataset[,tclvalue(Env$l.var$facteur2)]))
    if ("noms.list"%in%names(Env$l.fr4)) {
	tkconfigure(Env$l.fr4$noms.list,state="normal")
	tkdelete(Env$l.fr4$noms.list,0,"end")
	for (i in 1:length(Env$l.var$noms2)) {tkinsert(Env$l.fr4$noms.list,"end",Env$l.var$noms2[i])}
	tkconfigure(Env$l.fr4$colbarres.wdg,bg=Env$l.var$couleur1B[1])
	tkconfigure(Env$l.fr4$colbordures.wdg,bg=Env$l.var$col.borduresB[1])
	for (i in 1:9) {tkconfigure(Env$l.fr4$l.hachures[[i]],borderwidth=0)}
	tkconfigure(Env$l.fr4$l.hachures[[Env$l.var$hachuresB[1]]],borderwidth=2)
	tkconfigure(Env$l.fr4$stack.lab,foreground="black")
	tkconfigure(Env$l.fr4$stack.wdg,state="normal")
	tkdeselect(Env$l.fr4$stack.wdg)
    }
    active.legende()
    if ("noms.list"%in%names(Env$l.fr6)) {
	tkdelete(Env$l.fr6$noms.list,0,"end")
	for (i in 1:length(Env$l.var$noms2)) {tkinsert(Env$l.fr6$noms.list,"end",Env$l.var$noms2[i])}
    }
  } else {
    tclvalue(Env$l.var$plusieurs)<-0
    tclvalue(Env$l.var$couleur1A)<-"grey"
    tclvalue(Env$l.var$col.borduresA)<-"black"
    tclvalue(Env$l.var$hachuresA)<-"1"
    if ("noms.list"%in%names(Env$l.fr4)) {
	tkconfigure(Env$l.fr4$noms.list,state="normal")
	tkdelete(Env$l.fr4$noms.list,0,"end")
	tkconfigure(Env$l.fr4$noms.list,state="disabled")
	tkconfigure(Env$l.fr4$colbarres.wdg,bg=tclvalue(Env$l.var$couleur1A))
	tkconfigure(Env$l.fr4$colbordures.wdg,bg=tclvalue(Env$l.var$col.borduresA))
	for (i in 1:9) {tkconfigure(Env$l.fr4$l.hachures[[i]],borderwidth=0)}
	tkconfigure(Env$l.fr4$l.hachures[[as.numeric(tclvalue(Env$l.var$hachuresA))]],borderwidth=2)
	tkconfigure(Env$l.fr4$stack.lab,foreground="grey")
	tkdeselect(Env$l.fr4$stack.wdg)
	tkconfigure(Env$l.fr4$stack.wdg,state="disabled")
    }
    active.legende()
  }
  active.erreur()
}

barres.prop<-function() {
  tkconfigure(Env$l.fr1$titre1,foreground="grey")
  tkconfigure(Env$l.fr1$moyvar.lab,foreground="grey")
  tkconfigure(Env$l.fr1$moyvar.wdg,state="disabled")
  tkconfigure(Env$l.fr1$moyfac1.lab,foreground="grey")
  tkconfigure(Env$l.fr1$moyfac1.wdg,state="disabled")
  tkconfigure(Env$l.fr1$moyfac2.lab,foreground="grey")
  tkconfigure(Env$l.fr1$moyfac2.wdg,state="disabled")
  tkconfigure(Env$l.fr1$titre2,foreground="black")
  tkconfigure(Env$l.fr1$propvar.lab,foreground="black")
  tkconfigure(Env$l.fr1$propvar.wdg,state="readonly")
  tkconfigure(Env$l.fr1$propnivx.lab,foreground="black")
  tkconfigure(Env$l.fr1$propnivx.list,state="normal")
  tkconfigure(Env$l.fr1$propfac.lab,foreground="black")
  tkconfigure(Env$l.fr1$propfac.wdg,state="readonly")
  tclvalue(Env$l.var$stack)<-0
  if (nchar(Env$l.var$nomsprop)[1]>1) {
    tclvalue(Env$l.var$plusieurs)<-1
    Env$l.var$couleur1B<-grey.colors(length(Env$l.var$nomsprop))
    Env$l.var$col.borduresB<-rep("black",length(Env$l.var$nomsprop)) 
    Env$l.var$hachuresB<-rep(1,length(Env$l.var$nomsprop))
  } else {
    tclvalue(Env$l.var$plusieurs)<-0
    tclvalue(Env$l.var$couleur1A)<-"grey"
    tclvalue(Env$l.var$col.borduresA)<-"black"
    tclvalue(Env$l.var$hachuresA)<-"1"
  }
  if ("noms.list"%in%names(Env$l.fr3)) {
    tkdelete(Env$l.fr3$noms.list,0,"end")
    for (i in 1:length(Env$l.var$nomsprop.fac)) {tkinsert(Env$l.fr3$noms.list,"end",Env$l.var$nomsprop.fac[i])}
    tkdelete(Env$l.fr3$noms.wdg,0,"end")
  }
  if ("noms.list"%in%names(Env$l.fr4)) {
    tkconfigure(Env$l.fr4$noms.list,state="normal")
    tkdelete(Env$l.fr4$noms.list,0,"end")
    if (tclvalue(Env$l.var$plusieurs)==1) {
	for (i in 1:length(Env$l.var$nomsprop)) {tkinsert(Env$l.fr4$noms.list,"end",Env$l.var$nomsprop[i])}
	tkconfigure(Env$l.fr4$colbarres.wdg,bg=Env$l.var$couleur1B[1])
	tkconfigure(Env$l.fr4$colbordures.wdg,bg=Env$l.var$col.borduresB[1])
	for (i in 1:9) {tkconfigure(Env$l.fr4$l.hachures[[i]],borderwidth=0)}
	tkconfigure(Env$l.fr4$l.hachures[[Env$l.var$hachuresB[1]]],borderwidth=2)
	tkconfigure(Env$l.fr4$stack.lab,foreground="black")
	tkconfigure(Env$l.fr4$stack.wdg,state="normal")
	tkdeselect(Env$l.fr4$stack.wdg)
    } else {
	tkconfigure(Env$l.fr4$noms.list,state="disabled")
	tkconfigure(Env$l.fr4$colbarres.wdg,bg=tclvalue(Env$l.var$couleur1A))
	tkconfigure(Env$l.fr4$colbordures.wdg,bg=tclvalue(Env$l.var$col.borduresA))
	for (i in 1:9) {tkconfigure(Env$l.fr4$l.hachures[[i]],borderwidth=0)}
	tkconfigure(Env$l.fr4$l.hachures[[as.numeric(tclvalue(Env$l.var$hachuresA))]],borderwidth=2)
	tkconfigure(Env$l.fr4$stack.lab,foreground="grey")
	tkdeselect(Env$l.fr4$stack.wdg)
	tkconfigure(Env$l.fr4$stack.wdg,state="disabled")
    }
  }
  if ("type.wdg"%in%names(Env$l.fr5)) {
    tclvalue(Env$l.var$erreur)<-""
    tkconfigure(Env$l.fr5$type.wdg,values=Env$voc[c(95,97,98),1])
  }
  if ("noms.list"%in%names(Env$l.fr6)) {
    if (tclvalue(Env$l.var$plusieurs)==1) {
	tkconfigure(Env$l.fr6$noms.list,state="normal")
	tkdelete(Env$l.fr6$noms.list,0,"end")
	for (i in 1:length(Env$l.var$nomsprop)) {tkinsert(Env$l.fr6$noms.list,"end",Env$l.var$nomsprop[i])}
    } else {
	tkconfigure(Env$l.fr6$noms.list,state="normal")
	tkdelete(Env$l.fr6$noms.list,0,"end")
    }
  }
  active.legende()
  active.erreur()
}


#-------------------------------------------------
# Camembert - couleur des barres et hachures,
# activation/désactivation de la légende
#-------------------------------------------------

col.parts<-function() {
  if (nchar(tclvalue(tkcurselection(Env$l.fr4$noms.list)))>0) {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=Env$l.var$couleur1B[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1],title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	Env$l.var$couleur1B[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1]<-temp
	tkconfigure(Env$l.fr4$colparts.wdg,bg=Env$l.var$couleur1B[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1])
    }
  } else {
    msg(text=Env$voc[25,1],type="error")
  }
}

hachures2<-function(num) {
  if (nchar(tclvalue(tkcurselection(Env$l.fr4$noms.list)))>0) {
    Env$l.var$hachuresB[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1]<-num
    for (i in 1:9) {tkconfigure(Env$l.fr4$l.hachures[[i]],borderwidth=0)}
    tkconfigure(Env$l.fr4$l.hachures[[num]],borderwidth=2)
  } else {
    msg(text=Env$voc[25,1],type="error")
  }
}

active.legende2<-function() {
  if (tclvalue(Env$l.var$cam.lien)==1) {
    if ("legende.lab"%in%names(Env$l.fr5)) {
	tkconfigure(Env$l.fr5$legende.lab,foreground="grey")
	tkconfigure(Env$l.fr5$legende.wdg,state="disabled")
	tkconfigure(Env$l.fr5$titre.lab,foreground="grey")
	tkconfigure(Env$l.fr5$titre.wdg,state="disabled")
	tkconfigure(Env$l.fr5$position.lab,foreground="grey")
	tkconfigure(Env$l.fr5$position.wdg,state="disabled")
    }
  } else {
    if ("legende.lab"%in%names(Env$l.fr5)) {
	tkconfigure(Env$l.fr5$legende.lab,foreground="black")
	tkconfigure(Env$l.fr5$legende.wdg,state="normal")
	tkconfigure(Env$l.fr5$titre.lab,foreground="black")
	tkconfigure(Env$l.fr5$titre.wdg,state="normal")
	tkconfigure(Env$l.fr5$position.lab,foreground="black")
	tkconfigure(Env$l.fr5$position.wdg,state="readonly")
    }
  }
}


#-------------------------------------------------
# Courbe - type moyennes ou proportions
#-------------------------------------------------

courbe.moy<-function() {
  tkconfigure(Env$l.fr1$titre1,foreground="black")
  tkconfigure(Env$l.fr1$moyvarX.lab,foreground="black")
  tkconfigure(Env$l.fr1$moyvarX.wdg,state="readonly")
  tkconfigure(Env$l.fr1$moyvarY.lab,foreground="black")
  tkconfigure(Env$l.fr1$moyvarY.wdg,state="readonly")
  tkconfigure(Env$l.fr1$titre2,foreground="grey")
  tkconfigure(Env$l.fr1$propvarX.lab,foreground="grey")
  tkconfigure(Env$l.fr1$propvarX.wdg,state="disabled")
  tkconfigure(Env$l.fr1$propvarY.lab,foreground="grey")
  tkconfigure(Env$l.fr1$propvarY.wdg,state="disabled")
  tkconfigure(Env$l.fr1$propvarY.niv.lab,foreground="grey")
  tkconfigure(Env$l.fr1$propvarY.niv.wdg,state="disabled")
  if ("type.wdg"%in%names(Env$l.fr5)) {
    tclvalue(Env$l.var$erreur)<-""
    tkconfigure(Env$l.fr5$type.wdg,values=Env$voc[c(95:98),1])
  }
}

courbe.prop<-function() {
  tkconfigure(Env$l.fr1$titre1,foreground="grey")
  tkconfigure(Env$l.fr1$moyvarX.lab,foreground="grey")
  tkconfigure(Env$l.fr1$moyvarX.wdg,state="disabled")
  tkconfigure(Env$l.fr1$moyvarY.lab,foreground="grey")
  tkconfigure(Env$l.fr1$moyvarY.wdg,state="disabled")
  tkconfigure(Env$l.fr1$titre2,foreground="black")
  tkconfigure(Env$l.fr1$propvarX.lab,foreground="black")
  tkconfigure(Env$l.fr1$propvarX.wdg,state="readonly")
  tkconfigure(Env$l.fr1$propvarY.lab,foreground="black")
  tkconfigure(Env$l.fr1$propvarY.wdg,state="readonly")
  tkconfigure(Env$l.fr1$propvarY.niv.lab,foreground="black")
  tkconfigure(Env$l.fr1$propvarY.niv.wdg,state="readonly")
  if ("type.wdg"%in%names(Env$l.fr5)) {
    tclvalue(Env$l.var$erreur)<-""
    tkconfigure(Env$l.fr5$type.wdg,values=Env$voc[c(95,97,98),1])
  }
}


#-------------------------------------------------
# Courbe et nuage - symboles et couleur des points
#-------------------------------------------------

symboles<-function(num) {
  if (tclvalue(Env$l.var$plusieurs)==0) {
    tclvalue(Env$l.var$symboleA)<-as.character(num)
    for (i in 1:8) {tkconfigure(Env$l.fr4$l.symboles[[i]],borderwidth=0)}
    tkconfigure(Env$l.fr4$l.symboles[[num]],borderwidth=2)
  } else {
    if (nchar(tclvalue(tkcurselection(Env$l.fr4$noms.list)))>0) {
	Env$l.var$symboleB[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1]<-num
	for (i in 1:8) {tkconfigure(Env$l.fr4$l.symboles[[i]],borderwidth=0)}
	tkconfigure(Env$l.fr4$l.symboles[[num]],borderwidth=2)
    } else {
	msg(text=Env$voc[25,1],type="error")
    }
  }
}

col.symboles<-function() {
  if (tclvalue(Env$l.var$plusieurs)==0) {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$couleur2A),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$couleur2A)<-temp
	tkconfigure(Env$l.fr4$col.wdg,bg=tclvalue(Env$l.var$couleur2A))
    }
  } else {
    if (nchar(tclvalue(tkcurselection(Env$l.fr4$noms.list)))>0) {
	temp<-tclvalue(tcl("tk_chooseColor",initialcolor=Env$l.var$couleur2B[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1],title=Env$voc[64,1]))
	if (nchar(temp)>0) {
	  Env$l.var$couleur2B[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1]<-temp
	  tkconfigure(Env$l.fr4$col.wdg,bg=Env$l.var$couleur2B[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1])
	}
    } else {
	msg(text=Env$voc[25,1],type="error")
    }
  }
}


#-------------------------------------------------
# Frame 1
#-------------------------------------------------

fr1.close<-function() {
  Env$l.frames$Fr1.status<-0
  tkconfigure(Env$l.wdg$but.lab1,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_bas.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr1)) {tkdestroy(Env$l.fr1[[i]])}
  Env$l.fr1<-list()
  Env$l.fr1$vide<-tklabel(Env$l.frames$Fr1,text="",font=Env$police2)
  tkgrid(Env$l.fr1$vide)
}

fr1.openD<-function() {
  Env$l.frames$Fr1.status<-1
  tkconfigure(Env$l.wdg$but.lab1,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr1)) {tkdestroy(Env$l.fr1[[i]])}
  Env$l.fr1<-list()
  Env$l.fr1$titre1<-tklabel(Env$l.frames$Fr1,text=Env$voc[2,1],font=Env$police3)
  Env$l.fr1$ext.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[3,1],font=Env$police)
  Env$l.fr1$ext.wdg<-ttkcombobox(Env$l.frames$Fr1,values=c("txt","csv"),textvariable=Env$l.var$extension,font=Env$police,state="readonly")
  tkbind(Env$l.fr1$ext.wdg,"<<ComboboxSelected>>",function() {
    if (tclvalue(Env$l.var$extension)=="txt") {tclvalue(Env$l.var$sepcol)<-Env$voc[5,1]} else
	if (tclvalue(Env$l.var$extension)=="csv") {tclvalue(Env$l.var$sepcol)<-Env$voc[6,1]}
  })
  Env$l.fr1$col.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[4,1],font=Env$police)
  Env$l.fr1$col.wdg<-ttkcombobox(Env$l.frames$Fr1,values=Env$voc[5:7,1],textvariable=Env$l.var$sepcol,font=Env$police,state="readonly")
  Env$l.fr1$dec.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[8,1],font=Env$police)
  Env$l.fr1$dec.wdg<-ttkcombobox(Env$l.frames$Fr1,values=Env$voc[c(9,6),1],textvariable=Env$l.var$sepdec,font=Env$police,state="readonly")
  Env$l.fr1$na.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[10,1],font=Env$police)
  Env$l.fr1$na.wdg<-tkentry(Env$l.frames$Fr1,width=4,textvariable=Env$l.var$na,font=Env$police)
  Env$l.fr1$hdr.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[11,1],font=Env$police)
  Env$l.fr1$hdr.wdg<-tkcheckbutton(Env$l.frames$Fr1,variable=Env$l.var$header)
  Env$l.fr1$but1<-tkbutton(Env$l.frames$Fr1,text=Env$voc[13,1],font=Env$police,width=16,command=data.load1)
  Env$l.fr1$titre2<-tklabel(Env$l.frames$Fr1,text=Env$voc[12,1],font=Env$police3)
  Env$l.fr1$obj.list<-tklistbox(Env$l.frames$Fr1,height=7,font=Env$police,selectmode="single",yscrollcommand=function(...) tkset(Env$l.fr1$obj.scroll,...))
  Env$l.fr1$obj.scroll<-tkscrollbar(Env$l.frames$Fr1,repeatinterval=5,command=function(...) tkyview(Env$l.fr1$obj.list,...))
  tables<-NULL
  if (length(ls(.GlobalEnv))>0) {
    for (i in 1:length(ls(.GlobalEnv))) {
	if (is.data.frame(get(ls(.GlobalEnv)[i]))) {tables<-c(tables,ls(.GlobalEnv)[i])}
    }
  }
  tkdelete(Env$l.fr1$obj.list,0,"end")
  if (!is.null(tables)) {
    for (i in 1:length(tables)) {tkinsert(Env$l.fr1$obj.list,"end",tables[i])}
  } else {}
  tkbind(Env$l.fr1$obj.list,"<Enter>",function() {msg(text=Env$voc[14,1],type="info")})
  tkbind(Env$l.fr1$obj.list,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr1$but2<-tkbutton(Env$l.frames$Fr1,text=Env$voc[13,1],font=Env$police,width=16,command=data.load2)
  Env$l.fr1$espace.ver1<-tklabel(Env$l.frames$Fr1,text="",font=Env$police2)
  Env$l.fr1$espace.ver2<-tklabel(Env$l.frames$Fr1,text="",font=Env$police2)
  Env$l.fr1$espace.hor<-tklabel(Env$l.frames$Fr1,text="                                                  ",font=Env$police)
  tkgrid(Env$l.fr1$titre1,row=0,column=0,columnspan=2)
  tkgrid(Env$l.fr1$espace.ver1,row=1,column=0)
  tkgrid(Env$l.fr1$ext.lab,row=2,column=0,sticky="e")
  tkgrid(Env$l.fr1$ext.wdg,row=2,column=1,sticky="w")
  tkgrid(Env$l.fr1$col.lab,row=3,column=0,sticky="e")
  tkgrid(Env$l.fr1$col.wdg,row=3,column=1,sticky="w")
  tkgrid(Env$l.fr1$dec.lab,row=4,column=0,sticky="e")
  tkgrid(Env$l.fr1$dec.wdg,row=4,column=1,sticky="w")
  tkgrid(Env$l.fr1$na.lab,row=5,column=0,sticky="e")
  tkgrid(Env$l.fr1$na.wdg,row=5,column=1,sticky="w")
  tkgrid(Env$l.fr1$hdr.lab,row=6,column=0,sticky="e")
  tkgrid(Env$l.fr1$hdr.wdg,row=6,column=1,sticky="w")
  tkgrid(Env$l.fr1$espace.ver2,row=7,column=0)
  tkgrid(Env$l.fr1$but1,row=8,column=0,columnspan=2)
  tkgrid(Env$l.fr1$espace.hor,row=0,column=3)
  tkgrid(Env$l.fr1$titre2,row=0,column=4)
  tkgrid(Env$l.fr1$obj.list,Env$l.fr1$obj.scroll,row=2,column=4,rowspan=5);tkgrid.configure(Env$l.fr1$obj.scroll,sticky="ens")
  tkgrid(Env$l.fr1$but2,row=8,column=4)
}

fr1.openH<-function() {
  Env$l.frames$Fr1.status<-1
  tkconfigure(Env$l.wdg$but.lab1,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr1)) {tkdestroy(Env$l.fr1[[i]])}
  Env$l.fr1<-list()
  Env$l.fr1$variable.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[36,1],font=Env$police)
  Env$l.fr1$variable.wdg<-ttkcombobox(Env$l.frames$Fr1,values=Env$l.var$var.num,textvariable=Env$l.var$variable,font=Env$police,state="readonly")
  Env$l.fr1$facteur.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[37,1],font=Env$police)
  Env$l.fr1$facteur.wdg<-ttkcombobox(Env$l.frames$Fr1,values=c(Env$voc[82,1],Env$l.var$var.fact),textvariable=Env$l.var$facteur1,font=Env$police,state="readonly")
  tkbind(Env$l.fr1$facteur.wdg,"<<ComboboxSelected>>",function() {
    if (nchar(tclvalue(Env$l.var$facteur1))>0 & tclvalue(Env$l.var$facteur1)!=Env$voc[82,1]) {
	tkconfigure(Env$l.fr1$niveau.lab,foreground="black")
	tkconfigure(Env$l.fr1$niveau.wdg,state="readonly")
	tkdelete(Env$l.fr1$niveau.wdg,0,"end")
	tkconfigure(Env$l.fr1$niveau.wdg,values=levels(Env$dataset[,tclvalue(Env$l.var$facteur1)]))
	tclvalue(Env$l.var$niveau)<-levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[1]
    } else {
	tkconfigure(Env$l.fr1$niveau.lab,foreground="grey")
	tkdelete(Env$l.fr1$niveau.wdg,0,"end")
	tclvalue(Env$l.var$niveau)<-""
	tkconfigure(Env$l.fr1$niveau.wdg,state="disabled")
    }
  })
  Env$l.fr1$niveau.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[38,1],font=Env$police,foreground=ifelse(nchar(tclvalue(Env$l.var$facteur1))>0 & tclvalue(Env$l.var$facteur1)!=Env$voc[82,1],"black","grey"))
  Env$l.fr1$niveau.wdg<-ttkcombobox(Env$l.frames$Fr1,values="",textvariable=Env$l.var$niveau,font=Env$police,state=ifelse(nchar(tclvalue(Env$l.var$facteur1))>0 & tclvalue(Env$l.var$facteur1)!=Env$voc[82,1],"readonly","disabled"))
  Env$l.fr1$sysinfo.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[248,1],font=Env$police)
  Env$l.fr1$sysinfo.wdg<-tkcheckbutton(Env$l.frames$Fr1,variable=Env$l.var$sysinfo)
  Env$l.fr1$type.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[39,1],font=Env$police)
  Env$l.fr1$type.wdg<-ttkcombobox(Env$l.frames$Fr1,values=Env$voc[40:42,1],textvariable=Env$l.var$hist.type,font=Env$police,state="readonly")
  if (nchar(tclvalue(Env$l.var$facteur1))>0 & tclvalue(Env$l.var$facteur1)!=Env$voc[82,1]) {
    tkconfigure(Env$l.fr1$niveau.wdg,values=levels(Env$dataset[,tclvalue(Env$l.var$facteur1)]))
  }
  tkbind(Env$l.fr1$type.wdg,"<Enter>",function() {msg(text=Env$voc[114,1],type="warning")})
  tkbind(Env$l.fr1$type.wdg,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.fr1$type.wdg,"<<ComboboxSelected>>",function() {
    if ("hor.liminf.wdg"%in%names(Env$l.fr3)) {
	if (tclvalue(Env$l.var$hist.type)==Env$voc[40,1]) {
	  tkconfigure(Env$l.fr3$hor.liminf.wdg,state="disabled")
	  tkconfigure(Env$l.fr3$hor.limsup.wdg,state="disabled")
	} else {
	  tkconfigure(Env$l.fr3$hor.liminf.wdg,state="normal")
	  tkconfigure(Env$l.fr3$hor.limsup.wdg,state="normal")
	}
    }
    if ("tracer.lab"%in%names(Env$l.fr5)) {
	if (!tclvalue(Env$l.var$hist.type)==Env$voc[42,1]) {
	  tkconfigure(Env$l.fr5$tracer.lab,foreground="grey")
	  tkconfigure(Env$l.fr5$tracer.wdg,state="disabled")
	  tkconfigure(Env$l.fr5$col.lab,foreground="grey")
	  tkconfigure(Env$l.fr5$col.wdg,bg="grey")
	  tkconfigure(Env$l.fr5$trait.lab,foreground="grey")
	  tkconfigure(Env$l.fr5$trait.wdg,state="disabled")
	  tkconfigure(Env$l.fr5$epaisseur.lab,foreground="grey")
	  tkconfigure(Env$l.fr5$epaisseur.wdg,state="disabled")
	} else {
	  tkconfigure(Env$l.fr5$tracer.lab,foreground="black")
	  tkconfigure(Env$l.fr5$tracer.wdg,state="normal")
	  tkconfigure(Env$l.fr5$col.lab,foreground="black")
	  tkconfigure(Env$l.fr5$col.wdg,bg=tclvalue(Env$l.var$couleur2A))
	  tkconfigure(Env$l.fr5$trait.lab,foreground="black")
	  tkconfigure(Env$l.fr5$trait.wdg,state="readonly")
	  tkconfigure(Env$l.fr5$epaisseur.lab,foreground="black")
	  tkconfigure(Env$l.fr5$epaisseur.wdg,state="normal")
	}
    }
  })
  Env$l.fr1$encadre.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[43,1],font=Env$police)
  Env$l.fr1$encadre.wdg<-tkcheckbutton(Env$l.frames$Fr1,variable=Env$l.var$encadre)
  Env$l.fr1$espace.hor<-tklabel(Env$l.frames$Fr1,text="                                   ",font=Env$police)
  tkgrid(Env$l.fr1$variable.lab,row=0,column=0,sticky="e")
  tkgrid(Env$l.fr1$variable.wdg,row=0,column=1,sticky="w")
  tkgrid(Env$l.fr1$facteur.lab,row=1,column=0,sticky="e")
  tkgrid(Env$l.fr1$facteur.wdg,row=1,column=1,sticky="w")
  tkgrid(Env$l.fr1$niveau.lab,row=2,column=0,sticky="e")
  tkgrid(Env$l.fr1$niveau.wdg,row=2,column=1,sticky="w")
  tkgrid(Env$l.fr1$espace.hor,row=0,column=2)
  tkgrid(Env$l.fr1$type.lab,row=0,column=3,sticky="e")
  tkgrid(Env$l.fr1$type.wdg,row=0,column=4,sticky="w")
  tkgrid(Env$l.fr1$encadre.lab,row=1,column=3,sticky="e")
  tkgrid(Env$l.fr1$encadre.wdg,row=1,column=4,sticky="w")
  tkgrid(Env$l.fr1$sysinfo.lab,row=2,column=3,sticky="e")
  tkgrid(Env$l.fr1$sysinfo.wdg,row=2,column=4,sticky="w")
}

fr1.openM<-function() {
  Env$l.frames$Fr1.status<-1
  tkconfigure(Env$l.wdg$but.lab1,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr1)) {tkdestroy(Env$l.fr1[[i]])}
  Env$l.fr1<-list()
  Env$l.fr1$variable.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[36,1],font=Env$police)
  Env$l.fr1$variable.wdg<-ttkcombobox(Env$l.frames$Fr1,values=Env$l.var$var.num,textvariable=Env$l.var$variable,font=Env$police,state="readonly")
  Env$l.fr1$facteur1.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[86,1],font=Env$police)
  Env$l.fr1$facteur1.wdg<-ttkcombobox(Env$l.frames$Fr1,values=Env$l.var$var.fact,textvariable=Env$l.var$facteur1,font=Env$police,state="readonly")
  tkbind(Env$l.fr1$facteur1.wdg,"<<ComboboxSelected>>",function() {
    if (nchar(tclvalue(Env$l.var$facteur2))>0 & tclvalue(Env$l.var$facteur2)!=Env$voc[82,1]) {
	Env$l.var$facteur.interaction<-interaction(Env$dataset[,tclvalue(Env$l.var$facteur2)],Env$dataset[,tclvalue(Env$l.var$facteur1)])
    } else {
	Env$l.var$facteur.interaction<-""
    }
    Env$l.var$noms1<-levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])
    if ("noms.list"%in%names(Env$l.fr3)) {
	tkdelete(Env$l.fr3$noms.list,0,"end")
	for (i in 1:length(Env$l.var$noms1)) {tkinsert(Env$l.fr3$noms.list,"end",Env$l.var$noms1[i])}
    }
  })
  Env$l.fr1$facteur2.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[87,1],font=Env$police)
  Env$l.fr1$facteur2.wdg<-ttkcombobox(Env$l.frames$Fr1,values=c(Env$voc[82,1],Env$l.var$var.fact),textvariable=Env$l.var$facteur2,font=Env$police,state="readonly")
  tkbind(Env$l.fr1$facteur2.wdg,"<<ComboboxSelected>>",function() {
    if (nchar(tclvalue(Env$l.var$facteur2))>0 & tclvalue(Env$l.var$facteur2)!=Env$voc[82,1]) {
	Env$l.var$facteur.interaction<-interaction(Env$dataset[,tclvalue(Env$l.var$facteur2)],Env$dataset[,tclvalue(Env$l.var$facteur1)])
	Env$l.var$noms2<-levels(Env$dataset[,tclvalue(Env$l.var$facteur2)])
	tclvalue(Env$l.var$plusieurs)<-1
	Env$l.var$couleur1B<-grey.colors(nlevels(Env$dataset[,tclvalue(Env$l.var$facteur2)]))
	Env$l.var$col.borduresB<-rep("black",nlevels(Env$dataset[,tclvalue(Env$l.var$facteur2)]))
	if ("noms.list"%in%names(Env$l.fr4)) {
	  tkconfigure(Env$l.fr4$noms.list,state="normal")
	  tkdelete(Env$l.fr4$noms.list,0,"end")
	  for (i in 1:length(Env$l.var$noms2)) {tkinsert(Env$l.fr4$noms.list,"end",Env$l.var$noms2[i])}
	  tkconfigure(Env$l.fr4$colboites.wdg,bg=Env$l.var$couleur1B[1])
	  tkconfigure(Env$l.fr4$colbordures.wdg,bg=Env$l.var$col.borduresB[1])
	}
	active.legende()
	if ("noms.list"%in%names(Env$l.fr6)) {
	  tkdelete(Env$l.fr6$noms.list,0,"end")
	  for (i in 1:length(Env$l.var$noms2)) {tkinsert(Env$l.fr6$noms.list,"end",Env$l.var$noms2[i])}
	}
    } else {
	Env$l.var$facteur.interaction<-""
	Env$l.var$noms2<-""
	tclvalue(Env$l.var$plusieurs)<-0
	tclvalue(Env$l.var$couleur1A)<-"grey"
	tclvalue(Env$l.var$col.borduresA)<-"black"
	if ("noms.list"%in%names(Env$l.fr4)) {
	  tkconfigure(Env$l.fr4$noms.list,state="normal")
	  tkdelete(Env$l.fr4$noms.list,0,"end")
	  tkconfigure(Env$l.fr4$noms.list,state="disabled")
	  tkconfigure(Env$l.fr4$colboites.wdg,bg=tclvalue(Env$l.var$couleur1A))
	  tkconfigure(Env$l.fr4$colbordures.wdg,bg=tclvalue(Env$l.var$col.borduresA))
	}
	active.legende()
	if ("noms.list"%in%names(Env$l.fr6)) {
	  tkdelete(Env$l.fr6$noms.list,0,"end")
	}
    }
  })
  Env$l.fr1$sysinfo.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[248,1],font=Env$police)
  Env$l.fr1$sysinfo.wdg<-tkcheckbutton(Env$l.frames$Fr1,variable=Env$l.var$sysinfo)
  Env$l.fr1$orientation.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[66,1],font=Env$police)
  Env$l.fr1$orientation.wdg<-ttkcombobox(Env$l.frames$Fr1,values=Env$voc[67:68,1],textvariable=Env$l.var$box.orient,font=Env$police,state="readonly")
  Env$l.fr1$encadre.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[43,1],font=Env$police)
  Env$l.fr1$encadre.wdg<-tkcheckbutton(Env$l.frames$Fr1,variable=Env$l.var$encadre)
  Env$l.fr1$espace.hor<-tklabel(Env$l.frames$Fr1,text="                                   ",font=Env$police)
  tkgrid(Env$l.fr1$variable.lab,row=0,column=0,sticky="e")
  tkgrid(Env$l.fr1$variable.wdg,row=0,column=1,sticky="w")
  tkgrid(Env$l.fr1$facteur1.lab,row=1,column=0,sticky="e")
  tkgrid(Env$l.fr1$facteur1.wdg,row=1,column=1,sticky="w")
  tkgrid(Env$l.fr1$facteur2.lab,row=2,column=0,sticky="e")
  tkgrid(Env$l.fr1$facteur2.wdg,row=2,column=1,sticky="w")
  tkgrid(Env$l.fr1$espace.hor,row=0,column=2)
  tkgrid(Env$l.fr1$orientation.lab,row=0,column=3,sticky="e")
  tkgrid(Env$l.fr1$orientation.wdg,row=0,column=4,sticky="w")
  tkgrid(Env$l.fr1$encadre.lab,row=1,column=3,sticky="e")
  tkgrid(Env$l.fr1$encadre.wdg,row=1,column=4,sticky="w")
  tkgrid(Env$l.fr1$sysinfo.lab,row=2,column=3,sticky="e")
  tkgrid(Env$l.fr1$sysinfo.wdg,row=2,column=4,sticky="w")
}

fr1.openB<-function() {
  Env$l.frames$Fr1.status<-1
  tkconfigure(Env$l.wdg$but.lab1,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr1)) {tkdestroy(Env$l.fr1[[i]])}
  Env$l.fr1<-list()
  Env$l.fr1$titre1<-tklabel(Env$l.frames$Fr1,text=Env$voc[84,1],font=Env$police3)
  Env$l.fr1$rb1<-tkradiobutton(Env$l.frames$Fr1,variable=Env$l.var$moyprop,value="moy",command=barres.moy)
  Env$l.fr1$moyvar.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[36,1],font=Env$police)
  Env$l.fr1$moyvar.wdg<-ttkcombobox(Env$l.frames$Fr1,values=Env$l.var$var.num,textvariable=Env$l.var$variable,font=Env$police,state="readonly")
  Env$l.fr1$moyfac1.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[86,1],font=Env$police)
  Env$l.fr1$moyfac1.wdg<-ttkcombobox(Env$l.frames$Fr1,values=Env$l.var$var.fact,textvariable=Env$l.var$facteur1,font=Env$police,state="readonly")
  tkbind(Env$l.fr1$moyfac1.wdg,"<<ComboboxSelected>>",function() {
    if (nchar(tclvalue(Env$l.var$facteur1))>0) {
	Env$l.var$noms1<-levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])
	if ("noms.list"%in%names(Env$l.fr3)) {
	  tkdelete(Env$l.fr3$noms.list,0,"end")
	  for (i in 1:length(Env$l.var$noms1)) {tkinsert(Env$l.fr3$noms.list,"end",Env$l.var$noms1[i])}
	}
    }
  })
  Env$l.fr1$moyfac2.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[87,1],font=Env$police)
  Env$l.fr1$moyfac2.wdg<-ttkcombobox(Env$l.frames$Fr1,values=c(Env$voc[82,1],Env$l.var$var.fact),textvariable=Env$l.var$facteur2,font=Env$police,state="readonly")
  tkbind(Env$l.fr1$moyfac2.wdg,"<<ComboboxSelected>>",function() {
    tclvalue(Env$l.var$stack)<-0
    if (nchar(tclvalue(Env$l.var$facteur2))>0 & tclvalue(Env$l.var$facteur2)!=Env$voc[82,1]) {
	Env$l.var$noms2<-levels(Env$dataset[,tclvalue(Env$l.var$facteur2)])
	tclvalue(Env$l.var$plusieurs)<-1
	Env$l.var$couleur1B<-grey.colors(nlevels(Env$dataset[,tclvalue(Env$l.var$facteur2)]))
	Env$l.var$col.borduresB<-rep("black",nlevels(Env$dataset[,tclvalue(Env$l.var$facteur2)]))
	Env$l.var$hachuresB<-rep(1,nlevels(Env$dataset[,tclvalue(Env$l.var$facteur2)]))
	if ("noms.list"%in%names(Env$l.fr4)) {
	  tkconfigure(Env$l.fr4$noms.list,state="normal")
	  tkdelete(Env$l.fr4$noms.list,0,"end")
	  for (i in 1:length(Env$l.var$noms2)) {tkinsert(Env$l.fr4$noms.list,"end",Env$l.var$noms2[i])}
	  tkconfigure(Env$l.fr4$colbarres.wdg,bg=Env$l.var$couleur1B[1])
	  tkconfigure(Env$l.fr4$colbordures.wdg,bg=Env$l.var$col.borduresB[1])
	  for (i in 1:9) {tkconfigure(Env$l.fr4$l.hachures[[i]],borderwidth=0)}
	  tkconfigure(Env$l.fr4$l.hachures[[Env$l.var$hachuresB[1]]],borderwidth=2)
	  tkconfigure(Env$l.fr4$stack.lab,foreground="black")
	  tkconfigure(Env$l.fr4$stack.wdg,state="normal")
	  tkdeselect(Env$l.fr4$stack.wdg)
	}
	active.legende()
	if ("noms.list"%in%names(Env$l.fr6)) {
	  tkdelete(Env$l.fr6$noms.list,0,"end")
	  for (i in 1:length(Env$l.var$noms2)) {tkinsert(Env$l.fr6$noms.list,"end",Env$l.var$noms2[i])}
	}
    } else {
	Env$l.var$noms2<-""
	tclvalue(Env$l.var$plusieurs)<-0
	tclvalue(Env$l.var$couleur1A)<-"grey"
	tclvalue(Env$l.var$col.borduresA)<-"black"
	tclvalue(Env$l.var$hachuresA)<-"1"
	if ("noms.list"%in%names(Env$l.fr4)) {
	  tkconfigure(Env$l.fr4$noms.list,state="normal")
	  tkdelete(Env$l.fr4$noms.list,0,"end")
	  tkconfigure(Env$l.fr4$noms.list,state="disabled")
	  tkconfigure(Env$l.fr4$colbarres.wdg,bg=tclvalue(Env$l.var$couleur1A))
	  tkconfigure(Env$l.fr4$colbordures.wdg,bg=tclvalue(Env$l.var$col.borduresA))
	  for (i in 1:9) {tkconfigure(Env$l.fr4$l.hachures[[i]],borderwidth=0)}
	  tkconfigure(Env$l.fr4$l.hachures[[as.numeric(tclvalue(Env$l.var$hachuresA))]],borderwidth=2)
	  tkconfigure(Env$l.fr4$stack.lab,foreground="grey")
	  tkdeselect(Env$l.fr4$stack.wdg)
	  tkconfigure(Env$l.fr4$stack.wdg,state="disabled")
	}
	active.legende()
	if ("noms.list"%in%names(Env$l.fr6)) {
	  tkdelete(Env$l.fr6$noms.list,0,"end")
	}
    }
    active.erreur()
  })
  Env$l.fr1$nobar.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[264,1],font=Env$police)
  Env$l.fr1$nobar.wdg<-tkcheckbutton(Env$l.frames$Fr1,variable=Env$l.var$nobar)
  Env$l.fr1$encadre.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[43,1],font=Env$police)
  Env$l.fr1$encadre.wdg<-tkcheckbutton(Env$l.frames$Fr1,variable=Env$l.var$encadre)
  Env$l.fr1$sysinfo.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[248,1],font=Env$police)
  Env$l.fr1$sysinfo.wdg<-tkcheckbutton(Env$l.frames$Fr1,variable=Env$l.var$sysinfo)
  Env$l.fr1$titre2<-tklabel(Env$l.frames$Fr1,text=Env$voc[85,1],font=Env$police3)
  Env$l.fr1$rb2<-tkradiobutton(Env$l.frames$Fr1,variable=Env$l.var$moyprop,value="prop",command=barres.prop)
  Env$l.fr1$propvar.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[36,1],font=Env$police)
  Env$l.fr1$propvar.wdg<-ttkcombobox(Env$l.frames$Fr1,values=Env$l.var$var.fact,textvariable=Env$l.var$proportions,font=Env$police,state="readonly")
  tkbind(Env$l.fr1$propvar.wdg,"<<ComboboxSelected>>",function() {
    if (nchar(tclvalue(Env$l.var$proportions))>0) {
	tclvalue(Env$l.var$plusieurs)<-0
	tclvalue(Env$l.var$prop.niveaux)<-"0"
	Env$l.var$nomsprop<-""
	tkdelete(Env$l.fr1$propnivx.list,0,"end")
	for (i in 1:nlevels(Env$dataset[,tclvalue(Env$l.var$proportions)])) {tkinsert(Env$l.fr1$propnivx.list,"end",levels(Env$dataset[,tclvalue(Env$l.var$proportions)])[i])}
	tkselection.set(Env$l.fr1$propnivx.list,tclvalue(Env$l.var$prop.niveaux))
	if ("noms.list"%in%names(Env$l.fr4)) {
	  tkdelete(Env$l.fr4$noms.list,0,"end")
	  tkconfigure(Env$l.fr4$noms.list,state="disabled")
	}
	active.legende()
    }
  })
  Env$l.fr1$propnivx.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[88,1],font=Env$police)
  Env$l.fr1$propnivx.list<-tklistbox(Env$l.frames$Fr1,height=6,font=Env$police,selectmode="multiple",yscrollcommand=function(...) tkset(Env$l.fr1$propnivx.scroll,...))
  Env$l.fr1$propnivx.scroll<-tkscrollbar(Env$l.frames$Fr1,repeatinterval=5,command=function(...) tkyview(Env$l.fr1$propnivx.list,...))
  tkbind(Env$l.fr1$propnivx.list,"<Enter>",function() {if(tclvalue(Env$l.var$moyprop)=="prop") {msg(text=Env$voc[142,1],type="info")}})
  tkbind(Env$l.fr1$propnivx.list,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.fr1$propnivx.list,"<ButtonRelease-1>",function() {
    tclvalue(Env$l.var$stack)<-0
    if (nchar(tclvalue(tkcurselection(Env$l.fr1$propnivx.list)))>2) {
	tclvalue(Env$l.var$plusieurs)<-1
	tclvalue(Env$l.var$prop.niveaux)<-tclvalue(tkcurselection(Env$l.fr1$propnivx.list))
	num<-as.numeric(strsplit(tclvalue(tkcurselection(Env$l.fr1$propnivx.list)),split=" ")[[1]])+1
	Env$l.var$couleur1B<-grey.colors(length(num))
	Env$l.var$col.borduresB<-rep("black",length(num))
	Env$l.var$hachuresB<-rep(1,length(num))
	Env$l.var$nomsprop<-levels(Env$dataset[,tclvalue(Env$l.var$proportions)])[num]
	if ("noms.list"%in%names(Env$l.fr4)) {
	  tkconfigure(Env$l.fr4$noms.list,state="normal")
	  tkdelete(Env$l.fr4$noms.list,0,"end")
	  for (i in 1:length(Env$l.var$nomsprop)) {tkinsert(Env$l.fr4$noms.list,"end",Env$l.var$nomsprop[i])}
	  tkconfigure(Env$l.fr4$colbarres.wdg,bg=Env$l.var$couleur1B[1])
	  tkconfigure(Env$l.fr4$colbordures.wdg,bg=Env$l.var$col.borduresB[1])
	  for (i in 1:9) {tkconfigure(Env$l.fr4$l.hachures[[i]],borderwidth=0)}
	  tkconfigure(Env$l.fr4$l.hachures[[Env$l.var$hachuresB[1]]],borderwidth=2)
	  tkconfigure(Env$l.fr4$stack.lab,foreground="black")
	  tkconfigure(Env$l.fr4$stack.wdg,state="normal")
	  tkdeselect(Env$l.fr4$stack.wdg)
	}
	active.legende()
	if ("noms.list"%in%names(Env$l.fr6)) {
	  tkconfigure(Env$l.fr6$noms.list,state="normal")
	  tkdelete(Env$l.fr6$noms.list,0,"end")
	  for (i in 1:length(Env$l.var$nomsprop)) {tkinsert(Env$l.fr6$noms.list,"end",Env$l.var$nomsprop[i])}
	}
    } else {
	tclvalue(Env$l.var$plusieurs)<-0
	tclvalue(Env$l.var$prop.niveaux)<-tclvalue(tkcurselection(Env$l.fr1$propnivx.list))
	tclvalue(Env$l.var$couleur1A)<-"grey"
	tclvalue(Env$l.var$col.borduresA)<-"black"
	tclvalue(Env$l.var$hachuresA)<-"1"
	Env$l.var$nomsprop<-""
	if ("noms.list"%in%names(Env$l.fr4)) {
	  tkconfigure(Env$l.fr4$noms.list,state="normal")
	  tkdelete(Env$l.fr4$noms.list,0,"end")
	  tkconfigure(Env$l.fr4$noms.list,state="disabled")
	  tkconfigure(Env$l.fr4$colbarres.wdg,bg=tclvalue(Env$l.var$couleur1A))
	  tkconfigure(Env$l.fr4$colbordures.wdg,bg=tclvalue(Env$l.var$col.borduresA))
	  for (i in 1:9) {tkconfigure(Env$l.fr4$l.hachures[[i]],borderwidth=0)}
	  tkconfigure(Env$l.fr4$l.hachures[[as.numeric(tclvalue(Env$l.var$hachuresA))]],borderwidth=2)
	  tkconfigure(Env$l.fr4$stack.lab,foreground="grey")
	  tkdeselect(Env$l.fr4$stack.wdg)
	  tkconfigure(Env$l.fr4$stack.wdg,state="disabled")
	}
	active.legende()
	if ("noms.list"%in%names(Env$l.fr6)) {
	  tkconfigure(Env$l.fr6$noms.list,state="normal")
	  tkdelete(Env$l.fr6$noms.list,0,"end")
	  tkconfigure(Env$l.fr6$noms.list,state="disabled")
	}
    }
    active.erreur()
  })
  Env$l.fr1$propfac.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[65,1],font=Env$police)
  Env$l.fr1$propfac.wdg<-ttkcombobox(Env$l.frames$Fr1,values=Env$l.var$var.fact,textvariable=Env$l.var$facteurprop,font=Env$police,state="readonly")
  if (tclvalue(Env$l.var$moyprop)=="moy") {
    barres.moy()
  } else {
    barres.prop()
  }
  tkbind(Env$l.fr1$propfac.wdg,"<<ComboboxSelected>>",function() {
    if (nchar(tclvalue(Env$l.var$facteurprop))>0) {
	Env$l.var$nomsprop.fac<-levels(Env$dataset[,tclvalue(Env$l.var$facteurprop)])
	if ("noms.list"%in%names(Env$l.fr3)) {
	  tkdelete(Env$l.fr3$noms.list,0,"end")
	  for (i in 1:length(Env$l.var$nomsprop.fac)) {tkinsert(Env$l.fr3$noms.list,"end",Env$l.var$nomsprop.fac[i])}
	}
    } else {
	if ("noms.list"%in%names(Env$l.fr3)) {
	  tkdelete(Env$l.fr3$noms.list,0,"end")
	}
    }
    if ("noms.wdg"%in%names(Env$l.fr3)) {tkdelete(Env$l.fr3$noms.wdg,0,"end")}
  })
  Env$l.fr1$espace.ver<-tklabel(Env$l.frames$Fr1,text="",font=Env$police2)
  Env$l.fr1$espace.hor<-tklabel(Env$l.frames$Fr1,text="                    ",font=Env$police)
  tkgrid(Env$l.fr1$titre1,row=0,column=0,sticky="e")
  tkgrid(Env$l.fr1$rb1,row=0,column=1,sticky="w")
  tkgrid(Env$l.fr1$espace.ver,row=1,column=0)
  tkgrid(Env$l.fr1$moyvar.lab,row=2,column=0,sticky="e")
  tkgrid(Env$l.fr1$moyvar.wdg,row=2,column=1,sticky="w")
  tkgrid(Env$l.fr1$moyfac1.lab,row=3,column=0,sticky="e")
  tkgrid(Env$l.fr1$moyfac1.wdg,row=3,column=1,sticky="w")
  tkgrid(Env$l.fr1$moyfac2.lab,row=4,column=0,sticky="e")
  tkgrid(Env$l.fr1$moyfac2.wdg,row=4,column=1,sticky="w")
  tkgrid(Env$l.fr1$nobar.lab,row=5,column=0,sticky="e")
  tkgrid(Env$l.fr1$nobar.wdg,row=5,column=1,sticky="w")
  tkgrid(Env$l.fr1$encadre.lab,row=6,column=0,sticky="e")
  tkgrid(Env$l.fr1$encadre.wdg,row=6,column=1,sticky="w")
  tkgrid(Env$l.fr1$sysinfo.lab,row=7,column=0,sticky="e")
  tkgrid(Env$l.fr1$sysinfo.wdg,row=7,column=1,sticky="w")
  tkgrid(Env$l.fr1$espace.hor,row=0,column=2)
  tkgrid(Env$l.fr1$rb2,row=0,column=3,sticky="e")
  tkgrid(Env$l.fr1$titre2,row=0,column=4,sticky="w")
  tkgrid(Env$l.fr1$propvar.lab,row=2,column=3,sticky="e")
  tkgrid(Env$l.fr1$propvar.wdg,row=2,column=4,sticky="w")
  tkgrid(Env$l.fr1$propnivx.lab,row=3,column=3,sticky="e")
  tkgrid(Env$l.fr1$propnivx.list,Env$l.fr1$propnivx.scroll,row=3,column=4,rowspan=4,sticky="w");tkgrid.configure(Env$l.fr1$propnivx.scroll,sticky="ens")
  tkgrid(Env$l.fr1$propfac.lab,row=7,column=3,sticky="e")
  tkgrid(Env$l.fr1$propfac.wdg,row=7,column=4,sticky="w")
}

fr1.openCa<-function() {
  Env$l.frames$Fr1.status<-1
  tkconfigure(Env$l.wdg$but.lab1,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr1)) {tkdestroy(Env$l.fr1[[i]])}
  Env$l.fr1<-list()
  Env$l.fr1$var.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[36,1],font=Env$police)
  Env$l.fr1$var.wdg<-ttkcombobox(Env$l.frames$Fr1,values=c(Env$l.var$var.num,Env$l.var$var.fact),textvariable=Env$l.var$variable,font=Env$police,state="readonly")
  tkbind(Env$l.fr1$var.wdg,"<<ComboboxSelected>>",function() {
    if (is.factor(Env$dataset[,tclvalue(Env$l.var$variable)])) {
	Env$l.var$noms1<-levels(Env$dataset[,tclvalue(Env$l.var$variable)])
    } else {
	Env$l.var$noms1<-levels(factor(Env$dataset[,tclvalue(Env$l.var$variable)]))
    }
    tclvalue(Env$l.var$plusieurs)<-1
    Env$l.var$nomsparts<-Env$l.var$noms1
    Env$l.var$couleur1B<-grey.colors(length(Env$l.var$nomsparts))
    Env$l.var$hachuresB<-rep(1,length(Env$l.var$nomsparts))
    tclvalue(Env$l.var$parts.niveaux)<-paste(0:(length(Env$l.var$noms1)-1),collapse=" ")
    tkdelete(Env$l.fr1$parts.list,0,"end")
    for (i in 1:length(Env$l.var$noms1)) {tkinsert(Env$l.fr1$parts.list,"end",Env$l.var$noms1[i])}
    if ("noms.list"%in%names(Env$l.fr3)) {
	tkdelete(Env$l.fr3$noms.list,0,"end")
	for (i in 1:length(Env$l.var$nomsparts)) {tkinsert(Env$l.fr3$noms.list,"end",Env$l.var$nomsparts[i])}
	tkdelete(Env$l.fr3$noms.wdg,0,"end")
    }
    if ("noms.list"%in%names(Env$l.fr4)) {
	tkdelete(Env$l.fr4$noms.list,0,"end")
	for (i in 1:length(Env$l.var$nomsparts)) {tkinsert(Env$l.fr4$noms.list,"end",Env$l.var$nomsparts[i])}
	tkconfigure(Env$l.fr4$colparts.wdg,bg=Env$l.var$couleur1B[1])
	for (i in 1:9) {tkconfigure(Env$l.fr4$l.hachures[[i]],borderwidth=0)}
	tkconfigure(Env$l.fr4$l.hachures[[Env$l.var$hachuresB[1]]],borderwidth=2)
    }
  })
  Env$l.fr1$sysinfo.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[248,1],font=Env$police)
  Env$l.fr1$sysinfo.wdg<-tkcheckbutton(Env$l.frames$Fr1,variable=Env$l.var$sysinfo)
  Env$l.fr1$parts.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[115,1],font=Env$police)
  Env$l.fr1$parts.list<-tklistbox(Env$l.frames$Fr1,height=7,font=Env$police,selectmode="multiple",yscrollcommand=function(...) tkset(Env$l.fr1$parts.scroll,...))
  Env$l.fr1$parts.scroll<-tkscrollbar(Env$l.frames$Fr1,repeatinterval=5,command=function(...) tkyview(Env$l.fr1$parts.list,...))
  tkbind(Env$l.fr1$parts.list,"<Enter>",function() {msg(text=Env$voc[159,1],type="info")})
  tkbind(Env$l.fr1$parts.list,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.fr1$parts.list,"<ButtonRelease-1>",function() {
    tclvalue(Env$l.var$parts.niveaux)<-tclvalue(tkcurselection(Env$l.fr1$parts.list))
    if (nchar(tclvalue(tkcurselection(Env$l.fr1$parts.list)))>2) {
	tclvalue(Env$l.var$plusieurs)<-1
    } else {
	tclvalue(Env$l.var$plusieurs)<-0
    }
    Env$l.var$nomsparts<-Env$l.var$noms1[as.numeric(strsplit(tclvalue(Env$l.var$parts.niveaux),split=" ")[[1]])+1]
    Env$l.var$couleur1B<-grey.colors(length(Env$l.var$nomsparts))
    Env$l.var$hachuresB<-rep(1,length(Env$l.var$nomsparts))
    if ("noms.list"%in%names(Env$l.fr3)) {
	tkdelete(Env$l.fr3$noms.list,0,"end")
	for (i in 1:length(Env$l.var$nomsparts)) {tkinsert(Env$l.fr3$noms.list,"end",Env$l.var$nomsparts[i])}
	tkdelete(Env$l.fr3$noms.wdg,0,"end")
    }
    if ("noms.list"%in%names(Env$l.fr4)) {
	tkdelete(Env$l.fr4$noms.list,0,"end")
	for (i in 1:length(Env$l.var$nomsparts)) {tkinsert(Env$l.fr4$noms.list,"end",Env$l.var$nomsparts[i])}
	tkconfigure(Env$l.fr4$colparts.wdg,bg=Env$l.var$couleur1B[1])
	for (i in 1:9) {tkconfigure(Env$l.fr4$l.hachures[[i]],borderwidth=0)}
	tkconfigure(Env$l.fr4$l.hachures[[Env$l.var$hachuresB[1]]],borderwidth=2)
    }
  })
  Env$l.fr1$espace.hor<-tklabel(Env$l.frames$Fr1,text="                                        ",font=Env$police)
  tkgrid(Env$l.fr1$var.lab,row=0,column=0,sticky="e")
  tkgrid(Env$l.fr1$var.wdg,row=0,column=1,sticky="w")
  tkgrid(Env$l.fr1$sysinfo.lab,row=1,column=0,sticky="e")
  tkgrid(Env$l.fr1$sysinfo.wdg,row=1,column=1,sticky="w")
  tkgrid(Env$l.fr1$espace.hor,row=0,column=2)
  tkgrid(Env$l.fr1$parts.lab,row=0,column=3,sticky="e")
  tkgrid(Env$l.fr1$parts.list,Env$l.fr1$parts.scroll,row=0,column=4,rowspan=7,sticky="w");tkgrid.configure(Env$l.fr1$parts.scroll,sticky="ens")
}

fr1.openCo<-function() {
  Env$l.frames$Fr1.status<-1
  tkconfigure(Env$l.wdg$but.lab1,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr1)) {tkdestroy(Env$l.fr1[[i]])}
  Env$l.fr1<-list()
  Env$l.fr1$titre1<-tklabel(Env$l.frames$Fr1,text=Env$voc[84,1],font=Env$police3)
  Env$l.fr1$rb1<-tkradiobutton(Env$l.frames$Fr1,variable=Env$l.var$moyprop,value="moy",command=courbe.moy)
  Env$l.fr1$moyvarX.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[127,1],font=Env$police)
  Env$l.fr1$moyvarX.wdg<-ttkcombobox(Env$l.frames$Fr1,width=15,values=Env$l.var$var.num,textvariable=Env$l.var$varX,font=Env$police,state="readonly")
  Env$l.fr1$moyvarY.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[128,1],font=Env$police)
  Env$l.fr1$moyvarY.wdg<-ttkcombobox(Env$l.frames$Fr1,width=15,values=Env$l.var$var.num,textvariable=Env$l.var$varY,font=Env$police,state="readonly")
  Env$l.fr1$sysinfo.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[248,1],font=Env$police)
  Env$l.fr1$sysinfo.wdg<-tkcheckbutton(Env$l.frames$Fr1,variable=Env$l.var$sysinfo)
  Env$l.fr1$fact.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[37,1],font=Env$police)
  Env$l.fr1$fact.wdg<-ttkcombobox(Env$l.frames$Fr1,width=15,values=c(Env$voc[82,1],Env$l.var$var.fact),textvariable=Env$l.var$facteur1,font=Env$police,state="readonly")
  tkbind(Env$l.fr1$fact.wdg,"<<ComboboxSelected>>",function() {
    if (nchar(tclvalue(Env$l.var$facteur1))>0) {
	tclvalue(Env$l.var$plusieurs)<-0
	tclvalue(Env$l.var$symboleA)<-"1"
	tclvalue(Env$l.var$couleur2A)<-"black"
	tclvalue(Env$l.var$taille.ptsA)<-"1"
	tclvalue(Env$l.var$type.courbeA)<-""
	tclvalue(Env$l.var$trait1)<-""
	tclvalue(Env$l.var$epaisseur1)<-"1"
	tkdelete(Env$l.fr1$fact.list,0,"end")
	if ("noms.list"%in%names(Env$l.fr4)) {
	  tkdelete(Env$l.fr4$noms.list,0,"end")
	  tkconfigure(Env$l.fr4$noms.list,state="disabled")
	  for (i in 1:8) {tkconfigure(Env$l.fr4$l.symboles[[i]],borderwidth=0)}
	  tkconfigure(Env$l.fr4$l.symboles[[as.numeric(tclvalue(Env$l.var$symboleA))]],borderwidth=2)
	  tkconfigure(Env$l.fr4$col.wdg,bg=tclvalue(Env$l.var$couleur2A))
	}
	active.legende()
	if (tclvalue(Env$l.var$facteur1)!=Env$voc[82,1]) {
	  for (i in 1:nlevels(Env$dataset[,tclvalue(Env$l.var$facteur1)])) {tkinsert(Env$l.fr1$fact.list,"end",levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[i])}
	  tclvalue(Env$l.var$niveau)<-"0"
	  tkselection.set(Env$l.fr1$fact.list,tclvalue(Env$l.var$niveau))
	}
    }
  })
  Env$l.fr1$fact.lab2<-tklabel(Env$l.frames$Fr1,text=Env$voc[130,1],font=Env$police)
  Env$l.fr1$fact.list<-tklistbox(Env$l.frames$Fr1,height=4,width=15,font=Env$police,selectmode="multiple",yscrollcommand=function(...) tkset(Env$l.fr1$fact.scroll,...))
  Env$l.fr1$fact.scroll<-tkscrollbar(Env$l.frames$Fr1,repeatinterval=4,command=function(...) tkyview(Env$l.fr1$fact.list,...))
  tkbind(Env$l.fr1$fact.list,"<Enter>",function() {msg(text=Env$voc[139,1],type="info")})
  tkbind(Env$l.fr1$fact.list,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.fr1$fact.list,"<ButtonRelease-1>",function() {
    tclvalue(Env$l.var$niveau)<-tclvalue(tkcurselection(Env$l.fr1$fact.list))
    if (length(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]])>1) {
	tclvalue(Env$l.var$plusieurs)<-1
	Env$l.var$symboleB<-rep(1,length(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]]))
	Env$l.var$couleur2B<-grey.colors(length(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]]))
	Env$l.var$taille.ptsB<-rep(1,length(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]]))
	tclvalue(Env$l.var$taille.ptsA)<-as.character(Env$l.var$taille.ptsB[1])
	Env$l.var$type.courbeB<-rep("",length(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]]))
	tclvalue(Env$l.var$type.courbeA)<-Env$l.var$type.courbeB[1]
	Env$l.var$trait2<-rep("",length(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]]))
	tclvalue(Env$l.var$trait1)<-Env$l.var$trait2[1]
	Env$l.var$epaisseur2<-rep(1,length(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]]))
	tclvalue(Env$l.var$epaisseur1)<-as.character(Env$l.var$epaisseur2[1])
	Env$l.var$noms1<-levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]])+1]
	if ("noms.list"%in%names(Env$l.fr4)) {
	  tkconfigure(Env$l.fr4$noms.list,state="normal")
	  tkdelete(Env$l.fr4$noms.list,0,"end")
	  for (i in 1:length(Env$l.var$noms1)) {tkinsert(Env$l.fr4$noms.list,"end",Env$l.var$noms1[i])}
	  for (i in 1:8) {tkconfigure(Env$l.fr4$l.symboles[[i]],borderwidth=0)}
	  tkconfigure(Env$l.fr4$l.symboles[[Env$l.var$symboleB[1]]],borderwidth=2)
	  tkconfigure(Env$l.fr4$col.wdg,bg=Env$l.var$couleur2B[1])
	}
	active.legende()
	if ("noms.list"%in%names(Env$l.fr6)) {
	  tkconfigure(Env$l.fr6$noms.list,state="normal")
	  tkdelete(Env$l.fr6$noms.list,0,"end")
	  for (i in 1:length(Env$l.var$noms1)) {tkinsert(Env$l.fr6$noms.list,"end",Env$l.var$noms1[i])}
	  tkdelete(Env$l.fr6$noms.wdg,0,"end")
	}
    } else {
	tclvalue(Env$l.var$plusieurs)<-0
	tclvalue(Env$l.var$symboleA)<-"1"
	tclvalue(Env$l.var$couleur2A)<-"black"
	tclvalue(Env$l.var$taille.ptsA)<-"1"
	tclvalue(Env$l.var$type.courbeA)<-""
	tclvalue(Env$l.var$trait1)<-""
	tclvalue(Env$l.var$epaisseur1)<-"1"
	if ("noms.list"%in%names(Env$l.fr4)) {
	  tkdelete(Env$l.fr4$noms.list,0,"end")
	  tkconfigure(Env$l.fr4$noms.list,state="disabled")
	  for (i in 1:8) {tkconfigure(Env$l.fr4$l.symboles[[i]],borderwidth=0)}
	  tkconfigure(Env$l.fr4$l.symboles[[as.numeric(tclvalue(Env$l.var$symboleA))]],borderwidth=2)
	  tkconfigure(Env$l.fr4$col.wdg,bg=tclvalue(Env$l.var$couleur2A))
	}
	if ("noms.list"%in%names(Env$l.fr6)) {
	  tkdelete(Env$l.fr6$noms.list,0,"end")
	  tkdelete(Env$l.fr6$noms.wdg,0,"end")
	}
	active.legende()
    }
  })
  Env$l.fr1$encadre.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[43,1],font=Env$police)
  Env$l.fr1$encadre.wdg<-tkcheckbutton(Env$l.frames$Fr1,variable=Env$l.var$encadre)
  Env$l.fr1$titre2<-tklabel(Env$l.frames$Fr1,text=Env$voc[85,1],font=Env$police3)
  Env$l.fr1$rb2<-tkradiobutton(Env$l.frames$Fr1,variable=Env$l.var$moyprop,value="prop",command=courbe.prop)
  Env$l.fr1$propvarX.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[127,1],font=Env$police)
  Env$l.fr1$propvarX.wdg<-ttkcombobox(Env$l.frames$Fr1,width=15,values=Env$l.var$var.num,textvariable=Env$l.var$varX.prop,font=Env$police,state="readonly")
  Env$l.fr1$propvarY.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[128,1],font=Env$police)
  Env$l.fr1$propvarY.wdg<-ttkcombobox(Env$l.frames$Fr1,width=15,values=Env$l.var$var.fact,textvariable=Env$l.var$proportions,font=Env$police,state="readonly")
  tkbind(Env$l.fr1$propvarY.wdg,"<<ComboboxSelected>>",function() {
    tclvalue(Env$l.var$prop.niveaux)<-""
    tkconfigure(Env$l.fr1$propvarY.niv.wdg,values=levels(Env$dataset[,tclvalue(Env$l.var$proportions)]))
  })
  Env$l.fr1$propvarY.niv.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[129,1],font=Env$police)
  Env$l.fr1$propvarY.niv.wdg<-ttkcombobox(Env$l.frames$Fr1,width=15,values="",textvariable=Env$l.var$prop.niveaux,font=Env$police,state="readonly")
  tkbind(Env$l.fr1$propvarY.niv.wdg,"<Enter>",function() {if(tclvalue(Env$l.var$moyprop)=="prop") {msg(text=Env$voc[140,1],type="warning")}})
  tkbind(Env$l.fr1$propvarY.niv.wdg,"<Leave>",function() {msg(text="",type="info")})
  if (tclvalue(Env$l.var$moyprop)=="moy") {courbe.moy()} else {courbe.prop()}
  Env$l.fr1$espace.ver<-tklabel(Env$l.frames$Fr1,text="",font=Env$police2)
  Env$l.fr1$espace.hor1<-tklabel(Env$l.frames$Fr1,text="                    ",font=Env$police)
  Env$l.fr1$espace.hor2<-tklabel(Env$l.frames$Fr1,text="                    ",font=Env$police)
  tkgrid(Env$l.fr1$titre1,row=0,column=0,sticky="e")
  tkgrid(Env$l.fr1$rb1,row=0,column=1,sticky="w")
  tkgrid(Env$l.fr1$espace.ver,row=1,column=0)
  tkgrid(Env$l.fr1$moyvarX.lab,row=2,column=0,sticky="e")
  tkgrid(Env$l.fr1$moyvarX.wdg,row=2,column=1,sticky="w")
  tkgrid(Env$l.fr1$moyvarY.lab,row=3,column=0,sticky="e")
  tkgrid(Env$l.fr1$moyvarY.wdg,row=3,column=1,sticky="w")
  tkgrid(Env$l.fr1$espace.hor1,row=0,column=2)
  tkgrid(Env$l.fr1$fact.lab,row=2,column=3,sticky="e")
  tkgrid(Env$l.fr1$fact.wdg,row=2,column=4,sticky="w")
  tkgrid(Env$l.fr1$fact.lab2,row=3,column=3,sticky="e")
  tkgrid(Env$l.fr1$fact.list,Env$l.fr1$fact.scroll,row=3,column=4,rowspan=4,sticky="w");tkgrid.configure(Env$l.fr1$fact.scroll,sticky="ens")
  tkgrid(Env$l.fr1$encadre.lab,row=7,column=3,sticky="e")
  tkgrid(Env$l.fr1$encadre.wdg,row=7,column=4,sticky="w")
  tkgrid(Env$l.fr1$sysinfo.lab,row=8,column=3,sticky="e")
  tkgrid(Env$l.fr1$sysinfo.wdg,row=8,column=4,sticky="w")
  tkgrid(Env$l.fr1$espace.hor2,row=0,column=5)
  tkgrid(Env$l.fr1$rb2,row=0,column=6,sticky="e")
  tkgrid(Env$l.fr1$titre2,row=0,column=7,sticky="w")
  tkgrid(Env$l.fr1$propvarX.lab,row=2,column=6,sticky="e")
  tkgrid(Env$l.fr1$propvarX.wdg,row=2,column=7,sticky="w")
  tkgrid(Env$l.fr1$propvarY.lab,row=3,column=6,sticky="e")
  tkgrid(Env$l.fr1$propvarY.wdg,row=3,column=7,sticky="w")
  tkgrid(Env$l.fr1$propvarY.niv.lab,row=4,column=6,sticky="e")
  tkgrid(Env$l.fr1$propvarY.niv.wdg,row=4,column=7,sticky="w")
}

fr1.openN<-function() {
  Env$l.frames$Fr1.status<-1
  tkconfigure(Env$l.wdg$but.lab1,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr1)) {tkdestroy(Env$l.fr1[[i]])}
  Env$l.fr1<-list()
  Env$l.fr1$varX.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[127,1],font=Env$police)
  Env$l.fr1$varX.wdg<-ttkcombobox(Env$l.frames$Fr1,values=Env$l.var$var.num,textvariable=Env$l.var$varX,font=Env$police,state="readonly")
  Env$l.fr1$varY.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[128,1],font=Env$police)
  Env$l.fr1$varY.wdg<-ttkcombobox(Env$l.frames$Fr1,values=Env$l.var$var.num,textvariable=Env$l.var$varY,font=Env$police,state="readonly")
  Env$l.fr1$encadre.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[43,1],font=Env$police)
  Env$l.fr1$encadre.wdg<-tkcheckbutton(Env$l.frames$Fr1,variable=Env$l.var$encadre)
  Env$l.fr1$sysinfo.lab<-tklabel(Env$l.frames$Fr1,text=Env$voc[248,1],font=Env$police)
  Env$l.fr1$sysinfo.wdg<-tkcheckbutton(Env$l.frames$Fr1,variable=Env$l.var$sysinfo)
  Env$l.fr1$fact.lab1<-tklabel(Env$l.frames$Fr1,text=Env$voc[37,1],font=Env$police)
  Env$l.fr1$fact.wdg<-ttkcombobox(Env$l.frames$Fr1,values=c(Env$voc[82,1],Env$l.var$var.fact),textvariable=Env$l.var$facteur1,font=Env$police,state="readonly")
  tkbind(Env$l.fr1$fact.wdg,"<<ComboboxSelected>>",function() {
    if (nchar(tclvalue(Env$l.var$facteur1))>0) {
	tclvalue(Env$l.var$plusieurs)<-0
	tclvalue(Env$l.var$symboleA)<-"1"
	tclvalue(Env$l.var$couleur2A)<-"black"
	tclvalue(Env$l.var$taille.ptsA)<-"1"
	tclvalue(Env$l.var$droiteA)<-""
	tclvalue(Env$l.var$intervalA)<-""
	tclvalue(Env$l.var$trait1)<-""
	tclvalue(Env$l.var$epaisseur1)<-"1"
	tkdelete(Env$l.fr1$fact.list,0,"end")
	if ("noms.list"%in%names(Env$l.fr4)) {
	  tkdelete(Env$l.fr4$noms.list,0,"end")
	  tkconfigure(Env$l.fr4$noms.list,state="disabled")
	  for (i in 1:8) {tkconfigure(Env$l.fr4$l.symboles[[i]],borderwidth=0)}
	  tkconfigure(Env$l.fr4$l.symboles[[as.numeric(tclvalue(Env$l.var$symboleA))]],borderwidth=2)
	  tkconfigure(Env$l.fr4$col.wdg,bg=tclvalue(Env$l.var$couleur2A))
	}
	if ("noms.list"%in%names(Env$l.fr5)) {
	  tkdelete(Env$l.fr5$noms.list,0,"end")
	  tkconfigure(Env$l.fr5$noms.list,state="disabled")
	}
	active.legende()
	if (tclvalue(Env$l.var$facteur1)!=Env$voc[82,1]) {
	  for (i in 1:nlevels(Env$dataset[,tclvalue(Env$l.var$facteur1)])) {tkinsert(Env$l.fr1$fact.list,"end",levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[i])}
	  tclvalue(Env$l.var$niveau)<-"0"
	  tkselection.set(Env$l.fr1$fact.list,tclvalue(Env$l.var$niveau))
	}
    }
  })
  Env$l.fr1$fact.lab2<-tklabel(Env$l.frames$Fr1,text=Env$voc[130,1],font=Env$police)
  Env$l.fr1$fact.list<-tklistbox(Env$l.frames$Fr1,height=5,font=Env$police,selectmode="multiple",yscrollcommand=function(...) tkset(Env$l.fr1$fact.scroll,...))
  Env$l.fr1$fact.scroll<-tkscrollbar(Env$l.frames$Fr1,repeatinterval=4,command=function(...) tkyview(Env$l.fr1$fact.list,...))
  tkbind(Env$l.fr1$fact.list,"<Enter>",function() {msg(text=Env$voc[143,1],type="info")})
  tkbind(Env$l.fr1$fact.list,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.fr1$fact.list,"<ButtonRelease-1>",function() {
    tclvalue(Env$l.var$niveau)<-tclvalue(tkcurselection(Env$l.fr1$fact.list))
    if (length(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]])>1) {
	tclvalue(Env$l.var$plusieurs)<-1
	Env$l.var$symboleB<-rep(1,length(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]]))
	Env$l.var$couleur2B<-grey.colors(length(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]]))
	Env$l.var$taille.ptsB<-rep(1,length(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]]))
	tclvalue(Env$l.var$taille.ptsA)<-as.character(Env$l.var$taille.ptsB[1])
	Env$l.var$droiteB<-rep("",length(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]]))
	tclvalue(Env$l.var$droiteA)<-Env$l.var$droiteB[1]
	Env$l.var$intervalB<-rep("",length(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]]))
	tclvalue(Env$l.var$intervalA)<-Env$l.var$intervalB[1]
	Env$l.var$trait2<-rep("",length(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]]))
	tclvalue(Env$l.var$trait1)<-Env$l.var$trait2[1]
	Env$l.var$epaisseur2<-rep(1,length(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]]))
	tclvalue(Env$l.var$epaisseur1)<-as.character(Env$l.var$epaisseur2[1])
	Env$l.var$noms1<-levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]])+1]
	if ("noms.list"%in%names(Env$l.fr4)) {
	  tkconfigure(Env$l.fr4$noms.list,state="normal")
	  tkdelete(Env$l.fr4$noms.list,0,"end")
	  for (i in 1:length(Env$l.var$noms1)) {tkinsert(Env$l.fr4$noms.list,"end",Env$l.var$noms1[i])}
	  for (i in 1:8) {tkconfigure(Env$l.fr4$l.symboles[[i]],borderwidth=0)}
	  tkconfigure(Env$l.fr4$l.symboles[[Env$l.var$symboleB[1]]],borderwidth=2)
	  tkconfigure(Env$l.fr4$col.wdg,bg=Env$l.var$couleur2B[1])
	}
	if ("noms.list"%in%names(Env$l.fr5)) {
	  tkconfigure(Env$l.fr5$noms.list,state="normal")
	  tkdelete(Env$l.fr5$noms.list,0,"end")
	  for (i in 1:length(Env$l.var$noms1)) {tkinsert(Env$l.fr5$noms.list,"end",Env$l.var$noms1[i])}
	}
	active.legende()
	if ("noms.list"%in%names(Env$l.fr6)) {
	  tkconfigure(Env$l.fr6$noms.list,state="normal")
	  tkdelete(Env$l.fr6$noms.list,0,"end")
	  for (i in 1:length(Env$l.var$noms1)) {tkinsert(Env$l.fr6$noms.list,"end",Env$l.var$noms1[i])}
	  tkdelete(Env$l.fr6$noms.wdg,0,"end")
	}
    } else {
	tclvalue(Env$l.var$plusieurs)<-0
	tclvalue(Env$l.var$symboleA)<-"1"
	tclvalue(Env$l.var$couleur2A)<-"black"
	tclvalue(Env$l.var$taille.ptsA)<-"1"
	tclvalue(Env$l.var$droiteA)<-""
	tclvalue(Env$l.var$trait1)<-""
	tclvalue(Env$l.var$epaisseur1)<-"1"
	if ("noms.list"%in%names(Env$l.fr4)) {
	  tkdelete(Env$l.fr4$noms.list,0,"end")
	  tkconfigure(Env$l.fr4$noms.list,state="disabled")
	  for (i in 1:8) {tkconfigure(Env$l.fr4$l.symboles[[i]],borderwidth=0)}
	  tkconfigure(Env$l.fr4$l.symboles[[as.numeric(tclvalue(Env$l.var$symboleA))]],borderwidth=2)
	  tkconfigure(Env$l.fr4$col.wdg,bg=tclvalue(Env$l.var$couleur2A))
	}
	if ("noms.list"%in%names(Env$l.fr5)) {
	  tkdelete(Env$l.fr5$noms.list,0,"end")
	  tkconfigure(Env$l.fr5$noms.list,state="disabled")
	}
	if ("noms.list"%in%names(Env$l.fr6)) {
	  tkdelete(Env$l.fr6$noms.list,0,"end")
	  tkdelete(Env$l.fr6$noms.wdg,0,"end")
	}
	active.legende()
    }
  })
  Env$l.fr1$espace.hor<-tklabel(Env$l.frames$Fr1,text="                                        ",font=Env$police)
  tkgrid(Env$l.fr1$varX.lab,row=0,column=0,sticky="e")
  tkgrid(Env$l.fr1$varX.wdg,row=0,column=1,sticky="w")
  tkgrid(Env$l.fr1$varY.lab,row=1,column=0,sticky="e")
  tkgrid(Env$l.fr1$varY.wdg,row=1,column=1,sticky="w")
  tkgrid(Env$l.fr1$encadre.lab,row=2,column=0,sticky="e")
  tkgrid(Env$l.fr1$encadre.wdg,row=2,column=1,sticky="w")
  tkgrid(Env$l.fr1$espace.hor,row=0,column=2)
  tkgrid(Env$l.fr1$fact.lab1,row=0,column=3,sticky="e")
  tkgrid(Env$l.fr1$fact.wdg,row=0,column=4,sticky="w")
  tkgrid(Env$l.fr1$fact.lab2,row=1,column=3,sticky="e")
  tkgrid(Env$l.fr1$fact.list,Env$l.fr1$fact.scroll,row=1,column=4,rowspan=4,sticky="w");tkgrid.configure(Env$l.fr1$fact.scroll,sticky="ens")
  tkgrid(Env$l.fr1$sysinfo.lab,row=4,column=0,sticky="e")
  tkgrid(Env$l.fr1$sysinfo.wdg,row=4,column=1,sticky="w")
}


#-------------------------------------------------
# Frame 2
#-------------------------------------------------

fr2.close<-function() {
  Env$l.frames$Fr2.status<-0
  tkconfigure(Env$l.wdg$but.lab2,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_bas.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr2)) {tkdestroy(Env$l.fr2[[i]])}
  Env$l.fr2<-list()
  Env$l.fr2$vide<-tklabel(Env$l.frames$Fr2,text="",font=Env$police2)
  tkgrid(Env$l.fr2$vide)
}

fr2.openD<-function() {
  Env$l.frames$Fr2.status<-1
  tkconfigure(Env$l.wdg$but.lab2,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr2)) {tkdestroy(Env$l.fr2[[i]])}
  Env$l.fr2<-list()
  Env$l.fr2$titre1<-tklabel(Env$l.frames$Fr2,text=Env$voc[15,1],font=Env$police3)
  Env$l.fr2$var.list<-tklistbox(Env$l.frames$Fr2,height=7,font=Env$police,selectmode="single",yscrollcommand=function(...) tkset(Env$l.fr2$var.scroll,...))
  Env$l.fr2$var.scroll<-tkscrollbar(Env$l.frames$Fr2,repeatinterval=5,command=function(...) tkyview(Env$l.fr2$var.list,...))
  tkbind(Env$l.fr2$var.list,"<ButtonRelease-1>",function() {
    if(!is.null(Env$dataset)) {
	tkconfigure(Env$l.fr2$type.wdg,text=class(Env$dataset[,as.numeric(tclvalue(tkcurselection(Env$l.fr2$var.list)))+1]))
	tkconfigure(Env$l.fr2$resume.wdg,state="normal")
	tkdelete(Env$l.fr2$resume.wdg,"0.0","end")
	resume<-cbind(summary(Env$dataset)[,as.numeric(tclvalue(tkcurselection(Env$l.fr2$var.list)))+1])
	for (i in 1:nrow(resume)) {
	  if(!is.na(resume[i])) {tkinsert(Env$l.fr2$resume.wdg,"end",paste(resume[i],"\n",sep=""))}
	}
	tkconfigure(Env$l.fr2$resume.wdg,state="disabled")
    }
  })
  Env$l.fr2$type.lab<-tklabel(Env$l.frames$Fr2,text=Env$voc[16,1],font=Env$police)
  Env$l.fr2$type.wdg<-tklabel(Env$l.frames$Fr2,text="",font=Env$police)
  Env$l.fr2$titre2<-tklabel(Env$l.frames$Fr2,text=Env$voc[17,1],font=Env$police3)
  Env$l.fr2$resume.wdg<-tktext(Env$l.frames$Fr2,width=16,height=8,font=Env$police4,state="disabled")
  tkdelete(Env$l.fr2$var.list,0,"end")
  if(!is.null(Env$dataset)) {
    for (i in 1:ncol(Env$dataset)) {tkinsert(Env$l.fr2$var.list,"end",colnames(Env$dataset)[i])}
  }
  tkconfigure(Env$l.fr2$type.wdg,text="")
  tkconfigure(Env$l.fr2$resume.wdg,state="normal")
  tkdelete(Env$l.fr2$resume.wdg,"0.0","end")
  tkconfigure(Env$l.fr2$resume.wdg,state="disabled")
  Env$l.fr2$espace.ver<-tklabel(Env$l.frames$Fr2,text="",font=Env$police2)
  Env$l.fr2$espace.hor1<-tklabel(Env$l.frames$Fr2,text="                                                  ",font=Env$police2)
  Env$l.fr2$espace.hor2<-tklabel(Env$l.frames$Fr2,text="                                                  ",font=Env$police2)
  tkgrid(Env$l.fr2$titre1,row=0,column=0)
  tkgrid(Env$l.fr2$espace.ver,row=1,column=0)
  tkgrid(Env$l.fr2$var.list,Env$l.fr2$var.scroll,row=2,column=0,rowspan=5);tkgrid.configure(Env$l.fr2$var.scroll,sticky="ens")
  tkgrid(Env$l.fr2$espace.hor1,row=2,column=2)
  tkgrid(Env$l.fr2$type.lab,row=3,column=3,sticky="e")
  tkgrid(Env$l.fr2$type.wdg,row=3,column=4,sticky="w")
  tkgrid(Env$l.fr2$espace.hor2,row=2,column=5)
  tkgrid(Env$l.fr2$titre2,row=0,column=6)
  tkgrid(Env$l.fr2$resume.wdg,row=2,column=6,rowspan=5)
}

fr2.opengraphe<-function() {
  Env$l.frames$Fr2.status<-1
  tkconfigure(Env$l.wdg$but.lab2,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr2)) {tkdestroy(Env$l.fr2[[i]])}
  Env$l.fr2<-list()
  Env$l.fr2$titre.lab<-tklabel(Env$l.frames$Fr2,text=Env$voc[44,1],font=Env$police)
  Env$l.fr2$titre.wdg<-tkentry(Env$l.frames$Fr2,width=30,font=Env$police,textvariable=Env$l.var$titre)
  Env$l.fr2$col.lab<-tklabel(Env$l.frames$Fr2,text=Env$voc[45,1],font=Env$police)
  Env$l.fr2$col.wdg<-tkbutton(Env$l.frames$Fr2,text="Aa",font=Env$police3,foreground=tclvalue(Env$l.var$titre.col),activeforeground=tclvalue(Env$l.var$titre.col),command=function() {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$titre.col),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$titre.col)<-temp
	tkconfigure(Env$l.fr2$col.wdg,foreground=temp,activeforeground=temp)
    }
  })
  Env$l.fr2$taille.lab<-tklabel(Env$l.frames$Fr2,text=Env$voc[46,1],font=Env$police)
  Env$l.fr2$taille.wdg<-tkscale(Env$l.frames$Fr2,from=0.5,to=4,showvalue=TRUE,font=Env$police,variable=Env$l.var$titre.taille,resolution=0.1,orient="horizontal")
  Env$l.fr2$soustitre.lab<-tklabel(Env$l.frames$Fr2,text=Env$voc[247,1],font=Env$police)
  Env$l.fr2$soustitre.wdg<-tkentry(Env$l.frames$Fr2,width=30,font=Env$police,textvariable=Env$l.var$soustitre)
  Env$l.fr2$espace.hor1<-tklabel(Env$l.frames$Fr2,text="                              ",font=Env$police)
  Env$l.fr2$espace.hor2<-tklabel(Env$l.frames$Fr2,text="                                        ",font=Env$police)
  tkgrid(Env$l.fr2$titre.lab,row=0,column=0,sticky="e")
  tkgrid(Env$l.fr2$titre.wdg,row=0,column=1,sticky="w")
  tkgrid(Env$l.fr2$espace.hor1,row=0,column=2)
  tkgrid(Env$l.fr2$col.lab,row=0,column=3,sticky="e")
  tkgrid(Env$l.fr2$col.wdg,row=0,column=4,sticky="w")
  tkgrid(Env$l.fr2$espace.hor2,row=0,column=5)
  tkgrid(Env$l.fr2$taille.lab,row=0,column=6,sticky="e")
  tkgrid(Env$l.fr2$taille.wdg,row=0,column=7,sticky="w")
  tkgrid(Env$l.fr2$soustitre.lab,row=1,column=0,sticky="e")
  tkgrid(Env$l.fr2$soustitre.wdg,row=1,column=1,sticky="w")
}


#-------------------------------------------------
# Frame 3
#-------------------------------------------------

fr3.close<-function() {
  Env$l.frames$Fr3.status<-0
  tkconfigure(Env$l.wdg$but.lab3,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_bas.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr3)) {tkdestroy(Env$l.fr3[[i]])}
  Env$l.fr3<-list()
  Env$l.fr3$vide<-tklabel(Env$l.frames$Fr3,text="",font=Env$police2)
  tkgrid(Env$l.fr3$vide)
}

fr3.openD<-function() {
  Env$l.frames$Fr3.status<-1
  tkconfigure(Env$l.wdg$but.lab3,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr3)) {tkdestroy(Env$l.fr3[[i]])}
  Env$l.fr3<-list()
  Env$l.fr3$titre1<-tklabel(Env$l.frames$Fr3,text=Env$voc[15,1],font=Env$police3)
  Env$l.fr3$var.list<-tklistbox(Env$l.frames$Fr3,height=7,font=Env$police,selectmode="single",yscrollcommand=function(...) tkset(Env$l.fr3$var.scroll,...))
  Env$l.fr3$var.scroll<-tkscrollbar(Env$l.frames$Fr3,repeatinterval=5,command=function(...) tkyview(Env$l.fr3$var.list,...))
  tkbind(Env$l.fr3$var.list,"<ButtonRelease-1>",function() {
    if(!is.null(Env$dataset)) {
	tkdelete(Env$l.fr3$nom.wdg,0,"end")
	tkinsert(Env$l.fr3$nom.wdg,"end",names(Env$dataset[as.numeric(tclvalue(tkcurselection(Env$l.fr3$var.list)))+1]))
    }
  })
  Env$l.fr3$titre2<-tklabel(Env$l.frames$Fr3,text=Env$voc[22,1],font=Env$police3)
  Env$l.fr3$nom.wdg<-tkentry(Env$l.frames$Fr3,width=20,font=Env$police)
  tkbind(Env$l.fr3$nom.wdg,"<Enter>",function() {msg(text=Env$voc[26,1],type="info")})
  tkbind(Env$l.fr3$nom.wdg,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.fr3$nom.wdg,"<ButtonRelease-1>",function() {
    if(!is.null(Env$dataset)) {
	tkdelete(Env$l.fr3$nom.wdg,0,"end")
    }
  })
  tkbind(Env$l.fr3$nom.wdg,"<Return>",rename.variable)
  if(!is.null(Env$dataset)) {
    for (i in 1:ncol(Env$dataset)) {tkinsert(Env$l.fr3$var.list,"end",colnames(Env$dataset)[i])}
  }
  Env$l.fr3$espace.ver<-tklabel(Env$l.frames$Fr3,text="",font=Env$police2)
  Env$l.fr3$espace.hor<-tklabel(Env$l.frames$Fr3,text="                                                       ",font=Env$police)
  tkgrid(Env$l.fr3$titre1,row=0,column=0)
  tkgrid(Env$l.fr3$espace.ver)
  tkgrid(Env$l.fr3$var.list,Env$l.fr3$var.scroll,row=2,column=0,rowspan=5);tkgrid.configure(Env$l.fr3$var.scroll,sticky="ens")
  tkgrid(Env$l.fr3$espace.hor,row=2,column=1)
  tkgrid(Env$l.fr3$titre2,row=0,column=2)
  tkgrid(Env$l.fr3$nom.wdg,row=2,column=2)
}

fr3.openH<-function() {
  Env$l.frames$Fr3.status<-1
  tkconfigure(Env$l.wdg$but.lab3,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr3)) {tkdestroy(Env$l.fr3[[i]])}
  Env$l.fr3<-list()
  Env$l.fr3$grad.lab1<-tklabel(Env$l.frames$Fr3,text=Env$voc[47,1],font=Env$police)
  Env$l.fr3$grad.col<-tkbutton(Env$l.frames$Fr3,text="Aa",font=Env$police3,foreground=tclvalue(Env$l.var$graduations.col),activeforeground=tclvalue(Env$l.var$graduations.col),command=function() {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$graduations.col),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$graduations.col)<-temp
	tkconfigure(Env$l.fr3$grad.col,foreground=temp,activeforeground=temp)
    }
  })
  Env$l.fr3$grad.lab2<-tklabel(Env$l.frames$Fr3,text=Env$voc[48,1],font=Env$police)
  Env$l.fr3$grad.taille<-tkscale(Env$l.frames$Fr3,from=0.5,to=4,showvalue=TRUE,font=Env$police,variable=Env$l.var$graduations.taille,resolution=0.1,orient="horizontal")
  Env$l.fr3$grad.lab3<-tklabel(Env$l.frames$Fr3,text=Env$voc[245,1],font=Env$police)
  Env$l.fr3$grad.orient<-ttkcombobox(Env$l.frames$Fr3,width=16,values=c(Env$voc[246,1],Env$voc[67,1]),textvariable=Env$l.var$graduations.orient,font=Env$police,state="readonly")
  Env$l.fr3$leg.lab1<-tklabel(Env$l.frames$Fr3,text=Env$voc[49,1],font=Env$police)
  Env$l.fr3$leg.col<-tkbutton(Env$l.frames$Fr3,text="Aa",font=Env$police3,foreground=tclvalue(Env$l.var$legendes.col),activeforeground=tclvalue(Env$l.var$legendes.col),command=function() {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$legendes.col),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$legendes.col)<-temp
	tkconfigure(Env$l.fr3$leg.col,foreground=temp,activeforeground=temp)
    }
  })
  Env$l.fr3$leg.lab2<-tklabel(Env$l.frames$Fr3,text=Env$voc[50,1],font=Env$police)
  Env$l.fr3$leg.taille<-tkscale(Env$l.frames$Fr3,from=0.5,to=4,showvalue=TRUE,font=Env$police,variable=Env$l.var$legendes.taille,resolution=0.1,orient="horizontal")
  Env$l.fr3$titre1<-tklabel(Env$l.frames$Fr3,text=Env$voc[51,1],font=Env$police3)
  Env$l.fr3$hor.titre.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[44,1],font=Env$police)
  Env$l.fr3$hor.titre.wdg<-tkentry(Env$l.frames$Fr3,width=20,font=Env$police,textvariable=Env$l.var$titre.axehor)
  Env$l.fr3$hor.liminf.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[53,1],font=Env$police)
  Env$l.fr3$hor.liminf.wdg<-tkentry(Env$l.frames$Fr3,width=5,font=Env$police,textvariable=Env$l.var$liminf.axehor,state=ifelse(tclvalue(Env$l.var$hist.type)==Env$voc[40,1],"disabled","normal"))
  tkbind(Env$l.fr3$hor.liminf.wdg,"<Enter>",function() {msg(text=Env$voc[155,1],type="warning")})
  tkbind(Env$l.fr3$hor.liminf.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr3$hor.limsup.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[54,1],font=Env$police)
  Env$l.fr3$hor.limsup.wdg<-tkentry(Env$l.frames$Fr3,width=5,font=Env$police,textvariable=Env$l.var$limsup.axehor,state=ifelse(tclvalue(Env$l.var$hist.type)==Env$voc[40,1],"disabled","normal"))
  tkbind(Env$l.fr3$hor.limsup.wdg,"<Enter>",function() {msg(text=Env$voc[155,1],type="warning")})
  tkbind(Env$l.fr3$hor.limsup.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr3$titre2<-tklabel(Env$l.frames$Fr3,text=Env$voc[52,1],font=Env$police3)
  Env$l.fr3$ver.titre.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[44,1],font=Env$police)
  Env$l.fr3$ver.titre.wdg<-tkentry(Env$l.frames$Fr3,width=20,font=Env$police,textvariable=Env$l.var$titre.axever)
  Env$l.fr3$ver.limsup.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[54,1],font=Env$police)
  Env$l.fr3$ver.limsup.wdg<-tkentry(Env$l.frames$Fr3,width=5,font=Env$police,textvariable=Env$l.var$limsup.axever)
  tkbind(Env$l.fr3$ver.limsup.wdg,"<Enter>",function() {msg(text=Env$voc[155,1],type="warning")})
  tkbind(Env$l.fr3$ver.limsup.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr3$espace.ver1<-tklabel(Env$l.frames$Fr3,text="",font=Env$police2)
  Env$l.fr3$espace.ver2<-tklabel(Env$l.frames$Fr3,text="",font=Env$police2)
  Env$l.fr3$espace.hor<-tklabel(Env$l.frames$Fr3,text="                                   ",font=Env$police)
  tkgrid(Env$l.fr3$grad.lab1,row=0,column=0,sticky="e")
  tkgrid(Env$l.fr3$grad.col,row=0,column=1,sticky="w")
  tkgrid(Env$l.fr3$espace.hor,row=0,column=2)
  tkgrid(Env$l.fr3$grad.lab2,row=0,column=3,sticky="e")
  tkgrid(Env$l.fr3$grad.taille,row=0,column=4,sticky="w")
  tkgrid(Env$l.fr3$grad.lab3,row=1,column=0,sticky="e")
  tkgrid(Env$l.fr3$grad.orient,row=1,column=1,sticky="w")
  tkgrid(Env$l.fr3$leg.lab1,row=2,column=0,sticky="e")
  tkgrid(Env$l.fr3$leg.col,row=2,column=1,sticky="w")
  tkgrid(Env$l.fr3$leg.lab2,row=2,column=3,sticky="e")
  tkgrid(Env$l.fr3$leg.taille,row=2,column=4,sticky="w")
  tkgrid(Env$l.fr3$espace.ver1)
  tkgrid(Env$l.fr3$titre1,row=4,column=0,columnspan=2)
  tkgrid(Env$l.fr3$espace.ver2)
  tkgrid(Env$l.fr3$hor.titre.lab,row=6,column=0,sticky="e")
  tkgrid(Env$l.fr3$hor.titre.wdg,row=6,column=1,sticky="w")
  tkgrid(Env$l.fr3$hor.liminf.lab,row=7,column=0,sticky="e")
  tkgrid(Env$l.fr3$hor.liminf.wdg,row=7,column=1,sticky="w")
  tkgrid(Env$l.fr3$hor.limsup.lab,row=8,column=0,sticky="e")
  tkgrid(Env$l.fr3$hor.limsup.wdg,row=8,column=1,sticky="w")
  tkgrid(Env$l.fr3$titre2,row=4,column=4,columnspan=2)
  tkgrid(Env$l.fr3$ver.titre.lab,row=6,column=3,sticky="e")
  tkgrid(Env$l.fr3$ver.titre.wdg,row=6,column=4,sticky="w")
  tkgrid(Env$l.fr3$ver.limsup.lab,row=7,column=3,sticky="e")
  tkgrid(Env$l.fr3$ver.limsup.wdg,row=7,column=4,sticky="w")
}

fr3.openM<-function() {
  Env$l.frames$Fr3.status<-1
  tkconfigure(Env$l.wdg$but.lab3,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr3)) {tkdestroy(Env$l.fr3[[i]])}
  Env$l.fr3<-list()
  Env$l.fr3$grad.lab1<-tklabel(Env$l.frames$Fr3,text=Env$voc[47,1],font=Env$police)
  Env$l.fr3$grad.col<-tkbutton(Env$l.frames$Fr3,text="Aa",font=Env$police3,foreground=tclvalue(Env$l.var$graduations.col),activeforeground=tclvalue(Env$l.var$graduations.col),command=function() {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$graduations.col),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$graduations.col)<-temp
	tkconfigure(Env$l.fr3$grad.col,foreground=temp,activeforeground=temp)
    }
  })
  Env$l.fr3$grad.lab2<-tklabel(Env$l.frames$Fr3,text=Env$voc[48,1],font=Env$police)
  Env$l.fr3$grad.taille<-tkscale(Env$l.frames$Fr3,from=0.5,to=4,showvalue=TRUE,font=Env$police,variable=Env$l.var$graduations.taille,resolution=0.1,orient="horizontal")
  Env$l.fr3$grad.lab3<-tklabel(Env$l.frames$Fr3,text=Env$voc[245,1],font=Env$police)
  Env$l.fr3$grad.orient<-ttkcombobox(Env$l.frames$Fr3,width=16,values=c(Env$voc[246,1],Env$voc[67,1]),textvariable=Env$l.var$graduations.orient,font=Env$police,state="readonly")
  Env$l.fr3$leg.lab1<-tklabel(Env$l.frames$Fr3,text=Env$voc[49,1],font=Env$police)
  Env$l.fr3$leg.col<-tkbutton(Env$l.frames$Fr3,text="Aa",font=Env$police3,foreground=tclvalue(Env$l.var$legendes.col),activeforeground=tclvalue(Env$l.var$legendes.col),command=function() {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$legendes.col),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$legendes.col)<-temp
	tkconfigure(Env$l.fr3$leg.col,foreground=temp,activeforeground=temp)
    }
  })
  Env$l.fr3$leg.lab2<-tklabel(Env$l.frames$Fr3,text=Env$voc[50,1],font=Env$police)
  Env$l.fr3$leg.taille<-tkscale(Env$l.frames$Fr3,from=0.5,to=4,showvalue=TRUE,font=Env$police,variable=Env$l.var$legendes.taille,resolution=0.1,orient="horizontal")
  Env$l.fr3$titre1<-tklabel(Env$l.frames$Fr3,text=Env$voc[69,1],font=Env$police3)
  Env$l.fr3$noms.titre.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[44,1],font=Env$police)
  Env$l.fr3$noms.titre.wdg<-tkentry(Env$l.frames$Fr3,width=20,font=Env$police,textvariable=Env$l.var$titre.axenoms)
  Env$l.fr3$noms.lab1<-tklabel(Env$l.frames$Fr3,text=Env$voc[71,1],font=Env$police)
  Env$l.fr3$noms.list<-tklistbox(Env$l.frames$Fr3,height=4,font=Env$police,selectmode="single",yscrollcommand=function(...) tkset(Env$l.fr3$noms.scroll,...))
  Env$l.fr3$noms.scroll<-tkscrollbar(Env$l.frames$Fr3,repeatinterval=4,command=function(...) tkyview(Env$l.fr3$noms.list,...))
  for (i in 1:length(Env$l.var$noms1)) tkinsert(Env$l.fr3$noms.list,"end",Env$l.var$noms1[i])
  tkbind(Env$l.fr3$noms.list,"<ButtonRelease-1>",function() {
    tkdelete(Env$l.fr3$noms.wdg,0,"end")
    tkinsert(Env$l.fr3$noms.wdg,"end",Env$l.var$noms1[as.numeric(tclvalue(tkcurselection(Env$l.fr3$noms.list)))+1])
  })
  Env$l.fr3$noms.lab2<-tklabel(Env$l.frames$Fr3,text=Env$voc[83,1],font=Env$police)
  Env$l.fr3$noms.wdg<-tkentry(Env$l.frames$Fr3,width=20,font=Env$police)
  tkbind(Env$l.fr3$noms.wdg,"<Enter>",function() {msg(text=Env$voc[26,1],type="info")})
  tkbind(Env$l.fr3$noms.wdg,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.fr3$noms.wdg,"<ButtonRelease-1>",function() {tkdelete(Env$l.fr3$noms.wdg,0,"end")})
  tkbind(Env$l.fr3$noms.wdg,"<Return>",function() {rename.noms1(value.list=tclvalue(tkcurselection(Env$l.fr3$noms.list)),value.nom=tclvalue(tkget(Env$l.fr3$noms.wdg)))})
  Env$l.fr3$titre2<-tklabel(Env$l.frames$Fr3,text=Env$voc[70,1],font=Env$police3)
  Env$l.fr3$valeurs.titre.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[44,1],font=Env$police)
  Env$l.fr3$valeurs.titre.wdg<-tkentry(Env$l.frames$Fr3,width=20,font=Env$police,textvariable=Env$l.var$titre.axevaleurs)
  Env$l.fr3$valeurs.liminf.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[53,1],font=Env$police)
  Env$l.fr3$valeurs.liminf.wdg<-tkentry(Env$l.frames$Fr3,width=5,font=Env$police,textvariable=Env$l.var$liminf.axevaleurs)
  tkbind(Env$l.fr3$valeurs.liminf.wdg,"<Enter>",function() {msg(text=Env$voc[155,1],type="warning")})
  tkbind(Env$l.fr3$valeurs.liminf.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr3$valeurs.limsup.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[54,1],font=Env$police)
  Env$l.fr3$valeurs.limsup.wdg<-tkentry(Env$l.frames$Fr3,width=5,font=Env$police,textvariable=Env$l.var$limsup.axevaleurs)
  tkbind(Env$l.fr3$valeurs.limsup.wdg,"<Enter>",function() {msg(text=Env$voc[155,1],type="warning")})
  tkbind(Env$l.fr3$valeurs.limsup.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr3$valeurs.log.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[73,1],font=Env$police)
  Env$l.fr3$valeurs.log.wdg<-tkcheckbutton(Env$l.frames$Fr3,variable=Env$l.var$log.axevaleurs)
  Env$l.fr3$espace.ver1<-tklabel(Env$l.frames$Fr3,text="",font=Env$police2)
  Env$l.fr3$espace.ver2<-tklabel(Env$l.frames$Fr3,text="",font=Env$police2)
  Env$l.fr3$espace.hor<-tklabel(Env$l.frames$Fr3,text="                                   ",font=Env$police)
  tkgrid(Env$l.fr3$grad.lab1,row=0,column=0,sticky="e")
  tkgrid(Env$l.fr3$grad.col,row=0,column=1,sticky="w")
  tkgrid(Env$l.fr3$espace.hor,row=0,column=2)
  tkgrid(Env$l.fr3$grad.lab2,row=0,column=3,sticky="e")
  tkgrid(Env$l.fr3$grad.taille,row=0,column=4,sticky="w")
  tkgrid(Env$l.fr3$grad.lab3,row=1,column=0,sticky="e")
  tkgrid(Env$l.fr3$grad.orient,row=1,column=1,sticky="w")
  tkgrid(Env$l.fr3$leg.lab1,row=2,column=0,sticky="e")
  tkgrid(Env$l.fr3$leg.col,row=2,column=1,sticky="w")
  tkgrid(Env$l.fr3$leg.lab2,row=2,column=3,sticky="e")
  tkgrid(Env$l.fr3$leg.taille,row=2,column=4,sticky="w")
  tkgrid(Env$l.fr3$espace.ver1)
  tkgrid(Env$l.fr3$titre1,row=4,column=0,columnspan=2)
  tkgrid(Env$l.fr3$espace.ver2)
  tkgrid(Env$l.fr3$noms.titre.lab,row=6,column=0,sticky="e")
  tkgrid(Env$l.fr3$noms.titre.wdg,row=6,column=1,sticky="w")
  tkgrid(Env$l.fr3$noms.lab1,row=7,column=0,sticky="e")
  tkgrid(Env$l.fr3$noms.list,Env$l.fr3$noms.scroll,row=7,column=1,rowspan=4);tkgrid.configure(Env$l.fr3$noms.scroll,sticky="ens")
  tkgrid(Env$l.fr3$noms.lab2,row=11,column=0,sticky="e")
  tkgrid(Env$l.fr3$noms.wdg,row=11,column=1,sticky="w")
  tkgrid(Env$l.fr3$titre2,row=4,column=3,columnspan=2)
  tkgrid(Env$l.fr3$valeurs.titre.lab,row=6,column=3,sticky="e")
  tkgrid(Env$l.fr3$valeurs.titre.wdg,row=6,column=4,sticky="w")
  tkgrid(Env$l.fr3$valeurs.liminf.lab,row=7,column=3,sticky="e")
  tkgrid(Env$l.fr3$valeurs.liminf.wdg,row=7,column=4,sticky="w")
  tkgrid(Env$l.fr3$valeurs.limsup.lab,row=8,column=3,sticky="e")
  tkgrid(Env$l.fr3$valeurs.limsup.wdg,row=8,column=4,sticky="w")
  tkgrid(Env$l.fr3$valeurs.log.lab,row=9,column=3,sticky="e")
  tkgrid(Env$l.fr3$valeurs.log.wdg,row=9,column=4,sticky="w")
}

fr3.openB<-function() {
  Env$l.frames$Fr3.status<-1
  tkconfigure(Env$l.wdg$but.lab3,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr3)) {tkdestroy(Env$l.fr3[[i]])}
  Env$l.fr3<-list()
  Env$l.fr3$grad.lab1<-tklabel(Env$l.frames$Fr3,text=Env$voc[47,1],font=Env$police)
  Env$l.fr3$grad.col<-tkbutton(Env$l.frames$Fr3,text="Aa",font=Env$police3,foreground=tclvalue(Env$l.var$graduations.col),activeforeground=tclvalue(Env$l.var$graduations.col),command=function() {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$graduations.col),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$graduations.col)<-temp
	tkconfigure(Env$l.fr3$grad.col,foreground=temp,activeforeground=temp)
    }
  })
  Env$l.fr3$grad.lab2<-tklabel(Env$l.frames$Fr3,text=Env$voc[48,1],font=Env$police)
  Env$l.fr3$grad.taille<-tkscale(Env$l.frames$Fr3,from=0.5,to=4,showvalue=TRUE,font=Env$police,variable=Env$l.var$graduations.taille,resolution=0.1,orient="horizontal")
  Env$l.fr3$grad.lab3<-tklabel(Env$l.frames$Fr3,text=Env$voc[245,1],font=Env$police)
  Env$l.fr3$grad.orient<-ttkcombobox(Env$l.frames$Fr3,width=16,values=c(Env$voc[246,1],Env$voc[67,1]),textvariable=Env$l.var$graduations.orient,font=Env$police,state="readonly")
  Env$l.fr3$leg.lab1<-tklabel(Env$l.frames$Fr3,text=Env$voc[49,1],font=Env$police)
  Env$l.fr3$leg.col<-tkbutton(Env$l.frames$Fr3,text="Aa",font=Env$police3,foreground=tclvalue(Env$l.var$legendes.col),activeforeground=tclvalue(Env$l.var$legendes.col),command=function() {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$legendes.col),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$legendes.col)<-temp
	tkconfigure(Env$l.fr3$leg.col,foreground=temp,activeforeground=temp)
    }
  })
  Env$l.fr3$leg.lab2<-tklabel(Env$l.frames$Fr3,text=Env$voc[50,1],font=Env$police)
  Env$l.fr3$leg.taille<-tkscale(Env$l.frames$Fr3,from=0.5,to=4,showvalue=TRUE,font=Env$police,variable=Env$l.var$legendes.taille,resolution=0.1,orient="horizontal")
  Env$l.fr3$titre1<-tklabel(Env$l.frames$Fr3,text=Env$voc[51,1],font=Env$police3)
  Env$l.fr3$noms.titre.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[44,1],font=Env$police)
  Env$l.fr3$noms.titre.wdg<-tkentry(Env$l.frames$Fr3,width=20,font=Env$police,textvariable=Env$l.var$titre.axehor)
  Env$l.fr3$noms.lab1<-tklabel(Env$l.frames$Fr3,text=Env$voc[89,1],font=Env$police)
  Env$l.fr3$noms.list<-tklistbox(Env$l.frames$Fr3,height=4,font=Env$police,selectmode="single",yscrollcommand=function(...) tkset(Env$l.fr3$noms.scroll,...))
  Env$l.fr3$noms.scroll<-tkscrollbar(Env$l.frames$Fr3,repeatinterval=4,command=function(...) tkyview(Env$l.fr3$noms.list,...))
  if (tclvalue(Env$l.var$moyprop)=="moy") {
    for (i in 1:length(Env$l.var$noms1)) {tkinsert(Env$l.fr3$noms.list,"end",Env$l.var$noms1[i])}
  } else {
    for (i in 1:length(Env$l.var$nomsprop.fac)) {tkinsert(Env$l.fr3$noms.list,"end",Env$l.var$nomsprop.fac[i])}
  }
  tkbind(Env$l.fr3$noms.list,"<ButtonRelease-1>",function() {
    tkdelete(Env$l.fr3$noms.wdg,0,"end")
    if (tclvalue(Env$l.var$moyprop)=="moy") {
	tkinsert(Env$l.fr3$noms.wdg,"end",Env$l.var$noms1[as.numeric(tclvalue(tkcurselection(Env$l.fr3$noms.list)))+1])
    } else {
	tkinsert(Env$l.fr3$noms.wdg,"end",Env$l.var$nomsprop.fac[as.numeric(tclvalue(tkcurselection(Env$l.fr3$noms.list)))+1])
    }
  })
  Env$l.fr3$noms.lab2<-tklabel(Env$l.frames$Fr3,text=Env$voc[90,1],font=Env$police)
  Env$l.fr3$noms.wdg<-tkentry(Env$l.frames$Fr3,width=20,font=Env$police)
  tkbind(Env$l.fr3$noms.wdg,"<Enter>",function() {msg(text=Env$voc[26,1],type="info")})
  tkbind(Env$l.fr3$noms.wdg,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.fr3$noms.wdg,"<ButtonRelease-1>",function() {tkdelete(Env$l.fr3$noms.wdg,0,"end")})
  tkbind(Env$l.fr3$noms.wdg,"<Return>",function() {
    if (tclvalue(Env$l.var$moyprop)=="moy") {
	rename.noms1(value.list=tclvalue(tkcurselection(Env$l.fr3$noms.list)),value.nom=tclvalue(tkget(Env$l.fr3$noms.wdg)))
    } else {
	rename.nomsprop.fac(value.list=tclvalue(tkcurselection(Env$l.fr3$noms.list)),value.nom=tclvalue(tkget(Env$l.fr3$noms.wdg)))
    }
  })
  Env$l.fr3$titre2<-tklabel(Env$l.frames$Fr3,text=Env$voc[52,1],font=Env$police3)
  Env$l.fr3$valeurs.titre.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[44,1],font=Env$police)
  Env$l.fr3$valeurs.titre.wdg<-tkentry(Env$l.frames$Fr3,width=20,font=Env$police,textvariable=Env$l.var$titre.axever)
  Env$l.fr3$valeurs.liminf.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[53,1],font=Env$police)
  Env$l.fr3$valeurs.liminf.wdg<-tkentry(Env$l.frames$Fr3,width=5,font=Env$police,textvariable=Env$l.var$liminf.axever)
  tkbind(Env$l.fr3$valeurs.liminf.wdg,"<Enter>",function() {msg(text=Env$voc[155,1],type="warning")})
  tkbind(Env$l.fr3$valeurs.liminf.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr3$valeurs.limsup.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[54,1],font=Env$police)
  Env$l.fr3$valeurs.limsup.wdg<-tkentry(Env$l.frames$Fr3,width=5,font=Env$police,textvariable=Env$l.var$limsup.axever)
  tkbind(Env$l.fr3$valeurs.limsup.wdg,"<Enter>",function() {msg(text=Env$voc[155,1],type="warning")})
  tkbind(Env$l.fr3$valeurs.limsup.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr3$valeurs.log.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[73,1],font=Env$police)
  Env$l.fr3$valeurs.log.wdg<-tkcheckbutton(Env$l.frames$Fr3,variable=Env$l.var$log.axever)
  tkbind(Env$l.fr3$valeurs.log.wdg,"<Enter>",function() {msg(text=Env$voc[158,1],type="warning")})
  tkbind(Env$l.fr3$valeurs.log.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr3$espace.ver1<-tklabel(Env$l.frames$Fr3,text="",font=Env$police2)
  Env$l.fr3$espace.ver2<-tklabel(Env$l.frames$Fr3,text="",font=Env$police2)
  Env$l.fr3$espace.hor<-tklabel(Env$l.frames$Fr3,text="                                   ",font=Env$police)
  tkgrid(Env$l.fr3$grad.lab1,row=0,column=0,sticky="e")
  tkgrid(Env$l.fr3$grad.col,row=0,column=1,sticky="w")
  tkgrid(Env$l.fr3$espace.hor,row=0,column=2)
  tkgrid(Env$l.fr3$grad.lab2,row=0,column=3,sticky="e")
  tkgrid(Env$l.fr3$grad.taille,row=0,column=4,sticky="w")
  tkgrid(Env$l.fr3$grad.lab3,row=1,column=0,sticky="e")
  tkgrid(Env$l.fr3$grad.orient,row=1,column=1,sticky="w")
  tkgrid(Env$l.fr3$leg.lab1,row=2,column=0,sticky="e")
  tkgrid(Env$l.fr3$leg.col,row=2,column=1,sticky="w")
  tkgrid(Env$l.fr3$leg.lab2,row=2,column=3,sticky="e")
  tkgrid(Env$l.fr3$leg.taille,row=2,column=4,sticky="w")
  tkgrid(Env$l.fr3$espace.ver1)
  tkgrid(Env$l.fr3$titre1,row=4,column=0,columnspan=2)
  tkgrid(Env$l.fr3$espace.ver2)
  tkgrid(Env$l.fr3$noms.titre.lab,row=6,column=0,sticky="e")
  tkgrid(Env$l.fr3$noms.titre.wdg,row=6,column=1,sticky="w")
  tkgrid(Env$l.fr3$noms.lab1,row=7,column=0,sticky="e")
  tkgrid(Env$l.fr3$noms.list,Env$l.fr3$noms.scroll,row=7,column=1,rowspan=4);tkgrid.configure(Env$l.fr3$noms.scroll,sticky="ens")
  tkgrid(Env$l.fr3$noms.lab2,row=11,column=0,sticky="e")
  tkgrid(Env$l.fr3$noms.wdg,row=11,column=1,sticky="w")
  tkgrid(Env$l.fr3$titre2,row=4,column=3,columnspan=2)
  tkgrid(Env$l.fr3$valeurs.titre.lab,row=6,column=3,sticky="e")
  tkgrid(Env$l.fr3$valeurs.titre.wdg,row=6,column=4,sticky="w")
  tkgrid(Env$l.fr3$valeurs.liminf.lab,row=7,column=3,sticky="e")
  tkgrid(Env$l.fr3$valeurs.liminf.wdg,row=7,column=4,sticky="w")
  tkgrid(Env$l.fr3$valeurs.limsup.lab,row=8,column=3,sticky="e")
  tkgrid(Env$l.fr3$valeurs.limsup.wdg,row=8,column=4,sticky="w")
  tkgrid(Env$l.fr3$valeurs.log.lab,row=9,column=3,sticky="e")
  tkgrid(Env$l.fr3$valeurs.log.wdg,row=9,column=4,sticky="w")
}

fr3.openCa<-function() {
  Env$l.frames$Fr3.status<-1
  tkconfigure(Env$l.wdg$but.lab3,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr3)) {tkdestroy(Env$l.fr3[[i]])}
  Env$l.fr3<-list()
  Env$l.fr3$noms.lab1<-tklabel(Env$l.frames$Fr3,text=Env$voc[116,1],font=Env$police)
  Env$l.fr3$noms.list<-tklistbox(Env$l.frames$Fr3,height=7,font=Env$police,selectmode="single",yscrollcommand=function(...) tkset(Env$l.fr3$noms.scroll,...))
  Env$l.fr3$noms.scroll<-tkscrollbar(Env$l.frames$Fr3,repeatinterval=5,command=function(...) tkyview(Env$l.fr3$noms.list,...))
  for (i in 1:length(Env$l.var$nomsparts)) {tkinsert(Env$l.fr3$noms.list,"end",Env$l.var$nomsparts[i])}
  tkbind(Env$l.fr3$noms.list,"<ButtonRelease-1>",function() {
    tkdelete(Env$l.fr3$noms.wdg,0,"end")
    tkinsert(Env$l.fr3$noms.wdg,"end",Env$l.var$nomsparts[as.numeric(tclvalue(tkcurselection(Env$l.fr3$noms.list)))+1])
  })
  Env$l.fr3$noms.lab2<-tklabel(Env$l.frames$Fr3,text=Env$voc[117,1],font=Env$police)
  Env$l.fr3$noms.wdg<-tkentry(Env$l.frames$Fr3,width=20,font=Env$police)
  tkbind(Env$l.fr3$noms.wdg,"<Enter>",function() {msg(text=Env$voc[26,1],type="info")})
  tkbind(Env$l.fr3$noms.wdg,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.fr3$noms.wdg,"<ButtonRelease-1>",function() {tkdelete(Env$l.fr3$noms.wdg,0,"end")})
  tkbind(Env$l.fr3$noms.wdg,"<Return>",rename.nomsparts)
  Env$l.fr3$orient.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[118,1],font=Env$police)
  Env$l.fr3$orient.wdg<-ttkcombobox(Env$l.frames$Fr3,width=35,values=Env$voc[119:120,1],textvariable=Env$l.var$cam.orient,font=Env$police,state="readonly")
  Env$l.fr3$start.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[121,1],font=Env$police)
  Env$l.fr3$start.wdg<-tkscale(Env$l.frames$Fr3,font=Env$police,from=0,to=380,showvalue=TRUE,variable=Env$l.var$cam.start,orient="horizontal",length=125)
  Env$l.fr3$lien.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[122,1],font=Env$police)
  Env$l.fr3$lien.wdg<-tkcheckbutton(Env$l.frames$Fr3,variable=Env$l.var$cam.lien,command=active.legende2)
  tkbind(Env$l.fr3$lien.wdg,"<Enter>",function() {msg(text=Env$voc[125,1],type="warning")})
  tkbind(Env$l.fr3$lien.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr3$espace.hor<-tklabel(Env$l.frames$Fr3,text="                              ",font=Env$police)
  tkgrid(Env$l.fr3$noms.lab1,row=0,column=0,sticky="e")
  tkgrid(Env$l.fr3$noms.list,Env$l.fr3$noms.scroll,row=0,column=1,rowspan=7,sticky="w");tkgrid.configure(Env$l.fr3$noms.scroll,sticky="ens")
  tkgrid(Env$l.fr3$noms.lab2,row=7,column=0,sticky="e")
  tkgrid(Env$l.fr3$noms.wdg,row=7,column=1,sticky="w")
  tkgrid(Env$l.fr3$espace.hor,row=0,column=2)
  tkgrid(Env$l.fr3$orient.lab,row=0,column=3,sticky="e")
  tkgrid(Env$l.fr3$orient.wdg,row=0,column=4,sticky="w")
  tkgrid(Env$l.fr3$start.lab,row=1,column=3,sticky="e")
  tkgrid(Env$l.fr3$start.wdg,row=1,column=4,sticky="w")
  tkgrid(Env$l.fr3$lien.lab,row=2,column=3,sticky="e")
  tkgrid(Env$l.fr3$lien.wdg,row=2,column=4,sticky="w")
}

fr3.openCoN<-function() {
  Env$l.frames$Fr3.status<-1
  tkconfigure(Env$l.wdg$but.lab3,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr3)) {tkdestroy(Env$l.fr3[[i]])}
  Env$l.fr3<-list()
  Env$l.fr3$grad.lab1<-tklabel(Env$l.frames$Fr3,text=Env$voc[47,1],font=Env$police)
  Env$l.fr3$grad.col<-tkbutton(Env$l.frames$Fr3,text="Aa",font=Env$police3,foreground=tclvalue(Env$l.var$graduations.col),activeforeground=tclvalue(Env$l.var$graduations.col),command=function() {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$graduations.col),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$graduations.col)<-temp
	tkconfigure(Env$l.fr3$grad.col,foreground=temp,activeforeground=temp)
    }
  })
  Env$l.fr3$grad.lab2<-tklabel(Env$l.frames$Fr3,text=Env$voc[48,1],font=Env$police)
  Env$l.fr3$grad.taille<-tkscale(Env$l.frames$Fr3,from=0.5,to=4,showvalue=TRUE,font=Env$police,variable=Env$l.var$graduations.taille,resolution=0.1,orient="horizontal")
  Env$l.fr3$grad.lab3<-tklabel(Env$l.frames$Fr3,text=Env$voc[245,1],font=Env$police)
  Env$l.fr3$grad.orient<-ttkcombobox(Env$l.frames$Fr3,width=16,values=c(Env$voc[246,1],Env$voc[67,1]),textvariable=Env$l.var$graduations.orient,font=Env$police,state="readonly")
  Env$l.fr3$leg.lab1<-tklabel(Env$l.frames$Fr3,text=Env$voc[49,1],font=Env$police)
  Env$l.fr3$leg.col<-tkbutton(Env$l.frames$Fr3,text="Aa",font=Env$police3,foreground=tclvalue(Env$l.var$legendes.col),activeforeground=tclvalue(Env$l.var$legendes.col),command=function() {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$legendes.col),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$legendes.col)<-temp
	tkconfigure(Env$l.fr3$leg.col,foreground=temp,activeforeground=temp)
    }
  })
  Env$l.fr3$leg.lab2<-tklabel(Env$l.frames$Fr3,text=Env$voc[50,1],font=Env$police)
  Env$l.fr3$leg.taille<-tkscale(Env$l.frames$Fr3,from=0.5,to=4,showvalue=TRUE,font=Env$police,variable=Env$l.var$legendes.taille,resolution=0.1,orient="horizontal")
  Env$l.fr3$titre1<-tklabel(Env$l.frames$Fr3,text=Env$voc[51,1],font=Env$police3)
  Env$l.fr3$hor.titre.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[44,1],font=Env$police)
  Env$l.fr3$hor.titre.wdg<-tkentry(Env$l.frames$Fr3,width=20,font=Env$police,textvariable=Env$l.var$titre.axehor)
  Env$l.fr3$hor.liminf.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[53,1],font=Env$police)
  Env$l.fr3$hor.liminf.wdg<-tkentry(Env$l.frames$Fr3,width=5,font=Env$police,textvariable=Env$l.var$liminf.axehor)
  tkbind(Env$l.fr3$hor.liminf.wdg,"<Enter>",function() {msg(text=Env$voc[155,1],type="warning")})
  tkbind(Env$l.fr3$hor.liminf.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr3$hor.limsup.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[54,1],font=Env$police)
  Env$l.fr3$hor.limsup.wdg<-tkentry(Env$l.frames$Fr3,width=5,font=Env$police,textvariable=Env$l.var$limsup.axehor)
  tkbind(Env$l.fr3$hor.limsup.wdg,"<Enter>",function() {msg(text=Env$voc[155,1],type="warning")})
  tkbind(Env$l.fr3$hor.limsup.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr3$hor.log.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[73,1],font=Env$police)
  Env$l.fr3$hor.log.wdg<-tkcheckbutton(Env$l.frames$Fr3,variable=Env$l.var$log.axehor)
  Env$l.fr3$titre2<-tklabel(Env$l.frames$Fr3,text=Env$voc[52,1],font=Env$police3)
  Env$l.fr3$ver.titre.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[44,1],font=Env$police)
  Env$l.fr3$ver.titre.wdg<-tkentry(Env$l.frames$Fr3,width=20,font=Env$police,textvariable=Env$l.var$titre.axever)
  Env$l.fr3$ver.liminf.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[53,1],font=Env$police)
  Env$l.fr3$ver.liminf.wdg<-tkentry(Env$l.frames$Fr3,width=5,font=Env$police,textvariable=Env$l.var$liminf.axever)
  tkbind(Env$l.fr3$ver.liminf.wdg,"<Enter>",function() {msg(text=Env$voc[155,1],type="warning")})
  tkbind(Env$l.fr3$ver.liminf.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr3$ver.limsup.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[54,1],font=Env$police)
  Env$l.fr3$ver.limsup.wdg<-tkentry(Env$l.frames$Fr3,width=5,font=Env$police,textvariable=Env$l.var$limsup.axever)
  tkbind(Env$l.fr3$ver.limsup.wdg,"<Enter>",function() {msg(text=Env$voc[155,1],type="warning")})
  tkbind(Env$l.fr3$ver.limsup.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr3$ver.log.lab<-tklabel(Env$l.frames$Fr3,text=Env$voc[73,1],font=Env$police)
  Env$l.fr3$ver.log.wdg<-tkcheckbutton(Env$l.frames$Fr3,variable=Env$l.var$log.axever)
  Env$l.fr3$espace.ver1<-tklabel(Env$l.frames$Fr3,text="",font=Env$police2)
  Env$l.fr3$espace.ver2<-tklabel(Env$l.frames$Fr3,text="",font=Env$police2)
  Env$l.fr3$espace.hor<-tklabel(Env$l.frames$Fr3,text="                                   ",font=Env$police)
  tkgrid(Env$l.fr3$grad.lab1,row=0,column=0,sticky="e")
  tkgrid(Env$l.fr3$grad.col,row=0,column=1,sticky="w")
  tkgrid(Env$l.fr3$espace.hor,row=0,column=2)
  tkgrid(Env$l.fr3$grad.lab2,row=0,column=3,sticky="e")
  tkgrid(Env$l.fr3$grad.taille,row=0,column=4,sticky="w")
  tkgrid(Env$l.fr3$grad.lab3,row=1,column=0,sticky="e")
  tkgrid(Env$l.fr3$grad.orient,row=1,column=1,sticky="w")
  tkgrid(Env$l.fr3$leg.lab1,row=2,column=0,sticky="e")
  tkgrid(Env$l.fr3$leg.col,row=2,column=1,sticky="w")
  tkgrid(Env$l.fr3$leg.lab2,row=2,column=3,sticky="e")
  tkgrid(Env$l.fr3$leg.taille,row=2,column=4,sticky="w")
  tkgrid(Env$l.fr3$espace.ver1)
  tkgrid(Env$l.fr3$titre1,row=4,column=0,columnspan=2)
  tkgrid(Env$l.fr3$espace.ver2)
  tkgrid(Env$l.fr3$hor.titre.lab,row=6,column=0,sticky="e")
  tkgrid(Env$l.fr3$hor.titre.wdg,row=6,column=1,sticky="w")
  tkgrid(Env$l.fr3$hor.liminf.lab,row=7,column=0,sticky="e")
  tkgrid(Env$l.fr3$hor.liminf.wdg,row=7,column=1,sticky="w")
  tkgrid(Env$l.fr3$hor.limsup.lab,row=8,column=0,sticky="e")
  tkgrid(Env$l.fr3$hor.limsup.wdg,row=8,column=1,sticky="w")
  tkgrid(Env$l.fr3$hor.log.lab,row=9,column=0,sticky="e")
  tkgrid(Env$l.fr3$hor.log.wdg,row=9,column=1,sticky="w")
  tkgrid(Env$l.fr3$titre2,row=4,column=4,columnspan=2)
  tkgrid(Env$l.fr3$ver.titre.lab,row=6,column=3,sticky="e")
  tkgrid(Env$l.fr3$ver.titre.wdg,row=6,column=4,sticky="w")
  tkgrid(Env$l.fr3$ver.liminf.lab,row=7,column=3,sticky="e")
  tkgrid(Env$l.fr3$ver.liminf.wdg,row=7,column=4,sticky="w")
  tkgrid(Env$l.fr3$ver.limsup.lab,row=8,column=3,sticky="e")
  tkgrid(Env$l.fr3$ver.limsup.wdg,row=8,column=4,sticky="w")
  tkgrid(Env$l.fr3$ver.log.lab,row=9,column=3,sticky="e")
  tkgrid(Env$l.fr3$ver.log.wdg,row=9,column=4,sticky="w")
}


#-------------------------------------------------
# Frame 4
#-------------------------------------------------

fr4.close<-function() {
  Env$l.frames$Fr4.status<-0
  tkconfigure(Env$l.wdg$but.lab4,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_bas.gif",fsep=.Platform$file.sep)))
  if ("l.hachures"%in%names(Env$l.fr4)) {
    for (i in 1:9) {tkdestroy(Env$l.fr4$l.hachures[[i]])}
  }
  if ("l.symboles"%in%names(Env$l.fr4)) {
    for (i in 1:8) {tkdestroy(Env$l.fr4$l.symboles[[i]])}
  }
  for (i in 1:length(Env$l.fr4)) {
    if (names(Env$l.fr4)[i]!="l.hachures" & names(Env$l.fr4)[i]!="l.symboles") {
	tkdestroy(Env$l.fr4[[i]])
    }
  }
  Env$l.fr4<-list()
  Env$l.fr4$vide<-tklabel(Env$l.frames$Fr4,text="",font=Env$police2)
  tkgrid(Env$l.fr4$vide)
}

fr4.openD<-function() {
  Env$l.frames$Fr4.status<-1
  tkconfigure(Env$l.wdg$but.lab4,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr4)) {tkdestroy(Env$l.fr4[[i]])}
  Env$l.fr4<-list()
  type.var<-"unknown"
  Env$l.fr4$titre1<-tklabel(Env$l.frames$Fr4,text=Env$voc[15,1],font=Env$police3)
  Env$l.fr4$var.list<-tklistbox(Env$l.frames$Fr4,height=7,font=Env$police,selectmode="single",yscrollcommand=function(...) tkset(Env$l.fr3$var.scroll,...))
  Env$l.fr4$var.scroll<-tkscrollbar(Env$l.frames$Fr4,repeatinterval=5,command=function(...) tkyview(Env$l.fr4$var.list,...))
  tkbind(Env$l.fr4$var.list,"<Enter>",function() {msg(text=Env$voc[34,1],type="warning")})
  tkbind(Env$l.fr4$var.list,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.fr4$var.list,"<ButtonRelease-1>",function() {
    type.var<-class(Env$dataset[,as.numeric(tclvalue(tkcurselection(Env$l.fr4$var.list)))+1])
    if (type.var=="numeric" | type.var=="integer") {
	tkconfigure(Env$l.fr4$rb.noregroup,state="normal",text=paste(Env$voc[27,1],nlevels(as.factor(Env$dataset[,as.numeric(tclvalue(tkcurselection(Env$l.fr4$var.list)))+1])),Env$voc[28,1],sep=""))
	tkconfigure(Env$l.fr4$rb.regroup1,state="normal")
	tkconfigure(Env$l.fr4$but,state="normal")
    } else if (type.var=="character") {
	tkconfigure(Env$l.fr4$rb.noregroup,state="disabled")
	tkconfigure(Env$l.fr4$rb.regroup1,state="disabled")
	tkconfigure(Env$l.fr4$rb.regroup2,state="disabled")
	tkconfigure(Env$l.fr4$rb.regroup3,state="disabled")
	tkconfigure(Env$l.fr4$curs.wdg,state="disabled",foreground="grey")
	tkconfigure(Env$l.fr4$curs.lab,foreground="grey")
	tkconfigure(Env$l.fr4$but,state="normal")
    } else {
	tkconfigure(Env$l.fr4$rb.noregroup,state="disabled")
	tkconfigure(Env$l.fr4$rb.regroup1,state="disabled")
	tkconfigure(Env$l.fr4$rb.regroup2,state="disabled")
	tkconfigure(Env$l.fr4$rb.regroup3,state="disabled")
	tkconfigure(Env$l.fr4$curs.wdg,state="disabled",foreground="grey")
	tkconfigure(Env$l.fr4$curs.lab,foreground="grey")
	tkconfigure(Env$l.fr4$but,state="disabled")
    }
  })
  Env$l.fr4$rb.noregroup<-tkradiobutton(Env$l.frames$Fr4,font=Env$police,variable=Env$l.var$regroup1,value=0,text=paste(Env$voc[27,1],0,Env$voc[28,1],sep=""),command=function() {
    tkconfigure(Env$l.fr4$rb.regroup2,state="disabled")
    tkconfigure(Env$l.fr4$rb.regroup3,state="disabled")
    tkconfigure(Env$l.fr4$curs.wdg,state="disabled",foreground="grey")
    tkconfigure(Env$l.fr4$curs.lab,foreground="grey")
  },state="disabled")
  Env$l.fr4$rb.regroup1<-tkradiobutton(Env$l.frames$Fr4,font=Env$police,variable=Env$l.var$regroup1,value=1,text=Env$voc[29,1],command=function() {
    tkconfigure(Env$l.fr4$rb.regroup2,state="normal")
    tkconfigure(Env$l.fr4$rb.regroup3,state="normal")
    tkconfigure(Env$l.fr4$curs.wdg,state="normal",foreground="black")
    tkconfigure(Env$l.fr4$curs.lab,foreground="black")
  },state="disabled")
  Env$l.fr4$rb.regroup2<-tkradiobutton(Env$l.frames$Fr4,font=Env$police,variable=Env$l.var$regroup2,value="long",text=Env$voc[30,1],state="disabled")
  Env$l.fr4$rb.regroup3<-tkradiobutton(Env$l.frames$Fr4,font=Env$police,variable=Env$l.var$regroup2,value="eff",text=Env$voc[31,1],state="disabled")
  Env$l.fr4$curs.wdg<-tkscale(Env$l.frames$Fr4,from=2,to=20,showvalue=TRUE,font=Env$police,variable=Env$l.var$regroup3,resolution=1,orient="horizontal",state="disabled",foreground="grey")
  Env$l.fr4$curs.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[32,1],foreground="grey")
  Env$l.fr4$but<-tkbutton(Env$l.frames$Fr4,text=Env$voc[33,1],width=16,state="disabled",command=function() {convert.variable(type.var)})
  if(!is.null(Env$dataset)) {
    for (i in 1:ncol(Env$dataset)) {tkinsert(Env$l.fr4$var.list,"end",colnames(Env$dataset)[i])}
  }
  Env$l.fr4$espace.ver<-tklabel(Env$l.frames$Fr4,text="",font=Env$police2)
  Env$l.fr4$espace.hor1<-tklabel(Env$l.frames$Fr4,text="                    ",font=Env$police)
  Env$l.fr4$espace.hor2<-tklabel(Env$l.frames$Fr4,text="               ",font=Env$police)
  tkgrid(Env$l.fr4$titre1,row=0,column=0)
  tkgrid(Env$l.fr4$espace.ver)
  tkgrid(Env$l.fr4$var.list,Env$l.fr4$var.scroll,row=2,column=0,rowspan=5);tkgrid.configure(Env$l.fr4$var.scroll,sticky="ens")
  tkgrid(Env$l.fr4$espace.hor1,row=2,column=1)
  tkgrid(Env$l.fr4$rb.noregroup,row=2,column=2,columnspan=2,sticky="w")
  tkgrid(Env$l.fr4$rb.regroup1,row=3,column=2,columnspan=2,sticky="w")
  tkgrid(Env$l.fr4$rb.regroup2,row=4,column=3,sticky="w")
  tkgrid(Env$l.fr4$rb.regroup3,row=5,column=3,sticky="w")
  tkgrid(Env$l.fr4$curs.wdg,row=4,column=4)
  tkgrid(Env$l.fr4$curs.lab,row=5,column=4)
  tkgrid(Env$l.fr4$espace.hor2,row=2,column=5)
  tkgrid(Env$l.fr4$but,row=2,column=6,rowspan=2)
}

fr4.openH<-function() {
  Env$l.frames$Fr4.status<-1
  tkconfigure(Env$l.wdg$but.lab4,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr4)) {tkdestroy(Env$l.fr4[[i]])}
  Env$l.fr4<-list()
  Env$l.fr4$nombre.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[55,1],font=Env$police)
  Env$l.fr4$nombre.wdg<-tkentry(Env$l.frames$Fr4,width=5,font=Env$police,textvariable=Env$l.var$hist.barres)
  tkbind(Env$l.fr4$nombre.wdg,"<Enter>",function() {msg(text=Env$voc[156,1],type="info")})
  tkbind(Env$l.fr4$nombre.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr4$colbarres.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[56,1],font=Env$police)
  Env$l.fr4$colbarres.wdg<-tkcanvas(Env$l.frames$Fr4,width="25",height="20",bg=tclvalue(Env$l.var$couleur1A))
  tkbind(Env$l.fr4$colbarres.wdg,"<ButtonRelease-1>",function() {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$couleur1A),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$couleur1A)<-temp
	tkconfigure(Env$l.fr4$colbarres.wdg,bg=temp)
    }
  })
  Env$l.fr4$colbordures.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[57,1],font=Env$police)
  Env$l.fr4$colbordures.wdg<-tkcanvas(Env$l.frames$Fr4,width="25",height="20",bg=tclvalue(Env$l.var$col.borduresA))
  tkbind(Env$l.fr4$colbordures.wdg,"<ButtonRelease-1>",function() {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$col.borduresA),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$col.borduresA)<-temp
	tkconfigure(Env$l.fr4$colbordures.wdg,bg=temp)
    }
  })
  Env$l.fr4$espace.hor1<-tklabel(Env$l.frames$Fr4,text="                              ",font=Env$police)
  Env$l.fr4$espace.hor2<-tklabel(Env$l.frames$Fr4,text="                              ",font=Env$police)
  tkgrid(Env$l.fr4$nombre.lab,row=0,column=0,sticky="e")
  tkgrid(Env$l.fr4$nombre.wdg,row=0,column=1,sticky="w")
  tkgrid(Env$l.fr4$espace.hor1,row=0,column=2)
  tkgrid(Env$l.fr4$colbarres.lab,row=0,column=3,sticky="e")
  tkgrid(Env$l.fr4$colbarres.wdg,row=0,column=4,sticky="w")
  tkgrid(Env$l.fr4$espace.hor2,row=0,column=5)
  tkgrid(Env$l.fr4$colbordures.lab,row=0,column=6,sticky="e")
  tkgrid(Env$l.fr4$colbordures.wdg,row=0,column=7,sticky="w")
}

fr4.openM<-function() {
  Env$l.frames$Fr4.status<-1
  tkconfigure(Env$l.wdg$but.lab4,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr4)) {tkdestroy(Env$l.fr4[[i]])}
  Env$l.fr4<-list()
  Env$l.fr4$noms.list<-tklistbox(Env$l.frames$Fr4,height=7,font=Env$police,selectmode="single",yscrollcommand=function(...) tkset(Env$l.fr4$noms.scroll,...))
  Env$l.fr4$noms.scroll<-tkscrollbar(Env$l.frames$Fr4,repeatinterval=5,command=function(...) tkyview(Env$l.fr4$noms.list,...))
  for (i in 1:length(Env$l.var$noms2)) {tkinsert(Env$l.fr4$noms.list,"end",Env$l.var$noms2[i])}
  tkbind(Env$l.fr4$noms.list,"<Enter>",function() {if (tclvalue(Env$l.var$plusieurs)==1) {msg(text=Env$voc[243,1],type="info")}})
  tkbind(Env$l.fr4$noms.list,"<Leave>",function() {if (tclvalue(Env$l.var$plusieurs)==1) {msg(text="",type="info")}})
  tkbind(Env$l.fr4$noms.list,"<ButtonRelease-1>",function() {
    tkconfigure(Env$l.fr4$colboites.wdg,bg=Env$l.var$couleur1B[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1])
    tkconfigure(Env$l.fr4$colbordures.wdg,bg=Env$l.var$col.borduresB[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1])
  })
  Env$l.fr4$colboites.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[74,1],font=Env$police)
  Env$l.fr4$colboites.wdg<-tkcanvas(Env$l.frames$Fr4,width="25",height="20",bg=tclvalue(Env$l.var$couleur1A))
  tkbind(Env$l.fr4$colboites.wdg,"<ButtonRelease-1>",col.boites)
  Env$l.fr4$colbordures.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[57,1],font=Env$police)
  Env$l.fr4$colbordures.wdg<-tkcanvas(Env$l.frames$Fr4,width="25",height="20",bg=tclvalue(Env$l.var$col.borduresA))
  if (tclvalue(Env$l.var$plusieurs)==0) {
    tkconfigure(Env$l.fr4$noms.list,state="disabled")
    tkconfigure(Env$l.fr4$colboites.wdg,bg=tclvalue(Env$l.var$couleur1A))
    tkconfigure(Env$l.fr4$colbordures.wdg,bg=tclvalue(Env$l.var$col.borduresA))
  } else {
    tkconfigure(Env$l.fr4$noms.list,state="normal")
    tkconfigure(Env$l.fr4$colboites.wdg,bg=Env$l.var$couleur1B[1])
    tkconfigure(Env$l.fr4$colbordures.wdg,bg=Env$l.var$col.borduresB[1])
    tkselection.set(Env$l.fr4$noms.list,"0")
  }
  tkbind(Env$l.fr4$colbordures.wdg,"<ButtonRelease-1>",col.bordures)
  Env$l.fr4$IC.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[75,1],font=Env$police)
  Env$l.fr4$IC.wdg<-tkcheckbutton(Env$l.frames$Fr4,variable=Env$l.var$ICmediane)
  tkbind(Env$l.fr4$IC.wdg,"<Enter>",function() {msg(text=Env$voc[76,1],type="info")})
  tkbind(Env$l.fr4$IC.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr4$varwidth.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[257,1],font=Env$police)
  Env$l.fr4$varwidth.wdg<-tkcheckbutton(Env$l.frames$Fr4,variable=Env$l.var$varwidth)
  tkbind(Env$l.fr4$varwidth.wdg,"<Enter>",function() {msg(text=Env$voc[258,1],type="info")})
  tkbind(Env$l.fr4$varwidth.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr4$moy.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[255,1],font=Env$police)
  Env$l.fr4$moy.wdg<-tkcheckbutton(Env$l.frames$Fr4,variable=Env$l.var$boxmoy)
  tkbind(Env$l.fr4$moy.wdg,"<Enter>",function() {msg(text=Env$voc[256,1],type="info")})
  tkbind(Env$l.fr4$moy.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr4$espace.hor1<-tklabel(Env$l.frames$Fr4,text="                    ",font=Env$police)
  Env$l.fr4$espace.hor2<-tklabel(Env$l.frames$Fr4,text="                    ",font=Env$police)
  tkgrid(Env$l.fr4$noms.list,Env$l.fr4$noms.scroll,row=0,column=0,rowspan=4,sticky="e");tkgrid.configure(Env$l.fr4$noms.scroll,sticky="ens")
  tkgrid(Env$l.fr4$espace.hor1,row=0,column=1)
  tkgrid(Env$l.fr4$colboites.lab,row=0,column=2,sticky="e")
  tkgrid(Env$l.fr4$colboites.wdg,row=0,column=3,sticky="w")
  tkgrid(Env$l.fr4$colbordures.lab,row=1,column=2,sticky="e")
  tkgrid(Env$l.fr4$colbordures.wdg,row=1,column=3,sticky="w")
  tkgrid(Env$l.fr4$espace.hor2,row=0,column=4)
  tkgrid(Env$l.fr4$IC.lab,row=0,column=5,sticky="e")
  tkgrid(Env$l.fr4$IC.wdg,row=0,column=6,sticky="w")
  tkgrid(Env$l.fr4$varwidth.lab,row=1,column=5,sticky="e")
  tkgrid(Env$l.fr4$varwidth.wdg,row=1,column=6,sticky="w")
  tkgrid(Env$l.fr4$moy.lab,row=2,column=5,sticky="e")
  tkgrid(Env$l.fr4$moy.wdg,row=2,column=6,sticky="w")
}

fr4.openB<-function() {
  Env$l.frames$Fr4.status<-1
  tkconfigure(Env$l.wdg$but.lab4,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr4)) {tkdestroy(Env$l.fr4[[i]])}
  Env$l.fr4<-list()
  Env$l.fr4$noms.list<-tklistbox(Env$l.frames$Fr4,height=7,font=Env$police,selectmode="single",yscrollcommand=function(...) tkset(Env$l.fr4$noms.scroll,...))
  Env$l.fr4$noms.scroll<-tkscrollbar(Env$l.frames$Fr4,repeatinterval=5,command=function(...) tkyview(Env$l.fr4$noms.list,...))
  if (tclvalue(Env$l.var$moyprop)=="moy") {
    for (i in 1:length(Env$l.var$noms2)) {tkinsert(Env$l.fr4$noms.list,"end",Env$l.var$noms2[i])}
  } else {
    for (i in 1:length(Env$l.var$nomsprop)) {tkinsert(Env$l.fr4$noms.list,"end",Env$l.var$nomsprop[i])}
  }
  tkbind(Env$l.fr4$noms.list,"<Enter>",function() {if (tclvalue(Env$l.var$plusieurs)==1) {msg(text=Env$voc[112,1],type="info")}})
  tkbind(Env$l.fr4$noms.list,"<Leave>",function() {if (tclvalue(Env$l.var$plusieurs)==1) {msg(text="",type="info")}})
  tkbind(Env$l.fr4$noms.list,"<ButtonRelease-1>",function() {
    if (tclvalue(Env$l.var$plusieurs)==1) {
	tkconfigure(Env$l.fr4$colbarres.wdg,bg=Env$l.var$couleur1B[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1])
	tkconfigure(Env$l.fr4$colbordures.wdg,bg=Env$l.var$col.borduresB[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1])
	for (i in 1:9) {tkconfigure(Env$l.fr4$l.hachures[[i]],borderwidth=0)}
	tkconfigure(Env$l.fr4$l.hachures[[Env$l.var$hachuresB[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1]]],borderwidth=2)
    }
  })
  Env$l.fr4$colbarres.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[91,1],font=Env$police)
  Env$l.fr4$colbarres.wdg<-tkcanvas(Env$l.frames$Fr4,width="25",height="20",bg=tclvalue(Env$l.var$couleur1A))
  tkbind(Env$l.fr4$colbarres.wdg,"<ButtonRelease-1>",col.barres)
  Env$l.fr4$colbordures.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[57,1],font=Env$police)
  Env$l.fr4$colbordures.wdg<-tkcanvas(Env$l.frames$Fr4,width="25",height="20",bg=tclvalue(Env$l.var$col.borduresA))
  tkbind(Env$l.fr4$colbordures.wdg,"<ButtonRelease-1>",col.bordures)
  Env$l.fr4$stack.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[93,1],font=Env$police)
  Env$l.fr4$stack.wdg<-tkcheckbutton(Env$l.frames$Fr4,variable=Env$l.var$stack,command=function() {if (tclvalue(Env$l.var$plusieurs)==1) {active.erreur()}})
  tkbind(Env$l.fr4$stack.wdg,"<Enter>",function() {msg(text=Env$voc[113,1],type="warning")})
  tkbind(Env$l.fr4$stack.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr4$hachures.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[92,1],font=Env$police)
  Env$l.fr4$l.hachures<-list()
  for (i in 1:9) {
    Env$l.fr4$l.hachures[[i]]<-tklabel(Env$l.frames$Fr4,height=35,width=35,font=Env$police,relief="groove",borderwidth=0,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",paste("Hachures",i,".gif",sep=""),fsep=.Platform$file.sep)))
  }
  if (tclvalue(Env$l.var$plusieurs)==0) {
    tkconfigure(Env$l.fr4$noms.list,state="disabled")
    tkconfigure(Env$l.fr4$colbarres.wdg,bg=tclvalue(Env$l.var$couleur1A))
    tkconfigure(Env$l.fr4$colbordures.wdg,bg=tclvalue(Env$l.var$col.borduresA))
    tkconfigure(Env$l.fr4$stack.lab,foreground="grey")
    tkconfigure(Env$l.fr4$stack.wdg,state="disabled")
    tkconfigure(Env$l.fr4$l.hachures[[as.numeric(tclvalue(Env$l.var$hachuresA))]],borderwidth=2)
  } else {
    tkconfigure(Env$l.fr4$noms.list,state="normal")
    tkconfigure(Env$l.fr4$colbarres.wdg,bg=Env$l.var$couleur1B[1])
    tkconfigure(Env$l.fr4$colbordures.wdg,bg=Env$l.var$col.borduresB[1])
    tkconfigure(Env$l.fr4$stack.lab,foreground="black")
    tkconfigure(Env$l.fr4$stack.wdg,state="normal")
    tkconfigure(Env$l.fr4$l.hachures[[Env$l.var$hachuresB[1]]],borderwidth=2)
    tkselection.set(Env$l.fr4$noms.list,"0")
  }
  tkbind(Env$l.fr4$l.hachures[[1]],"<ButtonRelease-1>",function() {hachures(num=1)})
  tkbind(Env$l.fr4$l.hachures[[2]],"<ButtonRelease-1>",function() {hachures(num=2)})
  tkbind(Env$l.fr4$l.hachures[[3]],"<ButtonRelease-1>",function() {hachures(num=3)})
  tkbind(Env$l.fr4$l.hachures[[4]],"<ButtonRelease-1>",function() {hachures(num=4)})
  tkbind(Env$l.fr4$l.hachures[[5]],"<ButtonRelease-1>",function() {hachures(num=5)})
  tkbind(Env$l.fr4$l.hachures[[6]],"<ButtonRelease-1>",function() {hachures(num=6)})
  tkbind(Env$l.fr4$l.hachures[[7]],"<ButtonRelease-1>",function() {hachures(num=7)})
  tkbind(Env$l.fr4$l.hachures[[8]],"<ButtonRelease-1>",function() {hachures(num=8)})
  tkbind(Env$l.fr4$l.hachures[[9]],"<ButtonRelease-1>",function() {hachures(num=9)})
  Env$l.fr4$espace.hor1<-tklabel(Env$l.frames$Fr4,text="                         ",font=Env$police)
  Env$l.fr4$espace.hor2<-tklabel(Env$l.frames$Fr4,text="                              ",font=Env$police)
  tkgrid(Env$l.fr4$noms.list,Env$l.fr4$noms.scroll,row=0,column=0,rowspan=4,sticky="e");tkgrid.configure(Env$l.fr4$noms.scroll,sticky="ens")
  tkgrid(Env$l.fr4$espace.hor1,row=0,column=1)
  tkgrid(Env$l.fr4$colbarres.lab,row=0,column=2,sticky="e")
  tkgrid(Env$l.fr4$colbarres.wdg,row=0,column=3,sticky="w")
  tkgrid(Env$l.fr4$colbordures.lab,row=1,column=2,sticky="e")
  tkgrid(Env$l.fr4$colbordures.wdg,row=1,column=3,sticky="w")
  tkgrid(Env$l.fr4$stack.lab,row=2,column=2,sticky="e")
  tkgrid(Env$l.fr4$stack.wdg,row=2,column=3,sticky="w")
  tkgrid(Env$l.fr4$espace.hor2,row=0,column=4)
  tkgrid(Env$l.fr4$hachures.lab,row=0,column=5,sticky="e")
  tkgrid(Env$l.fr4$l.hachures[[1]],row=0,column=6)
  tkgrid(Env$l.fr4$l.hachures[[2]],row=0,column=7)
  tkgrid(Env$l.fr4$l.hachures[[3]],row=0,column=8)
  tkgrid(Env$l.fr4$l.hachures[[4]],row=1,column=6)
  tkgrid(Env$l.fr4$l.hachures[[5]],row=1,column=7)
  tkgrid(Env$l.fr4$l.hachures[[6]],row=1,column=8)
  tkgrid(Env$l.fr4$l.hachures[[7]],row=2,column=6)
  tkgrid(Env$l.fr4$l.hachures[[8]],row=2,column=7)
  tkgrid(Env$l.fr4$l.hachures[[9]],row=2,column=8)
}

fr4.openCa<-function() {
  Env$l.frames$Fr4.status<-1
  tkconfigure(Env$l.wdg$but.lab4,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr4)) {tkdestroy(Env$l.fr4[[i]])}
  Env$l.fr4<-list()
  Env$l.fr4$noms.list<-tklistbox(Env$l.frames$Fr4,height=7,font=Env$police,selectmode="single",yscrollcommand=function(...) tkset(Env$l.fr4$noms.scroll,...))
  Env$l.fr4$noms.scroll<-tkscrollbar(Env$l.frames$Fr4,repeatinterval=5,command=function(...) tkyview(Env$l.fr4$noms.list,...))
  for (i in 1:length(Env$l.var$nomsparts)) {tkinsert(Env$l.fr4$noms.list,"end",Env$l.var$nomsparts[i])}
  tkselection.set(Env$l.fr4$noms.list,"0") 
  tkbind(Env$l.fr4$noms.list,"<Enter>",function() {msg(text=Env$voc[123,1],type="info")})
  tkbind(Env$l.fr4$noms.list,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.fr4$noms.list,"<ButtonRelease-1>",function() {
    tkconfigure(Env$l.fr4$colparts.wdg,bg=Env$l.var$couleur1B[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1])
    for (i in 1:9) {tkconfigure(Env$l.fr4$l.hachures[[i]],borderwidth=0)}
    tkconfigure(Env$l.fr4$l.hachures[[Env$l.var$hachuresB[as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1]]],borderwidth=2)
  })
  Env$l.fr4$colparts.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[124,1],font=Env$police)
  Env$l.fr4$colparts.wdg<-tkcanvas(Env$l.frames$Fr4,width="25",height="20",bg=Env$l.var$couleur1B[1])
  tkbind(Env$l.fr4$colparts.wdg,"<ButtonRelease-1>",col.parts)
  Env$l.fr4$colbordures.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[57,1],font=Env$police)
  Env$l.fr4$colbordures.wdg<-tkcanvas(Env$l.frames$Fr4,width="25",height="20",bg=tclvalue(Env$l.var$col.borduresA))
  tkbind(Env$l.fr4$colbordures.wdg,"<Enter>",function() {msg(text=Env$voc[126,1],type="info")})
  tkbind(Env$l.fr4$colbordures.wdg,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.fr4$colbordures.wdg,"<ButtonRelease-1>",function() {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$col.borduresA),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$col.borduresA)<-temp
	tkconfigure(Env$l.fr4$colbordures.wdg,bg=tclvalue(Env$l.var$col.borduresA))
    }
  })
  Env$l.fr4$hachures.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[92,1],font=Env$police)
  Env$l.fr4$l.hachures<-list()
  for (i in 1:9) {
    Env$l.fr4$l.hachures[[i]]<-tklabel(Env$l.frames$Fr4,height=35,width=35,font=Env$police,relief="groove",borderwidth=0,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",paste("Hachures",i,".gif",sep=""),fsep=.Platform$file.sep)))
  }
  tkconfigure(Env$l.fr4$l.hachures[[Env$l.var$hachuresB[1]]],borderwidth=2)
  tkbind(Env$l.fr4$l.hachures[[1]],"<ButtonRelease-1>",function() {hachures2(num=1)})
  tkbind(Env$l.fr4$l.hachures[[2]],"<ButtonRelease-1>",function() {hachures2(num=2)})
  tkbind(Env$l.fr4$l.hachures[[3]],"<ButtonRelease-1>",function() {hachures2(num=3)})
  tkbind(Env$l.fr4$l.hachures[[4]],"<ButtonRelease-1>",function() {hachures2(num=4)})
  tkbind(Env$l.fr4$l.hachures[[5]],"<ButtonRelease-1>",function() {hachures2(num=5)})
  tkbind(Env$l.fr4$l.hachures[[6]],"<ButtonRelease-1>",function() {hachures2(num=6)})
  tkbind(Env$l.fr4$l.hachures[[7]],"<ButtonRelease-1>",function() {hachures2(num=7)})
  tkbind(Env$l.fr4$l.hachures[[8]],"<ButtonRelease-1>",function() {hachures2(num=8)})
  tkbind(Env$l.fr4$l.hachures[[9]],"<ButtonRelease-1>",function() {hachures2(num=9)})
  Env$l.fr4$espace.hor1<-tklabel(Env$l.frames$Fr4,text="                         ",font=Env$police)
  Env$l.fr4$espace.hor2<-tklabel(Env$l.frames$Fr4,text="                              ",font=Env$police)
  tkgrid(Env$l.fr4$noms.list,Env$l.fr4$noms.scroll,row=0,column=0,rowspan=7,sticky="e");tkgrid.configure(Env$l.fr4$noms.scroll,sticky="ens")
  tkgrid(Env$l.fr4$espace.hor1,row=0,column=1)
  tkgrid(Env$l.fr4$colparts.lab,row=0,column=2,sticky="e")
  tkgrid(Env$l.fr4$colparts.wdg,row=0,column=3,sticky="w")
  tkgrid(Env$l.fr4$colbordures.lab,row=1,column=2,sticky="e")
  tkgrid(Env$l.fr4$colbordures.wdg,row=1,column=3,sticky="w")
  tkgrid(Env$l.fr4$espace.hor2,row=0,column=4)
  tkgrid(Env$l.fr4$hachures.lab,row=0,column=5,sticky="e")
  tkgrid(Env$l.fr4$l.hachures[[1]],row=0,column=6)
  tkgrid(Env$l.fr4$l.hachures[[2]],row=0,column=7)
  tkgrid(Env$l.fr4$l.hachures[[3]],row=0,column=8)
  tkgrid(Env$l.fr4$l.hachures[[4]],row=1,column=6)
  tkgrid(Env$l.fr4$l.hachures[[5]],row=1,column=7)
  tkgrid(Env$l.fr4$l.hachures[[6]],row=1,column=8)
  tkgrid(Env$l.fr4$l.hachures[[7]],row=2,column=6)
  tkgrid(Env$l.fr4$l.hachures[[8]],row=2,column=7)
  tkgrid(Env$l.fr4$l.hachures[[9]],row=2,column=8)
}

fr4.openCo<-function() {
  Env$l.frames$Fr4.status<-1
  tkconfigure(Env$l.wdg$but.lab4,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr4)) {tkdestroy(Env$l.fr4[[i]])}
  Env$l.fr4<-list()
  Env$l.fr4$noms.list<-tklistbox(Env$l.frames$Fr4,height=10,font=Env$police,selectmode="single",yscrollcommand=function(...) tkset(Env$l.fr4$noms.scroll,...))
  Env$l.fr4$noms.scroll<-tkscrollbar(Env$l.frames$Fr4,repeatinterval=5,command=function(...) tkyview(Env$l.fr4$noms.list,...))
  tkbind(Env$l.fr4$noms.list,"<Enter>",function() {if (tclvalue(Env$l.var$plusieurs)==1) {msg(text=Env$voc[141,1],type="info")}})
  tkbind(Env$l.fr4$noms.list,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.fr4$noms.list,"<ButtonRelease-1>",function() {
    if (tclvalue(Env$l.var$plusieurs)==1) {
	Env$l.var$select<-as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1
	for (i in 1:8) {tkconfigure(Env$l.fr4$l.symboles[[i]],borderwidth=0)}
	tkconfigure(Env$l.fr4$l.symboles[[Env$l.var$symboleB[Env$l.var$select]]],borderwidth=2)
	tkconfigure(Env$l.fr4$col.wdg,bg=Env$l.var$couleur2B[Env$l.var$select])
	tclvalue(Env$l.var$taille.ptsA)<-as.character(Env$l.var$taille.ptsB[Env$l.var$select])
	tclvalue(Env$l.var$type.courbeA)<-as.character(Env$l.var$type.courbeB[Env$l.var$select])
	tclvalue(Env$l.var$trait1)<-as.character(Env$l.var$trait2[Env$l.var$select])
	tclvalue(Env$l.var$epaisseur1)<-as.character(Env$l.var$epaisseur2[Env$l.var$select])
    }
  })
  Env$l.fr4$symboles.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[138,1],font=Env$police)
  Env$l.fr4$l.symboles<-list()
  for (i in 1:8) {
    Env$l.fr4$l.symboles[[i]]<-tklabel(Env$l.frames$Fr4,height=35,width=35,font=Env$police,relief="groove",borderwidth=0,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",paste("Symbole",i,".gif",sep=""),fsep=.Platform$file.sep)))
  }
  tkbind(Env$l.fr4$l.symboles[[1]],"<ButtonRelease-1>",function() {symboles(num=1)})
  tkbind(Env$l.fr4$l.symboles[[2]],"<ButtonRelease-1>",function() {symboles(num=2)})
  tkbind(Env$l.fr4$l.symboles[[3]],"<ButtonRelease-1>",function() {symboles(num=3)})
  tkbind(Env$l.fr4$l.symboles[[4]],"<ButtonRelease-1>",function() {symboles(num=4)})
  tkbind(Env$l.fr4$l.symboles[[5]],"<ButtonRelease-1>",function() {symboles(num=5)})
  tkbind(Env$l.fr4$l.symboles[[6]],"<ButtonRelease-1>",function() {symboles(num=6)})
  tkbind(Env$l.fr4$l.symboles[[7]],"<ButtonRelease-1>",function() {symboles(num=7)})
  tkbind(Env$l.fr4$l.symboles[[8]],"<ButtonRelease-1>",function() {symboles(num=8)})
  Env$l.fr4$col.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[45,1],font=Env$police)
  Env$l.fr4$col.wdg<-tkcanvas(Env$l.frames$Fr4,width="25",height="20")
  tkbind(Env$l.fr4$col.wdg,"<ButtonRelease-1>",col.symboles)
  Env$l.fr4$taille.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[131,1],font=Env$police)
  Env$l.fr4$taille.wdg<-tkscale(Env$l.frames$Fr4,showvalue=TRUE,from=0.5,to=3,resolution=0.1,font=Env$police,variable=Env$l.var$taille.ptsA,orient="horizontal",command=function(...) {
    if (tclvalue(Env$l.var$plusieurs)==1) {
	Env$l.var$taille.ptsB[Env$l.var$select]<-as.numeric(tclvalue(Env$l.var$taille.ptsA))
	tkselection.set(Env$l.fr4$noms.list,as.character(Env$l.var$select-1))
    }
  })
  Env$l.fr4$type.courbe.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[132,1],font=Env$police)
  Env$l.fr4$type.courbe.wdg<-ttkcombobox(Env$l.frames$Fr4,values=Env$voc[133:137,1],textvariable=Env$l.var$type.courbeA,state="readonly",font=Env$police)
  tkbind(Env$l.fr4$type.courbe.wdg,"<<ComboboxSelected>>",function() {
    if (tclvalue(Env$l.var$plusieurs)==1) {
	Env$l.var$type.courbeB[Env$l.var$select]<-tclvalue(Env$l.var$type.courbeA)
	tkselection.set(Env$l.fr4$noms.list,as.character(Env$l.var$select-1))
    }
  })
  Env$l.fr4$type.trait.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[59,1],font=Env$police)
  Env$l.fr4$type.trait.wdg<-ttkcombobox(Env$l.frames$Fr4,values=Env$voc[60:62,1],textvariable=Env$l.var$trait1,state="readonly",font=Env$police)
  tkbind(Env$l.fr4$type.trait.wdg,"<<ComboboxSelected>>",function() {
    if (tclvalue(Env$l.var$plusieurs)==1) {
	Env$l.var$trait2[Env$l.var$select]<-tclvalue(Env$l.var$trait1)
	tkselection.set(Env$l.fr4$noms.list,as.character(Env$l.var$select-1))
    }
  })
  Env$l.fr4$epaisseur.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[63,1],font=Env$police)
  Env$l.fr4$epaisseur.wdg<-tkscale(Env$l.frames$Fr4,showvalue=TRUE,from=1,to=4,resolution=1,font=Env$police,variable=Env$l.var$epaisseur1,orient="horizontal",command=function(...) {
    if (tclvalue(Env$l.var$plusieurs)==1) {
	Env$l.var$epaisseur2[Env$l.var$select]<-as.numeric(tclvalue(Env$l.var$epaisseur1))
	tkselection.set(Env$l.fr4$noms.list,as.character(Env$l.var$select-1))
    }
  })
  if (tclvalue(Env$l.var$plusieurs)==1) {
    Env$l.var$select<-1
    for (i in 1:length(Env$l.var$noms1)) {tkinsert(Env$l.fr4$noms.list,"end",Env$l.var$noms1[i])}
    tkselection.set(Env$l.fr4$noms.list,"0")
    tkconfigure(Env$l.fr4$l.symboles[[Env$l.var$symboleB[1]]],borderwidth=2)
    tkconfigure(Env$l.fr4$col.wdg,bg=Env$l.var$couleur2B[1])
    tclvalue(Env$l.var$taille.ptsA)<-as.character(Env$l.var$taille.ptsB[1])
    tclvalue(Env$l.var$type.courbeA)<-as.character(Env$l.var$type.courbeB[1])
    tclvalue(Env$l.var$trait1)<-as.character(Env$l.var$trait2[1])
    tclvalue(Env$l.var$epaisseur1)<-as.character(Env$l.var$epaisseur2[1])
  } else {
    tkconfigure(Env$l.fr4$noms.list,state="disabled")
    tkconfigure(Env$l.fr4$l.symboles[[as.numeric(tclvalue(Env$l.var$symboleA))]],borderwidth=2)
    tkconfigure(Env$l.fr4$col.wdg,bg=tclvalue(Env$l.var$couleur2A))
  }
  Env$l.fr4$espace.hor1<-tklabel(Env$l.frames$Fr4,text="          ",font=Env$police)
  Env$l.fr4$espace.hor2<-tklabel(Env$l.frames$Fr4,text="                    ",font=Env$police)
  tkgrid(Env$l.fr4$noms.list,Env$l.fr4$noms.scroll,row=0,column=0,rowspan=10,sticky="w");tkgrid.configure(Env$l.fr4$noms.scroll,sticky="ens")
  tkgrid(Env$l.fr4$espace.hor1,row=0,column=1)
  tkgrid(Env$l.fr4$symboles.lab,row=0,column=2,sticky="e")
  tkgrid(Env$l.fr4$l.symboles[[1]],row=0,column=3)
  tkgrid(Env$l.fr4$l.symboles[[2]],row=0,column=4)
  tkgrid(Env$l.fr4$l.symboles[[3]],row=0,column=5)
  tkgrid(Env$l.fr4$l.symboles[[4]],row=0,column=6)
  tkgrid(Env$l.fr4$l.symboles[[5]],row=1,column=3)
  tkgrid(Env$l.fr4$l.symboles[[6]],row=1,column=4)
  tkgrid(Env$l.fr4$l.symboles[[7]],row=1,column=5)
  tkgrid(Env$l.fr4$l.symboles[[8]],row=1,column=6)
  tkgrid(Env$l.fr4$col.lab,row=2,column=2,sticky="e")
  tkgrid(Env$l.fr4$col.wdg,row=2,column=3,columnspan=4,sticky="w")
  tkgrid(Env$l.fr4$taille.lab,row=3,column=2,sticky="e")
  tkgrid(Env$l.fr4$taille.wdg,row=3,column=3,columnspan=4,sticky="w")
  tkgrid(Env$l.fr4$espace.hor2,row=0,column=7)
  tkgrid(Env$l.fr4$type.courbe.lab,row=0,column=8,sticky="e")
  tkgrid(Env$l.fr4$type.courbe.wdg,row=0,column=9,sticky="w")
  tkgrid(Env$l.fr4$type.trait.lab,row=1,column=8,sticky="e")
  tkgrid(Env$l.fr4$type.trait.wdg,row=1,column=9,sticky="w")
  tkgrid(Env$l.fr4$epaisseur.lab,row=2,column=8,sticky="e")
  tkgrid(Env$l.fr4$epaisseur.wdg,row=2,column=9,sticky="w")
}

fr4.openN<-function() {
  Env$l.frames$Fr4.status<-1
  tkconfigure(Env$l.wdg$but.lab4,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr4)) {tkdestroy(Env$l.fr4[[i]])}
  Env$l.fr4<-list()
  Env$l.fr4$noms.list<-tklistbox(Env$l.frames$Fr4,height=7,font=Env$police,selectmode="single",yscrollcommand=function(...) tkset(Env$l.fr4$noms.scroll,...))
  Env$l.fr4$noms.scroll<-tkscrollbar(Env$l.frames$Fr4,repeatinterval=5,command=function(...) tkyview(Env$l.fr4$noms.list,...))
  tkbind(Env$l.fr4$noms.list,"<Enter>",function() {if (tclvalue(Env$l.var$plusieurs)==1) {msg(text=Env$voc[141,1],type="info")}})
  tkbind(Env$l.fr4$noms.list,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.fr4$noms.list,"<ButtonRelease-1>",function() {
    if (tclvalue(Env$l.var$plusieurs)==1) {
	Env$l.var$select<-as.numeric(tclvalue(tkcurselection(Env$l.fr4$noms.list)))+1
	for (i in 1:8) {tkconfigure(Env$l.fr4$l.symboles[[i]],borderwidth=0)}
	tkconfigure(Env$l.fr4$l.symboles[[Env$l.var$symboleB[Env$l.var$select]]],borderwidth=2)
	tkconfigure(Env$l.fr4$col.wdg,bg=Env$l.var$couleur2B[Env$l.var$select])
	tclvalue(Env$l.var$taille.ptsA)<-as.character(Env$l.var$taille.ptsB[Env$l.var$select])
    }
  })
  Env$l.fr4$symboles.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[138,1],font=Env$police)
  Env$l.fr4$l.symboles<-list()
  for (i in 1:8) {
    Env$l.fr4$l.symboles[[i]]<-tklabel(Env$l.frames$Fr4,height=35,width=35,font=Env$police,relief="groove",borderwidth=0,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",paste("Symbole",i,".gif",sep=""),fsep=.Platform$file.sep)))
  }
  tkbind(Env$l.fr4$l.symboles[[1]],"<ButtonRelease-1>",function() {symboles(num=1)})
  tkbind(Env$l.fr4$l.symboles[[2]],"<ButtonRelease-1>",function() {symboles(num=2)})
  tkbind(Env$l.fr4$l.symboles[[3]],"<ButtonRelease-1>",function() {symboles(num=3)})
  tkbind(Env$l.fr4$l.symboles[[4]],"<ButtonRelease-1>",function() {symboles(num=4)})
  tkbind(Env$l.fr4$l.symboles[[5]],"<ButtonRelease-1>",function() {symboles(num=5)})
  tkbind(Env$l.fr4$l.symboles[[6]],"<ButtonRelease-1>",function() {symboles(num=6)})
  tkbind(Env$l.fr4$l.symboles[[7]],"<ButtonRelease-1>",function() {symboles(num=7)})
  tkbind(Env$l.fr4$l.symboles[[8]],"<ButtonRelease-1>",function() {symboles(num=8)})
  Env$l.fr4$col.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[45,1],font=Env$police)
  Env$l.fr4$col.wdg<-tkcanvas(Env$l.frames$Fr4,width="25",height="20")
  tkbind(Env$l.fr4$col.wdg,"<ButtonRelease-1>",col.symboles)
  Env$l.fr4$taille.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[131,1],font=Env$police)
  Env$l.fr4$taille.wdg<-tkscale(Env$l.frames$Fr4,showvalue=TRUE,from=0.5,to=3,resolution=0.1,font=Env$police,variable=Env$l.var$taille.ptsA,orient="horizontal",command=function(...) {
    if (tclvalue(Env$l.var$plusieurs)==1) {
	Env$l.var$taille.ptsB[Env$l.var$select]<-as.numeric(tclvalue(Env$l.var$taille.ptsA))
	tkselection.set(Env$l.fr4$noms.list,as.character(Env$l.var$select-1))
    }
  })
  if (tclvalue(Env$l.var$plusieurs)==1) {
    Env$l.var$select<-1
    for (i in 1:length(Env$l.var$noms1)) {tkinsert(Env$l.fr4$noms.list,"end",Env$l.var$noms1[i])}
    tkselection.set(Env$l.fr4$noms.list,"0")
    tkconfigure(Env$l.fr4$l.symboles[[Env$l.var$symboleB[1]]],borderwidth=2)
    tkconfigure(Env$l.fr4$col.wdg,bg=Env$l.var$couleur2B[1])
    tclvalue(Env$l.var$taille.ptsA)<-as.character(Env$l.var$taille.ptsB[1])
  } else {
    tkconfigure(Env$l.fr4$noms.list,state="disabled")
    tkconfigure(Env$l.fr4$l.symboles[[as.numeric(tclvalue(Env$l.var$symboleA))]],borderwidth=2)
    tkconfigure(Env$l.fr4$col.wdg,bg=tclvalue(Env$l.var$couleur2A))
  }
  Env$l.fr4$ptlab.lab<-tklabel(Env$l.frames$Fr4,text=Env$voc[242,1],font=Env$police)
  Env$l.fr4$ptlab.wdg<-tkcheckbutton(Env$l.frames$Fr4,variable=Env$l.var$ptlab)
  Env$l.fr4$espace.hor1<-tklabel(Env$l.frames$Fr4,text="                    ",font=Env$police)
  Env$l.fr4$espace.hor2<-tklabel(Env$l.frames$Fr4,text="                    ",font=Env$police)
  tkgrid(Env$l.fr4$noms.list,Env$l.fr4$noms.scroll,row=0,column=0,rowspan=3,sticky="w");tkgrid.configure(Env$l.fr4$noms.scroll,sticky="ens")
  tkgrid(Env$l.fr4$espace.hor1,row=0,column=1)
  tkgrid(Env$l.fr4$symboles.lab,row=0,column=2,sticky="e")
  tkgrid(Env$l.fr4$l.symboles[[1]],row=0,column=3)
  tkgrid(Env$l.fr4$l.symboles[[2]],row=0,column=4)
  tkgrid(Env$l.fr4$l.symboles[[3]],row=0,column=5)
  tkgrid(Env$l.fr4$l.symboles[[4]],row=0,column=6)
  tkgrid(Env$l.fr4$l.symboles[[5]],row=1,column=3)
  tkgrid(Env$l.fr4$l.symboles[[6]],row=1,column=4)
  tkgrid(Env$l.fr4$l.symboles[[7]],row=1,column=5)
  tkgrid(Env$l.fr4$l.symboles[[8]],row=1,column=6)
  tkgrid(Env$l.fr4$col.lab,row=2,column=2,sticky="e")
  tkgrid(Env$l.fr4$col.wdg,row=2,column=3,columnspan=4,sticky="w")
  tkgrid(Env$l.fr4$espace.hor2,row=0,column=7)
  tkgrid(Env$l.fr4$taille.lab,row=0,column=8,sticky="e")
  tkgrid(Env$l.fr4$taille.wdg,row=0,column=9,sticky="w")
  tkgrid(Env$l.fr4$ptlab.lab,row=1,column=8,sticky="e")
  tkgrid(Env$l.fr4$ptlab.wdg,row=1,column=9,sticky="w")
}


#-------------------------------------------------
# Frame 5
#-------------------------------------------------

fr5.close<-function() {
  Env$l.frames$Fr5.status<-0
  tkconfigure(Env$l.wdg$but.lab5,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_bas.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr5)) {tkdestroy(Env$l.fr5[[i]])}
  Env$l.fr5<-list()
  Env$l.fr5$vide<-tklabel(Env$l.frames$Fr5,text="",font=Env$police2)
  tkgrid(Env$l.fr5$vide)
}

fr5.openD<-function() {
  Env$l.frames$Fr5.status<-1
  tkconfigure(Env$l.wdg$but.lab5,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr5)) {tkdestroy(Env$l.fr5[[i]])}
  Env$l.fr5<-list()
  Env$l.var$levels.temp <- NULL
  Env$l.fr5$fact.lab <- tklabel(Env$l.frames$Fr5,text=Env$voc[249,1],font=Env$police)
  Env$l.fr5$fact.wdg <- ttkcombobox(Env$l.frames$Fr5,values=Env$l.var$var.fact,textvariable=Env$l.var$facteur1,state="readonly")
  tkbind(Env$l.fr5$fact.wdg,"<<ComboboxSelected>>",function() {
    tkdelete(Env$l.fr5$liste.actual,0,"end")
    for (i in 1:nlevels(Env$dataset[,tclvalue(Env$l.var$facteur1)])) {
	tkinsert(Env$l.fr5$liste.actual,"end",levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[i])
    }
    tkdelete(Env$l.fr5$liste.new,0,"end")
  })
  Env$l.fr5$titre1 <- tklabel(Env$l.frames$Fr5,text=Env$voc[250,1],font=Env$police3)
  Env$l.fr5$liste.actual <- tklistbox(Env$l.frames$Fr5,selectmode="single",height=5)
  Env$l.fr5$titre2 <- tklabel(Env$l.frames$Fr5,text=Env$voc[251,1],font=Env$police3)
  Env$l.fr5$liste.new <- tklistbox(Env$l.frames$Fr5,selectmode="single",height=5)
  Env$l.fr5$but1 <- tkbutton(Env$l.frames$Fr5,text=">",width=5,command=function() {
    if (tclvalue(tkcurselection(Env$l.fr5$liste.actual))!="") {
	selection <- as.numeric(tclvalue(tkcurselection(Env$l.fr5$liste.actual)))+1
	tkinsert(Env$l.fr5$liste.new,"end",levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[selection])
	Env$l.var$levels.temp <- c(Env$l.var$levels.temp,levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[selection])
    }
  })
  Env$l.fr5$but2 <- tkbutton(Env$l.frames$Fr5,text="<",width=5,command=function() {
    if (tclvalue(tkcurselection(Env$l.fr5$liste.new))!="") {
	selection <- tclvalue(tkcurselection(Env$l.fr5$liste.new))
	tkdelete(Env$l.fr5$liste.new,selection)
	Env$l.var$levels.temp <- Env$l.var$levels.temp[-(as.numeric(selection)+1)]
    }
  })
  Env$l.fr5$but3 <- tkbutton(Env$l.frames$Fr5,text=Env$voc[252,1],width=16,command=function() {
    if (length(Env$l.var$levels.temp)==nlevels(Env$dataset[,tclvalue(Env$l.var$facteur1)])) {
	Env$dataset[,tclvalue(Env$l.var$facteur1)] <- factor(Env$dataset[,tclvalue(Env$l.var$facteur1)],levels=Env$l.var$levels.temp)
	msg(text=Env$voc[254,1],type="info")
    } else {
	msg(text=Env$voc[253,1],type="error")
    }
  })
  Env$l.fr5$espace.hor1 <- tklabel(Env$l.frames$Fr5,text="          ",font=Env$police)
  Env$l.fr5$espace.hor2 <- tklabel(Env$l.frames$Fr5,text="          ",font=Env$police)
  tkgrid(Env$l.fr5$fact.lab,row=1,column=0,sticky="ne")
  tkgrid(Env$l.fr5$fact.wdg,row=1,column=1,sticky="nw")
  tkgrid(Env$l.fr5$espace.hor1,row=1,column=2)
  tkgrid(Env$l.fr5$titre1,row=0,column=3)
  tkgrid(Env$l.fr5$liste.actual,row=1,column=3,rowspan=2)
  tkgrid(Env$l.fr5$but1,row=1,column=4)
  tkgrid(Env$l.fr5$but2,row=2,column=4)
  tkgrid(Env$l.fr5$titre2,row=0,column=5)
  tkgrid(Env$l.fr5$liste.new,row=1,column=5,rowspan=2)
  tkgrid(Env$l.fr5$espace.hor2,row=1,column=6)
  tkgrid(Env$l.fr5$but3,row=1,column=7,sticky="n")
}

fr5.openH<-function() {
  Env$l.frames$Fr5.status<-1
  tkconfigure(Env$l.wdg$but.lab5,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr5)) {tkdestroy(Env$l.fr5[[i]])}
  Env$l.fr5<-list()
  Env$l.fr5$tracer.lab<-tklabel(Env$l.frames$Fr5,text=Env$voc[58,1],font=Env$police,foreground=ifelse(tclvalue(Env$l.var$hist.type)==Env$voc[42,1],"black","grey"))
  Env$l.fr5$tracer.wdg<-tkcheckbutton(Env$l.frames$Fr5,variable=Env$l.var$hist.dens,state=ifelse(tclvalue(Env$l.var$hist.type)==Env$voc[42,1],"normal","disabled"))
  Env$l.fr5$col.lab<-tklabel(Env$l.frames$Fr5,text=Env$voc[45,1],font=Env$police,foreground=ifelse(tclvalue(Env$l.var$hist.type)==Env$voc[42,1],"black","grey"))
  Env$l.fr5$col.wdg<-tkcanvas(Env$l.frames$Fr5,width="25",height="20",bg=ifelse(tclvalue(Env$l.var$hist.type)==Env$voc[42,1],tclvalue(Env$l.var$couleur2A),"grey"))
  tkbind(Env$l.fr5$col.wdg,"<ButtonRelease-1>",function() {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$couleur2A),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$couleur2A)<-temp
	tkconfigure(Env$l.fr5$col.wdg,bg=temp)
    }
  })
  Env$l.fr5$trait.lab<-tklabel(Env$l.frames$Fr5,text=Env$voc[59,1],font=Env$police,foreground=ifelse(tclvalue(Env$l.var$hist.type)==Env$voc[42,1],"black","grey"))
  Env$l.fr5$trait.wdg<-ttkcombobox(Env$l.frames$Fr5,font=Env$police,values=Env$voc[60:62,1],textvariable=Env$l.var$trait1,state=ifelse(tclvalue(Env$l.var$hist.type)==Env$voc[42,1],"readonly","disabled"),width=15)
  Env$l.fr5$epaisseur.lab<-tklabel(Env$l.frames$Fr5,text=Env$voc[63,1],font=Env$police,foreground=ifelse(tclvalue(Env$l.var$hist.type)==Env$voc[42,1],"black","grey"))
  Env$l.fr5$epaisseur.wdg<-tkscale(Env$l.frames$Fr5,from=1,to=5,showvalue=TRUE,font=Env$police,variable=Env$l.var$epaisseur1,resolution=1,orient="horizontal",state=ifelse(tclvalue(Env$l.var$hist.type)==Env$voc[42,1],"normal","disabled"))
  Env$l.fr5$espace.hor1<-tklabel(Env$l.frames$Fr5,text="               ",font=Env$police)
  Env$l.fr5$espace.hor2<-tklabel(Env$l.frames$Fr5,text="               ",font=Env$police)
  Env$l.fr5$espace.hor3<-tklabel(Env$l.frames$Fr5,text="               ",font=Env$police)
  tkgrid(Env$l.fr5$tracer.lab,row=0,column=0,sticky="e")
  tkgrid(Env$l.fr5$tracer.wdg,row=0,column=1,sticky="w")
  tkgrid(Env$l.fr5$espace.hor1,row=0,column=2)
  tkgrid(Env$l.fr5$col.lab,row=0,column=3,sticky="e")
  tkgrid(Env$l.fr5$col.wdg,row=0,column=4,sticky="w")
  tkgrid(Env$l.fr5$espace.hor2,row=0,column=5)
  tkgrid(Env$l.fr5$trait.lab,row=0,column=6,sticky="e")
  tkgrid(Env$l.fr5$trait.wdg,row=0,column=7,sticky="w")
  tkgrid(Env$l.fr5$espace.hor3,row=0,column=8)
  tkgrid(Env$l.fr5$epaisseur.lab,row=0,column=9,sticky="e")
  tkgrid(Env$l.fr5$epaisseur.wdg,row=0,column=10,sticky="w")
}

fr5.openM<-function() {
  Env$l.frames$Fr5.status<-1
  tkconfigure(Env$l.wdg$but.lab5,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr5)) {tkdestroy(Env$l.fr5[[i]])}
  Env$l.fr5<-list()
  Env$l.fr5$taille.lab<-tklabel(Env$l.frames$Fr5,text=Env$voc[46,1],font=Env$police)
  Env$l.fr5$taille.wdg<-tkscale(Env$l.frames$Fr5,from=0,to=3,showvalue=TRUE,font=Env$police,variable=Env$l.var$lg.moustaches,resolution=0.1,orient="horizontal")
  tkbind(Env$l.fr5$taille.wdg,"<Enter>",function() {msg(text=Env$voc[77,1],type="info")})
  tkbind(Env$l.fr5$taille.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr5$trait.lab<-tklabel(Env$l.frames$Fr5,text=Env$voc[59,1],font=Env$police)
  Env$l.fr5$trait.wdg<-ttkcombobox(Env$l.frames$Fr5,font=Env$police,values=Env$voc[60:62,1],textvariable=Env$l.var$trait1,state="readonly",width=15)
  Env$l.fr5$espace.hor<-tklabel(Env$l.frames$Fr5,text="                                   ",font=Env$police)
  tkgrid(Env$l.fr5$taille.lab,row=0,column=0,sticky="e")
  tkgrid(Env$l.fr5$taille.wdg,row=0,column=1,sticky="w")
  tkgrid(Env$l.fr5$espace.hor,row=0,column=2)
  tkgrid(Env$l.fr5$trait.lab,row=0,column=3,sticky="e")
  tkgrid(Env$l.fr5$trait.wdg,row=0,column=4,sticky="w")
}

fr5.openB<-function() {
  Env$l.frames$Fr5.status<-1
  tkconfigure(Env$l.wdg$but.lab5,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr5)) {tkdestroy(Env$l.fr5[[i]])}
  Env$l.fr5<-list()
  Env$l.fr5$type.lab<-tklabel(Env$l.frames$Fr5,text=Env$voc[94,1],font=Env$police)
  Env$l.fr5$type.wdg<-ttkcombobox(Env$l.frames$Fr5,values="",textvariable=Env$l.var$erreur,font=Env$police,state="readonly")
  if (tclvalue(Env$l.var$moyprop)=="moy") {
    tkconfigure(Env$l.fr5$type.wdg,values=Env$voc[c(95:98),1])
  } else {
    tkconfigure(Env$l.fr5$type.wdg,values=Env$voc[c(95,97,98),1])
  }
  Env$l.fr5$col.lab<-tklabel(Env$l.frames$Fr5,text=Env$voc[45,1],font=Env$police)
  Env$l.fr5$col.wdg<-tkcanvas(Env$l.frames$Fr5,width="25",height="20",bg=tclvalue(Env$l.var$couleur2A))
  tkbind(Env$l.fr5$col.wdg,"<ButtonRelease-1>",function() {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$couleur2A),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$couleur2A)<-temp
	tkconfigure(Env$l.fr5$col.wdg,bg=temp)
    }
  })
  active.erreur()
  Env$l.fr5$espace.hor<-tklabel(Env$l.frames$Fr5,text="                                                  ",font=Env$police)
  tkgrid(Env$l.fr5$type.lab,row=0,column=0,sticky="e")
  tkgrid(Env$l.fr5$type.wdg,row=0,column=1,sticky="w")
  tkgrid(Env$l.fr5$espace.hor,row=0,column=2)
  tkgrid(Env$l.fr5$col.lab,row=0,column=3)
  tkgrid(Env$l.fr5$col.wdg,row=0,column=4)
}

fr5.openCa<-function() {
  Env$l.frames$Fr5.status<-1
  tkconfigure(Env$l.wdg$but.lab5,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr5)) {tkdestroy(Env$l.fr5[[i]])}
  Env$l.fr5<-list()
  Env$l.fr5$legende.lab<-tklabel(Env$l.frames$Fr5,text=Env$voc[99,1],font=Env$police)
  Env$l.fr5$legende.wdg<-tkcheckbutton(Env$l.frames$Fr5,variable=Env$l.var$legende)
  Env$l.fr5$position.lab<-tklabel(Env$l.frames$Fr5,text=Env$voc[100,1],font=Env$police)
  Env$l.fr5$position.wdg<-ttkcombobox(Env$l.frames$Fr5,font=Env$police,values=Env$voc[101:109,1],textvariable=Env$l.var$legende.pos,state="readonly",width=15)
  Env$l.fr5$titre.lab<-tklabel(Env$l.frames$Fr5,text=Env$voc[44,1],font=Env$police)
  Env$l.fr5$titre.wdg<-tkentry(Env$l.frames$Fr5,width=40,font=Env$police,textvariable=Env$l.var$legende.titre)
  active.legende2()
  Env$l.fr5$espace.hor1<-tklabel(Env$l.frames$Fr5,text="               ",font=Env$police)
  Env$l.fr5$espace.hor2<-tklabel(Env$l.frames$Fr5,text="               ",font=Env$police)
  tkgrid(Env$l.fr5$legende.lab,row=0,column=0,sticky="e")
  tkgrid(Env$l.fr5$legende.wdg,row=0,column=1,sticky="w")
  tkgrid(Env$l.fr5$espace.hor1,row=0,column=2)
  tkgrid(Env$l.fr5$titre.lab,row=0,column=3,sticky="e")
  tkgrid(Env$l.fr5$titre.wdg,row=0,column=4,sticky="w")
  tkgrid(Env$l.fr5$espace.hor2,row=0,column=5)
  tkgrid(Env$l.fr5$position.lab,row=0,column=6,sticky="e")
  tkgrid(Env$l.fr5$position.wdg,row=0,column=7,sticky="w")
}

fr5.openCo<-function() {
  Env$l.frames$Fr5.status<-1
  tkconfigure(Env$l.wdg$but.lab5,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr5)) {tkdestroy(Env$l.fr5[[i]])}
  Env$l.fr5<-list()
  Env$l.fr5$type.lab<-tklabel(Env$l.frames$Fr5,text=Env$voc[94,1],font=Env$police)
  Env$l.fr5$type.wdg<-ttkcombobox(Env$l.frames$Fr5,values="",textvariable=Env$l.var$erreur,font=Env$police,state="readonly")
  if (tclvalue(Env$l.var$moyprop)=="moy") {
    tkconfigure(Env$l.fr5$type.wdg,values=Env$voc[c(95:98),1])
  } else {
    tkconfigure(Env$l.fr5$type.wdg,values=Env$voc[c(95,97,98),1])
  }
  tkgrid(Env$l.fr5$type.lab,row=0,column=0,sticky="e")
  tkgrid(Env$l.fr5$type.wdg,row=0,column=1,sticky="w")
}

fr5.openN<-function() {
  Env$l.frames$Fr5.status<-1
  tkconfigure(Env$l.wdg$but.lab5,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr5)) {tkdestroy(Env$l.fr5[[i]])}
  Env$l.fr5<-list()
  Env$l.fr5$noms.list<-tklistbox(Env$l.frames$Fr5,height=7,font=Env$police,selectmode="single",yscrollcommand=function(...) tkset(Env$l.fr5$noms.scroll,...))
  Env$l.fr5$noms.scroll<-tkscrollbar(Env$l.frames$Fr5,repeatinterval=5,command=function(...) tkyview(Env$l.fr5$noms.list,...))
  tkbind(Env$l.fr5$noms.list,"<Enter>",function() {if (tclvalue(Env$l.var$plusieurs)==1) {msg(text=Env$voc[141,1],type="info")}})
  tkbind(Env$l.fr5$noms.list,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.fr5$noms.list,"<ButtonRelease-1>",function() {
    if (tclvalue(Env$l.var$plusieurs)==1) {
	Env$l.var$select<-as.numeric(tclvalue(tkcurselection(Env$l.fr5$noms.list)))+1
	tclvalue(Env$l.var$droiteA)<-as.character(Env$l.var$droiteB[Env$l.var$select])
	tclvalue(Env$l.var$intervalA)<-as.character(Env$l.var$intervalB[Env$l.var$select])
	tclvalue(Env$l.var$trait1)<-as.character(Env$l.var$trait2[Env$l.var$select])
	tclvalue(Env$l.var$epaisseur1)<-as.character(Env$l.var$epaisseur2[Env$l.var$select])
    }
  })
  Env$l.fr5$droite.lab<-tklabel(Env$l.frames$Fr5,text=Env$voc[144,1],font=Env$police)
  Env$l.fr5$droite.wdg<-ttkcombobox(Env$l.frames$Fr5,width=43,values=Env$voc[c(95,145:148),1],textvariable=Env$l.var$droiteA,state="readonly",font=Env$police)
  tkbind(Env$l.fr5$droite.wdg,"<<ComboboxSelected>>",function() {
    if (tclvalue(Env$l.var$plusieurs)==1) {
	Env$l.var$droiteB[Env$l.var$select]<-tclvalue(Env$l.var$droiteA)
	tkselection.set(Env$l.fr5$noms.list,as.character(Env$l.var$select-1))
    }
  })
  Env$l.fr5$interval.lab<-tklabel(Env$l.frames$Fr5,text=Env$voc[259,1],font=Env$police)
  Env$l.fr5$interval.wdg<-ttkcombobox(Env$l.frames$Fr5,width=43,values=Env$voc[260:263,1],textvariable=Env$l.var$intervalA,state="readonly",font=Env$police)
  tkbind(Env$l.fr5$interval.wdg,"<<ComboboxSelected>>",function() {
    if (tclvalue(Env$l.var$plusieurs)==1) {
	Env$l.var$intervalB[Env$l.var$select]<-tclvalue(Env$l.var$intervalA)
	tkselection.set(Env$l.fr5$noms.list,as.character(Env$l.var$select-1))
    }
  })
  Env$l.fr5$type.trait.lab<-tklabel(Env$l.frames$Fr5,text=Env$voc[59,1],font=Env$police)
  Env$l.fr5$type.trait.wdg<-ttkcombobox(Env$l.frames$Fr5,values=Env$voc[60:62,1],textvariable=Env$l.var$trait1,state="readonly",font=Env$police)
  tkbind(Env$l.fr5$type.trait.wdg,"<<ComboboxSelected>>",function() {
    if (tclvalue(Env$l.var$plusieurs)==1) {
	Env$l.var$trait2[Env$l.var$select]<-tclvalue(Env$l.var$trait1)
	tkselection.set(Env$l.fr5$noms.list,as.character(Env$l.var$select-1))
    }
  })
  Env$l.fr5$epaisseur.lab<-tklabel(Env$l.frames$Fr5,text=Env$voc[63,1],font=Env$police)
  Env$l.fr5$epaisseur.wdg<-tkscale(Env$l.frames$Fr5,showvalue=TRUE,from=1,to=4,resolution=1,font=Env$police,variable=Env$l.var$epaisseur1,orient="horizontal",command=function(...) {
    if (tclvalue(Env$l.var$plusieurs)==1) {
	Env$l.var$epaisseur2[Env$l.var$select]<-as.numeric(tclvalue(Env$l.var$epaisseur1))
	tkselection.set(Env$l.fr5$noms.list,as.character(Env$l.var$select-1))
    }
  })
  if (tclvalue(Env$l.var$plusieurs)==1) {
    Env$l.var$select<-1
    for (i in 1:length(Env$l.var$noms1)) {tkinsert(Env$l.fr5$noms.list,"end",Env$l.var$noms1[i])}
    tkselection.set(Env$l.fr5$noms.list,"0")
    tclvalue(Env$l.var$droiteA)<-as.character(Env$l.var$droiteB[1])
    tclvalue(Env$l.var$trait1)<-as.character(Env$l.var$trait2[1])
    tclvalue(Env$l.var$epaisseur1)<-as.character(Env$l.var$epaisseur2[1])
  } else {
    tkconfigure(Env$l.fr5$noms.list,state="disabled")
  }
  Env$l.fr5$espace.hor<-tklabel(Env$l.frames$Fr5,text="                    ",font=Env$police)
  tkgrid(Env$l.fr5$noms.list,Env$l.fr5$noms.scroll,row=0,column=0,rowspan=4,sticky="w");tkgrid.configure(Env$l.fr5$noms.scroll,sticky="ens")
  tkgrid(Env$l.fr5$espace.hor,row=0,column=1)
  tkgrid(Env$l.fr5$droite.lab,row=0,column=2,sticky="e")
  tkgrid(Env$l.fr5$droite.wdg,row=0,column=3,sticky="w")
  tkgrid(Env$l.fr5$interval.lab,row=1,column=2,sticky="e")
  tkgrid(Env$l.fr5$interval.wdg,row=1,column=3,sticky="w")
  tkgrid(Env$l.fr5$type.trait.lab,row=2,column=2,sticky="e")
  tkgrid(Env$l.fr5$type.trait.wdg,row=2,column=3,sticky="w")
  tkgrid(Env$l.fr5$epaisseur.lab,row=3,column=2,sticky="e")
  tkgrid(Env$l.fr5$epaisseur.wdg,row=3,column=3,sticky="w")
}


#-------------------------------------------------
# Frame 6
#-------------------------------------------------

fr6.close<-function() {
  Env$l.frames$Fr6.status<-0
  tkconfigure(Env$l.wdg$but.lab6,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_bas.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr6)) {tkdestroy(Env$l.fr6[[i]])}
  Env$l.fr6<-list()
  Env$l.fr6$vide<-tklabel(Env$l.frames$Fr6,text="",font=Env$police2)
  tkgrid(Env$l.fr6$vide)
}

fr6.openM<-function() {
  Env$l.frames$Fr6.status<-1
  tkconfigure(Env$l.wdg$but.lab6,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr6)) {tkdestroy(Env$l.fr6[[i]])}
  Env$l.fr6<-list()
  Env$l.fr6$legende.lab<-tklabel(Env$l.frames$Fr6,text=Env$voc[99,1],font=Env$police)
  Env$l.fr6$legende.wdg<-tkcheckbutton(Env$l.frames$Fr6,variable=Env$l.var$legende)
  Env$l.fr6$titre.lab<-tklabel(Env$l.frames$Fr6,text=Env$voc[44,1],font=Env$police)
  Env$l.fr6$titre.wdg<-tkentry(Env$l.frames$Fr6,textvariable=Env$l.var$legende.titre,font=Env$police,width=40)
  Env$l.fr6$position.lab<-tklabel(Env$l.frames$Fr6,text=Env$voc[100,1],font=Env$police)
  Env$l.fr6$position.wdg<-ttkcombobox(Env$l.frames$Fr6,font=Env$police,values=Env$voc[101:109,1],textvariable=Env$l.var$legende.pos,state="readonly",width=15)
  Env$l.fr6$noms.lab1<-tklabel(Env$l.frames$Fr6,text=Env$voc[110,1],font=Env$police)
  Env$l.fr6$noms.list<-tklistbox(Env$l.frames$Fr6,height=4,font=Env$police,selectmode="single",yscrollcommand=function(...) tkset(Env$l.fr6$noms.scroll,...))
  Env$l.fr6$noms.scroll<-tkscrollbar(Env$l.frames$Fr6,repeatinterval=4,command=function(...) tkyview(Env$l.fr6$noms.list,...))
  tkbind(Env$l.fr6$noms.list,"<ButtonRelease-1>",function() {
    tkdelete(Env$l.fr6$noms.wdg,0,"end")
	tkinsert(Env$l.fr6$noms.wdg,"end",Env$l.var$noms2[as.numeric(tclvalue(tkcurselection(Env$l.fr6$noms.list)))+1])
  })
  Env$l.fr6$noms.lab2<-tklabel(Env$l.frames$Fr6,text=Env$voc[111,1],font=Env$police)
  Env$l.fr6$noms.wdg<-tkentry(Env$l.frames$Fr6,width=20,font=Env$police)
  tkbind(Env$l.fr6$noms.wdg,"<Enter>",function() {msg(text=Env$voc[26,1],type="info")})
  tkbind(Env$l.fr6$noms.wdg,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.fr6$noms.wdg,"<ButtonRelease-1>",function() {tkdelete(Env$l.fr6$noms.wdg,0,"end")})
  tkbind(Env$l.fr6$noms.wdg,"<Return>",function() {
    rename.legende3(value.list=tclvalue(tkcurselection(Env$l.fr6$noms.list)),value.nom=tclvalue(tkget(Env$l.fr6$noms.wdg)))
  })
  active.legende()
  for (i in 1:length(Env$l.var$noms2)) tkinsert(Env$l.fr6$noms.list,"end",Env$l.var$noms2[i])
  Env$l.fr6$espace.hor<-tklabel(Env$l.frames$Fr6,text="                    ",font=Env$police)
  tkgrid(Env$l.fr6$legende.lab,row=0,column=0,sticky="e")
  tkgrid(Env$l.fr6$legende.wdg,row=0,column=1,sticky="w")
  tkgrid(Env$l.fr6$titre.lab,row=1,column=0,sticky="e")
  tkgrid(Env$l.fr6$titre.wdg,row=1,column=1,sticky="w")
  tkgrid(Env$l.fr6$position.lab,row=2,column=0,sticky="e")
  tkgrid(Env$l.fr6$position.wdg,row=2,column=1,sticky="w")
  tkgrid(Env$l.fr6$espace.hor,row=0,column=2)
  tkgrid(Env$l.fr6$noms.lab1,row=0,column=3,sticky="e")
  tkgrid(Env$l.fr6$noms.list,Env$l.fr6$noms.scroll,row=0,column=4,rowspan=4,sticky="w");tkgrid.configure(Env$l.fr6$noms.scroll,sticky="ens")
  tkgrid(Env$l.fr6$noms.lab2,row=4,column=3,sticky="e")
  tkgrid(Env$l.fr6$noms.wdg,row=4,column=4,sticky="w")
}

fr6.openB<-function() {
  Env$l.frames$Fr6.status<-1
  tkconfigure(Env$l.wdg$but.lab6,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr6)) {tkdestroy(Env$l.fr6[[i]])}
  Env$l.fr6<-list()
  Env$l.fr6$legende.lab<-tklabel(Env$l.frames$Fr6,text=Env$voc[99,1],font=Env$police)
  Env$l.fr6$legende.wdg<-tkcheckbutton(Env$l.frames$Fr6,variable=Env$l.var$legende)
  Env$l.fr6$titre.lab<-tklabel(Env$l.frames$Fr6,text=Env$voc[44,1],font=Env$police)
  Env$l.fr6$titre.wdg<-tkentry(Env$l.frames$Fr6,textvariable=Env$l.var$legende.titre,font=Env$police,width=40)
  Env$l.fr6$position.lab<-tklabel(Env$l.frames$Fr6,text=Env$voc[100,1],font=Env$police)
  Env$l.fr6$position.wdg<-ttkcombobox(Env$l.frames$Fr6,font=Env$police,values=Env$voc[101:109,1],textvariable=Env$l.var$legende.pos,state="readonly",width=15)
  Env$l.fr6$noms.lab1<-tklabel(Env$l.frames$Fr6,text=Env$voc[110,1],font=Env$police)
  Env$l.fr6$noms.list<-tklistbox(Env$l.frames$Fr6,height=4,font=Env$police,selectmode="single",yscrollcommand=function(...) tkset(Env$l.fr6$noms.scroll,...))
  Env$l.fr6$noms.scroll<-tkscrollbar(Env$l.frames$Fr6,repeatinterval=4,command=function(...) tkyview(Env$l.fr6$noms.list,...))
  tkbind(Env$l.fr6$noms.list,"<ButtonRelease-1>",function() {
    tkdelete(Env$l.fr6$noms.wdg,0,"end")
    if (tclvalue(Env$l.var$moyprop)=="moy") {
	tkinsert(Env$l.fr6$noms.wdg,"end",Env$l.var$noms2[as.numeric(tclvalue(tkcurselection(Env$l.fr6$noms.list)))+1])
    } else {
	tkinsert(Env$l.fr6$noms.wdg,"end",Env$l.var$nomsprop[as.numeric(tclvalue(tkcurselection(Env$l.fr6$noms.list)))+1])
    }
  })
  Env$l.fr6$noms.lab2<-tklabel(Env$l.frames$Fr6,text=Env$voc[111,1],font=Env$police)
  Env$l.fr6$noms.wdg<-tkentry(Env$l.frames$Fr6,width=20,font=Env$police)
  tkbind(Env$l.fr6$noms.wdg,"<Enter>",function() {msg(text=Env$voc[26,1],type="info")})
  tkbind(Env$l.fr6$noms.wdg,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.fr6$noms.wdg,"<ButtonRelease-1>",function() {tkdelete(Env$l.fr6$noms.wdg,0,"end")})
  tkbind(Env$l.fr6$noms.wdg,"<Return>",function() {
    rename.legende(value.list=tclvalue(tkcurselection(Env$l.fr6$noms.list)),value.nom=tclvalue(tkget(Env$l.fr6$noms.wdg)))
  })
  active.legende()
  if (tclvalue(Env$l.var$moyprop)=="moy") {
    for (i in 1:length(Env$l.var$noms2)) tkinsert(Env$l.fr6$noms.list,"end",Env$l.var$noms2[i])
  } else {
    for (i in 1:length(Env$l.var$nomsprop)) tkinsert(Env$l.fr6$noms.list,"end",Env$l.var$nomsprop[i])
  }
  Env$l.fr6$espace.hor<-tklabel(Env$l.frames$Fr6,text="                    ",font=Env$police)
  tkgrid(Env$l.fr6$legende.lab,row=0,column=0,sticky="e")
  tkgrid(Env$l.fr6$legende.wdg,row=0,column=1,sticky="w")
  tkgrid(Env$l.fr6$titre.lab,row=1,column=0,sticky="e")
  tkgrid(Env$l.fr6$titre.wdg,row=1,column=1,sticky="w")
  tkgrid(Env$l.fr6$position.lab,row=2,column=0,sticky="e")
  tkgrid(Env$l.fr6$position.wdg,row=2,column=1,sticky="w")
  tkgrid(Env$l.fr6$espace.hor,row=0,column=2)
  tkgrid(Env$l.fr6$noms.lab1,row=0,column=3,sticky="e")
  tkgrid(Env$l.fr6$noms.list,Env$l.fr6$noms.scroll,row=0,column=4,rowspan=4,sticky="w");tkgrid.configure(Env$l.fr6$noms.scroll,sticky="ens")
  tkgrid(Env$l.fr6$noms.lab2,row=4,column=3,sticky="e")
  tkgrid(Env$l.fr6$noms.wdg,row=4,column=4,sticky="w")
}

fr6.openCo<-function() {
  Env$l.frames$Fr6.status<-1
  tkconfigure(Env$l.wdg$but.lab6,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr6)) {tkdestroy(Env$l.fr6[[i]])}
  Env$l.fr6<-list()
  Env$l.fr6$legende.lab<-tklabel(Env$l.frames$Fr6,text=Env$voc[99,1],font=Env$police)
  Env$l.fr6$legende.wdg<-tkcheckbutton(Env$l.frames$Fr6,variable=Env$l.var$legende)
  Env$l.fr6$titre.lab<-tklabel(Env$l.frames$Fr6,text=Env$voc[44,1],font=Env$police)
  Env$l.fr6$titre.wdg<-tkentry(Env$l.frames$Fr6,textvariable=Env$l.var$legende.titre,font=Env$police,width=40)
  Env$l.fr6$position.lab<-tklabel(Env$l.frames$Fr6,text=Env$voc[100,1],font=Env$police)
  Env$l.fr6$position.wdg<-ttkcombobox(Env$l.frames$Fr6,font=Env$police,values=Env$voc[101:109,1],textvariable=Env$l.var$legende.pos,state="readonly",width=15)
  Env$l.fr6$noms.lab1<-tklabel(Env$l.frames$Fr6,text=Env$voc[110,1],font=Env$police)
  Env$l.fr6$noms.list<-tklistbox(Env$l.frames$Fr6,height=4,font=Env$police,selectmode="single",yscrollcommand=function(...) tkset(Env$l.fr6$noms.scroll,...))
  Env$l.fr6$noms.scroll<-tkscrollbar(Env$l.frames$Fr6,repeatinterval=4,command=function(...) tkyview(Env$l.fr6$noms.list,...))
  if (tclvalue(Env$l.var$plusieurs)==1) {
	for (i in 1:length(Env$l.var$noms1)) {tkinsert(Env$l.fr6$noms.list,"end",Env$l.var$noms1[i])}
  }
  tkbind(Env$l.fr6$noms.list,"<ButtonRelease-1>",function() {
    tkdelete(Env$l.fr6$noms.wdg,0,"end")
    tkinsert(Env$l.fr6$noms.wdg,"end",Env$l.var$noms1[as.numeric(tclvalue(tkcurselection(Env$l.fr6$noms.list)))+1])
  })
  Env$l.fr6$noms.lab2<-tklabel(Env$l.frames$Fr6,text=Env$voc[111,1],font=Env$police)
  Env$l.fr6$noms.wdg<-tkentry(Env$l.frames$Fr6,width=20,font=Env$police)
  tkbind(Env$l.fr6$noms.wdg,"<Enter>",function() {msg(text=Env$voc[26,1],type="info")})
  tkbind(Env$l.fr6$noms.wdg,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.fr6$noms.wdg,"<ButtonRelease-1>",function() {tkdelete(Env$l.fr6$noms.wdg,0,"end")})
  tkbind(Env$l.fr6$noms.wdg,"<Return>",function() {
    rename.legende2(value.list=tclvalue(tkcurselection(Env$l.fr6$noms.list)),value.nom=tclvalue(tkget(Env$l.fr6$noms.wdg)))
  })
  active.legende()
  Env$l.fr6$espace.hor<-tklabel(Env$l.frames$Fr6,text="                    ",font=Env$police)
  tkgrid(Env$l.fr6$legende.lab,row=0,column=0,sticky="e")
  tkgrid(Env$l.fr6$legende.wdg,row=0,column=1,sticky="w")
  tkgrid(Env$l.fr6$titre.lab,row=1,column=0,sticky="e")
  tkgrid(Env$l.fr6$titre.wdg,row=1,column=1,sticky="w")
  tkgrid(Env$l.fr6$position.lab,row=2,column=0,sticky="e")
  tkgrid(Env$l.fr6$position.wdg,row=2,column=1,sticky="w")
  tkgrid(Env$l.fr6$espace.hor,row=0,column=2)
  tkgrid(Env$l.fr6$noms.lab1,row=0,column=3,sticky="e")
  tkgrid(Env$l.fr6$noms.list,Env$l.fr6$noms.scroll,row=0,column=4,rowspan=4,sticky="w");tkgrid.configure(Env$l.fr6$noms.scroll,sticky="ens")
  tkgrid(Env$l.fr6$noms.lab2,row=4,column=3,sticky="e")
  tkgrid(Env$l.fr6$noms.wdg,row=4,column=4,sticky="w")
}

fr6.openN<-function() {
  Env$l.frames$Fr6.status<-1
  tkconfigure(Env$l.wdg$but.lab6,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)))
  for (i in 1:length(Env$l.fr6)) {tkdestroy(Env$l.fr6[[i]])}
  Env$l.fr6<-list()
  Env$l.fr6$legende.lab<-tklabel(Env$l.frames$Fr6,text=Env$voc[99,1],font=Env$police)
  Env$l.fr6$legende.wdg<-tkcheckbutton(Env$l.frames$Fr6,variable=Env$l.var$legende)
  Env$l.fr6$titre.lab<-tklabel(Env$l.frames$Fr6,text=Env$voc[44,1],font=Env$police)
  Env$l.fr6$titre.wdg<-tkentry(Env$l.frames$Fr6,textvariable=Env$l.var$legende.titre,font=Env$police,width=40)
  Env$l.fr6$position.lab<-tklabel(Env$l.frames$Fr6,text=Env$voc[100,1],font=Env$police)
  Env$l.fr6$position.wdg<-ttkcombobox(Env$l.frames$Fr6,font=Env$police,values=Env$voc[101:109,1],textvariable=Env$l.var$legende.pos,state="readonly",width=15)
  Env$l.fr6$noms.lab1<-tklabel(Env$l.frames$Fr6,text=Env$voc[110,1],font=Env$police)
  Env$l.fr6$noms.list<-tklistbox(Env$l.frames$Fr6,height=4,font=Env$police,selectmode="single",yscrollcommand=function(...) tkset(Env$l.fr6$noms.scroll,...))
  Env$l.fr6$noms.scroll<-tkscrollbar(Env$l.frames$Fr6,repeatinterval=4,command=function(...) tkyview(Env$l.fr6$noms.list,...))
  if (tclvalue(Env$l.var$plusieurs)==1) {
	for (i in 1:length(Env$l.var$noms1)) {tkinsert(Env$l.fr6$noms.list,"end",Env$l.var$noms1[i])}
  }
  tkbind(Env$l.fr6$noms.list,"<ButtonRelease-1>",function() {
    tkdelete(Env$l.fr6$noms.wdg,0,"end")
    tkinsert(Env$l.fr6$noms.wdg,"end",Env$l.var$noms1[as.numeric(tclvalue(tkcurselection(Env$l.fr6$noms.list)))+1])
  })
  Env$l.fr6$noms.lab2<-tklabel(Env$l.frames$Fr6,text=Env$voc[111,1],font=Env$police)
  Env$l.fr6$noms.wdg<-tkentry(Env$l.frames$Fr6,width=20,font=Env$police)
  tkbind(Env$l.fr6$noms.wdg,"<Enter>",function() {msg(text=Env$voc[26,1],type="info")})
  tkbind(Env$l.fr6$noms.wdg,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.fr6$noms.wdg,"<ButtonRelease-1>",function() {tkdelete(Env$l.fr6$noms.wdg,0,"end")})
  tkbind(Env$l.fr6$noms.wdg,"<Return>",function() {
    rename.legende2(value.list=tclvalue(tkcurselection(Env$l.fr6$noms.list)),value.nom=tclvalue(tkget(Env$l.fr6$noms.wdg)))
  })
  active.legende()
  Env$l.fr6$espace.hor<-tklabel(Env$l.frames$Fr6,text="                    ",font=Env$police)
  tkgrid(Env$l.fr6$legende.lab,row=0,column=0,sticky="e")
  tkgrid(Env$l.fr6$legende.wdg,row=0,column=1,sticky="w")
  tkgrid(Env$l.fr6$titre.lab,row=1,column=0,sticky="e")
  tkgrid(Env$l.fr6$titre.wdg,row=1,column=1,sticky="w")
  tkgrid(Env$l.fr6$position.lab,row=2,column=0,sticky="e")
  tkgrid(Env$l.fr6$position.wdg,row=2,column=1,sticky="w")
  tkgrid(Env$l.fr6$espace.hor,row=0,column=2)
  tkgrid(Env$l.fr6$noms.lab1,row=0,column=3,sticky="e")
  tkgrid(Env$l.fr6$noms.list,Env$l.fr6$noms.scroll,row=0,column=4,rowspan=4,sticky="w");tkgrid.configure(Env$l.fr6$noms.scroll,sticky="ens")
  tkgrid(Env$l.fr6$noms.lab2,row=4,column=3,sticky="e")
  tkgrid(Env$l.fr6$noms.wdg,row=4,column=4,sticky="w")
}


#-------------------------------------------------
# Remise à zéro des variables de la Frame 7
#-------------------------------------------------

reinit.fr7<-function() {
  tclvalue(Env$l.var$nw.col)<-"white"
  tclvalue(Env$l.var$nw.lignes)<-"1"
  tclvalue(Env$l.var$nw.colonnes)<-"1"
  tclvalue(Env$l.var$add.param1)<-""
  tclvalue(Env$l.var$add.trait)<-""
  tclvalue(Env$l.var$add.epaisseur1)<-"1"
  tclvalue(Env$l.var$add.col1)<-"black"
  tclvalue(Env$l.var$add.param2)<-""
  tclvalue(Env$l.var$add.param3)<-""
  tclvalue(Env$l.var$add.distrib)<-"norm"
  tclvalue(Env$l.var$add.txt)<-""
  tclvalue(Env$l.var$add.epaisseur2)<-"1"
  tclvalue(Env$l.var$add.col2)<-"black"
  tclvalue(Env$l.var$fen.num)<-""
  tclvalue(Env$l.var$fen.larg)<-"1250"
  tclvalue(Env$l.var$fen.type)<-"jpg"
  tclvalue(Env$l.var$fen.res)<-"150"
}


#-------------------------------------------------
# Fermer la Frame 7
#-------------------------------------------------

fr7.close<-function() {
  for (i in 1:length(Env$l.fr7)) {tkdestroy(Env$l.fr7[[i]])}
  Env$l.fr7<-list()
  Env$l.fr7$vide<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  tkgrid(Env$l.fr7$vide)
  tkconfigure(Env$l.frames$Fr7,borderwidth=0)
  reinit.fr7()
}


#-------------------------------------------------
# Ouvrir une nouvelle fenêtre
#-------------------------------------------------

new.window<-function() {
  fr7.close()
  for (i in 1:length(Env$l.fr7)) {tkdestroy(Env$l.fr7[[i]])}
  Env$l.fr7<-list()
  tkconfigure(Env$l.frames$Fr7,borderwidth=3)
  Env$l.fr7$titre<-tklabel(Env$l.frames$Fr7,text=Env$voc[149,1],font=Env$police3)
  Env$l.fr7$col.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[150,1],font=Env$police)
  Env$l.fr7$col.wdg<-tkcanvas(Env$l.frames$Fr7,width="40",height="25",bg=tclvalue(Env$l.var$nw.col))
  tkbind(Env$l.fr7$col.wdg,"<ButtonRelease-1>",function() {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$nw.col),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$nw.col)<-temp
	tkconfigure(Env$l.fr7$col.wdg,bg=tclvalue(Env$l.var$nw.col))
    }
  })
  Env$l.fr7$lignes.wdg<-tkscale(Env$l.frames$Fr7,from=1,to=4,showvalue=TRUE,font=Env$police5,variable=Env$l.var$nw.lignes,length=200,resolution=1,orient="vertical",command=function(...) {
    tkconfigure(Env$l.fr7$fenetre,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",paste("Fenetre",tclvalue(Env$l.var$nw.lignes),"-",tclvalue(Env$l.var$nw.colonnes),".gif",sep=""),fsep=.Platform$file.sep)))
  })
  Env$l.fr7$colonnes.wdg<-tkscale(Env$l.frames$Fr7,from=1,to=4,showvalue=TRUE,font=Env$police5,variable=Env$l.var$nw.colonnes,length=200,resolution=1,orient="horizontal",command=function(...) {
    tkconfigure(Env$l.fr7$fenetre,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",paste("Fenetre",tclvalue(Env$l.var$nw.lignes),"-",tclvalue(Env$l.var$nw.colonnes),".gif",sep=""),fsep=.Platform$file.sep)))
  })
  Env$l.fr7$fenetre<-tklabel(Env$l.frames$Fr7,height=200,width=200,font=Env$police,borderwidth=0,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fenetre1-1.gif",fsep=.Platform$file.sep)))
  Env$l.fr7$ok<-tkbutton(Env$l.frames$Fr7,width=16,text=Env$voc[151,1],font=Env$police,command=function() {
    dimensions<-if (tclvalue(Env$l.var$nw.lignes)=="1" & tclvalue(Env$l.var$nw.colonnes)=="1") {c(7,7)} else
      if (tclvalue(Env$l.var$nw.lignes)=="2" & tclvalue(Env$l.var$nw.colonnes)=="1") {c(6,12)} else
      if (tclvalue(Env$l.var$nw.lignes)=="3" & tclvalue(Env$l.var$nw.colonnes)=="1") {c(5,15)} else
      if (tclvalue(Env$l.var$nw.lignes)=="4" & tclvalue(Env$l.var$nw.colonnes)=="1") {c(4,16)} else
      if (tclvalue(Env$l.var$nw.lignes)=="1" & tclvalue(Env$l.var$nw.colonnes)=="2") {c(12,6)} else
      if (tclvalue(Env$l.var$nw.lignes)=="2" & tclvalue(Env$l.var$nw.colonnes)=="2") {c(10,10)} else
      if (tclvalue(Env$l.var$nw.lignes)=="3" & tclvalue(Env$l.var$nw.colonnes)=="2") {c(10,15)} else
      if (tclvalue(Env$l.var$nw.lignes)=="4" & tclvalue(Env$l.var$nw.colonnes)=="2") {c(8,16)} else
      if (tclvalue(Env$l.var$nw.lignes)=="1" & tclvalue(Env$l.var$nw.colonnes)=="3") {c(15,5)} else
      if (tclvalue(Env$l.var$nw.lignes)=="2" & tclvalue(Env$l.var$nw.colonnes)=="3") {c(15,10)} else
      if (tclvalue(Env$l.var$nw.lignes)=="3" & tclvalue(Env$l.var$nw.colonnes)=="3") {c(12,12)} else
      if (tclvalue(Env$l.var$nw.lignes)=="4" & tclvalue(Env$l.var$nw.colonnes)=="3") {c(12,16)} else
      if (tclvalue(Env$l.var$nw.lignes)=="1" & tclvalue(Env$l.var$nw.colonnes)=="4") {c(16,4)} else
      if (tclvalue(Env$l.var$nw.lignes)=="2" & tclvalue(Env$l.var$nw.colonnes)=="4") {c(16,8)} else
      if (tclvalue(Env$l.var$nw.lignes)=="3" & tclvalue(Env$l.var$nw.colonnes)=="4") {c(16,12)} else
      if (tclvalue(Env$l.var$nw.lignes)=="4" & tclvalue(Env$l.var$nw.colonnes)=="4") {c(16,16)}
    dev.new(width=dimensions[1],height=dimensions[2])
    par(bg=tclvalue(Env$l.var$nw.col),mfrow=c(as.numeric(tclvalue(Env$l.var$nw.lignes)),as.numeric(tclvalue(Env$l.var$nw.colonnes))),mar=c(5,6,4,2))
  })
  Env$l.fr7$fermer<-tkbutton(Env$l.frames$Fr7,width=16,text=Env$voc[152,1],font=Env$police,command=fr7.close)
  Env$l.fr7$espace.ver1<-tklabel(Env$l.frames$Fr7,text="",font=Env$police5)
  Env$l.fr7$espace.ver2<-tklabel(Env$l.frames$Fr7,text="",font=Env$police5)
  Env$l.fr7$espace.ver3<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver4<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver5<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.hor1<-tklabel(Env$l.frames$Fr7,text="     ",font=Env$police3)
  Env$l.fr7$espace.hor2<-tklabel(Env$l.frames$Fr7,text="     ",font=Env$police3)
  tkgrid(Env$l.fr7$espace.hor1,row=0,column=0)
  tkgrid(Env$l.fr7$titre,row=0,column=1,columnspan=2)
  tkgrid(Env$l.fr7$espace.hor2,row=0,column=3)
  tkgrid(Env$l.fr7$colonnes.wdg,column=2)
  tkgrid(Env$l.fr7$espace.ver1)
  tkgrid(Env$l.fr7$lignes.wdg,row=3,column=1)
  tkgrid(Env$l.fr7$fenetre,row=3,column=2)
  tkgrid(Env$l.fr7$espace.ver2)
  tkgrid(Env$l.fr7$col.lab,row=5,column=1,sticky="e")
  tkgrid(Env$l.fr7$col.wdg,row=5,column=2,sticky="w")
  tkgrid(Env$l.fr7$espace.ver3)
  tkgrid(Env$l.fr7$ok,column=1,columnspan=2)
  tkgrid(Env$l.fr7$espace.ver4)
  tkgrid(Env$l.fr7$fermer,column=1,columnspan=2)
  tkgrid(Env$l.fr7$espace.ver5)
}


#-------------------------------------------------
# Ajouter une droite horizontale sur le graphe
#-------------------------------------------------

horizontal<-function() {
  fr7.close()
  for (i in 1:length(Env$l.fr7)) {tkdestroy(Env$l.fr7[[i]])}
  Env$l.fr7<-list()
  tkconfigure(Env$l.frames$Fr7,borderwidth=3)
  Env$l.fr7$ord.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[164,1],font=Env$police)
  Env$l.fr7$ord.wdg<-tkentry(Env$l.frames$Fr7,width=8,textvariable=Env$l.var$add.param1,font=Env$police)
  tkbind(Env$l.fr7$ord.wdg,"<Enter>",function() {msg(text=Env$voc[155,1],type="warning")})
  tkbind(Env$l.fr7$ord.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr7$ord.but<-tkbutton(Env$l.frames$Fr7,text=Env$voc[162,1],width=20,font=Env$police,command=function() {
    if (dev.cur()>1) {
	ord.val<-round(locator(n=1)$y[1],4)
	tclvalue(Env$l.var$add.param1)<-ord.val
    } else {
	msg(text=Env$voc[163,1],type="error")
    }
  })
  Env$l.fr7$trait.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[59,1],font=Env$police)
  Env$l.fr7$trait.wdg<-ttkcombobox(Env$l.frames$Fr7,font=Env$police,values=c(Env$voc[60:62,1]),textvariable=Env$l.var$add.trait,state="readonly")
  Env$l.fr7$epaisseur.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[63,1],font=Env$police)
  Env$l.fr7$epaisseur.wdg<-tkscale(Env$l.frames$Fr7,from=1,to=5,showvalue=TRUE,font=Env$police,variable=Env$l.var$add.epaisseur1,resolution=1,orient="horizontal")
  Env$l.fr7$col.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[45,1],font=Env$police)
  Env$l.fr7$col.wdg<-tkcanvas(Env$l.frames$Fr7,width="40",height="25",bg=tclvalue(Env$l.var$add.col1))
  tkbind(Env$l.fr7$col.wdg,"<ButtonRelease-1>",function() {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$add.col1),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$add.col1)<-temp
	tkconfigure(Env$l.fr7$col.wdg,bg=tclvalue(Env$l.var$add.col1))
    }
  })
  Env$l.fr7$tracer<-tkbutton(Env$l.frames$Fr7,text=Env$voc[72,1],font=Env$police,font=Env$police,width=16,command=function() {
    if (dev.cur()>1) {
      abline(h=as.numeric(tclvalue(Env$l.var$add.param1)),lty=type.trait(type=tclvalue(Env$l.var$add.trait)),
	  lwd=as.numeric(tclvalue(Env$l.var$add.epaisseur1)),col=tclvalue(Env$l.var$add.col1))
	if (Env$l.code$save==TRUE) {
	  sink(file=file.path(Env$l.code$folder,paste(paste("GrapheR",paste(strsplit(as.character(Sys.Date()),split="-")[[1]],collapse="."),
	    sep="-"),".R",sep=""),fsep=.Platform$file.sep),append=TRUE)
	  cat("# Added: horizontal line\n\n")
	  texte<-paste("abline(h=",tclvalue(Env$l.var$add.param1),sep="")
	  if (tclvalue(Env$l.var$add.col1)!="black" & tclvalue(Env$l.var$add.col1)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$add.col1),"\"",sep="")}
	  texte<-paste(texte,", lty=",type.trait(type=tclvalue(Env$l.var$add.trait)),sep="")
	  texte<-paste(texte,", lwd=",tclvalue(Env$l.var$add.epaisseur1),")\n\n",sep="")
	  cat(texte)
	  sink(NULL)
	}
    } else {
	msg(text=Env$voc[163,1],type="error")
    }
  })
  Env$l.fr7$fermer<-tkbutton(Env$l.frames$Fr7,width=16,text=Env$voc[152,1],font=Env$police,command=fr7.close)
  Env$l.fr7$espace.ver1<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver2<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver3<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.hor1<-tklabel(Env$l.frames$Fr7,text="     ",font=Env$police3)
  Env$l.fr7$espace.hor2<-tklabel(Env$l.frames$Fr7,text="     ",font=Env$police3)
  Env$l.fr7$espace.hor3<-tklabel(Env$l.frames$Fr7,text="     ",font=Env$police3)
  tkgrid(Env$l.fr7$espace.ver1)
  tkgrid(Env$l.fr7$espace.hor1,row=1,column=0)
  tkgrid(Env$l.fr7$ord.lab,row=1,column=1,sticky="e")
  tkgrid(Env$l.fr7$ord.wdg,row=1,column=2)
  tkgrid(Env$l.fr7$espace.hor2,row=1,column=3)
  tkgrid(Env$l.fr7$ord.but,row=1,column=4)
  tkgrid(Env$l.fr7$espace.hor3,row=1,column=5)
  tkgrid(Env$l.fr7$trait.lab,row=2,column=1,sticky="e")
  tkgrid(Env$l.fr7$trait.wdg,row=2,column=2,columnspan=3,sticky="w")
  tkgrid(Env$l.fr7$epaisseur.lab,row=3,column=1,sticky="e")
  tkgrid(Env$l.fr7$epaisseur.wdg,row=3,column=2,columnspan=3,sticky="w")
  tkgrid(Env$l.fr7$col.lab,row=4,column=1,sticky="e")
  tkgrid(Env$l.fr7$col.wdg,row=4,column=2,columnspan=3,sticky="w")
  tkgrid(Env$l.fr7$espace.ver2)
  tkgrid(Env$l.fr7$tracer,row=6,column=1,columnspan=2,sticky="e")
  tkgrid(Env$l.fr7$fermer,row=6,column=3,columnspan=2)
  tkgrid(Env$l.fr7$espace.ver3)
}


#-------------------------------------------------
# Ajouter une droite verticale sur le graphe
#-------------------------------------------------

vertical<-function() {
  fr7.close()
  for (i in 1:length(Env$l.fr7)) {tkdestroy(Env$l.fr7[[i]])}
  Env$l.fr7<-list()
  tkconfigure(Env$l.frames$Fr7,borderwidth=3)
  Env$l.fr7$abs.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[161,1],font=Env$police)
  Env$l.fr7$abs.wdg<-tkentry(Env$l.frames$Fr7,width=8,textvariable=Env$l.var$add.param1,font=Env$police)
  tkbind(Env$l.fr7$abs.wdg,"<Enter>",function() {msg(text=Env$voc[155,1],type="warning")})
  tkbind(Env$l.fr7$abs.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr7$abs.but<-tkbutton(Env$l.frames$Fr7,text=Env$voc[162,1],width=20,font=Env$police,command=function() {
    if (dev.cur()>1) {
	abs.val<-round(locator(n=1)$x[1],4)
	tclvalue(Env$l.var$add.param1)<-abs.val
    } else {
	msg(text=Env$voc[163,1],type="error")
    }
  })
  Env$l.fr7$trait.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[59,1],font=Env$police)
  Env$l.fr7$trait.wdg<-ttkcombobox(Env$l.frames$Fr7,font=Env$police,values=c(Env$voc[60:62,1]),textvariable=Env$l.var$add.trait,state="readonly")
  Env$l.fr7$epaisseur.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[63,1],font=Env$police)
  Env$l.fr7$epaisseur.wdg<-tkscale(Env$l.frames$Fr7,from=1,to=5,showvalue=TRUE,font=Env$police,variable=Env$l.var$add.epaisseur1,resolution=1,orient="horizontal")
  Env$l.fr7$col.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[45,1],font=Env$police)
  Env$l.fr7$col.wdg<-tkcanvas(Env$l.frames$Fr7,width="40",height="25",bg=tclvalue(Env$l.var$add.col1))
  tkbind(Env$l.fr7$col.wdg,"<ButtonRelease-1>",function() {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$add.col1),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$add.col1)<-temp
	tkconfigure(Env$l.fr7$col.wdg,bg=tclvalue(Env$l.var$add.col1))
    }
  })
  Env$l.fr7$tracer<-tkbutton(Env$l.frames$Fr7,text=Env$voc[72,1],font=Env$police,font=Env$police,width=16,command=function() {
    if (dev.cur()>1) {
      abline(v=as.numeric(tclvalue(Env$l.var$add.param1)),lty=type.trait(type=tclvalue(Env$l.var$add.trait)),
	  lwd=as.numeric(tclvalue(Env$l.var$add.epaisseur1)),col=tclvalue(Env$l.var$add.col1))
	if (Env$l.code$save==TRUE) {
	  sink(file=file.path(Env$l.code$folder,paste(paste("GrapheR",paste(strsplit(as.character(Sys.Date()),split="-")[[1]],collapse="."),
	    sep="-"),".R",sep=""),fsep=.Platform$file.sep),append=TRUE)
	  cat("# Added: vertical line\n\n")
	  texte<-paste("abline(v=",tclvalue(Env$l.var$add.param1),sep="")
	  if (tclvalue(Env$l.var$add.col1)!="black" & tclvalue(Env$l.var$add.col1)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$add.col1),"\"",sep="")}
	  texte<-paste(texte,", lty=",type.trait(type=tclvalue(Env$l.var$add.trait)),sep="")
	  texte<-paste(texte,", lwd=",tclvalue(Env$l.var$add.epaisseur1),")\n\n",sep="")
	  cat(texte)
	  sink(NULL)
	}
    } else {
	msg(text=Env$voc[163,1],type="error")
    }
  })
  Env$l.fr7$fermer<-tkbutton(Env$l.frames$Fr7,width=16,text=Env$voc[152,1],font=Env$police,command=fr7.close)
  Env$l.fr7$espace.ver1<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver2<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver3<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.hor1<-tklabel(Env$l.frames$Fr7,text="     ",font=Env$police3)
  Env$l.fr7$espace.hor2<-tklabel(Env$l.frames$Fr7,text="     ",font=Env$police3)
  Env$l.fr7$espace.hor3<-tklabel(Env$l.frames$Fr7,text="     ",font=Env$police3)
  tkgrid(Env$l.fr7$espace.ver1)
  tkgrid(Env$l.fr7$espace.hor1,row=1,column=0)
  tkgrid(Env$l.fr7$abs.lab,row=1,column=1,sticky="e")
  tkgrid(Env$l.fr7$abs.wdg,row=1,column=2)
  tkgrid(Env$l.fr7$espace.hor2,row=1,column=3)
  tkgrid(Env$l.fr7$abs.but,row=1,column=4)
  tkgrid(Env$l.fr7$espace.hor3,row=1,column=5)
  tkgrid(Env$l.fr7$trait.lab,row=2,column=1,sticky="e")
  tkgrid(Env$l.fr7$trait.wdg,row=2,column=2,columnspan=3,sticky="w")
  tkgrid(Env$l.fr7$epaisseur.lab,row=3,column=1,sticky="e")
  tkgrid(Env$l.fr7$epaisseur.wdg,row=3,column=2,columnspan=3,sticky="w")
  tkgrid(Env$l.fr7$col.lab,row=4,column=1,sticky="e")
  tkgrid(Env$l.fr7$col.wdg,row=4,column=2,columnspan=3,sticky="w")
  tkgrid(Env$l.fr7$espace.ver2)
  tkgrid(Env$l.fr7$tracer,row=6,column=1,columnspan=2,sticky="e")
  tkgrid(Env$l.fr7$fermer,row=6,column=3,columnspan=2)
  tkgrid(Env$l.fr7$espace.ver3)
}


#--------------------------------------------------------------------------
# Ajouter une droite quelconque sur le graphe
#--------------------------------------------------------------------------

affine<-function() {
  fr7.close()
  for (i in 1:length(Env$l.fr7)) {tkdestroy(Env$l.fr7[[i]])}
  Env$l.fr7<-list()
  tkconfigure(Env$l.frames$Fr7,borderwidth=3)
  Env$l.fr7$equation.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[165,1],font=Env$police)
  Env$l.fr7$coeffa.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[166,1],font=Env$police)
  Env$l.fr7$coeffa.wdg<-tkentry(Env$l.frames$Fr7,width=8,textvariable=Env$l.var$add.param1,font=Env$police)
  tkbind(Env$l.fr7$coeffa.wdg,"<Enter>",function() {msg(text=Env$voc[155,1],type="warning")})
  tkbind(Env$l.fr7$coeffa.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr7$coeffb.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[167,1],font=Env$police)
  Env$l.fr7$coeffb.wdg<-tkentry(Env$l.frames$Fr7,width=8,textvariable=Env$l.var$add.param2,font=Env$police)
  tkbind(Env$l.fr7$coeffb.wdg,"<Enter>",function() {msg(text=Env$voc[155,1],type="warning")})
  tkbind(Env$l.fr7$coeffb.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr7$coeffs.but<-tkbutton(Env$l.frames$Fr7,text=Env$voc[168,1],width=30,font=Env$police,command=function() {
    if (dev.cur()>1) {
	coord<-locator(n=2)
	pt1x<-coord$x[1]
	pt1y<-coord$y[1]
	pt2x<-coord$x[2]
	pt2y<-coord$y[2]
	a<-round(((pt2y-pt1y)/(pt2x-pt1x)),4)
	b<-round(pt1y-a*pt1x,4)
	tclvalue(Env$l.var$add.param1)<-a
	tclvalue(Env$l.var$add.param2)<-b
    } else {
	msg(text=Env$voc[163,1],type="error")
    }
  })
  Env$l.fr7$trait.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[59,1],font=Env$police)
  Env$l.fr7$trait.wdg<-ttkcombobox(Env$l.frames$Fr7,font=Env$police,values=c(Env$voc[60:62,1]),textvariable=Env$l.var$add.trait,state="readonly")
  Env$l.fr7$epaisseur.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[63,1],font=Env$police)
  Env$l.fr7$epaisseur.wdg<-tkscale(Env$l.frames$Fr7,from=1,to=5,showvalue=TRUE,font=Env$police,variable=Env$l.var$add.epaisseur1,resolution=1,orient="horizontal")
  Env$l.fr7$col.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[45,1],font=Env$police)
  Env$l.fr7$col.wdg<-tkcanvas(Env$l.frames$Fr7,width="40",height="25",bg=tclvalue(Env$l.var$add.col1))
  tkbind(Env$l.fr7$col.wdg,"<ButtonRelease-1>",function() {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$add.col1),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$add.col1)<-temp
	tkconfigure(Env$l.fr7$col.wdg,bg=tclvalue(Env$l.var$add.col1))
    }
  })
  Env$l.fr7$tracer<-tkbutton(Env$l.frames$Fr7,text=Env$voc[72,1],font=Env$police,font=Env$police,width=16,command=function() {
    if (dev.cur()>1) {
      abline(as.numeric(tclvalue(Env$l.var$add.param2)),as.numeric(tclvalue(Env$l.var$add.param1)),
	  lty=type.trait(type=tclvalue(Env$l.var$add.trait)),lwd=as.numeric(tclvalue(Env$l.var$add.epaisseur1)),
	  col=tclvalue(Env$l.var$add.col1),untf=ifelse(graphe.log()=="",FALSE,TRUE))
	if (Env$l.code$save==TRUE) {
	  sink(file=file.path(Env$l.code$folder,paste(paste("GrapheR",paste(strsplit(as.character(Sys.Date()),split="-")[[1]],collapse="."),
	    sep="-"),".R",sep=""),fsep=.Platform$file.sep),append=TRUE)
	  cat("# Added: straight line\n\n")
	  texte<-paste("abline(",tclvalue(Env$l.var$add.param2),", ",tclvalue(Env$l.var$add.param1),sep="")
	  if (tclvalue(Env$l.var$add.col1)!="black" & tclvalue(Env$l.var$add.col1)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$add.col1),"\"",sep="")}
	  texte<-paste(texte,", lty=",type.trait(type=tclvalue(Env$l.var$add.trait)),sep="")
	  texte<-paste(texte,", lwd=",tclvalue(Env$l.var$add.epaisseur1),sep="")
	  texte<-paste(texte,", untf=",ifelse(graphe.log()=="","FALSE","TRUE"),")\n\n",sep="")
	  cat(texte)
	  sink(NULL)
	}
    } else {
	msg(text=Env$voc[163,1],type="error")
    }
  })
  Env$l.fr7$fermer<-tkbutton(Env$l.frames$Fr7,width=16,text=Env$voc[152,1],font=Env$police,command=fr7.close)
  Env$l.fr7$espace.ver1<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver2<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver3<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.hor1<-tklabel(Env$l.frames$Fr7,text="     ",font=Env$police3)
  Env$l.fr7$espace.hor2<-tklabel(Env$l.frames$Fr7,text="     ",font=Env$police3)
  tkgrid(Env$l.fr7$espace.ver1)
  tkgrid(Env$l.fr7$espace.hor1,row=1,column=0)
  tkgrid(Env$l.fr7$equation.lab,row=1,column=1,sticky="e")
  tkgrid(Env$l.fr7$coeffa.lab,row=1,column=2)
  tkgrid(Env$l.fr7$coeffa.wdg,row=1,column=3)
  tkgrid(Env$l.fr7$coeffb.lab,row=1,column=4)
  tkgrid(Env$l.fr7$coeffb.wdg,row=1,column=5)
  tkgrid(Env$l.fr7$espace.hor2,row=1,column=6)
  tkgrid(Env$l.fr7$coeffs.but,row=2,column=2,columnspan=4)
  tkgrid(Env$l.fr7$trait.lab,row=3,column=1,sticky="e")
  tkgrid(Env$l.fr7$trait.wdg,row=3,column=2,columnspan=4,sticky="w")
  tkgrid(Env$l.fr7$epaisseur.lab,row=4,column=1,sticky="e")
  tkgrid(Env$l.fr7$epaisseur.wdg,row=4,column=2,columnspan=4,sticky="w")
  tkgrid(Env$l.fr7$col.lab,row=5,column=1,sticky="e")
  tkgrid(Env$l.fr7$col.wdg,row=5,column=2,columnspan=4,sticky="w")
  tkgrid(Env$l.fr7$espace.ver2)
  tkgrid(Env$l.fr7$tracer,row=7,column=1,columnspan=3)
  tkgrid(Env$l.fr7$fermer,row=7,column=3,columnspan=3,sticky="e")
  tkgrid(Env$l.fr7$espace.ver3)
}


#--------------------------------------------------------------------------
# Ajouter une courbe de distribution théorique sur le graphe
#--------------------------------------------------------------------------

distrib<-function() {
  fr7.close()
  for (i in 1:length(Env$l.fr7)) {tkdestroy(Env$l.fr7[[i]])}
  Env$l.fr7<-list()
  tkconfigure(Env$l.frames$Fr7,borderwidth=3)
  type.distrib<-function(dist,par1,par2,par3,act2,act3) {
    tclvalue(Env$l.var$add.param1)<-""
    tclvalue(Env$l.var$add.param2)<-""
    tclvalue(Env$l.var$add.param3)<-""
    tkconfigure(Env$l.fr7$distrib.lab,text=dist)
    tkconfigure(Env$l.fr7$param1.lab,text=par1)
    tkconfigure(Env$l.fr7$param2.lab,text=par2)
    tkconfigure(Env$l.fr7$param3.lab,text=par3)
    tkconfigure(Env$l.fr7$param1.wdg,state="normal")
    tkconfigure(Env$l.fr7$param2.wdg,state=ifelse(act2==TRUE,"normal","disabled"))
    tkconfigure(Env$l.fr7$param3.wdg,state=ifelse(act3==TRUE,"normal","disabled"))
  }
  Env$l.fr7$titre.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[169,1],font=Env$police3)
  Env$l.fr7$rb1<-tkradiobutton(Env$l.frames$Fr7,font=Env$police,variable=Env$l.var$add.distrib,value="norm",text=Env$voc[170,1],
    command=function() {type.distrib(dist="N (\U03BC,\U03C3)",par1=paste(Env$voc[189,1],"\U03BC",Env$voc[200,1],sep=""),
    par2=paste(Env$voc[190,1],"\U03C3",Env$voc[200,1],sep=""),par3="",act2=TRUE,act3=FALSE)})
  Env$l.fr7$rb2<-tkradiobutton(Env$l.frames$Fr7,font=Env$police,variable=Env$l.var$add.distrib,value="binom",text=Env$voc[171,1],
    command=function() {type.distrib(dist="B (n,p)",par1=Env$voc[183,1],par2=Env$voc[184,1],par3="",act2=TRUE,act3=FALSE)})
  Env$l.fr7$rb3<-tkradiobutton(Env$l.frames$Fr7,font=Env$police,variable=Env$l.var$add.distrib,value="gamma",text=Env$voc[172,1],
    command=function() {type.distrib(dist="G (\U03B1,\U03B2)",par1=paste(Env$voc[191,1],"\U03B1",Env$voc[200,1],sep=""),
    par2=paste(Env$voc[192,1],"\U03B2",Env$voc[200,1],sep=""),par3="",act2=TRUE,act3=FALSE)})
  Env$l.fr7$rb4<-tkradiobutton(Env$l.frames$Fr7,font=Env$police,variable=Env$l.var$add.distrib,value="poiss",text=Env$voc[173,1],
    command=function() {type.distrib(dist="P (n.p)",par1=Env$voc[183,1],par2=Env$voc[185,1],par3="",act2=TRUE,act3=FALSE)})
  Env$l.fr7$rb5<-tkradiobutton(Env$l.frames$Fr7,font=Env$police,variable=Env$l.var$add.distrib,value="expo",text=Env$voc[174,1],
    command=function() {type.distrib(dist="exp (\U03BB)",par1=paste(Env$voc[193,1],"\U03BB",Env$voc[200,1],sep=""),
    par2="",par3="",act2=FALSE,act3=FALSE)})
  Env$l.fr7$rb6<-tkradiobutton(Env$l.frames$Fr7,font=Env$police,variable=Env$l.var$add.distrib,value="nbinom",text=Env$voc[175,1],
    command=function() {type.distrib(dist="BN (k,p)",par1=Env$voc[186,1],par2=Env$voc[184,1],par3="",act2=TRUE,act3=FALSE)})
  Env$l.fr7$rb7<-tkradiobutton(Env$l.frames$Fr7,font=Env$police,variable=Env$l.var$add.distrib,value="chi",
    text=if (Env$lang=="fr") {paste(Env$voc[176,1],"\U03C7\U00B2",sep="")} else
    if (Env$lang=="en") {paste("\U03C7\U00B2",Env$voc[176,1],sep="")} else
    if (Env$lang=="es") {paste("\U03C7\U00B2",Env$voc[176,1],sep="")} else
    if (Env$lang=="de") {paste("\U03C7\U00B2",Env$voc[176,1],sep="")},
    command=function() {type.distrib(dist="\U03C7\U00B2 (\U03BD)",par1=paste(Env$voc[194,1],"\U03BD",Env$voc[200,1],sep=""),
    par2="",par3="",act2=FALSE,act3=FALSE)})
  Env$l.fr7$rb8<-tkradiobutton(Env$l.frames$Fr7,font=Env$police,variable=Env$l.var$add.distrib,value="geom",text=Env$voc[177,1],
    command=function() {type.distrib(dist="G (p)",par1=Env$voc[184,1],par2="",par3="",act2=FALSE,act3=FALSE)})
  Env$l.fr7$rb9<-tkradiobutton(Env$l.frames$Fr7,font=Env$police,variable=Env$l.var$add.distrib,value="fish",text=Env$voc[178,1],
    command=function() {type.distrib(dist=paste("F (\U03BD","1,\U03BD","2)",sep=""),par1=paste(Env$voc[195,1],"\U03BD","1",Env$voc[200,1],sep=""),
    par2=paste(Env$voc[196,1],"\U03BD","1",Env$voc[200,1],sep=""),par3="",act2=TRUE,act3=FALSE)})
  Env$l.fr7$rb10<-tkradiobutton(Env$l.frames$Fr7,font=Env$police,variable=Env$l.var$add.distrib,value="hyper",text=Env$voc[179,1],
    command=function() {type.distrib(dist="H (N,n,p)",par1=Env$voc[187,1],par2=Env$voc[188,1],par3=Env$voc[184,1],
    act2=TRUE,act3=TRUE)})
  Env$l.fr7$rb11<-tkradiobutton(Env$l.frames$Fr7,font=Env$police,variable=Env$l.var$add.distrib,value="stud",text=Env$voc[180,1],
    command=function() {type.distrib(dist="t (\U03BD)",par1=paste(Env$voc[194,1],"\U03BD",Env$voc[200,1],sep=""),
    par2="",par3="",act2=FALSE,act3=FALSE)})
  Env$l.fr7$rb12<-tkradiobutton(Env$l.frames$Fr7,font=Env$police,variable=Env$l.var$add.distrib,value="mann",text=Env$voc[181,1],
    command=function() {type.distrib(dist="U (n1,n2)",par1=Env$voc[197,1],par2=Env$voc[198,1],par3="",act2=TRUE,act3=FALSE)})
  Env$l.fr7$rb13<-tkradiobutton(Env$l.frames$Fr7,font=Env$police,variable=Env$l.var$add.distrib,value="wilcox",text=Env$voc[182,1],
    command=function() {type.distrib(dist="V (n)",par1=Env$voc[199,1],par2="",par3="",act2=FALSE,act3=FALSE)})
  Env$l.fr7$distrib.lab<-tklabel(Env$l.frames$Fr7,text="N (\U03BC,\U03C3)",font=Env$police6)
  Env$l.fr7$param1.lab<-tklabel(Env$l.frames$Fr7,text=paste(Env$voc[189,1],"\U03BC",Env$voc[200,1],sep=""),font=Env$police)
  Env$l.fr7$param1.wdg<-tkentry(Env$l.frames$Fr7,width=5,font=Env$police,textvariable=Env$l.var$add.param1)
  tkbind(Env$l.fr7$param1.wdg,"<Enter>",function() {msg(text=Env$voc[155,1],type="warning")})
  tkbind(Env$l.fr7$param1.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr7$param2.lab<-tklabel(Env$l.frames$Fr7,text=paste(Env$voc[190,1],"\U03C3",Env$voc[200,1],sep=""),font=Env$police)
  Env$l.fr7$param2.wdg<-tkentry(Env$l.frames$Fr7,width=5,font=Env$police,textvariable=Env$l.var$add.param2)
  tkbind(Env$l.fr7$param2.wdg,"<Enter>",function() {msg(text=Env$voc[155,1],type="warning")})
  tkbind(Env$l.fr7$param2.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr7$param3.lab<-tklabel(Env$l.frames$Fr7,text="",font=Env$police)
  Env$l.fr7$param3.wdg<-tkentry(Env$l.frames$Fr7,width=5,font=Env$police,textvariable=Env$l.var$add.param3,state="disabled")
  tkbind(Env$l.fr7$param3.wdg,"<Enter>",function() {msg(text=Env$voc[155,1],type="warning")})
  tkbind(Env$l.fr7$param3.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr7$trait.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[59,1],font=Env$police)
  Env$l.fr7$trait.wdg<-ttkcombobox(Env$l.frames$Fr7,font=Env$police,values=c(Env$voc[60:62,1]),textvariable=Env$l.var$add.trait,state="readonly")
  Env$l.fr7$epaisseur.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[63,1],font=Env$police)
  Env$l.fr7$epaisseur.wdg<-tkscale(Env$l.frames$Fr7,from=1,to=5,showvalue=TRUE,font=Env$police,variable=Env$l.var$add.epaisseur1,resolution=1,orient="horizontal")
  Env$l.fr7$col.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[45,1],font=Env$police)
  Env$l.fr7$col.wdg<-tkcanvas(Env$l.frames$Fr7,width="40",height="25",bg=tclvalue(Env$l.var$add.col1))
  tkbind(Env$l.fr7$col.wdg,"<ButtonRelease-1>",function() {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$add.col1),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$add.col1)<-temp
	tkconfigure(Env$l.fr7$col.wdg,bg=tclvalue(Env$l.var$add.col1))
    }
  })
  Env$l.fr7$tracer<-tkbutton(Env$l.frames$Fr7,width=16,text=Env$voc[72,1],font=Env$police,command=function() {
    if (dev.cur()>1) {
	if (!is.null(Env$l.var$add.seq)) {
	  variable<-""
	  sequence<-""
	  distrib<-""
	  if (nchar(tclvalue(Env$l.var$facteur1))>0 & tclvalue(Env$l.var$facteur1)!=Env$voc[82,1]) {
	    variable<-"variable"
	  } else {
	    variable<-tclvalue(Env$l.var$variable)
	  }
	  if (tclvalue(Env$l.var$add.distrib)%in%c("norm","gamma","expo","chi","fish","stud")) {
	    sequence<-paste("sequence <- seq(min(",variable,", na.rm=TRUE), max(",variable,", na.rm=TRUE), abs(max(",variable,", na.rm=TRUE) - min(",
		variable,", na.rm=TRUE)) / 1000)\n",sep="")
	  } else {
	    sequence<-paste("sequence <- floor(min(",variable,", na.rm=TRUE)) : ceiling(max(",variable,", na.rm=TRUE))\n",sep="")
	  }
	  if (tclvalue(Env$l.var$add.distrib)=="norm") {
	    lines(Env$l.var$add.seq,dnorm(Env$l.var$add.seq,as.numeric(tclvalue(Env$l.var$add.param1)),
		as.numeric(tclvalue(Env$l.var$add.param2))),lty=type.trait(type=tclvalue(Env$l.var$add.trait)),lwd=as.numeric(tclvalue(Env$l.var$add.epaisseur1)),
		col=tclvalue(Env$l.var$add.col1))
	    distrib<-paste("dnorm(sequence, mean=",tclvalue(Env$l.var$add.param1),", sd=",tclvalue(Env$l.var$add.param2),")",sep="")
	  } else if (tclvalue(Env$l.var$add.distrib)=="binom") {
	    lines(Env$l.var$add.seq2,dbinom(Env$l.var$add.seq2,as.numeric(tclvalue(Env$l.var$add.param1)),
		as.numeric(tclvalue(Env$l.var$add.param2))),lty=type.trait(type=tclvalue(Env$l.var$add.trait)),lwd=as.numeric(tclvalue(Env$l.var$add.epaisseur1)),
		col=tclvalue(Env$l.var$add.col1),type="o",pch=16)
	    distrib<-paste("dbinom(sequence, size=",tclvalue(Env$l.var$add.param1),", prob=",tclvalue(Env$l.var$add.param2),")",sep="")
	  } else if (tclvalue(Env$l.var$add.distrib)=="gamma") {
	    lines(Env$l.var$add.seq,dgamma(Env$l.var$add.seq,shape=as.numeric(tclvalue(Env$l.var$add.param1)),
		scale=as.numeric(tclvalue(Env$l.var$add.param2))),lty=type.trait(type=tclvalue(Env$l.var$add.trait)),lwd=as.numeric(tclvalue(Env$l.var$add.epaisseur1)),
		col=tclvalue(Env$l.var$add.col1))
	    distrib<-paste("dgamma(sequence, shape=",tclvalue(Env$l.var$add.param1),", rate=",tclvalue(Env$l.var$add.param2),")",sep="")
	  } else if (tclvalue(Env$l.var$add.distrib)=="poiss") {
	    lines(Env$l.var$add.seq2,dpois(Env$l.var$add.seq2,as.numeric(tclvalue(Env$l.var$add.param1))*
		as.numeric(tclvalue(Env$l.var$add.param2))),lty=type.trait(type=tclvalue(Env$l.var$add.trait)),lwd=as.numeric(tclvalue(Env$l.var$add.epaisseur1)),
		col=tclvalue(Env$l.var$add.col1),type="o",pch=16)
	    distrib<-paste("dpois(sequence, lambda=",tclvalue(Env$l.var$add.param1)," * ",tclvalue(Env$l.var$add.param2),")",sep="")
	  } else if (tclvalue(Env$l.var$add.distrib)=="expo") {
	    lines(Env$l.var$add.seq,dexp(Env$l.var$add.seq,as.numeric(tclvalue(Env$l.var$add.param1))),
		lty=type.trait(type=tclvalue(Env$l.var$add.trait)),lwd=as.numeric(tclvalue(Env$l.var$add.epaisseur1)),
		col=tclvalue(Env$l.var$add.col1))
	    distrib<-paste("dexp(sequence, rate=",tclvalue(Env$l.var$add.param1),")",sep="")
	  } else if (tclvalue(Env$l.var$add.distrib)=="nbinom") {
	    lines(Env$l.var$add.seq2,dnbinom(Env$l.var$add.seq2,as.numeric(tclvalue(Env$l.var$add.param1)),
		as.numeric(tclvalue(Env$l.var$add.param2))),lty=type.trait(type=tclvalue(Env$l.var$add.trait)),lwd=as.numeric(tclvalue(Env$l.var$add.epaisseur1)),
		col=tclvalue(Env$l.var$add.col1),type="o",pch=16)
	    distrib<-paste("dnbinom(sequence, size=",tclvalue(Env$l.var$add.param1),", prob=",tclvalue(Env$l.var$add.param2),")",sep="")
	  } else if (tclvalue(Env$l.var$add.distrib)=="chi") {
	    lines(Env$l.var$add.seq,dchisq(Env$l.var$add.seq,as.numeric(tclvalue(Env$l.var$add.param1))),
		lty=type.trait(type=tclvalue(Env$l.var$add.trait)),lwd=as.numeric(tclvalue(Env$l.var$add.epaisseur1)),
		col=tclvalue(Env$l.var$add.col1))
	    distrib<-paste("dchisq(sequence, df=",tclvalue(Env$l.var$add.param1),")",sep="")
	  } else if (tclvalue(Env$l.var$add.distrib)=="geom") {
	    lines(Env$l.var$add.seq2,dgeom(Env$l.var$add.seq2,as.numeric(tclvalue(Env$l.var$add.param1))),
		lty=type.trait(type=tclvalue(Env$l.var$add.trait)),lwd=as.numeric(tclvalue(Env$l.var$add.epaisseur1)),
		col=tclvalue(Env$l.var$add.col1),type="o",pch=16)
	    distrib<-paste("dgeom(sequence, prob=",tclvalue(Env$l.var$add.param1),")",sep="")
	  } else if (tclvalue(Env$l.var$add.distrib)=="fish") {
	    lines(Env$l.var$add.seq,df(Env$l.var$add.seq,as.numeric(tclvalue(Env$l.var$add.param1)),
		as.numeric(tclvalue(Env$l.var$add.param2))),lty=type.trait(type=tclvalue(Env$l.var$add.trait)),lwd=as.numeric(tclvalue(Env$l.var$add.epaisseur1)),
		col=tclvalue(Env$l.var$add.col1))
	    distrib<-paste("df(sequence, df1=",tclvalue(Env$l.var$add.param1),", df2=",tclvalue(Env$l.var$add.param2),")",sep="")
	  } else if (tclvalue(Env$l.var$add.distrib)=="hyper") {
	    lines(Env$l.var$add.seq2,dhyper(Env$l.var$add.seq2,round(as.numeric(tclvalue(Env$l.var$add.param3))*
		as.numeric(tclvalue(Env$l.var$add.param1)),0),round((1-as.numeric(tclvalue(Env$l.var$add.param3)))*as.numeric(tclvalue(Env$l.var$add.param1)),0),
		as.numeric(tclvalue(Env$l.var$add.param2))),lty=type.trait(type=tclvalue(Env$l.var$add.trait)),lwd=as.numeric(tclvalue(Env$l.var$add.epaisseur1)),
		col=tclvalue(Env$l.var$add.col1),type="o",pch=16)
	    distrib<-paste("dhyper(sequence, m=round(",tclvalue(Env$l.var$add.param3)," * ",tclvalue(Env$l.var$add.param1),", 0), n=round((1 - ",tclvalue(Env$l.var$add.param3),") * ",
		tclvalue(Env$l.var$add.param1),", 0), k=",tclvalue(Env$l.var$add.param2),")",sep="")
	  } else if (tclvalue(Env$l.var$add.distrib)=="stud") {
	    lines(Env$l.var$add.seq,dt(Env$l.var$add.seq,as.numeric(tclvalue(Env$l.var$add.param1))),
		lty=type.trait(type=tclvalue(Env$l.var$add.trait)),lwd=as.numeric(tclvalue(Env$l.var$add.epaisseur1)),
		col=tclvalue(Env$l.var$add.col1))
	    distrib<-paste("dt(sequence, df=",tclvalue(Env$l.var$add.param1),")",sep="")
	  } else if (tclvalue(Env$l.var$add.distrib)=="mann") {
	    lines(Env$l.var$add.seq2,dwilcox(Env$l.var$add.seq2,as.numeric(tclvalue(Env$l.var$add.param1)),
		as.numeric(tclvalue(Env$l.var$add.param2))),lty=type.trait(type=tclvalue(Env$l.var$add.trait)),lwd=as.numeric(tclvalue(Env$l.var$add.epaisseur1)),
		col=tclvalue(Env$l.var$add.col1),type="o",pch=16)
	    distrib<-paste("dwilcox(sequence, m=",tclvalue(Env$l.var$add.param1),", n=",tclvalue(Env$l.var$add.param2),")",sep="")
	  } else if (tclvalue(Env$l.var$add.distrib)=="wilcox") {
	    lines(Env$l.var$add.seq2,dsignrank(Env$l.var$add.seq2,as.numeric(tclvalue(Env$l.var$add.param1))),
		lty=type.trait(type=tclvalue(Env$l.var$add.trait)),lwd=as.numeric(tclvalue(Env$l.var$add.epaisseur1)),
		col=tclvalue(Env$l.var$add.col1),type="o",pch=16)
	    distrib<-paste("dsignrank(sequence, n=",tclvalue(Env$l.var$add.param1),")",sep="")
	  }
	  if (Env$l.code$save==TRUE) {
	    sink(file=file.path(Env$l.code$folder,paste(paste("GrapheR",paste(strsplit(as.character(Sys.Date()),split="-")[[1]],collapse="."),
		sep="-"),".R",sep=""),fsep=.Platform$file.sep),append=TRUE)
	    cat("# Added: theoretical distribution curve\n\n")
	    cat(sequence)
	    texte<-paste("lines(sequence, ",distrib,sep="")
	    if (tclvalue(Env$l.var$add.col1)!="black" & tclvalue(Env$l.var$add.col1)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$add.col1),"\"",sep="")}
	    texte<-paste(texte,", lty=",type.trait(type=tclvalue(Env$l.var$add.trait)),sep="")
	    texte<-paste(texte,", lwd=",tclvalue(Env$l.var$add.epaisseur1),sep="")
	    if (tclvalue(Env$l.var$add.distrib)%in%c("binom","pois","nbinom","geom","hyper","mann","wilcox")) {texte<-paste(texte,", type=\"o\", pch=16",sep="")}
	    cat(paste(texte,")\n\n",sep=""))
	    sink(NULL)
	  }
	} else {
	  msg(text=Env$voc[201,1],type="error")
	}
    } else {
	msg(text=Env$voc[163,1],type="error")
    }
  })
  Env$l.fr7$fermer<-tkbutton(Env$l.frames$Fr7,width=16,text=Env$voc[152,1],font=Env$police,command=fr7.close)
  Env$l.fr7$espace.ver1<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver2<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver3<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver4<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver5<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver6<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver7<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.hor1<-tklabel(Env$l.frames$Fr7,text="     ",font=Env$police3)
  Env$l.fr7$espace.hor2<-tklabel(Env$l.frames$Fr7,text="     ",font=Env$police3)
  tkgrid(Env$l.fr7$espace.hor1,row=1,column=0)
  tkgrid(Env$l.fr7$titre.lab,row=1,column=1,columnspan=3)
  tkgrid(Env$l.fr7$espace.hor2,row=1,column=4)
  tkgrid(Env$l.fr7$espace.ver1)
  tkgrid(Env$l.fr7$rb1,row=3,column=1,sticky="w")
  tkgrid(Env$l.fr7$rb2,row=3,column=2,sticky="w")
  tkgrid(Env$l.fr7$rb3,row=4,column=1,sticky="w")
  tkgrid(Env$l.fr7$rb4,row=4,column=2,sticky="w")
  tkgrid(Env$l.fr7$rb5,row=5,column=1,sticky="w")
  tkgrid(Env$l.fr7$rb6,row=5,column=2,sticky="w")
  tkgrid(Env$l.fr7$rb7,row=6,column=1,sticky="w")
  tkgrid(Env$l.fr7$rb8,row=6,column=2,sticky="w")
  tkgrid(Env$l.fr7$rb9,row=7,column=1,sticky="w")
  tkgrid(Env$l.fr7$rb10,row=7,column=2,sticky="w")
  tkgrid(Env$l.fr7$rb11,row=8,column=1,sticky="w")
  tkgrid(Env$l.fr7$espace.ver2)
  tkgrid(Env$l.fr7$rb12,row=10,column=1,columnspan=2,sticky="w")
  tkgrid(Env$l.fr7$rb13,row=11,column=1,columnspan=2,sticky="w")
  tkgrid(Env$l.fr7$espace.ver3)
  tkgrid(Env$l.fr7$distrib.lab,row=13,column=1,columnspan=2)
  tkgrid(Env$l.fr7$espace.ver4)
  tkgrid(Env$l.fr7$param1.lab,row=15,column=1,sticky="e")
  tkgrid(Env$l.fr7$param1.wdg,row=15,column=2,sticky="w")
  tkgrid(Env$l.fr7$param2.lab,row=16,column=1,sticky="e")
  tkgrid(Env$l.fr7$param2.wdg,row=16,column=2,sticky="w")
  tkgrid(Env$l.fr7$param3.lab,row=17,column=1,sticky="e")
  tkgrid(Env$l.fr7$param3.wdg,row=17,column=2,sticky="w")
  tkgrid(Env$l.fr7$espace.ver5)
  tkgrid(Env$l.fr7$trait.lab,row=19,column=1,sticky="e")
  tkgrid(Env$l.fr7$trait.wdg,row=19,column=2,sticky="w")
  tkgrid(Env$l.fr7$epaisseur.lab,row=20,column=1,sticky="e")
  tkgrid(Env$l.fr7$epaisseur.wdg,row=20,column=2,sticky="w")
  tkgrid(Env$l.fr7$col.lab,row=21,column=1,sticky="e")
  tkgrid(Env$l.fr7$col.wdg,row=21,column=2,sticky="w")
  tkgrid(Env$l.fr7$espace.ver6)
  tkgrid(Env$l.fr7$tracer,row=23,column=1)
  tkgrid(Env$l.fr7$fermer,row=23,column=2)
  tkgrid(Env$l.fr7$espace.ver7)
}


#---------------------------------------------------------
# Ajouter du texte sur le graphe
#---------------------------------------------------------

texte<-function() {
  fr7.close()
  for (i in 1:length(Env$l.fr7)) {tkdestroy(Env$l.fr7[[i]])}
  Env$l.fr7<-list()
  tkconfigure(Env$l.frames$Fr7,borderwidth=3)
  Env$l.fr7$txt.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[202,1],font=Env$police)
  Env$l.fr7$txt.wdg<-tkentry(Env$l.frames$Fr7,width=33,font=Env$police,textvariable=Env$l.var$add.txt)
  Env$l.fr7$coords.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[203,1],font=Env$police)
  Env$l.fr7$x.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[204,1],font=Env$police)
  Env$l.fr7$x.wdg<-tkentry(Env$l.frames$Fr7,width=8,font=Env$police,textvariable=Env$l.var$add.param1)
  tkbind(Env$l.fr7$x.wdg,"<Enter>",function() {msg(text=Env$voc[155,1],type="warning")})
  tkbind(Env$l.fr7$x.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr7$y.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[205,1],font=Env$police)
  Env$l.fr7$y.wdg<-tkentry(Env$l.frames$Fr7,width=8,font=Env$police,textvariable=Env$l.var$add.param2)
  tkbind(Env$l.fr7$y.wdg,"<Enter>",function() {msg(text=Env$voc[155,1],type="warning")})
  tkbind(Env$l.fr7$y.wdg,"<Leave>",function() {msg(text="",type="info")})
  Env$l.fr7$coords.but<-tkbutton(Env$l.frames$Fr7,text=Env$voc[162,1],width=20,font=Env$police,command=function() {
    if (dev.cur()>1) {
	coords<-locator(n=1)
	tclvalue(Env$l.var$add.param1)<-round(coords$x[1],4)
	tclvalue(Env$l.var$add.param2)<-round(coords$y[1],4)
    } else {
	msg(text=Env$voc[163,1],type="error")
    }
  })
  Env$l.fr7$taille.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[46,1],font=Env$police)
  Env$l.fr7$taille.wdg<-tkscale(Env$l.frames$Fr7,from=0.5,to=5,showvalue=TRUE,font=Env$police,variable=Env$l.var$add.epaisseur1,resolution=0.25,orient="horizontal")
  Env$l.fr7$col.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[45,1],font=Env$police)
  Env$l.fr7$col.wdg<-tkbutton(Env$l.frames$Fr7,text="Aa",font=Env$police3,foreground=tclvalue(Env$l.var$add.col1),activeforeground=tclvalue(Env$l.var$add.col1),command=function() {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$add.col1),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$add.col1)<-temp
	tkconfigure(Env$l.fr7$col.wdg,foreground=temp,activeforeground=temp)
    }
  })
  Env$l.fr7$tracer<-tkbutton(Env$l.frames$Fr7,text=Env$voc[72,1],font=Env$police,width=16,command=function() {
    if (dev.cur()>1) {
	text(as.numeric(tclvalue(Env$l.var$add.param1)),as.numeric(tclvalue(Env$l.var$add.param2)),
	  labels=tclvalue(Env$l.var$add.txt),cex=as.numeric(tclvalue(Env$l.var$add.epaisseur1)),
	  col=tclvalue(Env$l.var$add.col1))
	if (Env$l.code$save==TRUE) {
	  sink(file=file.path(Env$l.code$folder,paste(paste("GrapheR",paste(strsplit(as.character(Sys.Date()),split="-")[[1]],collapse="."),
	    sep="-"),".R",sep=""),fsep=.Platform$file.sep),append=TRUE)
	  cat("# Added: text\n\n")
	  texte<-paste("text(",tclvalue(Env$l.var$add.param1),", ",tclvalue(Env$l.var$add.param2),sep="")
	  texte<-paste(texte,", labels=\"",tclvalue(Env$l.var$add.txt),"\"",sep="")
	  if (tclvalue(Env$l.var$add.col1)!="black" & tclvalue(Env$l.var$add.col1)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$add.col1),"\"",sep="")}
	  texte<-paste(texte,", cex=",tclvalue(Env$l.var$add.epaisseur1),")\n\n",sep="")
	  cat(texte)
	  sink(NULL)
	}
    } else {
	msg(text=Env$voc[163,1],type="error")
    }
  })
  Env$l.fr7$fermer<-tkbutton(Env$l.frames$Fr7,text=Env$voc[152,1],font=Env$police,width=16,command=fr7.close)
  Env$l.fr7$espace.ver1<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver2<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver3<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.hor1<-tklabel(Env$l.frames$Fr7,text="     ",font=Env$police3)
  Env$l.fr7$espace.hor2<-tklabel(Env$l.frames$Fr7,text="     ",font=Env$police3)
  Env$l.fr7$espace.hor3<-tklabel(Env$l.frames$Fr7,text="   ",font=Env$police3)
  tkgrid(Env$l.fr7$espace.ver1)
  tkgrid(Env$l.fr7$espace.hor1,row=1,column=0)
  tkgrid(Env$l.fr7$txt.lab,row=1,column=1,sticky="e")
  tkgrid(Env$l.fr7$txt.wdg,row=1,column=2,columnspan=3,sticky="w")
  tkgrid(Env$l.fr7$espace.hor2,row=1,column=5)
  tkgrid(Env$l.fr7$coords.lab,row=2,column=1,sticky="e")
  tkgrid(Env$l.fr7$x.lab,row=3,column=1,sticky="e")
  tkgrid(Env$l.fr7$x.wdg,row=3,column=2)
  tkgrid(Env$l.fr7$y.lab,row=4,column=1,sticky="e")
  tkgrid(Env$l.fr7$y.wdg,row=4,column=2)
  tkgrid(Env$l.fr7$espace.hor3,row=3,column=3)
  tkgrid(Env$l.fr7$coords.but,row=3,column=4,rowspan=2)
  tkgrid(Env$l.fr7$taille.lab,row=5,column=1,sticky="e")
  tkgrid(Env$l.fr7$taille.wdg,row=5,column=2,columnspan=3,sticky="w")
  tkgrid(Env$l.fr7$col.lab,row=6,column=1,sticky="e")
  tkgrid(Env$l.fr7$col.wdg,row=6,column=2,columnspan=3,sticky="w")
  tkgrid(Env$l.fr7$espace.ver2)
  tkgrid(Env$l.fr7$tracer,row=8,column=1,columnspan=2,sticky="e")
  tkgrid(Env$l.fr7$fermer,row=8,column=4)
  tkgrid(Env$l.fr7$espace.ver3)
}


#---------------------------------------------------------
# Ajouter des p-values sur le graphe
#---------------------------------------------------------

pval<-function() {
  fr7.close()
  for (i in 1:length(Env$l.fr7)) {tkdestroy(Env$l.fr7[[i]])}
  Env$l.fr7<-list()
  tkconfigure(Env$l.frames$Fr7,borderwidth=3)
  Env$l.fr7$titre.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[210,1],font=Env$police3)
  Env$l.fr7$img1<-tklabel(Env$l.frames$Fr7,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Pvalue1.gif",fsep=.Platform$file.sep)))
  Env$l.fr7$img2<-tklabel(Env$l.frames$Fr7,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Pvalue2.gif",fsep=.Platform$file.sep)))
  Env$l.fr7$txt1.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[202,1],font=Env$police)
  Env$l.fr7$txt1.wdg<-tkentry(Env$l.frames$Fr7,width=10,textvariable=Env$l.var$add.param1,font=Env$police)
  Env$l.fr7$txt2.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[202,1],font=Env$police)
  Env$l.fr7$txt2.wdg<-tkentry(Env$l.frames$Fr7,width=10,textvariable=Env$l.var$add.param2,font=Env$police)
  Env$l.fr7$taille1.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[46,1],font=Env$police)
  Env$l.fr7$taille1.wdg<-tkscale(Env$l.frames$Fr7,from=0.5,to=3,showvalue=TRUE,font=Env$police,variable=Env$l.var$add.epaisseur1,resolution=0.25,orient="horizontal")
  Env$l.fr7$taille2.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[46,1],font=Env$police)
  Env$l.fr7$taille2.wdg<-tkscale(Env$l.frames$Fr7,from=0.5,to=3,showvalue=TRUE,font=Env$police,variable=Env$l.var$add.epaisseur2,resolution=0.25,orient="horizontal")
  Env$l.fr7$col1.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[45,1],font=Env$police)
  Env$l.fr7$col1.wdg<-tkbutton(Env$l.frames$Fr7,text="Aa",font=Env$police3,command=function() {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$add.col1),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$add.col1)<-temp
	tkconfigure(Env$l.fr7$col1.wdg,foreground=temp,activeforeground=temp)
    }
  })
  Env$l.fr7$col2.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[45,1],font=Env$police)
  Env$l.fr7$col2.wdg<-tkbutton(Env$l.frames$Fr7,text="Aa",font=Env$police3,command=function() {
    temp<-tclvalue(tcl("tk_chooseColor",initialcolor=tclvalue(Env$l.var$add.col2),title=Env$voc[64,1]))
    if (nchar(temp)>0) {
	tclvalue(Env$l.var$add.col2)<-temp
	tkconfigure(Env$l.fr7$col2.wdg,foreground=temp,activeforeground=temp)
    }
  })
  Env$l.fr7$expl1.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[206,1],font=Env$police)
  Env$l.fr7$expl2.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[207,1],font=Env$police)
  if (!is.null(Env$l.var$add.hauteurs)) {
    ecart1<-0.15*max(Env$l.var$add.hauteurs)
    ecart2<-0.08*max(Env$l.var$add.hauteurs)
  }
  Env$l.fr7$but1<-tkbutton(Env$l.frames$Fr7,text=Env$voc[208,1],font=Env$police,width=16,command=function(...) {
    if (dev.cur()>1) {
	if (!is.null(Env$l.var$add.abscisses)) {
	  coords<-locator(n=2)
	  diff1<-numeric(length(Env$l.var$add.abscisses))
	  diff2<-numeric(length(Env$l.var$add.abscisses))
	  for (i in 1:length(Env$l.var$add.abscisses)) {
	    diff1[i]<-abs(coords$x[1]-Env$l.var$add.abscisses[i])
	    diff2[i]<-abs(coords$x[2]-Env$l.var$add.abscisses[i])
	  }
	  x1<-which(diff1==min(diff1))
	  x2<-which(diff2==min(diff2))
	  segments(min(Env$l.var$add.abscisses[x1],Env$l.var$add.abscisses[x2]),Env$l.var$add.matrice[x2,x1]+ecart1,
	    max(Env$l.var$add.abscisses[x1],Env$l.var$add.abscisses[x2]),Env$l.var$add.matrice[x2,x1]+ecart1)
	  segments(Env$l.var$add.abscisses[x1],Env$l.var$add.hauteurs[x1]+ecart2,Env$l.var$add.abscisses[x1],
	    Env$l.var$add.matrice[x2,x1]+ecart1)
	  segments(Env$l.var$add.abscisses[x2],Env$l.var$add.hauteurs[x2]+ecart2,Env$l.var$add.abscisses[x2],
	    Env$l.var$add.matrice[x2,x1]+ecart1)
	  text(Env$l.var$add.abscisses[x1]+(Env$l.var$add.abscisses[x2]-Env$l.var$add.abscisses[x1])/2,
	    Env$l.var$add.matrice[x2,x1]+ecart1+ecart2,tclvalue(Env$l.var$add.param1),
	    cex=as.numeric(tclvalue(Env$l.var$add.epaisseur1)),col=tclvalue(Env$l.var$add.col1))
	  if (Env$l.code$save==TRUE) {
	    sink(file=file.path(Env$l.code$folder,paste(paste("GrapheR",paste(strsplit(as.character(Sys.Date()),split="-")[[1]],collapse="."),
		sep="-"),".R",sep=""),fsep=.Platform$file.sep),append=TRUE)
	    cat("# Added: p-value\n\n")
	    cat(paste("segments(",round(min(Env$l.var$add.abscisses[x1],Env$l.var$add.abscisses[x2]),2),", ",round(Env$l.var$add.matrice[x2,x1]+ecart1,2),
		", ",round(max(Env$l.var$add.abscisses[x1],Env$l.var$add.abscisses[x2]),2),", ",round(Env$l.var$add.matrice[x2,x1]+ecart1,2),")\n",sep=""))
	    cat(paste("segments(",round(Env$l.var$add.abscisses[x1],2),", ",round(Env$l.var$add.hauteurs[x1]+ecart2,2),", ",round(Env$l.var$add.abscisses[x1],2),", ",round(Env$l.var$add.matrice[x2,x1]+ecart1,2),")\n",sep=""))
	    cat(paste("segments(",round(Env$l.var$add.abscisses[x2],2),", ",round(Env$l.var$add.hauteurs[x2]+ecart2,2),", ",round(Env$l.var$add.abscisses[x2],2),", ",round(Env$l.var$add.matrice[x2,x1]+ecart1,2),")\n",sep=""))
	    texte<-paste("text(",round(Env$l.var$add.abscisses[x1]+(Env$l.var$add.abscisses[x2]-Env$l.var$add.abscisses[x1])/2,2),", ",round(Env$l.var$add.matrice[x2,x1]+ecart1+ecart2,2),sep="")
	    texte<-paste(texte,", labels=\"",tclvalue(Env$l.var$add.param1),"\"",sep="")
	    if (tclvalue(Env$l.var$add.col1)!="black" & tclvalue(Env$l.var$add.col1)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$add.col1),"\"",sep="")}
	    texte<-paste(texte,", cex=",tclvalue(Env$l.var$add.epaisseur1),")\n\n",sep="")
	    cat(texte)
	    sink(NULL)
	  }
	  Env$l.var$add.hauteurs[x1:x2]<-Env$l.var$add.matrice[x2,x1]+ecart1+ecart2
	  for (i in 1:length(Env$l.var$add.abscisses)) {
	    for (j in 1:length(Env$l.var$add.abscisses)) {
		Env$l.var$add.matrice[i,j]<-max(Env$l.var$add.hauteurs[i:j])
	    }
	  }
	} else {
	  msg(text=Env$voc[209,1],type="error")
	}
    } else {
	msg(text=Env$voc[163,1],type="error")
    }
  })
  Env$l.fr7$but2<-tkbutton(Env$l.frames$Fr7,text=Env$voc[208,1],font=Env$police,width=16,command=function(...) {
    if (dev.cur()>1) {
	if (!is.null(Env$l.var$add.abscisses)) {
	  coords<-locator(n=2)
	  diff1<-numeric(length(Env$l.var$add.abscisses))
	  diff2<-numeric(length(Env$l.var$add.abscisses))
	  for (i in 1:length(Env$l.var$add.abscisses)) {
	    diff1[i]<-abs(coords$x[1]-Env$l.var$add.abscisses[i])
	    diff2[i]<-abs(coords$x[2]-Env$l.var$add.abscisses[i])
	  }
	  x1<-which(diff1==min(diff1))
	  x2<-which(diff2==min(diff2))
	  segments(min(Env$l.var$add.abscisses[x1],Env$l.var$add.abscisses[x2])-0.5,Env$l.var$add.matrice[x2,x1]+ecart1,
	    max(Env$l.var$add.abscisses[x1],Env$l.var$add.abscisses[x2])+0.5,Env$l.var$add.matrice[x2,x1]+ecart1)
	  text(Env$l.var$add.abscisses[x1]+(Env$l.var$add.abscisses[x2]-Env$l.var$add.abscisses[x1])/2,
	    Env$l.var$add.matrice[x2,x1]+ecart1+ecart2,tclvalue(Env$l.var$add.param2),
	    cex=as.numeric(tclvalue(Env$l.var$add.epaisseur2)),col=tclvalue(Env$l.var$add.col2))
	  if (Env$l.code$save==TRUE) {
	    sink(file=file.path(Env$l.code$folder,paste(paste("GrapheR",paste(strsplit(as.character(Sys.Date()),split="-")[[1]],collapse="."),
		sep="-"),".R",sep=""),fsep=.Platform$file.sep),append=TRUE)
	    cat("# Added: p-value\n\n")
	    cat(paste("segments(",round(min(Env$l.var$add.abscisses[x1],Env$l.var$add.abscisses[x2])-0.5,2),", ",round(Env$l.var$add.matrice[x2,x1]+ecart1,2),
		", ",round(max(Env$l.var$add.abscisses[x1],Env$l.var$add.abscisses[x2])+0.5,2),", ",round(Env$l.var$add.matrice[x2,x1]+ecart1,2),")\n",sep=""))
	    texte<-paste("text(",round(Env$l.var$add.abscisses[x1]+(Env$l.var$add.abscisses[x2]-Env$l.var$add.abscisses[x1])/2,2),", ",round(Env$l.var$add.matrice[x2,x1]+ecart1+ecart2,2),sep="")
	    texte<-paste(texte,", labels=\"",tclvalue(Env$l.var$add.param2),"\"",sep="")
	    if (tclvalue(Env$l.var$add.col2)!="black" & tclvalue(Env$l.var$add.col2)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$add.col2),"\"",sep="")}
	    texte<-paste(texte,", cex=",tclvalue(Env$l.var$add.epaisseur2),")\n\n",sep="")
	    cat(texte)
	    sink(NULL)
	  }
	  Env$l.var$add.hauteurs[x1:x2]<-Env$l.var$add.matrice[x2,x1]+ecart1+ecart2
	  for (i in 1:length(Env$l.var$add.abscisses)) {
	    for (j in 1:length(Env$l.var$add.abscisses)) {
		Env$l.var$add.matrice[i,j]<-max(Env$l.var$add.hauteurs[i:j])
	    }
	  }
	} else {
	  msg(text=Env$voc[209,1],type="error")
	}
    } else {
	msg(text=Env$voc[163,1],type="error")
    }
  })
  Env$l.fr7$fermer<-tkbutton(Env$l.frames$Fr7,text=Env$voc[152,1],font=Env$police,width=16,command=fr7.close)
  Env$l.fr7$espace.ver1<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver2<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver3<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver4<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver5<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.hor1<-tklabel(Env$l.frames$Fr7,text="     ",font=Env$police3)
  Env$l.fr7$espace.hor2<-tklabel(Env$l.frames$Fr7,text="     ",font=Env$police3)
  Env$l.fr7$espace.hor3<-tklabel(Env$l.frames$Fr7,text="               ",font=Env$police3)
  tkgrid(Env$l.fr7$espace.hor1,row=0,column=0)
  tkgrid(Env$l.fr7$titre.lab,row=0,column=1,columnspan=5)
  tkgrid(Env$l.fr7$espace.hor2,row=0,column=6)
  tkgrid(Env$l.fr7$espace.ver1)
  tkgrid(Env$l.fr7$img1,row=2,column=1,columnspan=2)
  tkgrid(Env$l.fr7$espace.hor3,row=2,column=3)
  tkgrid(Env$l.fr7$img2,row=2,column=4,columnspan=2)
  tkgrid(Env$l.fr7$espace.ver2)
  tkgrid(Env$l.fr7$txt1.lab,row=4,column=1,sticky="e")
  tkgrid(Env$l.fr7$txt1.wdg,row=4,column=2,sticky="w")
  tkgrid(Env$l.fr7$txt2.lab,row=4,column=4,sticky="e")
  tkgrid(Env$l.fr7$txt2.wdg,row=4,column=5,sticky="w")
  tkgrid(Env$l.fr7$taille1.lab,row=5,column=1,sticky="e")
  tkgrid(Env$l.fr7$taille1.wdg,row=5,column=2,sticky="w")
  tkgrid(Env$l.fr7$taille2.lab,row=5,column=4,sticky="e")
  tkgrid(Env$l.fr7$taille2.wdg,row=5,column=5,sticky="w")
  tkgrid(Env$l.fr7$col1.lab,row=6,column=1,sticky="e")
  tkgrid(Env$l.fr7$col1.wdg,row=6,column=2,sticky="w")
  tkgrid(Env$l.fr7$col2.lab,row=6,column=4,sticky="e")
  tkgrid(Env$l.fr7$col2.wdg,row=6,column=5,sticky="w")
  tkgrid(Env$l.fr7$espace.ver3)
  tkgrid(Env$l.fr7$expl1,row=8,column=1,columnspan=2)
  tkgrid(Env$l.fr7$expl2,row=8,column=4,columnspan=2)
  tkgrid(Env$l.fr7$but1,row=9,column=1,columnspan=2)
  tkgrid(Env$l.fr7$but2,row=9,column=4,columnspan=2)
  tkgrid(Env$l.fr7$espace.ver4)
  tkgrid(Env$l.fr7$fermer,column=1,columnspan=5)
  tkgrid(Env$l.fr7$espace.ver5)
}


#--------------------------------------------------
# Enregistrer le graphe
#--------------------------------------------------

enregistrer<-function() {
  fr7.close()
  for (i in 1:length(Env$l.fr7)) {tkdestroy(Env$l.fr7[[i]])}
  Env$l.fr7<-list()
  tkconfigure(Env$l.frames$Fr7,borderwidth=3)
  Env$l.fr7$fenetre.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[211,1],font=Env$police)
  fenetres<-character(length(dev.list()))
  for (i in 1:length(dev.list())) {
    fenetres[i]<-paste(Env$voc[212,1],dev.list()[i],sep="_")
  }
  Env$l.fr7$fenetre.wdg<-ttkcombobox(Env$l.frames$Fr7,font=Env$police,values=fenetres,textvariable=Env$l.var$fen.num,state="readonly")
  Env$l.fr7$fichier.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[213,1],font=Env$police)
  formats <- if(Env$figurej) {c("jpg","png","tiff","FigureJ")} else {c("jpg","png","tiff")}
  Env$l.fr7$fichier.wdg<-ttkcombobox(Env$l.frames$Fr7,width=8,font=Env$police,values=formats,textvariable=Env$l.var$fen.type,state="readonly")
  Env$l.fr7$largeur.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[214,1],font=Env$police)
  Env$l.fr7$largeur.wdg<-tkscale(Env$l.frames$Fr7,from=400,to=5000,showvalue=TRUE,font=Env$police,variable=Env$l.var$fen.larg,resolution=50,orient="horizontal")
  Env$l.fr7$res.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[265,1],font=Env$police)
  Env$l.fr7$res.wdg<-ttkcombobox(Env$l.frames$Fr7,width=8,font=Env$police,values=c("72","150","300","600"),textvariable=Env$l.var$fen.res,state="readonly")
  tkbind(Env$l.fr7$res.wdg,"<<ComboboxSelected>>",function() {
    if (tclvalue(Env$l.var$fen.res)=="72") {
	tclvalue(Env$l.var$fen.larg)=="600"
	tkset(Env$l.fr7$largeur.wdg,"600")
    }
    if (tclvalue(Env$l.var$fen.res)=="150") {
	tclvalue(Env$l.var$fen.larg)=="1250"
	tkset(Env$l.fr7$largeur.wdg,"1250")
    }
    if (tclvalue(Env$l.var$fen.res)=="300") {
	tclvalue(Env$l.var$fen.larg)=="2500"
	tkset(Env$l.fr7$largeur.wdg,"2500")
    }
    if (tclvalue(Env$l.var$fen.res)=="600") {
	tclvalue(Env$l.var$fen.larg)=="5000"
	tkset(Env$l.fr7$largeur.wdg,"5000")
    }
  })
  Env$l.fr7$sauver<-tkbutton(Env$l.frames$Fr7,text=Env$voc[215,1],font=Env$police,width=16,command=function() {
    if (dev.cur()>1) {
	if (nchar(tclvalue(Env$l.var$fen.num))>0) {
	  dev.set(as.numeric(strsplit(tclvalue(Env$l.var$fen.num),split="_")[[1]][2]))
	} else {
	  dev.set(dev.list()[1])
	}
	if (tclvalue(Env$l.var$fen.type)=="jpg") {
	  file<-tclvalue(tkgetSaveFile(filetypes="{Jpg {.jpg}}"))
	  if (nchar(file)>0) {
	    dev.print(jpeg,filename=paste(strsplit(file,".jpg"),".jpg",sep=""),quality=100,units="px",
		width=as.numeric(tclvalue(Env$l.var$fen.larg)),res=as.numeric(tclvalue(Env$l.var$fen.res)))
	  }
	}
	if (tclvalue(Env$l.var$fen.type)=="png") {
	  file<-tclvalue(tkgetSaveFile(filetypes="{Png {.png}}"))
	  if (nchar(file)>0) {
	    dev.print(png,filename=paste(strsplit(file,".png"),".png",sep=""),units="px",
		width=as.numeric(tclvalue(Env$l.var$fen.larg)),res=as.numeric(tclvalue(Env$l.var$fen.res)))
	  }
	}
	if (tclvalue(Env$l.var$fen.type)=="tiff") {
	  file<-tclvalue(tkgetSaveFile(filetypes="{Tiff {.tiff}}"))
	  if (nchar(file)>0) {
	    dev.print(tiff,filename=paste(strsplit(file,".tiff"),".tiff",sep=""),units="px",
		width=as.numeric(tclvalue(Env$l.var$fen.larg)),compression="none",res=as.numeric(tclvalue(Env$l.var$fen.res)))
	  }
	}
	if (tclvalue(Env$l.var$fen.type)=="FigureJ") {
	  file <- if (grepl("apple",Sys.getenv("R_PLATFORM"))) {
	    paste(Sys.getenv("HOME"),"/Library/Preferences/IJ_Prefs.txt",sep="")
	  } else { 
	    paste(Sys.getenv("HOME"),".imagej","IJ_Prefs.txt",sep=.Platform$file.sep)
	  }
	  if ( file.exists (file) ) {
	    # read prefs file into dataframe with keys/values
	    ijPrefs <- read.table(file,header=FALSE,sep="=",col.names=c("Key","Value"))
	    # read values from figurej panel properties
	    panelTempDir <- paste("",with(ijPrefs,Value[Key==".figurej.tempDir"]),sep="")
	    panelFilename <- paste("",with(ijPrefs,Value[Key==".figurej.panelFilename"]),sep="")
	    panelWidth <- paste("",with(ijPrefs,Value[Key==".figurej.panelWidth"]),sep="")
	    panelHeight <- paste("",with(ijPrefs,Value[Key==".figurej.panelHeight"]),sep="")	
	    if ((panelTempDir!="")&&(panelFilename!="")&&(panelWidth!="")&&(panelHeight!="")) {
		# save tif at defined location and resolution, and the code used for this graph
		if (nchar(file)>0) {
		  dev.print(tiff,filename=paste(panelTempDir,panelFilename,sep=""),units="px",width=as.numeric(panelWidth),height=as.numeric(panelHeight),res=as.numeric(tclvalue(Env$l.var$fen.res)))
		  msg(text="Panel exported",type="info")
		  sink(file=file.path(paste(panelTempDir,paste(strsplit(panelFilename,".tif"),".R",sep=""),sep=""),fsep=.Platform$file.sep),append=FALSE)
		  code.graphtype()
		  cat("# Loading of the dataset\n\n")
		  cat(paste("dataset <- ",Env$loading,"\n",sep=""))
		  cat("attach(dataset)\n\n")
		  code.data()
		  code.graph()
		  cat("detach(dataset)\n\n")
		  sink(NULL)
		}
	    } else {
		msg(text=Env$voc[266,1],type="error")
	    }
	  } else {
	    msg(text=Env$voc[266,1],type="error")
	  }
	}
    } else {
	msg(text=Env$voc[163,1],type="error")
    }
  })
  Env$l.fr7$fermer<-tkbutton(Env$l.frames$Fr7,text=Env$voc[152,1],font=Env$police,width=16,command=fr7.close)
  Env$l.fr7$espace.ver1<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver2<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver3<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.hor1<-tklabel(Env$l.frames$Fr7,text="     ",font=Env$police3)
  Env$l.fr7$espace.hor2<-tklabel(Env$l.frames$Fr7,text="     ",font=Env$police3)
  tkgrid(Env$l.fr7$espace.ver1)
  tkgrid(Env$l.fr7$espace.hor1,row=1,column=0)
  tkgrid(Env$l.fr7$fenetre.lab,row=1,column=1,sticky="e")
  tkgrid(Env$l.fr7$fenetre.wdg,row=1,column=2,sticky="w")
  tkgrid(Env$l.fr7$espace.hor2,row=1,column=3)
  tkgrid(Env$l.fr7$fichier.lab,row=2,column=1,sticky="e")
  tkgrid(Env$l.fr7$fichier.wdg,row=2,column=2,sticky="w")
  tkgrid(Env$l.fr7$largeur.lab,row=3,column=1,sticky="se")
  tkgrid(Env$l.fr7$largeur.wdg,row=3,column=2,sticky="w")
  tkgrid(Env$l.fr7$res.lab,row=4,column=1,sticky="se")
  tkgrid(Env$l.fr7$res.wdg,row=4,column=2,sticky="w")
  tkgrid(Env$l.fr7$espace.ver2)
  tkgrid(Env$l.fr7$sauver,row=6,column=1)
  tkgrid(Env$l.fr7$fermer,row=6,column=2)
  tkgrid(Env$l.fr7$espace.ver3)
}


#--------------------------------------------------------------------------
# Changer le langage utilisateur
#--------------------------------------------------------------------------

language.change<-function() {
  fr7.close()
  for (i in 1:length(Env$l.fr7)) {tkdestroy(Env$l.fr7[[i]])}
  Env$l.fr7<-list()
  tkconfigure(Env$l.frames$Fr7,borderwidth=3)
  lang.var<-tclVar("")
  save.var<-tclVar("1")
  Env$l.fr7$lang.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[216,1])
  Env$l.fr7$lang.wdg<-ttkcombobox(Env$l.frames$Fr7,width=10,font=Env$police,values=Env$voc[c(217,218,239,241),1],textvariable=lang.var,state="readonly")
  Env$l.fr7$save.lab<-tklabel(Env$l.frames$Fr7,text=Env$voc[219,1])
  Env$l.fr7$save.wdg<-tkcheckbutton(Env$l.frames$Fr7,variable=save.var)
  Env$l.fr7$ok<-tkbutton(Env$l.frames$Fr7,text=Env$voc[151,1],font=Env$police,width=16,command=function() {
    if (nchar(tclvalue(lang.var))>0) {
	langue<-NULL
	if (tclvalue(lang.var)==Env$voc[217,1]) {langue<-"en"} else
	if (tclvalue(lang.var)==Env$voc[218,1]) {langue<-"fr"} else
	if (tclvalue(lang.var)==Env$voc[239,1]) {langue<-"es"} else
	if (tclvalue(lang.var)==Env$voc[241,1]) {langue<-"de"}
	if (tclvalue(save.var)==1) {
	  write(langue,file=file.path(path.package("GrapheR"),"lang","Language.txt",fsep=.Platform$file.sep))
	}
	fermer.GrapheR()
	run.GrapheR(lang=langue,path.to.save=Env$path.to.save,figurej=Env$figurej)
    } else {
	msg(text=Env$voc[220,1],type="error")
    }
  })
  Env$l.fr7$fermer<-tkbutton(Env$l.frames$Fr7,text=Env$voc[152,1],font=Env$police,width=16,command=fr7.close)
  Env$l.fr7$espace.ver1<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver2<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.ver3<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  Env$l.fr7$espace.hor1<-tklabel(Env$l.frames$Fr7,text="     ",font=Env$police3)
  Env$l.fr7$espace.hor2<-tklabel(Env$l.frames$Fr7,text="     ",font=Env$police3)
  tkgrid(Env$l.fr7$espace.ver1)
  tkgrid(Env$l.fr7$espace.hor1,row=1,column=0)
  tkgrid(Env$l.fr7$lang.lab,row=1,column=1,sticky="e")
  tkgrid(Env$l.fr7$lang.wdg,row=1,column=2,sticky="w")
  tkgrid(Env$l.fr7$espace.hor2,row=1,column=3)
  tkgrid(Env$l.fr7$save.lab,row=2,column=1,sticky="e")
  tkgrid(Env$l.fr7$save.wdg,row=2,column=2,sticky="w")
  tkgrid(Env$l.fr7$espace.ver2)
  tkgrid(Env$l.fr7$ok,row=4,column=1)
  tkgrid(Env$l.fr7$fermer,row=4,column=2)
  tkgrid(Env$l.fr7$espace.ver3)
}


#-------------------------------------------------
# Graphe - ajouter le titre
#-------------------------------------------------

graphe.titre<-function(type="",orient=NULL) {
  if (type=="moust") {
    titre.x<-if (orient=="ver") {
	tclvalue(Env$l.var$titre.axenoms)
    } else {
	tclvalue(Env$l.var$titre.axevaleurs)
    }
    titre.y<-if (orient=="ver") {
	tclvalue(Env$l.var$titre.axevaleurs)
    } else {
	tclvalue(Env$l.var$titre.axenoms)
    }
    if (nchar(tclvalue(Env$l.var$soustitre))!=0) {
	title(main=tclvalue(Env$l.var$titre),col.main=tclvalue(Env$l.var$titre.col),cex.main=as.numeric(tclvalue(Env$l.var$titre.taille)),line=2.2)
	title(main=tclvalue(Env$l.var$soustitre),col.main=tclvalue(Env$l.var$titre.col),cex.main=0.7*as.numeric(tclvalue(Env$l.var$titre.taille)),line=0.9)
	title(xlab=titre.x,ylab=titre.y,cex.lab=as.numeric(tclvalue(Env$l.var$legendes.taille)),col.lab=tclvalue(Env$l.var$legendes.col))
    } else {
	title(main=tclvalue(Env$l.var$titre),col.main=tclvalue(Env$l.var$titre.col),cex.main=as.numeric(tclvalue(Env$l.var$titre.taille)),
	  xlab=titre.x,ylab=titre.y,cex.lab=as.numeric(tclvalue(Env$l.var$legendes.taille)),col.lab=tclvalue(Env$l.var$legendes.col))
    }
  } else {
    if (nchar(tclvalue(Env$l.var$soustitre))!=0) {
	title(main=tclvalue(Env$l.var$titre),col.main=tclvalue(Env$l.var$titre.col),cex.main=as.numeric(tclvalue(Env$l.var$titre.taille)),line=2.2)
	title(main=tclvalue(Env$l.var$soustitre),col.main=tclvalue(Env$l.var$titre.col),cex.main=0.7*as.numeric(tclvalue(Env$l.var$titre.taille)),line=0.9)
	title(xlab=tclvalue(Env$l.var$titre.axehor),ylab=tclvalue(Env$l.var$titre.axever),cex.lab=as.numeric(tclvalue(Env$l.var$legendes.taille)),
	  col.lab=tclvalue(Env$l.var$legendes.col))
    } else {
	title(main=tclvalue(Env$l.var$titre),col.main=tclvalue(Env$l.var$titre.col),cex.main=as.numeric(tclvalue(Env$l.var$titre.taille)),
	  xlab=tclvalue(Env$l.var$titre.axehor),ylab=tclvalue(Env$l.var$titre.axever),cex.lab=as.numeric(tclvalue(Env$l.var$legendes.taille)),
	  col.lab=tclvalue(Env$l.var$legendes.col))
    }
  }
  if (tclvalue(Env$l.var$sysinfo)==1) {
    mtext(paste("R version: ",getRversion()," - GrapheR version: ",packageVersion("GrapheR")," - Date: ",paste(strsplit(as.character(Sys.Date()),split="-")[[1]],
	collapse="."),sep=""),side=4,cex=0.6)
  }
}


#-------------------------------------------------
# Graphe - ajouter les axes
#-------------------------------------------------

graphe.axes<-function(type="",mids=NULL,longueur=NULL,orient=NULL,ordonnee=NULL) {
  parametres<-par()
  par(col.axis=tclvalue(Env$l.var$graduations.col),cex.axis=as.numeric(tclvalue(Env$l.var$graduations.taille)),
    las=ifelse(tclvalue(Env$l.var$graduations.orient)==Env$voc[246,1],0,1))
  if (type=="hist.freq") {
    axis(1,labels=mids,at=(0.5:(longueur-0.5)))
  } else if (type=="moust") {
    if (orient=="ver") {
	if (tclvalue(Env$l.var$plusieurs)==0) {
	  axis(1,labels=Env$l.var$noms1,at=1:length(Env$l.var$noms1))
	} else {
	  n.lab <- length(Env$l.var$noms1)
	  n.tot <- nlevels(Env$l.var$facteur.interaction)
	  n <- n.tot/n.lab
	  deb <- ((n.tot+1)-n*(n.lab-1))/2
	  mtext(Env$l.var$noms1,side=1,line=1,at=seq(deb,deb+n*(n.lab-1),n))
	}
    } else {
	axis(1)
    }
  } else if (type=="bar") {
    if (tclvalue(Env$l.var$encadre)==0) {abline(h=ordonnee)}
    if (tclvalue(Env$l.var$nobar)==1) {
	if (nrow(t(Env$l.var$add.abscisses))==1) {
	  mtext(text=if(tclvalue(Env$l.var$moyprop)=="moy"){Env$l.var$noms1}else{Env$l.var$nomsprop.fac},at=Env$l.var$add.abscisses,side=1,line=1.1,
	    cex=as.numeric(tclvalue(Env$l.var$legendes.taille)),col=tclvalue(Env$l.var$legendes.col))
	} else {
	  mtext(text=if(tclvalue(Env$l.var$moyprop)=="moy"){Env$l.var$noms1}else{Env$l.var$nomsprop.fac},at=apply(Env$l.var$add.abscisses,2,mean),side=1,
	    line=1.1,cex=as.numeric(tclvalue(Env$l.var$legendes.taille)),col=tclvalue(Env$l.var$legendes.col))
	}
    }
  } else {
    axis(1)
  }
  if (type=="moust") {
    if (orient=="ver") {
	axis(2)
    } else {
	if (tclvalue(Env$l.var$plusieurs)==0) {
	  axis(2,labels=Env$l.var$noms1,at=1:length(Env$l.var$noms1))
	} else {
	  n.lab <- length(Env$l.var$noms1)
	  n.tot <- nlevels(Env$l.var$facteur.interaction)
	  n <- n.tot/n.lab
	  deb <- ((n.tot+1)-n*(n.lab-1))/2
	  mtext(Env$l.var$noms1,side=2,line=1,at=seq(deb,deb+n*(n.lab-1),n))
	}
    }
  } else {
    axis(2)
  }
  par<-parametres
  par(las=0)
}


#-------------------------------------------------
# Graphe - échelles log
#-------------------------------------------------

graphe.log<-function() {
  result<-if (tclvalue(Env$l.var$log.axehor)==1) {
    if (tclvalue(Env$l.var$log.axever)==1) {
	"xy"
    } else {
	"x"
    }
  } else {
    if (tclvalue(Env$l.var$log.axever)==1) {
	"y"
    } else {
	""
    }
  }
  return(result)
}


#-------------------------------------------------
# Graphe - calculer les barres d'erreur
# (sauf courbe type proportion)
#-------------------------------------------------

graphe.erreurs.calculer<-function(variable,facteur1,facteur2=NULL,valeurs=NULL,prop.nvx=NULL) {
  erreur.inf<-NULL
  erreur.sup<-NULL
  if (tclvalue(Env$l.var$moyprop)=="moy") {
    if (tclvalue(Env$l.var$plusieurs)==0) {
	if (tclvalue(Env$l.var$erreur)==Env$voc[96,1]) {
	  erreur<-tapply(variable,facteur1,function(x) sd(x,na.rm=TRUE))
	  erreur.inf<-erreur.sup<-erreur
	} else if (tclvalue(Env$l.var$erreur)==Env$voc[97,1]) {
	  erreur<-tapply(variable,facteur1,function(x) sd(x,na.rm=TRUE)/sqrt(length(na.omit(x))))
	  erreur.inf<-erreur.sup<-erreur
	} else if (tclvalue(Env$l.var$erreur)==Env$voc[98,1]) {
	  moyennes<-tapply(variable,facteur1,function(x) mean(x,na.rm=TRUE))
	  erreur.inf<-erreur.sup<-moyennes-tapply(variable,facteur1,function(x) t.test(x)$conf.int[1])
	} else {
	  erreur.inf<-erreur.sup<-rep(0,nlevels(facteur1))
	}
    } else {
	if (tclvalue(Env$l.var$erreur)==Env$voc[96,1]) {
	  erreur<-tapply(variable,list(facteur2,facteur1),function(x) sd(x,na.rm=TRUE))
	  erreur.inf<-erreur.sup<-erreur
	} else if (tclvalue(Env$l.var$erreur)==Env$voc[97,1]) {
	  erreur<-tapply(variable,list(facteur2,facteur1),function(x) sd(x,na.rm=TRUE)/sqrt(length(na.omit(x))))
	  erreur.inf<-erreur.sup<-erreur
	} else if (tclvalue(Env$l.var$erreur)==Env$voc[98,1]) {
	  moyennes<-tapply(variable,list(facteur2,facteur1),function(x) mean(x,na.rm=TRUE))
	  erreur.inf<-erreur.sup<-moyennes-tapply(variable,list(facteur2,facteur1),function(x) t.test(x)$conf.int[1])
	} else {
	  erreur.inf<-erreur.sup<-matrix(0,ncol=nlevels(facteur1),nrow=nlevels(facteur2))
	}
    }
  } else {
    if (tclvalue(Env$l.var$plusieurs)==0) {
	if (tclvalue(Env$l.var$erreur)==Env$voc[97,1]) {
	  erreur<-NULL
	  for (i in 1:nlevels(facteur1)) {
	    erreur<-c(erreur,sqrt((valeurs[i]*(1-valeurs[i]))/(length(na.omit(variable[facteur1==levels(facteur1)[i]]))-1)))
	  }
	  erreur.inf<-erreur.sup<-erreur
	} else if (tclvalue(Env$l.var$erreur)==Env$voc[98,1]) {
	  for (i in 1:nlevels(facteur1)) {
	    erreur.inf<-c(erreur.inf,valeurs[i]-binom.test(length(na.omit(variable[variable==levels(variable)[prop.nvx] & facteur1==levels(facteur1)[i]])),length(na.omit(variable[facteur1==levels(facteur1)[i]])))$conf.int[1])
	    erreur.sup<-c(erreur.sup,binom.test(length(na.omit(variable[variable==levels(variable)[prop.nvx] & facteur1==levels(facteur1)[i]])),length(na.omit(variable[facteur1==levels(facteur1)[i]])))$conf.int[2]-valeurs[i])
	  }
	} else {
	  erreur.inf<-erreur.sup<-rep(0,nlevels(facteur1))
	}
    } else {
	if (tclvalue(Env$l.var$erreur)==Env$voc[97,1]) {
	  erreur<-matrix(0,nrow=nlevels(variable),ncol=nlevels(facteur1))
	  for (i in 1:nlevels(facteur1)) {
	    for (j in 1:length(prop.nvx)) {
		erreur[j,i]<-sqrt((valeurs[j,i]*(1-valeurs[j,i]))/(length(na.omit(variable[facteur1==levels(facteur1)[i]]))-1))
	    }
	  }
	  erreur.inf<-erreur.sup<-erreur
	} else if (tclvalue(Env$l.var$erreur)==Env$voc[98,1]) {
	  erreur.inf<-matrix(0,nrow=length(prop.nvx),ncol=nlevels(facteur1))
	  erreur.sup<-matrix(0,nrow=length(prop.nvx),ncol=nlevels(facteur1))
	  for (i in 1:nlevels(facteur1)) {
	    for (j in 1:length(prop.nvx)) {
		erreur.inf[j,i]<-valeurs[j,i]-binom.test(length(na.omit(variable[variable==levels(variable)[prop.nvx[j]] & facteur1==levels(facteur1)[i]])),
		  length(na.omit(variable[facteur1==levels(facteur1)[i]])))$conf.int[1]
		erreur.sup[j,i]<-binom.test(length(na.omit(variable[variable==levels(variable)[prop.nvx[j]] & facteur1==levels(facteur1)[i]])),
		  length(na.omit(variable[facteur1==levels(facteur1)[i]])))$conf.int[2]-valeurs[j,i]
	    }
	  }
	} else {
	  erreur.inf<-erreur.sup<-matrix(0,nrow=nlevels(variable),ncol=nlevels(facteur1))
	}
    }
  }
  return(list(erreur.inf=erreur.inf,erreur.sup=erreur.sup))
}


#-------------------------------------------------
# Graphe - ajouter les barres d'erreur
#-------------------------------------------------

graphe.erreurs.tracer<-function(abscisses,valeurs,erreur.inf,erreur.sup,couleur,amplitude=NULL) {
  if (tclvalue(Env$l.var$erreur)%in%Env$voc[96:98,1]) {
    if (is.null(amplitude)) {
	arrows(abscisses,valeurs-erreur.inf,abscisses,valeurs+erreur.sup,col=couleur,code=3,angle=90,length=0.1)
    } else {
	arrows(abscisses,valeurs-erreur.inf,abscisses,valeurs+erreur.sup,col=couleur,code=3,angle=90,length=0.015*amplitude)
    }
    msg(text="",type="info")
  }
}


#-------------------------------------------------
# Graphe - ajouter la légende
#-------------------------------------------------

graphe.legende<-function(type,symboles=NULL,lignes=NULL,traits=NULL) {
  par(xpd=TRUE)
  position<-if (tclvalue(Env$l.var$legende.pos)==Env$voc[102,1]) {"top"} else
    if (tclvalue(Env$l.var$legende.pos)==Env$voc[103,1]) {"topright"} else
    if (tclvalue(Env$l.var$legende.pos)==Env$voc[104,1]) {"left"} else
    if (tclvalue(Env$l.var$legende.pos)==Env$voc[105,1]) {"center"} else
    if (tclvalue(Env$l.var$legende.pos)==Env$voc[106,1]) {"right"} else
    if (tclvalue(Env$l.var$legende.pos)==Env$voc[107,1]) {"bottomleft"} else
    if (tclvalue(Env$l.var$legende.pos)==Env$voc[108,1]) {"bottom"} else
    if (tclvalue(Env$l.var$legende.pos)==Env$voc[109,1]) {"bottomright"} else
    {"topleft"}
  if (type=="moust") {
    if (nchar(tclvalue(Env$l.var$legende.titre))>0) {
	legend(position,legend=Env$l.var$noms2,fill=Env$l.var$couleur1B,title=tclvalue(Env$l.var$legende.titre))
    } else {
	legend(position,legend=Env$l.var$noms2,fill=Env$l.var$couleur1B)
    }
  }
  if (type=="bar") {
    if (nchar(tclvalue(Env$l.var$legende.titre))>0) {
	if (tclvalue(Env$l.var$moyprop)=="moy") {
	  if (tclvalue(Env$l.var$nobar)==0) {
	    legend(position,legend=Env$l.var$noms2,fill=Env$l.var$couleur1B,title=tclvalue(Env$l.var$legende.titre))
	  } else {
	    legend(position,legend=Env$l.var$noms2,pch=16,pt.cex=1.7,col=Env$l.var$couleur1B,title=tclvalue(Env$l.var$legende.titre))
	  }
	} else {
	  if (tclvalue(Env$l.var$nobar)==0) {
	    legend(position,legend=Env$l.var$nomsprop,fill=Env$l.var$couleur1B,title=tclvalue(Env$l.var$legende.titre))
	  } else {
	    legend(position,legend=Env$l.var$nomsprop,pch=16,pt.cex=1.7,col=Env$l.var$couleur1B,title=tclvalue(Env$l.var$legende.titre))
	  }
	}
    } else {
	if (tclvalue(Env$l.var$moyprop)=="moy") {
	  if (tclvalue(Env$l.var$nobar)==0) {
	    legend(position,legend=Env$l.var$noms2,fill=Env$l.var$couleur1B)
	  } else {
	    legend(position,legend=Env$l.var$noms2,pch=16,pt.cex=1.7,col=Env$l.var$couleur1B)
	  }
	} else {
	  if (tclvalue(Env$l.var$nobar)==0) {
	    legend(position,legend=Env$l.var$nomsprop,fill=Env$l.var$couleur1B)
	  } else {
	    legend(position,legend=Env$l.var$nomsprop,pch=16,pt.cex=1.7,col=Env$l.var$couleur1B)
	  }
	}
    }
  }
  if (type=="cam") {
    if (nchar(tclvalue(Env$l.var$legende.titre))>0) {
	legend(position,legend=Env$l.var$nomsparts,fill=Env$l.var$couleur1B,title=tclvalue(Env$l.var$legende.titre))
    } else {
	legend(position,legend=Env$l.var$nomsparts,fill=Env$l.var$couleur1B)
    }
  }
  if (type=="courbe") {
    epaisseur<-Env$l.var$epaisseur2
    taille.pts<-Env$l.var$taille.ptsB
    if (any(lignes=="p")) {
	traits[which(lignes=="p")]<-NA
	epaisseur[which(lignes=="p")]<-NA
    }
    if (any(lignes%in%c("l","h"))) {
	symboles[which(lignes%in%c("l","h"))]<-NA
	taille.pts[which(lignes%in%c("l","h"))]<-NA
    }
    if (nchar(tclvalue(Env$l.var$legende.titre))>0) {
	legend(position,legend=Env$l.var$noms1,col=Env$l.var$couleur2B,pch=symboles,pt.cex=taille.pts,
	  lty=traits,lwd=epaisseur,title=tclvalue(Env$l.var$legende.titre))
    } else {
	legend(position,legend=Env$l.var$noms1,col=Env$l.var$couleur2B,pch=symboles,pt.cex=taille.pts,
	  lty=traits,lwd=epaisseur)
    }
  }
  if (type=="nuage") {
    if (nchar(tclvalue(Env$l.var$legende.titre))>0) {
	legend(position,legend=Env$l.var$noms1,col=Env$l.var$couleur2B,pch=symboles,pt.cex=Env$l.var$taille.ptsB,
	  title=tclvalue(Env$l.var$legende.titre))
    } else {
	legend(position,legend=Env$l.var$noms1,col=Env$l.var$couleur2B,pch=symboles,pt.cex=Env$l.var$taille.ptsB)
    }
  }
  par(xpd=FALSE)
}


#-------------------------------------------------
# Graphe - ajouter un cadre
#-------------------------------------------------

graphe.box<-function() {
  if (tclvalue(Env$l.var$encadre)==1) {box()}
}


#-------------------------------------------------
# Graphe - type de trait
#-------------------------------------------------

type.trait<-function(type) {
  result<-NULL
  for (i in 1:length(type)) {
    result<-c(result,if (type[i]==Env$voc[61,1]) {3} else
    if (type[i]==Env$voc[62,1]) {2} else
    {1})
  }
  return(result)
}


#-------------------------------------------------
# Graphe - type de trait
#-------------------------------------------------

type.ligne<-function(type) {
  result<-NULL
  for (i in 1:length(type)) {
    result<-c(result,if (type[i]==Env$voc[133,1]) {"p"} else
    if (type[i]==Env$voc[134,1]) {"l"} else
    if (type[i]==Env$voc[136,1]) {"o"} else
    if (type[i]==Env$voc[137,1]) {"h"} else
    {"b"})
  }
  return(result)
}


#-------------------------------------------------
# Graphe - hachures
#-------------------------------------------------

graphe.hachures<-function(num) {
  densite<-NULL
  angle<-NULL
  if (tclvalue(Env$l.var$plusieurs)==0) {
    if (num==1) {densite<-0;angle<-0} else
    if (num==2) {densite<-4;angle<-135} else
    if (num==3) {densite<-4;angle<-90} else
    if (num==4) {densite<-4;angle<-45} else
    if (num==5) {densite<-4;angle<-0} else
    if (num==6) {densite<-15;angle<-135} else
    if (num==7) {densite<-15;angle<-90} else
    if (num==8) {densite<-15;angle<-45} else
    if (num==9) {densite<-15;angle<-0}
  } else {
    for (i in 1:length(num)) {
	if (num[i]==1) {densite<-c(densite,0);angle<-c(angle,0)} else
	if (num[i]==2) {densite<-c(densite,4);angle<-c(angle,135)} else
	if (num[i]==3) {densite<-c(densite,4);angle<-c(angle,90)} else
	if (num[i]==4) {densite<-c(densite,4);angle<-c(angle,45)} else
	if (num[i]==5) {densite<-c(densite,4);angle<-c(angle,0)} else
	if (num[i]==6) {densite<-c(densite,15);angle<-c(angle,135)} else
	if (num[i]==7) {densite<-c(densite,15);angle<-c(angle,90)} else
	if (num[i]==8) {densite<-c(densite,15);angle<-c(angle,45)} else
	if (num[i]==9) {densite<-c(densite,15);angle<-c(angle,0)}
    }
  }
  return(list(densite=densite,angle=angle))
}


#-------------------------------------------------
# Graphe - symboles
#-------------------------------------------------

graphe.symboles<-function(num) {
  symboles<-NULL
  if (tclvalue(Env$l.var$plusieurs)==0) {
    if (num==1) {symboles<-1} else
    if (num==2) {symboles<-0} else
    if (num==3) {symboles<-2} else
    if (num==4) {symboles<-3} else
    if (num==5) {symboles<-16} else
    if (num==6) {symboles<-15} else
    if (num==7) {symboles<-17} else
    if (num==8) {symboles<-4}
  } else {
    for (i in 1:length(num)) {
	if (num[i]==1) {symboles<-c(symboles,1)} else
	if (num[i]==2) {symboles<-c(symboles,0)} else
	if (num[i]==3) {symboles<-c(symboles,2)} else
	if (num[i]==4) {symboles<-c(symboles,3)} else
	if (num[i]==5) {symboles<-c(symboles,16)} else
	if (num[i]==6) {symboles<-c(symboles,15)} else
	if (num[i]==7) {symboles<-c(symboles,17)} else
	if (num[i]==8) {symboles<-c(symboles,4)}
    }
  }
  return(symboles)
}


#-------------------------------------------------
# Tracer le graphe
#-------------------------------------------------

tracer<-function() {
  if(!is.null(Env$dataset)) {
    if (Env$l.var$ecran=="H") {
	if (nchar(tclvalue(Env$l.var$variable))>0 & nchar(tclvalue(Env$l.var$hist.type))>0) {
	  par(mar=c(5,5,4,2),bg="white")
	  tracer.hist()
	  return(TRUE)
	} else {
	  msg(text=Env$voc[154,1],type="error")
	  return(FALSE)
	}
    } else
    if (Env$l.var$ecran=="M") {
	if (nchar(tclvalue(Env$l.var$variable))>0 & nchar(tclvalue(Env$l.var$facteur1))>0) {
	  par(mar=c(5,5,4,2),bg="white")
	  tracer.moust()
	  return(TRUE)
	} else {
	  msg(text=Env$voc[154,1],type="error")
	  return(FALSE)
	}
    } else
    if (Env$l.var$ecran=="B") {
	if (tclvalue(Env$l.var$moyprop)=="moy") {
	  if (nchar(tclvalue(Env$l.var$variable))>0 & nchar(tclvalue(Env$l.var$facteur1))>0) {
	    if (tclvalue(Env$l.var$plusieurs)==0) {
		par(mar=c(5,5,4,2),bg="white")
		tracer.barres.moyun()
	    } else {
		if (!any(table(Env$dataset[,tclvalue(Env$l.var$facteur1)],Env$dataset[,tclvalue(Env$l.var$facteur2)])==0)) {
		  par(mar=c(5,5,4,2),bg="white")
		  tracer.barres.moyplusieurs()
		} else {
		  msg(text=Env$voc[244,1],type="error")
		  return(FALSE)
		}
	    }
	    return(TRUE)
	  } else {
	    msg(text=Env$voc[154,1],type="error")
	    return(FALSE)
	  }
	} else {
	  if (nchar(tclvalue(Env$l.var$proportions))>0 & nchar(tclvalue(Env$l.var$prop.niveaux))>0 & nchar(tclvalue(Env$l.var$facteurprop))>0) {
	    if (tclvalue(Env$l.var$plusieurs)==0) {
		par(mar=c(5,5,4,2),bg="white")
		tracer.barres.propun()
	    } else {
		par(mar=c(5,5,4,2),bg="white")
		tracer.barres.propplusieurs()
	    }
	    return(TRUE)
	  } else {
	    msg(text=Env$voc[154,1],type="error")
	    return(FALSE)
	  }
	}
    } else
    if (Env$l.var$ecran=="Ca") {
	if (nchar(tclvalue(Env$l.var$variable))>0 & nchar(tclvalue(Env$l.var$parts.niveaux))>0) {
	  par(mar=c(5,5,4,2),bg="white")
	  tracer.cam()
	  return(TRUE)
	} else {
	  msg(text=Env$voc[154,1],type="error")
	  return(FALSE)
	}
    } else
    if (Env$l.var$ecran=="Co") {
	if (tclvalue(Env$l.var$moyprop)=="moy") {
	  if (nchar(tclvalue(Env$l.var$varX))>0 & nchar(tclvalue(Env$l.var$varY))>0) {
	    if (tclvalue(Env$l.var$plusieurs)==0) {
		if (nchar(tclvalue(Env$l.var$facteur1))>0 & tclvalue(Env$l.var$facteur1)!=Env$voc[82,1]) {
		  if (nchar(tclvalue(Env$l.var$niveau))>0) {
		    par(mar=c(5,5,4,2),bg="white")
		    tracer.courbe.moyun()
		    return(TRUE)
		  } else {
		    msg(text=Env$voc[154,1],type="error")
		    return(FALSE)
		  }
		} else {
		  par(mar=c(5,5,4,2),bg="white")
		  tracer.courbe.moyun()
		  return(TRUE)
		}
	    } else {
		par(mar=c(5,5,4,2),bg="white")
		tracer.courbe.moyplusieurs()
		return(TRUE)
	    }
	  } else {
	    msg(text=Env$voc[154,1],type="error")
	    return(FALSE)
	  }
	} else {
	  if (nchar(tclvalue(Env$l.var$varX.prop))>0 & nchar(tclvalue(Env$l.var$proportions))>0 & nchar(tclvalue(Env$l.var$prop.niveaux))>0) {
	    if (tclvalue(Env$l.var$plusieurs)==0) {
		if (nchar(tclvalue(Env$l.var$facteur1))>0 & tclvalue(Env$l.var$facteur1)!=Env$voc[82,1]) {
		  if (nchar(tclvalue(Env$l.var$niveau))>0) {
		    par(mar=c(5,5,4,2),bg="white")
		    tracer.courbe.propun()
		    return(TRUE)
		  } else {
		    msg(text=Env$voc[154,1],type="error")
		    return(FALSE)
		  }
		} else {
		  par(mar=c(5,5,4,2),bg="white")
		  tracer.courbe.propun()
		  return(TRUE)
		}
	    } else {
		par(mar=c(5,5,4,2),bg="white")
		tracer.courbe.propplusieurs()
		return(TRUE)
	    }
	  } else {
	    msg(text=Env$voc[154,1],type="error")
	    return(FALSE)
	  }
	}
    } else
    if (Env$l.var$ecran=="N") {
	if (nchar(tclvalue(Env$l.var$varX))>0 & nchar(tclvalue(Env$l.var$varY))>0) {
	  if (tclvalue(Env$l.var$plusieurs)==0) {
	    if (nchar(tclvalue(Env$l.var$facteur1))>0 & tclvalue(Env$l.var$facteur1)!=Env$voc[82,1]) {
		if (nchar(tclvalue(Env$l.var$niveau))>0) {
		  par(mar=c(5,5,4,2),bg="white")
		  tracer.nuage.un()
		  return(TRUE)
		} else {
		  msg(text=Env$voc[154,1],type="error")
		  return(FALSE)
		}
	    } else {
		par(mar=c(5,5,4,2),bg="white")
		tracer.nuage.un()
		return(TRUE)
	    }
	  } else {
	    par(mar=c(5,5,4,2),bg="white")
	    tracer.nuage.plusieurs()
	    return(TRUE)
	  }
	} else {
	  msg(text=Env$voc[154,1],type="error")
	  return(FALSE)
	}
    } else {
	return(FALSE)
    }
  } else {
    msg(text=Env$voc[153,1],type="error")
    return(FALSE)
  }
}


#-------------------------------------------------
# Tracer l'histogramme - limites des axes
#-------------------------------------------------

tracer.hist.limites<-function(variable,type,frequence=NULL) {
  x.inf<-if (type=="freq") {
    0
  } else if (type=="eff" | type=="dens") {
    if (tclvalue(Env$l.var$liminf.axehor)=="Auto") {
	min(hist(variable,breaks=ifelse(tclvalue(Env$l.var$hist.barres)=="Auto","Sturges",as.numeric(tclvalue(Env$l.var$hist.barres))-1),plot=FALSE)$breaks)
    } else {
	as.numeric(tclvalue(Env$l.var$liminf.axehor))
    }
  }
  x.sup<-if (type=="freq") {
    length(hist(variable,breaks=ifelse(tclvalue(Env$l.var$hist.barres)=="Auto","Sturges",as.numeric(tclvalue(Env$l.var$hist.barres))-1),plot=FALSE)$counts)
  } else if (type=="eff" | type=="dens") {
    if (tclvalue(Env$l.var$limsup.axehor)=="Auto") {
	max(hist(variable,breaks=ifelse(tclvalue(Env$l.var$hist.barres)=="Auto","Sturges",as.numeric(tclvalue(Env$l.var$hist.barres))-1),plot=FALSE)$breaks)
    } else {
	as.numeric(tclvalue(Env$l.var$limsup.axehor))
    }
  }
  y.sup<-if (type=="freq") {
    if (tclvalue(Env$l.var$limsup.axever)=="Auto") {
	1.1*max(frequence)
    } else {
	as.numeric(tclvalue(Env$l.var$limsup.axever))
    }
  } else if (type=="eff") {
    if (tclvalue(Env$l.var$limsup.axever)=="Auto") {
	max(hist(variable,breaks=ifelse(tclvalue(Env$l.var$hist.barres)=="Auto","Sturges",as.numeric(tclvalue(Env$l.var$hist.barres))-1),plot=FALSE)$counts)
    } else {
	as.numeric(tclvalue(Env$l.var$limsup.axever))
    }
  } else if (type=="dens") {
    if (tclvalue(Env$l.var$limsup.axever)=="Auto") {
	max(hist(variable,breaks=ifelse(tclvalue(Env$l.var$hist.barres)=="Auto","Sturges",as.numeric(tclvalue(Env$l.var$hist.barres))-1),plot=FALSE)$density)
    } else {
	as.numeric(tclvalue(Env$l.var$limsup.axever))
    }
  }
  return(list(xinf=x.inf,xsup=x.sup,ysup=y.sup))
}


#-------------------------------------------------
# Tracer l'histogramme
#-------------------------------------------------

tracer.hist<-function() {
  variable<-if (nchar(tclvalue(Env$l.var$facteur1))>0 & tclvalue(Env$l.var$facteur1)!=Env$voc[82,1]) {
    Env$dataset[,tclvalue(Env$l.var$variable)][Env$dataset[,tclvalue(Env$l.var$facteur1)]==tclvalue(Env$l.var$niveau)]
  } else {
    Env$dataset[,tclvalue(Env$l.var$variable)]
  }
  if (tclvalue(Env$l.var$hist.type)==Env$voc[40,1]) {
    total<-sum(hist(variable,breaks=ifelse(tclvalue(Env$l.var$hist.barres)=="Auto","Sturges",as.numeric(tclvalue(Env$l.var$hist.barres))-1),
	plot=FALSE)$counts)
    frequence<-hist(variable,breaks=ifelse(tclvalue(Env$l.var$hist.barres)=="Auto","Sturges",as.numeric(tclvalue(Env$l.var$hist.barres))-1),
	plot=FALSE)$counts/total
    limites<-tracer.hist.limites(variable=variable,type="freq",frequence=frequence)
    Env$l.code$x.inf<-x.inf<-limites$xinf
    Env$l.code$x.sup<-x.sup<-limites$xsup
    Env$l.code$y.sup<-y.sup<-limites$ysup
    barplot(frequence,axes=FALSE,ann=FALSE,space=0,col=tclvalue(Env$l.var$couleur1A),border=tclvalue(Env$l.var$col.borduresA),
	xlim=c(x.inf,x.sup),ylim=c(0,y.sup))
    graphe.axes(type="hist.freq",mids=hist(variable,breaks=ifelse(tclvalue(Env$l.var$hist.barres)=="Auto","Sturges",as.numeric(tclvalue(Env$l.var$hist.barres))-1),
	plot=FALSE)$mids,longueur=length(frequence))
  } else if (tclvalue(Env$l.var$hist.type)==Env$voc[41,1]) {
    limites<-tracer.hist.limites(variable=variable,type="eff")
    Env$l.code$x.inf<-x.inf<-limites$xinf
    Env$l.code$x.sup<-x.sup<-limites$xsup
    Env$l.code$y.sup<-y.sup<-limites$ysup
    hist(variable,axes=FALSE,ann=FALSE,freq=TRUE,col=tclvalue(Env$l.var$couleur1A),border=tclvalue(Env$l.var$col.borduresA),
	breaks=ifelse(tclvalue(Env$l.var$hist.barres)=="Auto","Sturges",as.numeric(tclvalue(Env$l.var$hist.barres))-1),
	xlim=c(x.inf,x.sup),ylim=c(0,y.sup))
    graphe.axes()
  } else if (tclvalue(Env$l.var$hist.type)==Env$voc[42,1]) {
    Env$l.var$add.seq<-seq(min(variable,na.rm=TRUE),max(variable,na.rm=TRUE),abs(max(variable,na.rm=TRUE)-min(variable,na.rm=TRUE))/1000)
    Env$l.var$add.seq2<-floor(min(variable,na.rm=TRUE)):ceiling(max(variable,na.rm=TRUE))
    limites<-tracer.hist.limites(variable=variable,type="dens")
    Env$l.code$x.inf<-x.inf<-limites$xinf
    Env$l.code$x.sup<-x.sup<-limites$xsup
    Env$l.code$y.sup<-y.sup<-limites$ysup
    hist(variable,axes=FALSE,ann=FALSE,freq=FALSE,col=tclvalue(Env$l.var$couleur1A),border=tclvalue(Env$l.var$col.borduresA),
      breaks=ifelse(tclvalue(Env$l.var$hist.barres)=="Auto","Sturges",as.numeric(tclvalue(Env$l.var$hist.barres))-1),
	xlim=c(x.inf,x.sup),ylim=c(0,y.sup))
    graphe.axes()
    if (tclvalue(Env$l.var$hist.dens)=="1") {
	lines(density(na.omit(variable)),col=tclvalue(Env$l.var$couleur2A),lwd=as.numeric(tclvalue(Env$l.var$epaisseur1)),lty=type.trait(type=tclvalue(Env$l.var$trait1)))
    }
  }
  graphe.titre()
  graphe.box()
}


#-------------------------------------------------
# Tracer les boîtes à moustaches
#-------------------------------------------------

tracer.moust<-function() {
  variable<-Env$dataset[,tclvalue(Env$l.var$variable)]
  facteur<-if (nchar(tclvalue(Env$l.var$facteur2))>0) {
    if (tclvalue(Env$l.var$facteur2)!=Env$voc[82,1]) {
	Env$l.var$facteur.interaction
    } else {
	Env$dataset[,tclvalue(Env$l.var$facteur1)]
    }
  } else {
    Env$dataset[,tclvalue(Env$l.var$facteur1)]
  }
  orient<-ifelse (tclvalue(Env$l.var$box.orient)==Env$voc[67,1],"hor","ver")
  log.axes<-if (orient=="ver" & tclvalue(Env$l.var$log.axevaleurs)==1) {
    "y"
  } else if (orient=="hor" & tclvalue(Env$l.var$log.axevaleurs)==1) {
    "x"
  } else {""}
  Env$l.code$y.inf<-y.inf<-if (tclvalue(Env$l.var$liminf.axevaleurs)=="Auto") {
    min(na.omit(variable))
  } else {
    as.numeric(tclvalue(Env$l.var$liminf.axevaleurs))
  }
  Env$l.code$y.sup<-y.sup<-if (tclvalue(Env$l.var$limsup.axevaleurs)=="Auto") {
    max(na.omit(variable))
  } else {
    as.numeric(tclvalue(Env$l.var$limsup.axevaleurs))
  }
  boxplot(variable~facteur,axes=FALSE,ann=FALSE,horizontal=ifelse(orient=="hor",TRUE,FALSE),
    col=if(tclvalue(Env$l.var$plusieurs)==0){tclvalue(Env$l.var$couleur1A)}else{Env$l.var$couleur1B},
    boxcol=if(tclvalue(Env$l.var$plusieurs)==0){tclvalue(Env$l.var$col.borduresA)}else{Env$l.var$col.borduresB},
    medcol=if(tclvalue(Env$l.var$plusieurs)==0){tclvalue(Env$l.var$col.borduresA)}else{Env$l.var$col.borduresB},
    whiskcol=if(tclvalue(Env$l.var$plusieurs)==0){tclvalue(Env$l.var$col.borduresA)}else{Env$l.var$col.borduresB},
    staplecol=if(tclvalue(Env$l.var$plusieurs)==0){tclvalue(Env$l.var$col.borduresA)}else{Env$l.var$col.borduresB},
    outcol=if(tclvalue(Env$l.var$plusieurs)==0){tclvalue(Env$l.var$col.borduresA)}else{Env$l.var$col.borduresB},
    whisklty=type.trait(type=tclvalue(Env$l.var$trait1)),
    range=as.numeric(tclvalue(Env$l.var$lg.moustaches)),notch=as.logical(as.numeric(tclvalue(Env$l.var$ICmediane))),log=log.axes,
    ylim=c(y.inf,y.sup),names=Env$l.var$noms1,varwidth=ifelse(tclvalue(Env$l.var$varwidth)==1,TRUE,FALSE))
  if (tclvalue(Env$l.var$boxmoy)==1) {
    if (tclvalue(Env$l.var$box.orient)==Env$voc[67,1]) {
	points(tapply(variable,facteur,function(x) mean(x,na.rm=TRUE)),1:nlevels(facteur),cex=2,col=Env$l.var$col.borduresB,pch="+")
    } else {
	points(1:nlevels(facteur),tapply(variable,facteur,function(x) mean(x,na.rm=TRUE)),cex=2,col=Env$l.var$col.borduresB,pch="+")
    }
  }
  graphe.titre(type="moust",orient=orient)
  graphe.axes(type="moust",orient=orient)
  graphe.box()
  if (tclvalue(Env$l.var$plusieurs)==1 & tclvalue(Env$l.var$legende)==1) {
    graphe.legende(type="moust")
  }
}


#-------------------------------------------------
# Tracer les barres - limites des axes
#-------------------------------------------------

tracer.barres.limites<-function(valeurs,erreur.inf,erreur.sup) {
  y.inf<-if(tclvalue(Env$l.var$liminf.axever)=="Auto") {
    if (any(valeurs-erreur.inf<0)) {
	if (tclvalue(Env$l.var$plusieurs)==1) {
	  if (tclvalue(Env$l.var$stack)==1) {
	    1.2*min(colSums(valeurs))
	  } else {
	    1.2*min(valeurs-erreur.inf)
	  }
	} else {
	  1.2*min(valeurs-erreur.inf)
	}
    } else {
	if (tclvalue(Env$l.var$log.axever)==1) {
	  0.01
	} else {
	  0
	}
    }
  } else {
    as.numeric(tclvalue(Env$l.var$liminf.axever))
  }
  y.sup<-if(tclvalue(Env$l.var$limsup.axever)=="Auto") {
    if (any(valeurs+erreur.sup>0)) {
	if (tclvalue(Env$l.var$plusieurs)==1) {
	  if (tclvalue(Env$l.var$stack)==1) {
	    1.2*max(colSums(valeurs))
	  } else {
	    1.2*max(valeurs+erreur.sup)
	  }
	} else {
	  1.2*max(valeurs+erreur.sup)
	}
    } else {
	if (tclvalue(Env$l.var$log.axever)==1) {
	  -0.01
	} else {
	  0
	}
    }
  } else {
    as.numeric(tclvalue(Env$l.var$limsup.axever))
  }
  ordonnee<-if (y.inf>=0 & y.sup>0) {
    y.inf
  } else if (y.inf<0 & y.sup>0){
    0
  } else if (y.inf<0 & y.sup<=0) {
    y.sup
  }
  return(list(yinf=y.inf,ysup=y.sup,ordonnee=ordonnee))
}


#-------------------------------------------------
# Tracer les barres - type moyenne, plusieurs=0
#-------------------------------------------------

tracer.barres.moyun<-function() {
  variable<-Env$dataset[,tclvalue(Env$l.var$variable)]
  facteur<-Env$dataset[,tclvalue(Env$l.var$facteur1)]
  valeurs<-tapply(variable,facteur,function(x) mean(x,na.rm=TRUE))
  erreurs<-graphe.erreurs.calculer(variable=variable,facteur1=facteur)
  erreurs$erreur.inf[is.na(erreurs$erreur.inf)] <- 0
  erreurs$erreur.sup[is.na(erreurs$erreur.sup)] <- 0
  limites<-tracer.barres.limites(valeurs=valeurs,erreur.inf=erreurs$erreur.inf,erreur.sup=erreurs$erreur.sup)
  Env$l.var$add.hauteurs<-valeurs+erreurs$erreur.sup
  Env$l.code$y.inf<-y.inf<-limites$yinf
  Env$l.code$y.sup<-y.sup<-limites$ysup
  if (tclvalue(Env$l.var$nobar)==0) {
    Env$l.var$add.abscisses<-barplot(valeurs,axes=FALSE,ann=FALSE,col=tclvalue(Env$l.var$couleur1A),log=graphe.log(),
	border=tclvalue(Env$l.var$col.borduresA),ylim=c(y.inf,y.sup),names.arg=Env$l.var$noms1,xpd=FALSE)
  } else {
    Env$l.var$add.abscisses<-barplot(valeurs,log=graphe.log(),ylim=c(y.inf,y.sup),plot=FALSE)
    plot(Env$l.var$add.abscisses,valeurs,cex=1.7,pch=16,col=tclvalue(Env$l.var$couleur1A),xlim=c(min(Env$l.var$add.abscisses)-0.5,max(Env$l.var$add.abscisses)+0.5),
	ylim=c(y.inf,y.sup),axes=FALSE,ann=FALSE,log=graphe.log())
  }
  Env$l.var$add.matrice<-matrix(numeric(length(Env$l.var$add.abscisses)^2),nrow=length(Env$l.var$add.abscisses),
    dimnames=list(1:length(Env$l.var$add.abscisses),1:length(Env$l.var$add.abscisses)))
  for (i in 1:length(Env$l.var$add.abscisses)) {
    for (j in 1:length(Env$l.var$add.abscisses)) {
	Env$l.var$add.matrice[j,i]<-max(Env$l.var$add.hauteurs[i:j])
    }
  }
  if (tclvalue(Env$l.var$nobar)==0 & graphe.log()=="" & tclvalue(Env$l.var$hachuresA)!="1") {
    hachures<-graphe.hachures(num=as.numeric(tclvalue(Env$l.var$hachuresA)))
    barplot(valeurs,axes=FALSE,ann=FALSE,col=tclvalue(Env$l.var$col.borduresA),border=tclvalue(Env$l.var$col.borduresA),
	log=graphe.log(),ylim=c(y.inf,y.sup),density=hachures$densite,angle=hachures$angle,names.arg="",xpd=FALSE,add=TRUE)
  }
  if (nchar(tclvalue(Env$l.var$erreur))>0 & tclvalue(Env$l.var$erreur)!=Env$voc[95,1]) {
    graphe.erreurs.tracer(abscisses=Env$l.var$add.abscisses,valeurs=valeurs,erreur.inf=erreurs$erreur.inf,
	erreur.sup=erreurs$erreur.sup,couleur=tclvalue(Env$l.var$couleur2A))
    if (tclvalue(Env$l.var$nobar)==1) {
	points(Env$l.var$add.abscisses,valeurs,cex=1.7,pch=16,col=tclvalue(Env$l.var$couleur1A))
    }
  }
  graphe.titre()
  graphe.axes(type="bar",ordonnee=limites$ordonnee)
  graphe.box()
}


#-------------------------------------------------
# Tracer les barres - type moyenne, plusieurs=1
#-------------------------------------------------

tracer.barres.moyplusieurs<-function() {
  variable<-Env$dataset[,tclvalue(Env$l.var$variable)]
  facteur1<-Env$dataset[,tclvalue(Env$l.var$facteur1)]
  facteur2<-Env$dataset[,tclvalue(Env$l.var$facteur2)]
  valeurs<-tapply(variable,list(facteur2,facteur1),function(x) mean(x,na.rm=TRUE))
  erreurs<-graphe.erreurs.calculer(variable=variable,facteur1=facteur1,facteur2=facteur2)
  erreurs$erreur.inf[is.na(erreurs$erreur.inf)] <- 0
  erreurs$erreur.sup[is.na(erreurs$erreur.sup)] <- 0
  limites<-tracer.barres.limites(valeurs=valeurs,erreur.inf=erreurs$erreur.inf,erreur.sup=erreurs$erreur.sup)
  Env$l.var$add.hauteurs<-valeurs+erreurs$erreur.sup
  Env$l.code$y.inf<-y.inf<-limites$yinf
  Env$l.code$y.sup<-y.sup<-limites$ysup
  if (tclvalue(Env$l.var$nobar)==0) {
    Env$l.var$add.abscisses<-barplot(valeurs,axes=FALSE,ann=FALSE,col=Env$l.var$couleur1B,log=graphe.log(),
	border=Env$l.var$col.borduresB,ylim=c(y.inf,y.sup),names.arg=Env$l.var$noms1,
	beside=ifelse(tclvalue(Env$l.var$stack)==1,FALSE,TRUE),xpd=FALSE)
  } else {
    Env$l.var$add.abscisses<-barplot(valeurs,log=graphe.log(),ylim=c(y.inf,y.sup),beside=ifelse(tclvalue(Env$l.var$stack)==1,FALSE,TRUE),plot=FALSE)
    plot(if(tclvalue(Env$l.var$stack)==1){rep(Env$l.var$add.abscisses,each=length(Env$l.var$add.abscisses))}else{Env$l.var$add.abscisses},valeurs,cex=1.7,pch=16,
	col=if(tclvalue(Env$l.var$stack)==1){Env$l.var$couleur1B}else{rep(Env$l.var$couleur1B,ncol(Env$l.var$add.abscisses))},xlim=c(min(Env$l.var$add.abscisses)-0.5,
	max(Env$l.var$add.abscisses)+0.5),ylim=c(y.inf,y.sup),axes=FALSE,ann=FALSE,log=graphe.log())
  }
  Env$l.var$add.matrice<-matrix(numeric(length(Env$l.var$add.abscisses)^2),nrow=length(Env$l.var$add.abscisses),
    dimnames=list(1:length(Env$l.var$add.abscisses),1:length(Env$l.var$add.abscisses)))
  for (i in 1:length(Env$l.var$add.abscisses)) {
    for (j in 1:length(Env$l.var$add.abscisses)) {
	Env$l.var$add.matrice[j,i]<-max(Env$l.var$add.hauteurs[i:j])
    }
  }
  if (tclvalue(Env$l.var$nobar)==0 & graphe.log()=="" & any(Env$l.var$hachuresB!=1)) {
    hachures<-graphe.hachures(num=Env$l.var$hachuresB)
    barplot(valeurs,axes=FALSE,ann=FALSE,col=Env$l.var$col.borduresB,border=Env$l.var$col.borduresB,
	log=graphe.log(),ylim=c(y.inf,y.sup),density=hachures$densite,angle=hachures$angle,
	beside=ifelse(tclvalue(Env$l.var$stack)==1,FALSE,TRUE),names.arg=rep("",nlevels(facteur1)),xpd=FALSE,add=TRUE)
  }
  if (tclvalue(Env$l.var$stack)==0 & nchar(tclvalue(Env$l.var$erreur))>0 & tclvalue(Env$l.var$erreur)!=Env$voc[95,1]) {
    graphe.erreurs.tracer(abscisses=Env$l.var$add.abscisses,valeurs=valeurs,erreur.inf=erreurs$erreur.inf,
	erreur.sup=erreurs$erreur.sup,couleur=tclvalue(Env$l.var$couleur2A))
    if (tclvalue(Env$l.var$nobar)==1) {
	points(if(tclvalue(Env$l.var$stack)==1){rep(Env$l.var$add.abscisses,each=length(Env$l.var$add.abscisses))}else{Env$l.var$add.abscisses},valeurs,cex=1.7,pch=16,
	  col=if(tclvalue(Env$l.var$stack)==1){Env$l.var$couleur1B}else{rep(Env$l.var$couleur1B,ncol(Env$l.var$add.abscisses))})
    }
  }
  graphe.titre()
  graphe.axes(type="bar",ordonnee=limites$ordonnee)
  graphe.box()
  if (tclvalue(Env$l.var$legende)==1) {
    graphe.legende(type="bar")
  }
}


#-------------------------------------------------
# Tracer les barres - type proportion, plusieurs=0
#-------------------------------------------------

tracer.barres.propun<-function() {
  niveau<-as.numeric(tclvalue(Env$l.var$prop.niveaux))+1
  variable<-Env$dataset[,tclvalue(Env$l.var$proportions)]
  facteur<-Env$dataset[,tclvalue(Env$l.var$facteurprop)]
  valeurs<-matrix(0,nrow=nlevels(variable),ncol=nlevels(facteur))
  for (i in 1:nlevels(facteur)) {
    for (j in 1:nlevels(variable)) {
	valeurs[j,i]<-length(na.omit(variable[variable==levels(variable)[j] & facteur==levels(facteur)[i]]))/length(na.omit(variable[facteur==levels(facteur)[i]]))
    }
  }
  valeurs<-valeurs[niveau,]
  erreurs<-graphe.erreurs.calculer(variable=variable,facteur1=facteur,valeurs=valeurs,prop.nvx=niveau)
  erreurs$erreur.inf[is.na(erreurs$erreur.inf)] <- 0
  erreurs$erreur.sup[is.na(erreurs$erreur.sup)] <- 0
  limites<-tracer.barres.limites(valeurs=valeurs,erreur.inf=erreurs$erreur.inf,erreur.sup=erreurs$erreur.sup)
  Env$l.var$add.hauteurs<-valeurs+erreurs$erreur.sup
  Env$l.code$y.inf<-y.inf<-limites$yinf
  Env$l.code$y.sup<-y.sup<-limites$ysup
  if (tclvalue(Env$l.var$nobar)==0) {
    Env$l.var$add.abscisses<-barplot(valeurs,axes=FALSE,ann=FALSE,col=tclvalue(Env$l.var$couleur1A),log=graphe.log(),
	border=tclvalue(Env$l.var$col.borduresA),ylim=c(y.inf,y.sup),names.arg=Env$l.var$nomsprop.fac,xpd=FALSE)
  } else {
    Env$l.var$add.abscisses<-barplot(valeurs,log=graphe.log(),ylim=c(y.inf,y.sup),plot=FALSE)
    plot(Env$l.var$add.abscisses,valeurs,cex=1.7,pch=16,col=tclvalue(Env$l.var$couleur1A),xlim=c(min(Env$l.var$add.abscisses)-0.5,max(Env$l.var$add.abscisses)+0.5),
	ylim=c(y.inf,y.sup),axes=FALSE,ann=FALSE,log=graphe.log())
  }
  Env$l.var$add.matrice<-matrix(numeric(length(Env$l.var$add.abscisses)^2),nrow=length(Env$l.var$add.abscisses),
    dimnames=list(1:length(Env$l.var$add.abscisses),1:length(Env$l.var$add.abscisses)))
  for (i in 1:length(Env$l.var$add.abscisses)) {
    for (j in 1:length(Env$l.var$add.abscisses)) {
	Env$l.var$add.matrice[j,i]<-max(Env$l.var$add.hauteurs[i:j])
    }
  }
  if (tclvalue(Env$l.var$nobar)==0 & graphe.log()=="" & tclvalue(Env$l.var$hachuresA)!="1") {
    hachures<-graphe.hachures(num=as.numeric(tclvalue(Env$l.var$hachuresA)))
    barplot(valeurs,axes=FALSE,ann=FALSE,col=tclvalue(Env$l.var$col.borduresA),border=tclvalue(Env$l.var$col.borduresA),
	log=graphe.log(),ylim=c(y.inf,y.sup),density=hachures$densite,angle=hachures$angle,names.arg="",xpd=FALSE,add=TRUE)
  }
  if (nchar(tclvalue(Env$l.var$erreur))>0 & tclvalue(Env$l.var$erreur)!=Env$voc[95,1]) {
    graphe.erreurs.tracer(abscisses=Env$l.var$add.abscisses,valeurs=valeurs,erreur.inf=erreurs$erreur.inf,
	erreur.sup=erreurs$erreur.sup,couleur=tclvalue(Env$l.var$couleur2A))
    if (tclvalue(Env$l.var$nobar)==1) {
	points(Env$l.var$add.abscisses,valeurs,cex=1.7,pch=16,col=tclvalue(Env$l.var$couleur1A))
    }
  }
  graphe.titre()
  graphe.axes(type="bar",ordonnee=limites$ordonnee)
  graphe.box()
}


#-------------------------------------------------
# Tracer les barres - type proportion, plusieurs=1
#-------------------------------------------------

tracer.barres.propplusieurs<-function() {
  niveaux<-as.numeric(strsplit(tclvalue(Env$l.var$prop.niveaux),split=" ")[[1]])+1
  variable<-Env$dataset[,tclvalue(Env$l.var$proportions)]
  facteur<-Env$dataset[,tclvalue(Env$l.var$facteurprop)]
  valeurs<-matrix(0,nrow=nlevels(variable),ncol=nlevels(facteur))
  for (i in 1:nlevels(facteur)) {
    for (j in 1:nlevels(variable)) {
	valeurs[j,i]<-length(na.omit(variable[variable==levels(variable)[j] & facteur==levels(facteur)[i]]))/length(na.omit(variable[facteur==levels(facteur)[i]]))
    }
  }
  valeurs-valeurs[niveaux,]
  erreurs<-graphe.erreurs.calculer(variable=variable,facteur1=facteur,valeurs=valeurs,prop.nvx=niveaux)
  erreurs$erreur.inf[is.na(erreurs$erreur.inf)] <- 0
  erreurs$erreur.sup[is.na(erreurs$erreur.sup)] <- 0
  limites<-tracer.barres.limites(valeurs=valeurs,erreur.inf=erreurs$erreur.inf,erreur.sup=erreurs$erreur.sup)
  Env$l.var$add.hauteurs<-valeurs+erreurs$erreur.sup
  Env$l.code$y.inf<-y.inf<-limites$yinf
  Env$l.code$y.sup<-y.sup<-limites$ysup
  if (tclvalue(Env$l.var$nobar)==0) {
    Env$l.var$add.abscisses<-barplot(valeurs,axes=FALSE,ann=FALSE,col=Env$l.var$couleur1B,log=graphe.log(),
	border=Env$l.var$col.borduresB,ylim=c(y.inf,y.sup),names.arg=Env$l.var$nomsprop.fac,
	beside=ifelse(tclvalue(Env$l.var$stack)==1,FALSE,TRUE),xpd=FALSE)
  } else {
    Env$l.var$add.abscisses<-barplot(valeurs,log=graphe.log(),ylim=c(y.inf,y.sup),beside=ifelse(tclvalue(Env$l.var$stack)==1,FALSE,TRUE),
	plot=FALSE)
    plot(if(tclvalue(Env$l.var$stack)==1){rep(Env$l.var$add.abscisses,each=length(Env$l.var$add.abscisses))}else{Env$l.var$add.abscisses},valeurs,cex=1.7,pch=16,
	col=if(tclvalue(Env$l.var$stack)==1){Env$l.var$couleur1B}else{rep(Env$l.var$couleur1B,ncol(Env$l.var$add.abscisses))},xlim=c(min(Env$l.var$add.abscisses)-0.5,
	max(Env$l.var$add.abscisses)+0.5),ylim=c(y.inf,y.sup),axes=FALSE,ann=FALSE,log=graphe.log())
  }
  Env$l.var$add.matrice<-matrix(numeric(length(Env$l.var$add.abscisses)^2),nrow=length(Env$l.var$add.abscisses),
    dimnames=list(1:length(Env$l.var$add.abscisses),1:length(Env$l.var$add.abscisses)))
  for (i in 1:length(Env$l.var$add.abscisses)) {
    for (j in 1:length(Env$l.var$add.abscisses)) {
	Env$l.var$add.matrice[j,i]<-max(Env$l.var$add.hauteurs[i:j])
    }
  }
  if (tclvalue(Env$l.var$nobar)==0 & graphe.log()=="" & any(Env$l.var$hachuresB!=1)) {
    hachures<-graphe.hachures(num=Env$l.var$hachuresB)
    barplot(valeurs,axes=FALSE,ann=FALSE,col=Env$l.var$col.borduresB,border=Env$l.var$col.borduresB,
	log=graphe.log(),ylim=c(y.inf,y.sup),density=hachures$densite,angle=hachures$angle,
	beside=ifelse(tclvalue(Env$l.var$stack)==1,FALSE,TRUE),names.arg=rep("",nlevels(facteur)),xpd=FALSE,add=TRUE)
  }
  if (tclvalue(Env$l.var$stack)==0 & nchar(tclvalue(Env$l.var$erreur))>0 & tclvalue(Env$l.var$erreur)!=Env$voc[95,1]) {
    graphe.erreurs.tracer(abscisses=Env$l.var$add.abscisses,valeurs=valeurs,erreur.inf=erreurs$erreur.inf,
	erreur.sup=erreurs$erreur.sup,couleur=tclvalue(Env$l.var$couleur2A))
    if (tclvalue(Env$l.var$nobar)==1) {
	points(if(tclvalue(Env$l.var$stack)==1){rep(Env$l.var$add.abscisses,each=length(Env$l.var$add.abscisses))}else{Env$l.var$add.abscisses},valeurs,cex=1.7,pch=16,
	  col=if(tclvalue(Env$l.var$stack)==1){Env$l.var$couleur1B}else{rep(Env$l.var$couleur1B,ncol(Env$l.var$add.abscisses))})
    }
  }
  graphe.titre()
  graphe.axes(type="bar",ordonnee=limites$ordonnee)
  graphe.box()
  if (tclvalue(Env$l.var$legende)==1) {
    graphe.legende(type="bar")
  }
}


#-------------------------------------------------
# Tracer le camembert
#-------------------------------------------------

tracer.cam<-function() {
  niveaux<-as.numeric(strsplit(tclvalue(Env$l.var$parts.niveaux),split=" ")[[1]])+1
  prevariable1<-factor(Env$dataset[,tclvalue(Env$l.var$variable)])
  prevariable2<-prevariable1[prevariable1%in%levels(prevariable1)[niveaux]]
  variable<-summary(prevariable2)
  pie(variable,col=Env$l.var$couleur1B,border=tclvalue(Env$l.var$col.borduresA),
    labels=if (tclvalue(Env$l.var$cam.lien)==1) {Env$l.var$nomsparts} else {NA},
    clockwise=ifelse(tclvalue(Env$l.var$cam.orient)==Env$voc[120,1],TRUE,FALSE),
    init.angle=ifelse(tclvalue(Env$l.var$cam.orient)==Env$voc[120,1],as.numeric(tclvalue(Env$l.var$cam.start))+90,
    -1*as.numeric(tclvalue(Env$l.var$cam.start))+90))
  if (any(Env$l.var$hachuresB!=1)) {
    hachures<-graphe.hachures(num=Env$l.var$hachuresB)
    par(new=TRUE)
    pie(variable,col=tclvalue(Env$l.var$col.borduresA),labels=NA,density=hachures$densite,angle=hachures$angle,
	clockwise=ifelse(tclvalue(Env$l.var$cam.orient)==Env$voc[120,1],TRUE,FALSE),
	init.angle=ifelse(tclvalue(Env$l.var$cam.orient)==Env$voc[120,1],as.numeric(tclvalue(Env$l.var$cam.start))+90,
	-1*as.numeric(tclvalue(Env$l.var$cam.start))+90))
  }
  graphe.titre()
  if (tclvalue(Env$l.var$cam.lien)==0 & tclvalue(Env$l.var$legende)==1) {
    graphe.legende(type="cam")
  }
}


#-------------------------------------------------
# Graphe - calculer les barres d'erreur (courbe)
#-------------------------------------------------

graphe.erreurs.calculer2<-function(varX,varY,niveau=NULL,valeurs=NULL,facteur=NULL) {
  erreur.inf<-NULL
  erreur.sup<-NULL
  if (tclvalue(Env$l.var$plusieurs)==0) {
    if (tclvalue(Env$l.var$erreur)==Env$voc[97,1]) {
	erreur<-NULL
	for (i in 1:nlevels(varX)) {
	  erreur<-c(erreur,sqrt((valeurs[i]*(1-valeurs[i]))/(length(varY[varX==levels(varX)[i]])-1)))
	}
	erreur.inf<-erreur.sup<-erreur
    } else if (tclvalue(Env$l.var$erreur)==Env$voc[98,1]) {
	for (i in 1:nlevels(varX)) {
	  erreur.inf<-c(erreur.inf,valeurs[i]-binom.test(length(na.omit(varY[varY==niveau & varX==levels(varX)[i]])),length(na.omit(varY[varX==levels(varX)[i]])))$conf.int[1])
	  erreur.sup<-c(erreur.sup,binom.test(length(na.omit(varY[varY==niveau & varX==levels(varX)[i]])),length(na.omit(varY[varX==levels(varX)[i]])))$conf.int[2]-valeurs[i])
	}
    } else {
	erreur.inf<-erreur.sup<-rep(0,nlevels(varX))
    }
  } else {
    if (tclvalue(Env$l.var$erreur)==Env$voc[97,1]) {
	erreur<-matrix(0,nrow=nlevels(facteur),ncol=nlevels(varX))
	for (i in 1:nlevels(varX)) {
	  for (j in 1:nlevels(facteur)) {
	    erreur[j,i]<-sqrt((valeurs[j,i]*(1-valeurs[j,i]))/(length(varY[facteur==levels(facteur)[j] & varX==levels(varX)[i]])-1))
	  }
	}
	erreur.inf<-erreur.sup<-erreur
    } else if (tclvalue(Env$l.var$erreur)==Env$voc[98,1]) {
	erreur.inf<-matrix(0,nrow=nlevels(facteur),ncol=nlevels(varX))
	erreur.sup<-matrix(0,nrow=nlevels(facteur),ncol=nlevels(varX))
	for (i in 1:nlevels(varX)) {
	  for (j in 1:nlevels(facteur)) {
	    erreur.inf[j,i]<-valeurs[j,i]-binom.test(length(na.omit(varY[varY==niveau & varX==levels(varX)[i] & facteur==levels(facteur)[j]])),
		length(varY[varX==levels(varX)[i] & facteur==levels(facteur)[j]]))$conf.int[1]
	    erreur.sup[j,i]<-binom.test(length(na.omit(varY[varY==niveau & varX==levels(varX)[i] & facteur==levels(facteur)[j]])),
		length(varY[varX==levels(varX)[i] & facteur==levels(facteur)[j]]))$conf.int[2]-valeurs[j,i]
	  }
	}
    } else {
	erreur.inf<-erreur.sup<-matrix(0,nrow=nlevels(facteur),ncol=nlevels(varX))
    }
  }
  return(list(erreur.inf=erreur.inf,erreur.sup=erreur.sup))
}


#-------------------------------------------------
# Tracer la courbe - limites des axes
#-------------------------------------------------

tracer.courbe.limites<-function(varX,valeurs,erreur.inf,erreur.sup) {
  x.inf<-if(tclvalue(Env$l.var$liminf.axehor)=="Auto") {
    0.8*min(varX,na.rm=TRUE)
  } else {
    as.numeric(tclvalue(Env$l.var$liminf.axehor))
  }
  x.sup<-if(tclvalue(Env$l.var$limsup.axehor)=="Auto") {
    1.1*max(varX,na.rm=TRUE)
  } else {
    as.numeric(tclvalue(Env$l.var$limsup.axehor))
  }
  y.inf<-if(tclvalue(Env$l.var$liminf.axever)=="Auto") {
    0.8*min(valeurs-erreur.inf,na.rm=TRUE)
  } else {
    as.numeric(tclvalue(Env$l.var$liminf.axever))
  }
  y.sup<-if(tclvalue(Env$l.var$limsup.axever)=="Auto") {
    1.1*max(valeurs+erreur.sup,na.rm=TRUE)
  } else {
    as.numeric(tclvalue(Env$l.var$limsup.axever))
  }
  return(list(xinf=x.inf,xsup=x.sup,yinf=y.inf,ysup=y.sup))
}


#-------------------------------------------------
# Tracer la courbe - type moyenne, plusieurs=0
#-------------------------------------------------

tracer.courbe.moyun<-function() {
  if(nchar(tclvalue(Env$l.var$facteur1))>0 & tclvalue(Env$l.var$facteur1)!=Env$voc[82,1]) {
    varX<-Env$dataset[,tclvalue(Env$l.var$varX)][Env$dataset[,tclvalue(Env$l.var$facteur1)]==levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(tclvalue(Env$l.var$niveau))+1]]
    varY<-Env$dataset[,tclvalue(Env$l.var$varY)][Env$dataset[,tclvalue(Env$l.var$facteur1)]==levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(tclvalue(Env$l.var$niveau))+1]]
  } else {
    varX<-Env$dataset[,tclvalue(Env$l.var$varX)]
    varY<-Env$dataset[,tclvalue(Env$l.var$varY)]
  }
  valeurs<-tapply(varY,varX,function(x) mean(x,na.rm=TRUE))
  erreurs<-graphe.erreurs.calculer(variable=varY,facteur1=factor(varX))
  erreurs$erreur.inf[is.na(erreurs$erreur.inf)] <- 0
  erreurs$erreur.sup[is.na(erreurs$erreur.sup)] <- 0
  limites<-tracer.courbe.limites(varX=varX,valeurs=valeurs,erreur.inf=erreurs$erreur.inf,erreur.sup=erreurs$erreur.sup)
  Env$l.code$x.inf<-x.inf<-limites$xinf
  Env$l.code$x.sup<-x.sup<-limites$xsup
  Env$l.code$y.inf<-y.inf<-limites$yinf
  Env$l.code$y.sup<-y.sup<-limites$ysup
  symbole<-graphe.symboles(num=as.numeric(tclvalue(Env$l.var$symboleA)))
  plot(valeurs~as.numeric(names(valeurs)),axes=FALSE,ann=FALSE,xlim=c(x.inf,x.sup),ylim=c(y.inf,y.sup),log=graphe.log(),
    col=tclvalue(Env$l.var$couleur2A),pch=symbole,cex=as.numeric(tclvalue(Env$l.var$taille.ptsA)),
    type=type.ligne(type=tclvalue(Env$l.var$type.courbeA)),lty=type.trait(type=tclvalue(Env$l.var$trait1)),
    lwd=as.numeric(tclvalue(Env$l.var$epaisseur1)))
  if (nchar(tclvalue(Env$l.var$erreur))>0 & tclvalue(Env$l.var$erreur)!=Env$voc[95,1]) {
    graphe.erreurs.tracer(abscisses=as.numeric(names(valeurs)),valeurs=valeurs,erreur.inf=erreurs$erreur.inf,
	erreur.sup=erreurs$erreur.sup,couleur=tclvalue(Env$l.var$couleur2A),amplitude=x.sup-x.inf)
  }
  graphe.titre()
  graphe.axes()
  graphe.box()
}


#-------------------------------------------------
# Tracer la courbe - type moyenne, plusieurs=1
#-------------------------------------------------

tracer.courbe.moyplusieurs<-function() {
  varX<-Env$dataset[,tclvalue(Env$l.var$varX)][Env$dataset[,tclvalue(Env$l.var$facteur1)]%in%levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]])+1]]
  varY<-Env$dataset[,tclvalue(Env$l.var$varY)][Env$dataset[,tclvalue(Env$l.var$facteur1)]%in%levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]])+1]]
  facteur<-Env$dataset[,tclvalue(Env$l.var$facteur1)][Env$dataset[,tclvalue(Env$l.var$facteur1)]%in%levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]])+1]]
  valeurs<-tapply(varY,list(facteur,varX),function(x) mean(x,na.rm=TRUE))
  erreurs<-graphe.erreurs.calculer(variable=varY,facteur1=factor(varX),facteur2=facteur)
  erreurs$erreur.inf[is.na(erreurs$erreur.inf)] <- 0
  erreurs$erreur.sup[is.na(erreurs$erreur.sup)] <- 0
  limites<-tracer.courbe.limites(varX=varX,valeurs=valeurs,erreur.inf=erreurs$erreur.inf,erreur.sup=erreurs$erreur.sup)
  Env$l.code$x.inf<-x.inf<-limites$xinf
  Env$l.code$x.sup<-x.sup<-limites$xsup
  Env$l.code$y.inf<-y.inf<-limites$yinf
  Env$l.code$y.sup<-y.sup<-limites$ysup
  symboles<-graphe.symboles(num=Env$l.var$symboleB)
  lignes<-type.ligne(type=Env$l.var$type.courbeB)
  traits<-type.trait(type=Env$l.var$trait2)
  plot(valeurs[1,]~as.numeric(colnames(valeurs)),axes=FALSE,ann=FALSE,xlim=c(x.inf,x.sup),ylim=c(y.inf,y.sup),log=graphe.log(),
    col=Env$l.var$couleur2B[1],pch=symboles[1],cex=Env$l.var$taille.ptsB[1],type=lignes[1],lty=traits[1],
    lwd=Env$l.var$epaisseur2[1])
  for (i in 2:length(Env$l.var$noms1)) {
    lines(as.numeric(colnames(valeurs)),valeurs[i,],col=Env$l.var$couleur2B[i],pch=symboles[i],cex=Env$l.var$taille.ptsB[i],
    type=lignes[i],lty=traits[i],lwd=Env$l.var$epaisseur2[i])
  }
  if (nchar(tclvalue(Env$l.var$erreur))>0 & tclvalue(Env$l.var$erreur)!=Env$voc[95,1]) {
    for (i in 1:length(Env$l.var$noms1)) {
	graphe.erreurs.tracer(abscisses=as.numeric(colnames(valeurs)),valeurs=valeurs[i,],erreur.inf=erreurs$erreur.inf[i,],
	  erreur.sup=erreurs$erreur.sup[i,],couleur=Env$l.var$couleur2B[i],amplitude=x.sup-x.inf)
    }
  }
  if (tclvalue(Env$l.var$legende)==1) {
    graphe.legende(type="courbe",symboles=symboles,lignes=lignes,traits=traits)
  }
  graphe.titre()
  graphe.axes()
  graphe.box()
}


#-------------------------------------------------
# Tracer la courbe - type proportion, plusieurs=0
#-------------------------------------------------

tracer.courbe.propun<-function() {
  if(nchar(tclvalue(Env$l.var$facteur1))>0 & tclvalue(Env$l.var$facteur1)!=Env$voc[82,1]) {
    varX<-Env$dataset[,tclvalue(Env$l.var$varX.prop)][Env$dataset[,tclvalue(Env$l.var$facteur1)]==levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(tclvalue(Env$l.var$niveau))+1]]
    varY<-Env$dataset[,tclvalue(Env$l.var$proportions)][Env$dataset[,tclvalue(Env$l.var$facteur1)]==levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(tclvalue(Env$l.var$niveau))+1]]
  } else {
    varX<-Env$dataset[,tclvalue(Env$l.var$varX.prop)]
    varY<-Env$dataset[,tclvalue(Env$l.var$proportions)]
  }
  valeurs<-integer(nlevels(factor(varX)))
  for (i in 1:nlevels(factor(varX))) {
    valeurs[i]<-length(na.omit(varY[varY==tclvalue(Env$l.var$prop.niveaux) & factor(varX)==levels(factor(varX))[i]]))/length(na.omit(varY[factor(varX)==levels(factor(varX))[i]]))
  }
  erreurs<-graphe.erreurs.calculer2(varX=factor(varX),varY=varY,niveau=tclvalue(Env$l.var$prop.niveaux),valeurs=valeurs)
  erreurs$erreur.inf[is.na(erreurs$erreur.inf)] <- 0
  erreurs$erreur.sup[is.na(erreurs$erreur.sup)] <- 0
  limites<-tracer.courbe.limites(varX=varX,valeurs=valeurs,erreur.inf=erreurs$erreur.inf,erreur.sup=erreurs$erreur.sup)
  Env$l.code$x.inf<-x.inf<-limites$xinf
  Env$l.code$x.sup<-x.sup<-limites$xsup
  Env$l.code$y.inf<-y.inf<-limites$yinf
  Env$l.code$y.sup<-y.sup<-limites$ysup
  symbole<-graphe.symboles(num=as.numeric(tclvalue(Env$l.var$symboleA)))
  plot(valeurs~as.numeric(as.character(levels(factor(varX)))),axes=FALSE,ann=FALSE,xlim=c(x.inf,x.sup),ylim=c(y.inf,y.sup),log=graphe.log(),
    col=tclvalue(Env$l.var$couleur2A),pch=symbole,cex=as.numeric(tclvalue(Env$l.var$taille.ptsA)),
    type=type.ligne(type=tclvalue(Env$l.var$type.courbeA)),lty=type.trait(type=tclvalue(Env$l.var$trait1)),
    lwd=as.numeric(tclvalue(Env$l.var$epaisseur1)))
  if (nchar(tclvalue(Env$l.var$erreur))>0 & tclvalue(Env$l.var$erreur)!=Env$voc[95,1]) {
    graphe.erreurs.tracer(abscisses=as.numeric(as.character(levels(factor(varX)))),valeurs=valeurs,erreur.inf=erreurs$erreur.inf,
	erreur.sup=erreurs$erreur.sup,couleur=tclvalue(Env$l.var$couleur2A),amplitude=x.sup-x.inf)
  }
  graphe.titre()
  graphe.axes()
  graphe.box()
}


#-------------------------------------------------
# Tracer la courbe - type proportion, plusieurs=1
#-------------------------------------------------

tracer.courbe.propplusieurs<-function() {
  varX<-Env$dataset[,tclvalue(Env$l.var$varX.prop)][Env$dataset[,tclvalue(Env$l.var$facteur1)]%in%levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]])+1]]
  varY<-Env$dataset[,tclvalue(Env$l.var$proportions)][Env$dataset[,tclvalue(Env$l.var$facteur1)]%in%levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]])+1]]
  facteur<-Env$dataset[,tclvalue(Env$l.var$facteur1)][Env$dataset[,tclvalue(Env$l.var$facteur1)]%in%levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]])+1]]
  valeurs<-matrix(0,nrow=nlevels(facteur),ncol=nlevels(factor(varX)))
  for (i in 1:nlevels(factor(varX))) {
    for (j in 1:nlevels(facteur)) {
	valeurs[j,i]<-length(na.omit(varY[varY==tclvalue(Env$l.var$prop.niveaux) & varX==levels(factor(varX))[i] & facteur==levels(facteur)[j]]))/length(na.omit(varY[varX==levels(factor(varX))[i] & facteur==levels(facteur)[j]]))
    }
  }
  erreurs<-graphe.erreurs.calculer2(varX=factor(varX),varY=varY,niveau=tclvalue(Env$l.var$prop.niveaux),valeurs=valeurs,facteur=facteur)
  erreurs$erreur.inf[is.na(erreurs$erreur.inf)] <- 0
  erreurs$erreur.sup[is.na(erreurs$erreur.sup)] <- 0
  limites<-tracer.courbe.limites(varX=varX,valeurs=valeurs,erreur.inf=erreurs$erreur.inf,erreur.sup=erreurs$erreur.sup)
  Env$l.code$x.inf<-x.inf<-limites$xinf
  Env$l.code$x.sup<-x.sup<-limites$xsup
  Env$l.code$y.inf<-y.inf<-limites$yinf
  Env$l.code$y.sup<-y.sup<-limites$ysup
  symboles<-graphe.symboles(num=Env$l.var$symboleB)
  lignes<-type.ligne(type=Env$l.var$type.courbeB)
  traits<-type.trait(type=Env$l.var$trait2)
  plot(valeurs[1,]~as.numeric(as.character(levels(factor(varX)))),axes=FALSE,ann=FALSE,xlim=c(x.inf,x.sup),ylim=c(y.inf,y.sup),log=graphe.log(),
    col=Env$l.var$couleur2B[1],pch=symboles[1],cex=Env$l.var$taille.ptsB[1],type=lignes[1],lty=traits[1],
    lwd=Env$l.var$epaisseur2[1])
  for (i in 2:length(Env$l.var$noms1)) {
    lines(as.numeric(as.character(levels(factor(varX)))),valeurs[i,],col=Env$l.var$couleur2B[i],pch=symboles[i],cex=Env$l.var$taille.ptsB[i],
    type=lignes[i],lty=traits[i],lwd=Env$l.var$epaisseur2[i])
  }
  if (nchar(tclvalue(Env$l.var$erreur))>0 & tclvalue(Env$l.var$erreur)!=Env$voc[95,1]) {
    for (i in 1:length(Env$l.var$noms1)) {
	graphe.erreurs.tracer(abscisses=as.numeric(as.character(levels(factor(varX)))),valeurs=valeurs[i,],erreur.inf=erreurs$erreur.inf[i,],
	  erreur.sup=erreurs$erreur.sup[i,],couleur=Env$l.var$couleur2B[i],amplitude=x.sup-x.inf)
    }
  }
  if (tclvalue(Env$l.var$legende)==1) {
    graphe.legende(type="courbe",symboles=symboles,lignes=lignes,traits=traits)
  }
  graphe.titre()
  graphe.axes()
  graphe.box()
}


#-------------------------------------------------
# Tracer le nuage - limites des axes
#-------------------------------------------------

tracer.nuage.limites<-function(varX,varY) {
  x.inf<-if(tclvalue(Env$l.var$liminf.axehor)=="Auto") {
    0.8*min(varX,na.rm=TRUE)
  } else {
    as.numeric(tclvalue(Env$l.var$liminf.axehor))
  }
  x.sup<-if(tclvalue(Env$l.var$limsup.axehor)=="Auto") {
    1.1*max(varX,na.rm=TRUE)
  } else {
    as.numeric(tclvalue(Env$l.var$limsup.axehor))
  }
  y.inf<-if(tclvalue(Env$l.var$liminf.axever)=="Auto") {
    0.8*min(varY,na.rm=TRUE)
  } else {
    as.numeric(tclvalue(Env$l.var$liminf.axever))
  }
  y.sup<-if(tclvalue(Env$l.var$limsup.axever)=="Auto") {
    1.1*max(varY,na.rm=TRUE)
  } else {
    as.numeric(tclvalue(Env$l.var$limsup.axever))
  }
  return(list(xinf=x.inf,xsup=x.sup,yinf=y.inf,ysup=y.sup))
}


#-------------------------------------------------
# Tracer le nuage - plusieurs=0
#-------------------------------------------------

tracer.nuage.un<-function() {
  if(nchar(tclvalue(Env$l.var$facteur1))>0 & tclvalue(Env$l.var$facteur1)!=Env$voc[82,1]) {
    varX<-Env$dataset[,tclvalue(Env$l.var$varX)][Env$dataset[,tclvalue(Env$l.var$facteur1)]==levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(tclvalue(Env$l.var$niveau))+1]]
    varY<-Env$dataset[,tclvalue(Env$l.var$varY)][Env$dataset[,tclvalue(Env$l.var$facteur1)]==levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(tclvalue(Env$l.var$niveau))+1]]
  } else {
    varX<-Env$dataset[,tclvalue(Env$l.var$varX)]
    varY<-Env$dataset[,tclvalue(Env$l.var$varY)]
  }
  limites<-tracer.nuage.limites(varX=varX,varY=varY)
  Env$l.code$x.inf<-x.inf<-limites$xinf
  Env$l.code$x.sup<-x.sup<-limites$xsup
  Env$l.code$y.inf<-y.inf<-limites$yinf
  Env$l.code$y.sup<-y.sup<-limites$ysup
  symbole<-graphe.symboles(num=as.numeric(tclvalue(Env$l.var$symboleA)))
  plot(varY~varX,axes=FALSE,ann=FALSE,xlim=c(x.inf,x.sup),ylim=c(y.inf,y.sup),log=graphe.log(),
    col=tclvalue(Env$l.var$couleur2A),pch=symbole,cex=as.numeric(tclvalue(Env$l.var$taille.ptsA)))
  if (nchar(tclvalue(Env$l.var$droiteA))>0 & tclvalue(Env$l.var$droiteA)!=Env$voc[95,1]) {
    if (tclvalue(Env$l.var$droiteA)==Env$voc[145,1]) {
	model<-lm(varY~varX)
	abline(model$coefficients,col=tclvalue(Env$l.var$couleur2A),
	  lty=type.trait(type=tclvalue(Env$l.var$trait1)),lwd=as.numeric(tclvalue(Env$l.var$epaisseur1)))
	if (tclvalue(Env$l.var$intervalA)==Env$voc[261,1]) {
	  varX2<-seq(min(varX,na.rm=TRUE),max(varX,na.rm=TRUE),abs(max(varX,na.rm=TRUE)-min(varX,na.rm=TRUE))/1000)
	  varY2<-predict(model,list(varX=varX2),interval="confidence")
	  lines(varX2,varY2[,"lwr"],lty=2,col=tclvalue(Env$l.var$couleur2A))
	  lines(varX2,varY2[,"upr"],lty=2,col=tclvalue(Env$l.var$couleur2A))
	}
	if (tclvalue(Env$l.var$intervalA)==Env$voc[262,1]) {
	  varX2<-seq(min(varX,na.rm=TRUE),max(varX,na.rm=TRUE),abs(max(varX,na.rm=TRUE)-min(varX,na.rm=TRUE))/1000)
	  varY2<-predict(model,list(varX=varX2),interval="prediction")
	  lines(varX2,varY2[,"lwr"],lty=3,col=tclvalue(Env$l.var$couleur2A))
	  lines(varX2,varY2[,"upr"],lty=3,col=tclvalue(Env$l.var$couleur2A))
	}
	if (tclvalue(Env$l.var$intervalA)==Env$voc[263,1]) {
	  varX2<-seq(min(varX,na.rm=TRUE),max(varX,na.rm=TRUE),abs(max(varX,na.rm=TRUE)-min(varX,na.rm=TRUE))/1000)
	  varY2.a<-predict(model,list(varX=varX2),interval="confidence")
	  lines(varX2,varY2.a[,"lwr"],lty=2,col=tclvalue(Env$l.var$couleur2A))
	  lines(varX2,varY2.a[,"upr"],lty=2,col=tclvalue(Env$l.var$couleur2A))
	  varY2.b<-predict(model,list(varX=varX2),interval="prediction")
	  lines(varX2,varY2.b[,"lwr"],lty=3,col=tclvalue(Env$l.var$couleur2A))
	  lines(varX2,varY2.b[,"upr"],lty=3,col=tclvalue(Env$l.var$couleur2A))
	}
    } else
    if (tclvalue(Env$l.var$droiteA)==Env$voc[146,1]) {
	b<-sd(varY,na.rm=TRUE)/sd(varX,na.rm=TRUE)*sign(cov(varX,varY,use="complete.obs"))
	a<-mean(varY,na.rm=TRUE)-b*mean(varX,na.rm=TRUE)
	abline(a,b,col=tclvalue(Env$l.var$couleur2A),
	  lty=type.trait(type=tclvalue(Env$l.var$trait1)),lwd=as.numeric(tclvalue(Env$l.var$epaisseur1)))
    } else
    if (tclvalue(Env$l.var$droiteA)==Env$voc[147,1]) {
	model<-lm(varY~varX+I(varX^2))
	varX2<-seq(min(varX,na.rm=TRUE),max(varX,na.rm=TRUE),abs(max(varX,na.rm=TRUE)-min(varX,na.rm=TRUE))/1000)
	varY2<-predict(model,list(varX=varX2))
	lines(varX2,varY2,col=tclvalue(Env$l.var$couleur2A),lty=type.trait(type=tclvalue(Env$l.var$trait1)),
	  lwd=as.numeric(tclvalue(Env$l.var$epaisseur1)))
	if (tclvalue(Env$l.var$intervalA)==Env$voc[261,1]) {
	  varY3<-predict(model,list(varX=varX2),interval="confidence")
	  lines(varX2,varY3[,"lwr"],lty=2,col=tclvalue(Env$l.var$couleur2A))
	  lines(varX2,varY3[,"upr"],lty=2,col=tclvalue(Env$l.var$couleur2A))
	}
	if (tclvalue(Env$l.var$intervalA)==Env$voc[262,1]) {
	  varY3<-predict(model,list(varX=varX2),interval="prediction")
	  lines(varX2,varY3[,"lwr"],lty=3,col=tclvalue(Env$l.var$couleur2A))
	  lines(varX2,varY3[,"upr"],lty=3,col=tclvalue(Env$l.var$couleur2A))
	}
	if (tclvalue(Env$l.var$intervalA)==Env$voc[263,1]) {
	  varY3.a<-predict(model,list(varX=varX2),interval="confidence")
	  lines(varX2,varY3.a[,"lwr"],lty=2,col=tclvalue(Env$l.var$couleur2A))
	  lines(varX2,varY3.a[,"upr"],lty=2,col=tclvalue(Env$l.var$couleur2A))
	  varY3.b<-predict(model,list(varX=varX2),interval="prediction")
	  lines(varX2,varY3.b[,"lwr"],lty=3,col=tclvalue(Env$l.var$couleur2A))
	  lines(varX2,varY3.b[,"upr"],lty=3,col=tclvalue(Env$l.var$couleur2A))
	}
    } else
    if (tclvalue(Env$l.var$droiteA)==Env$voc[148,1]) {
	panel.smooth(varX,varY,pch=symbole,cex=as.numeric(tclvalue(Env$l.var$taille.ptsA)),
	  col=tclvalue(Env$l.var$couleur2A),col.smooth=tclvalue(Env$l.var$couleur2A),
	  lty=type.trait(type=tclvalue(Env$l.var$trait1)),lwd=as.numeric(tclvalue(Env$l.var$epaisseur1)))
    }
  }
  if (tclvalue(Env$l.var$ptlab)==1) {
    text(varX,varY,pos=3,offset=0.4,cex=0.65*as.numeric(tclvalue(Env$l.var$taille.ptsA)),col=tclvalue(Env$l.var$couleur2A))
  }
  graphe.titre()
  graphe.axes()
  graphe.box()
}


#-------------------------------------------------
# Tracer le nuage - plusieurs=1
#-------------------------------------------------

tracer.nuage.plusieurs<-function() {
  varX<-Env$dataset[,tclvalue(Env$l.var$varX)][Env$dataset[,tclvalue(Env$l.var$facteur1)]%in%levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]])+1]]
  varY<-Env$dataset[,tclvalue(Env$l.var$varY)][Env$dataset[,tclvalue(Env$l.var$facteur1)]%in%levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]])+1]]
  facteur<-Env$dataset[,tclvalue(Env$l.var$facteur1)][Env$dataset[,tclvalue(Env$l.var$facteur1)]%in%levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]])+1]]
  niveaux<-levels(facteur)[as.numeric(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]])+1]
  limites<-tracer.nuage.limites(varX=varX,varY=varY)
  Env$l.code$x.inf<-x.inf<-limites$xinf
  Env$l.code$x.sup<-x.sup<-limites$xsup
  Env$l.code$y.inf<-y.inf<-limites$yinf
  Env$l.code$y.sup<-y.sup<-limites$ysup
  symboles<-graphe.symboles(num=Env$l.var$symboleB)
  plot(varY[facteur==niveaux[1]]~varX[facteur==niveaux[1]],axes=FALSE,ann=FALSE,xlim=c(x.inf,x.sup),ylim=c(y.inf,y.sup),
    log=graphe.log(),col=Env$l.var$couleur2B[1],pch=symboles[1],cex=Env$l.var$taille.ptsB[1])
  for (i in 2:length(niveaux)) {
    points(varY[facteur==niveaux[i]]~varX[facteur==niveaux[i]],col=Env$l.var$couleur2B[i],pch=symboles[i],
	cex=Env$l.var$taille.ptsB[i])
  }
  if (any(nchar(Env$l.var$droiteB)>0 & Env$l.var$droiteB!=Env$voc[95,1])) {
    for (i in 1:length(niveaux)) {
	if (Env$l.var$droiteB[i]==Env$voc[145,1]) {
	  x<-varX[facteur==niveaux[i]]
	  y<-varY[facteur==niveaux[i]]
	  model<-lm(y~x)
	  abline(model$coefficients,col=Env$l.var$couleur2B[i],
	    lty=type.trait(type=Env$l.var$trait2[i]),lwd=Env$l.var$epaisseur2[i])
	  if (Env$l.var$intervalB[i]==Env$voc[261,1]) {
	    varX2<-seq(min(varX[facteur==niveaux[i]],na.rm=TRUE),max(varX[facteur==niveaux[i]],na.rm=TRUE),abs(max(varX[facteur==niveaux[i]],na.rm=TRUE)-min(varX[facteur==niveaux[i]],na.rm=TRUE))/1000)
	    varY2<-predict(model,list(x=varX2),interval="confidence")
	    lines(varX2,varY2[,"lwr"],lty=2,col=Env$l.var$couleur2B[i])
	    lines(varX2,varY2[,"upr"],lty=2,col=Env$l.var$couleur2B[i])
	  }
	  if (Env$l.var$intervalB[i]==Env$voc[262,1]) {
	    varX2<-seq(min(varX[facteur==niveaux[i]],na.rm=TRUE),max(varX[facteur==niveaux[i]],na.rm=TRUE),abs(max(varX[facteur==niveaux[i]],na.rm=TRUE)-min(varX[facteur==niveaux[i]],na.rm=TRUE))/1000)
	    varY2<-predict(model,list(x=varX2),interval="prediction")
	    lines(varX2,varY2[,"lwr"],lty=3,col=Env$l.var$couleur2B[i])
	    lines(varX2,varY2[,"upr"],lty=3,col=Env$l.var$couleur2B[i])
	  }
	  if (Env$l.var$intervalB[i]==Env$voc[263,1]) {
	    varX2<-seq(min(varX[facteur==niveaux[i]],na.rm=TRUE),max(varX[facteur==niveaux[i]],na.rm=TRUE),abs(max(varX[facteur==niveaux[i]],na.rm=TRUE)-min(varX[facteur==niveaux[i]],na.rm=TRUE))/1000)
	    varY2.a<-predict(model,list(x=varX2),interval="confidence")
	    lines(varX2,varY2.a[,"lwr"],lty=2,col=Env$l.var$couleur2B[i])
	    lines(varX2,varY2.a[,"upr"],lty=2,col=Env$l.var$couleur2B[i])
	    varY2.b<-predict(model,list(x=varX2),interval="prediction")
	    lines(varX2,varY2.b[,"lwr"],lty=3,col=Env$l.var$couleur2B[i])
	    lines(varX2,varY2.b[,"upr"],lty=3,col=Env$l.var$couleur2B[i])
	  }
	} else
	if (Env$l.var$droiteB[i]==Env$voc[146,1]) {
	  b<-sd(varY[facteur==niveaux[i]],na.rm=TRUE)/sd(varX[facteur==niveaux[i]],na.rm=TRUE)*sign(cov(varX[facteur==niveaux[i]],varY[facteur==niveaux[i]],use="complete.obs"))
	  a<-mean(varY[facteur==niveaux[i]],na.rm=TRUE)-b*mean(varX[facteur==niveaux[i]],na.rm=TRUE)
	  abline(a,b,col=Env$l.var$couleur2B[i],lty=type.trait(type=Env$l.var$trait2[i]),lwd=Env$l.var$epaisseur2[i])
	} else
	if (Env$l.var$droiteB[i]==Env$voc[147,1]) {
	  x<-varX[facteur==niveaux[i]]
	  y<-varY[facteur==niveaux[i]]
	  model<-lm(y~x+I(x^2))
	  varX2<-seq(min(varX[facteur==niveaux[i]],na.rm=TRUE),max(varX[facteur==niveaux[i]],na.rm=TRUE),abs(max(varX[facteur==niveaux[i]],na.rm=TRUE)-min(varX[facteur==niveaux[i]],na.rm=TRUE))/1000)
	  varY2<-predict(model,list(x=varX2))
	  lines(varX2,varY2,col=Env$l.var$couleur2B[i],lty=type.trait(type=Env$l.var$trait2[i]),
	    lwd=Env$l.var$epaisseur2[i])
	  if (Env$l.var$intervalB[i]==Env$voc[261,1]) {
	    varX2<-seq(min(varX[facteur==niveaux[i]],na.rm=TRUE),max(varX[facteur==niveaux[i]],na.rm=TRUE),abs(max(varX[facteur==niveaux[i]],na.rm=TRUE)-min(varX[facteur==niveaux[i]],na.rm=TRUE))/1000)
	    varY2<-predict(model,list(x=varX2),interval="confidence")
	    lines(varX2,varY2[,"lwr"],lty=2,col=Env$l.var$couleur2B[i])
	    lines(varX2,varY2[,"upr"],lty=2,col=Env$l.var$couleur2B[i])
	  }
	  if (Env$l.var$intervalB[i]==Env$voc[262,1]) {
	    varX2<-seq(min(varX[facteur==niveaux[i]],na.rm=TRUE),max(varX[facteur==niveaux[i]],na.rm=TRUE),abs(max(varX[facteur==niveaux[i]],na.rm=TRUE)-min(varX[facteur==niveaux[i]],na.rm=TRUE))/1000)
	    varY2<-predict(model,list(x=varX2),interval="prediction")
	    lines(varX2,varY2[,"lwr"],lty=3,col=Env$l.var$couleur2B[i])
	    lines(varX2,varY2[,"upr"],lty=3,col=Env$l.var$couleur2B[i])
	  }
	  if (Env$l.var$intervalB[i]==Env$voc[263,1]) {
	    varX2<-seq(min(varX[facteur==niveaux[i]],na.rm=TRUE),max(varX[facteur==niveaux[i]],na.rm=TRUE),abs(max(varX[facteur==niveaux[i]],na.rm=TRUE)-min(varX[facteur==niveaux[i]],na.rm=TRUE))/1000)
	    varY2.a<-predict(model,list(x=varX2),interval="confidence")
	    lines(varX2,varY2.a[,"lwr"],lty=2,col=Env$l.var$couleur2B[i])
	    lines(varX2,varY2.a[,"upr"],lty=2,col=Env$l.var$couleur2B[i])
	    varY2.b<-predict(model,list(x=varX2),interval="prediction")
	    lines(varX2,varY2.b[,"lwr"],lty=3,col=Env$l.var$couleur2B[i])
	    lines(varX2,varY2.b[,"upr"],lty=3,col=Env$l.var$couleur2B[i])
	  }
	} else
	if (Env$l.var$droiteB[i]==Env$voc[148,1]) {
	  panel.smooth(varX[facteur==niveaux[i]],varY[facteur==niveaux[i]],pch=symboles[i],cex=Env$l.var$taille.ptsB[i],
	    col=Env$l.var$couleur2B[i],col.smooth=Env$l.var$couleur2B[i],lty=type.trait(type=Env$l.var$trait2[i]),
	    lwd=Env$l.var$epaisseur2[i])
	}
    }
  }
  if (tclvalue(Env$l.var$ptlab)==1) {
    for (i in 1:length(niveaux)) {
	text(varX[facteur==niveaux[i]],varY[facteur==niveaux[i]],pos=3,offset=0.4,cex=0.65*Env$l.var$taille.ptsB[i],col=Env$l.var$couleur2B[i])
    }
  }
  if (tclvalue(Env$l.var$legende)==1) {
    graphe.legende(type="nuage",symboles=symboles)
  }
  graphe.titre()
  graphe.axes()
  graphe.box()
}


#-------------------------------------------------
# Clic sur un bouton de la barre de navigation
#-------------------------------------------------

navigation<-function(type) {
  if (type%in%c("data","hist","moust","barres","cam","courbe","nuage")) {
  fr2.close()
  fr3.close()
  fr4.close()
  fr5.close()
  fr6.close()
    if (type=="data") {
	Env$l.var$ecran<-"D"
	tkconfigure(Env$l.lab$lab1,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[1,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab2,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[2,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab3,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[3,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab4,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[4,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab5,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[18,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab6,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Lab0.gif",fsep=.Platform$file.sep)))
	tkconfigure(Env$l.wdg$but.lab6,state="disabled")
	fr1.openD()
    }
    if (type=="hist") {
	Env$l.var$ecran<-"H"
	msg(text="",type="info")
	tkconfigure(Env$l.lab$lab1,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[5,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab2,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[6,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab3,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[7,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab4,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[8,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab5,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[9,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab6,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Lab0.gif",fsep=.Platform$file.sep)))
	tkconfigure(Env$l.wdg$but.lab5,state="normal")
	tkconfigure(Env$l.wdg$but.lab6,state="disabled")
	reinit.variables()
	fr1.openH()
    }
    if (type=="moust") {
	Env$l.var$ecran<-"M"
	msg(text="",type="info")
	tkconfigure(Env$l.lab$lab1,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[5,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab2,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[6,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab3,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[7,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab4,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[10,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab5,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[11,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab6,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[13,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.wdg$but.lab5,state="normal")
	tkconfigure(Env$l.wdg$but.lab6,state="normal")
	reinit.variables()
	fr1.openM()
    }
    if (type=="barres") {
	Env$l.var$ecran<-"B"
	msg(text="",type="info")
	tkconfigure(Env$l.lab$lab1,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[5,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab2,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[6,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab3,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[7,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab4,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[8,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab5,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[12,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab6,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[13,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.wdg$but.lab5,state="normal")
	tkconfigure(Env$l.wdg$but.lab6,state="normal")
	reinit.variables()
	fr1.openB()
    }
    if (type=="cam") {
	Env$l.var$ecran<-"Ca"
	msg(text="",type="info")
	tkconfigure(Env$l.lab$lab1,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[5,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab2,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[6,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab3,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[14,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab4,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[15,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab5,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[13,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab6,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Lab0.gif",fsep=.Platform$file.sep)))
	tkconfigure(Env$l.wdg$but.lab5,state="normal")
	tkconfigure(Env$l.wdg$but.lab6,state="disabled")
	reinit.variables()
	fr1.openCa()
    }
    if (type=="courbe") {
	Env$l.var$ecran<-"Co"
	msg(text="",type="info")
	tkconfigure(Env$l.lab$lab1,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[5,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab2,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[6,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab3,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[7,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab4,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[16,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab5,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[12,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab6,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[13,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.wdg$but.lab5,state="normal")
	tkconfigure(Env$l.wdg$but.lab6,state="normal")
	reinit.variables()
	fr1.openCo()
    }
    if (type=="nuage") {
	Env$l.var$ecran<-"N"
	msg(text="",type="info")
	tkconfigure(Env$l.lab$lab1,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[5,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab2,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[6,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab3,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[7,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab4,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[17,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab5,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[19,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.lab$lab6,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[13,1],fsep=.Platform$file.sep)))
	tkconfigure(Env$l.wdg$but.lab5,state="normal")
	tkconfigure(Env$l.wdg$but.lab6,state="normal")
	reinit.variables()
	fr1.openN()
    }
  }
}


#-------------------------------------------------
# Ouverture de l'aide
#-------------------------------------------------

aide<-function() {
  if (Env$lang=="en") {browseURL(file.path(path.package("GrapheR"),"doc","manual_en.pdf",fsep=.Platform$file.sep))}
  if (Env$lang=="fr") {browseURL(file.path(path.package("GrapheR"),"doc","manual_fr.pdf",fsep=.Platform$file.sep))}
  if (Env$lang=="es") {browseURL(file.path(path.package("GrapheR"),"doc","manual_en.pdf",fsep=.Platform$file.sep))}
  if (Env$lang=="de") {browseURL(file.path(path.package("GrapheR"),"doc","manual_de.pdf",fsep=.Platform$file.sep))}
}


#-------------------------------------------------
# code - Demande si sauvegarde du code
#-------------------------------------------------

code.ask<-function() {
  question<-tkmessageBox(message=Env$voc[240,1],icon="info",type="yesno")
  if (tclvalue(question)=="yes") {
    Env$l.code$save<-TRUE
    Env$l.code$folder<-tk_choose.dir()
    code.open()
  } else {
    Env$l.code$save<-FALSE
  }
}


#-------------------------------------------------
# Code - Ouverture de la connexion avec
#   le fichier externe
#-------------------------------------------------

code.open<-function() {
  sink(file=file.path(Env$l.code$folder,paste(paste("GrapheR",paste(strsplit(as.character(Sys.Date()),split="-")[[1]],collapse="."),
    sep="-"),".R",sep=""),fsep=.Platform$file.sep),append=TRUE)
  if (Env$l.code$graphsnb==0) {
    cat("#----------------------------------------\n")
    cat(paste("# GrapheR - session of ",paste(strsplit(as.character(Sys.Date()),split="-")[[1]],collapse="."),"\n",sep=""))
    cat("#----------------------------------------")
  }
  code.graphtype()
  cat("# Loading of the dataset\n\n")
  cat(paste("dataset <- ",Env$loading,"\n",sep=""))
  cat("attach(dataset)\n\n")
  code.data()
  code.graph()
  cat("detach(dataset)\n\n")
  sink(NULL)
}


#-------------------------------------------------
# Code - Type de graphe
#-------------------------------------------------

code.graphtype<-function() {
  graph <- if (Env$l.var$ecran=="H") {"Histogram"} else if (Env$l.var$ecran=="M") {"Box plot"} else if (Env$l.var$ecran=="B") {"Bar plot"} else
    if (Env$l.var$ecran=="Ca") {"Pie chart"} else if (Env$l.var$ecran=="Co") {"Curve"} else if (Env$l.var$ecran=="N") {"Scatter plot"}
  if (Env$l.code$graphsnb>0) {cat("\n")}
  cat("\n#------------------------------\n")
  Env$l.code$graphsnb <- Env$l.code$graphsnb+1
  cat(paste("# GRAPH ",Env$l.code$graphsnb,": ",graph,sep=""))
  cat("\n#------------------------------\n\n")
}


#-------------------------------------------------
# Code - Création des données
#-------------------------------------------------

code.data<-function() {
  variable<-tclvalue(Env$l.var$variable)
  facteur1<-tclvalue(Env$l.var$facteur1)
  facteur2<-tclvalue(Env$l.var$facteur2)
  niveau<-tclvalue(Env$l.var$niveau)
  if (Env$l.var$ecran=="H") {
    if (nchar(facteur1)>0 & facteur1!=Env$voc[82,1]) {
	cat("# Preliminary data creation\n\n")
	cat(paste("variable <- ",variable,"[",facteur1,"==\"",niveau,"\"]\n\n",sep=""))
    }
    if (tclvalue(Env$l.var$hist.type)==Env$voc[40,1]) {
	if (nchar(facteur1)==0 | facteur1==Env$voc[82,1]) {
	  cat("# Preliminary data creation\n\n")
	}
	texte<-""
	if (nchar(facteur1)>0 & facteur1!=Env$voc[82,1]) {
	  texte<-"frequencies <- hist(variable"
	} else {
	  texte<-paste("frequencies <- hist(",variable,sep="")
	}
	if (tclvalue(Env$l.var$hist.barres)!="Auto") {texte<-paste(texte,", breaks=",as.numeric(tclvalue(Env$l.var$hist.barres))-1,sep="")}
	texte<-paste(texte,", plot=FALSE)$counts / length(na.omit(",sep="")
	if (nchar(facteur1)>0 & facteur1!=Env$voc[82,1]) {
	  texte<- paste(texte,"variable))\n\n",sep="")
	} else {
	  texte<-paste(texte,variable,"))\n\n",sep="")
	}
	if (nchar(facteur1)>0 & facteur1!=Env$voc[82,1]) {
	  texte<-paste(texte,"middles <- hist(variable",sep="")
	} else {
	  texte<-paste(texte,"middles <- hist(",variable,sep="")
	}
	if (tclvalue(Env$l.var$hist.barres)!="Auto") {texte<-paste(texte,", breaks=",as.numeric(tclvalue(Env$l.var$hist.barres))-1,sep="")}
	texte<-paste(texte,", plot=FALSE)$mids\n\n",sep="")
	cat(texte)
    }
  } else if (Env$l.var$ecran=="M") {
    if (nchar(facteur2)>0 & facteur2!=Env$voc[82,1]) {
	cat("# Preliminary data creation\n\n")
	cat(paste("interaction <- interaction(",facteur2,",",facteur1,")\n",sep=""))
	if (tclvalue(Env$l.var$boxmoy)==1) {
	  cat(paste("means <- tapply(",variable,", interaction, function(x) mean(x,na.rm=TRUE))\n",sep=""))
	}
	cat("\n")
    } else {
	if (tclvalue(Env$l.var$boxmoy)==1) {
	  cat("# Preliminary data creation\n\n")	
	  cat(paste("means <- tapply(",variable,", ",facteur1,", function(x) mean(x,na.rm=TRUE))\n\n",sep=""))
	}
    }
  } else if (Env$l.var$ecran=="B") {
    proportions<-tclvalue(Env$l.var$proportions)
    facteurprop<-tclvalue(Env$l.var$facteurprop)
    if (tclvalue(Env$l.var$moyprop)=="moy") {
	if (tclvalue(Env$l.var$plusieurs)==0) {
	  cat("# Preliminary data creation\n\n")
	  cat(paste("means <- tapply(",variable,",",facteur1,",function(x) mean(x,na.rm=TRUE))\n",sep=""))
	  if (tclvalue(Env$l.var$erreur)==Env$voc[96,1]) {
	    cat(paste("\nstd.dev <- tapply(",variable,",",facteur1,",function(x) sd(x,na.rm=TRUE))\n\n",sep=""))
	  } else if (tclvalue(Env$l.var$erreur)==Env$voc[97,1]) {
	    cat(paste("\nstd.err <- tapply(",variable,",",facteur1,",function(x) sd(x,na.rm=TRUE)/sqrt(length(na.omit(x))))\n\n",sep=""))
	  } else if (tclvalue(Env$l.var$erreur)==Env$voc[98,1]) {
	    cat(paste("\nci <- means - tapply(",variable,",",facteur1,",function(x) t.test(x)$conf.int[1])\n\n",sep=""))
	  } else {
	    cat("\n")
	  }
	} else {
	  cat("# Preliminary data creation\n\n")
	  cat(paste("means <- tapply(",variable,",list(",facteur2,",",facteur1,"),function(x) mean(x,na.rm=TRUE))\n",sep=""))
	  if (tclvalue(Env$l.var$erreur)==Env$voc[96,1]) {
	    cat(paste("\nstd.dev <- tapply(",variable,",list(",facteur2,",",facteur1,"),function(x) sd(x,na.rm=TRUE))\n\n",sep=""))
	  } else if (tclvalue(Env$l.var$erreur)==Env$voc[97,1]) {
	    cat(paste("\nstd.err <- tapply(",variable,",list(",facteur2,",",facteur1,"),function(x) sd(x,na.rm=TRUE)/sqrt(length(na.omit(x))))\n\n",sep=""))
	  } else if (tclvalue(Env$l.var$erreur)==Env$voc[98,1]) {
	    cat(paste("\nci <- means - tapply(",variable,",list(",facteur2,",",facteur1,"),function(x) t.test(x)$conf.int[1])\n\n",sep=""))
	  } else {
	    cat("\n")
	  }
	}
    } else {
	cat("# Preliminary data creation\n\n")
	cat(paste("prop.total <- matrix(0,nrow=nlevels(",proportions,"),ncol=nlevels(",facteurprop,"),dimnames=list(levels(",proportions,"),levels(",facteurprop,")))\n",sep=""))
	cat(paste("for (i in 1:nlevels(",facteurprop,")) {\n",sep=""))
	cat(paste("  for (j in 1:nlevels(",proportions,")) {\n",sep=""))
	cat(paste("    prop.total[j,i] <- length(na.omit(",proportions,"[",proportions,"==levels(",proportions,")[j] & ",facteurprop,"==levels(",facteurprop,")[i]]))/length(na.omit(",
	  proportions,"[",facteurprop,"==levels(",facteurprop,")[i]]))\n",sep=""))
	cat("  }\n}\n")
	if (tclvalue(Env$l.var$plusieurs)==0) {
	  niveaux<-as.numeric(tclvalue(Env$l.var$prop.niveaux))+1
	  cat(paste("proportions <- prop.total[",niveaux,",]\n",sep=""))
	  if (tclvalue(Env$l.var$erreur)==Env$voc[97,1]) {
	    cat("\nstd.err <- NULL\n")
	    cat(paste("for (i in 1:nlevels(",facteurprop,")) {\n",sep=""))
	    cat(paste("  std.err <- c(std.err,sqrt((proportions[i]*(1-proportions[i]))/(length(na.omit(",proportions,"[",facteurprop,"==levels(",facteurprop,")[i]]))-1)))\n",sep=""))
	    cat("}\n\n")
	  } else if (tclvalue(Env$l.var$erreur)==Env$voc[98,1]) {
	    cat("\nci.inf <- NULL\n")
	    cat("ci.sup <- NULL\n")
	    cat(paste("for (i in 1:nlevels(",facteurprop,")) {\n",sep=""))
	    cat(paste("  ci.inf <- c(ci.inf,proportions[i] - binom.test(length(na.omit(",proportions,"[",proportions,"==levels(",proportions,")[",niveaux,"] & ",facteurprop,
		"==levels(",facteurprop,")[i]])),length(na.omit(",proportions,"[",facteurprop,"==levels(",facteurprop,")[i]])))$conf.int[1])\n",sep=""))
	    cat(paste("  ci.sup <- c(ci.sup,binom.test(length(na.omit(",proportions,"[",proportions,"==levels(",proportions,")[",niveaux,"] & ",facteurprop,
		"==levels(",facteurprop,")[i]])),length(na.omit(",proportions,"[",facteurprop,"==levels(",facteurprop,")[i]])))$conf.int[2] - proportions[i])\n",sep=""))
	    cat("}\n\n")
	  } else {
	    cat("\n")
	  }
	} else {
  	  niveaux<-as.numeric(strsplit(tclvalue(Env$l.var$prop.niveaux),split=" ")[[1]])+1
	  cat(paste("proportions <- prop.total[c(",paste(niveaux,collapse=","),"),]\n",sep=""))
	  if (tclvalue(Env$l.var$erreur)==Env$voc[97,1]) {
	    cat("\nstd.err <- matrix(0,nrow=nrow(proportions),ncol=ncol(proportions),dimnames=list(rownames(proportions),colnames(proportions)))\n")
	    cat("for (i in 1:ncol(proportions)) {\n")
	    cat("  for (j in 1:nrow(proportions)) {\n")
	    cat(paste("    std.err[j,i] <- sqrt((proportions[j,i]*(1-proportions[j,i]))/(length(na.omit(",proportions,"[",facteurprop,"==levels(",facteurprop,")[i]]))-1))\n",sep=""))
	    cat("  }\n}\n\n")
	  } else if (tclvalue(Env$l.var$erreur)==Env$voc[98,1]) {
	    cat("\nci.inf <- matrix(0,nrow=nrow(proportions),ncol=ncol(proportions),dimnames=list(rownames(proportions),colnames(proportions)))\n")
	    cat("ci.sup <- matrix(0,nrow=nrow(proportions),ncol=ncol(proportions),dimnames=list(rownames(proportions),colnames(proportions)))\n")
	    cat("for (i in 1:ncol(proportions)) {\n")
	    cat(paste("  for (j in c(",paste(niveaux,collapse=","),")) {\n",sep=""))
	    cat(paste("    ci.inf[j,i] <- proportions[j,i] - binom.test(length(na.omit(",proportions,"[",proportions,"==levels(",proportions,")[j] & ",facteurprop,
		"==levels(",facteurprop,")[i]])),length(na.omit(",proportions,"[",facteurprop,"==levels(",facteurprop,")[i]])))$conf.int[1]\n",sep=""))
	    cat(paste("    ci.sup[j,i] <- binom.test(length(na.omit(",proportions,"[",proportions,"==levels(",proportions,")[j] & ",facteurprop,
		"==levels(",facteurprop,")[i]])),length(na.omit(",proportions,"[",facteurprop,"==levels(",facteurprop,")[i]])))$conf.int[2] - proportions[j,i]\n",sep=""))
	    cat("  }\n}\n\n")
	  } else {
	    cat("\n")
	  }
	}
    }
  } else if (Env$l.var$ecran=="Ca") {
    niveaux<-levels(factor(Env$dataset[,tclvalue(Env$l.var$variable)]))[as.numeric(strsplit(tclvalue(Env$l.var$parts.niveaux),split=" ")[[1]])+1]
    cat("# Preliminary data creation\n\n")
    cat(paste("prevariable1 <- factor(",variable,")\n",sep=""))
    cat(paste("prevariable2 <- droplevels(prevariable1[prevariable1 %in% c(\"",paste(niveaux,collapse="\",\""),"\")])\n",sep=""))
    cat("variable <- summary(prevariable2)\n\n")
  } else if (Env$l.var$ecran=="Co") {
    cat("# Preliminary data creation\n\n")
    if (tclvalue(Env$l.var$moyprop)=="moy") {
	varX<-tclvalue(Env$l.var$varX)
	varY<-tclvalue(Env$l.var$varY)
	if (tclvalue(Env$l.var$plusieurs)==0) {
	  if (nchar(tclvalue(Env$l.var$facteur1))>0 & tclvalue(Env$l.var$facteur1)!=Env$voc[82,1]) {
	    niveau<-levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(tclvalue(Env$l.var$niveau))+1]
	    cat(paste("varX <- ",varX,"[",facteur1,"==\"",niveau,"\"]\n",sep=""))
	    cat(paste("varY <- ",varY,"[",facteur1,"==\"",niveau,"\"]\n",sep=""))
	    cat("means <- tapply(varY,varX,function(x) mean(x,na.rm=TRUE))\n")
	    if (tclvalue(Env$l.var$erreur)==Env$voc[96,1]) {
		cat("\nstd.dev <- tapply(varY,varX,function(x) sd(x,na.rm=TRUE))\n\n")
	    } else if (tclvalue(Env$l.var$erreur)==Env$voc[97,1]) {
		cat("\nstd.err <- tapply(varY,varX,function(x) sd(x,na.rm=TRUE)/sqrt(length(na.omit(x))))\n\n")
	    } else if (tclvalue(Env$l.var$erreur)==Env$voc[98,1]) {
		cat("\nci <- means - tapply(varY,varX,function(x) t.test(x)$conf.int[1])\n\n")
	    } else {
		cat("\n")
	    }
	  } else {
	    cat(paste("means <- tapply(",varY,",",varX,",function(x) mean(x,na.rm=TRUE))\n",sep=""))
	    if (tclvalue(Env$l.var$erreur)==Env$voc[96,1]) {
		cat(paste("\nstd.dev <- tapply(",varY,",",varX,",function(x) sd(x,na.rm=TRUE))\n\n",sep=""))
	    } else if (tclvalue(Env$l.var$erreur)==Env$voc[97,1]) {
		cat(paste("\nstd.err <- tapply(",varY,",",varX,",function(x) sd(x,na.rm=TRUE)/sqrt(length(na.omit(x))))\n\n",sep=""))
	    } else if (tclvalue(Env$l.var$erreur)==Env$voc[98,1]) {
		cat(paste("\nci <- means - tapply(",varY,",",varX,",function(x) t.test(x)$conf.int[1])\n\n",sep=""))
	    } else {
		cat("\n")
	    }
	  }
	} else {
	  niveaux<-levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]])+1]
	  cat(paste("varX <- ",varX,"[",facteur1," %in% c(\"",paste(niveaux,collapse="\",\""),"\")]\n",sep=""))
	  cat(paste("varY <- ",varY,"[",facteur1," %in% c(\"",paste(niveaux,collapse="\",\""),"\")]\n",sep=""))
	  cat(paste("fact <- ",facteur1,"[",facteur1," %in% c(\"",paste(niveaux,collapse="\",\""),"\")]\n",sep=""))
	  cat("means <- tapply(varY,list(fact,varX),function(x) mean(x,na.rm=TRUE))\n")
	  if (tclvalue(Env$l.var$erreur)==Env$voc[96,1]) {
	    cat("\nstd.dev <- tapply(varY,list(fact,varX),function(x) sd(x,na.rm=TRUE))\n\n")
	  } else if (tclvalue(Env$l.var$erreur)==Env$voc[97,1]) {
	    cat("\nstd.err <- tapply(varY,list(fact,varX),function(x) sd(x,na.rm=TRUE)/sqrt(length(na.omit(x))))\n\n")
	  } else if (tclvalue(Env$l.var$erreur)==Env$voc[98,1]) {
	    cat("\nci <- means - tapply(varY,list(fact,varX),function(x) t.test(x)$conf.int[1])\n\n")
	  } else {
	    cat("\n")
	  }
	}
    } else {
	varX<-tclvalue(Env$l.var$varX.prop)
	varY<-tclvalue(Env$l.var$proportions)
	varY.niv <- tclvalue(Env$l.var$prop.niveaux)
	if (tclvalue(Env$l.var$plusieurs)==0) {
	  if (nchar(tclvalue(Env$l.var$facteur1))>0 & tclvalue(Env$l.var$facteur1)!=Env$voc[82,1]) {
	    niveau<-levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(tclvalue(Env$l.var$niveau))+1]
	    cat(paste("varX <- factor(",varX,"[",facteur1,"==\"",niveau,"\"])\n",sep=""))
	    cat(paste("varY <- ",varY,"[",facteur1,"==\"",niveau,"\"]\n",sep=""))
	    cat("proportions <- NULL\n")
	    cat("for (i in 1:nlevels(varX)) {\n")
	    cat(paste("  proportions <- c(proportions,length(na.omit(varY[varY==\"",varY.niv,"\" & varX==levels(varX)[i]]))/",
		"length(na.omit(varY[varX==levels(varX)[i]])))","\n",sep=""))
	    cat("}\n")
	    if (tclvalue(Env$l.var$erreur)==Env$voc[97,1]) {
		cat("\nstd.err <- NULL\n")
		cat("for (i in 1:nlevels(varX)) {\n")
		cat("  std.err <- c(std.err,sqrt((proportions[i]*(1-proportions[i]))/(length(na.omit(varY[varX==levels(varX)[i]]))-1)))\n")
		cat("}\n\n")
	    } else if (tclvalue(Env$l.var$erreur)==Env$voc[98,1]) {
		cat("\nci.inf <- NULL\n")
		cat("ci.sup <- NULL\n")
		cat("for (i in 1:nlevels(varX)) {\n")
		cat(paste("  ci.inf <- c(ci.inf,proportions[i] - binom.test(length(na.omit(varY[varY==\"",varY.niv,"\" & varX==levels(varX)[i]])",
		  "),length(na.omit(varY[varX==levels(varX)[i]])))$conf.int[1])\n",sep=""))
		cat(paste("  ci.sup <- c(ci.sup,binom.test(length(na.omit(varY[varY==\"",varY.niv,"\" & varX==levels(varX)[i]])",
		  "),length(na.omit(varY[varX==levels(varX)[i]])))$conf.int[2] - proportions[i])\n",sep=""))
		cat("}\n\n")
	    } else {
		cat("\n")
	    }
	  } else {
	    cat(paste("varX <- factor(",varX,")\n",sep=""))
	    cat("proportions <- NULL\n")
	    cat(paste("for (i in 1:nlevels(varX)) {\n",sep=""))
	    cat(paste("  proportions <- c(proportions,length(",varY,"[",varY,"==\"",varY.niv,"\" & varX==levels(varX)[i]])/",
		"length(na.omit(",varY,"[varX==levels(varX)[i]])))","\n",sep=""))
	    cat("}\n")
	    if (tclvalue(Env$l.var$erreur)==Env$voc[97,1]) {
		cat("\nstd.err <- NULL\n")
		cat("for (i in 1:nlevels(varX)) {\n")
		cat(paste("  std.err <- c(std.err,sqrt((proportions[i]*(1-proportions[i]))/(length(na.omit(",varY,"[varX==levels(varX)[i]]))-1)))\n",sep=""))
		cat("}\n\n")
	    } else if (tclvalue(Env$l.var$erreur)==Env$voc[98,1]) {
		cat("\nci.inf <- NULL\n")
		cat("ci.sup <- NULL\n")
		cat(paste("for (i in 1:nlevels(varX)) {\n",sep=""))
		cat(paste("  ci.inf <- c(ci.inf,proportions[i] - binom.test(length(na.omit(",varY,"[",varY,"==\"",varY.niv,"\" & varX==levels(varX)[i]])",
		  "),length(na.omit(",varY,"[varX==levels(varX)[i]])))$conf.int[1])\n",sep=""))
		cat(paste("  ci.sup <- c(ci.sup,binom.test(length(na.omit(",varY,"[",varY,"==\"",varY.niv,"\" & varX==levels(varX)[i]])",
		  "),length(na.omit(",varY,"[varX==levels(varX)[i]])))$conf.int[2] - proportions[i])\n",sep=""))
		cat("}\n\n")
	    } else {
		cat("\n")
	    }
	  }
	} else {
	  niveaux<-levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]])+1]
	  cat(paste("varX <- factor(",varX,"[",facteur1," %in% c(\"",paste(niveaux,collapse="\",\""),"\")])\n",sep=""))
	  cat(paste("varY <- ",varY,"[",facteur1," %in% c(\"",paste(niveaux,collapse="\",\""),"\")]\n",sep=""))
	  cat(paste("fact <- ",facteur1,"[",facteur1," %in% c(\"",paste(niveaux,collapse="\",\""),"\")]\n",sep=""))
	  cat("proportions <- matrix(0,nrow=nlevels(fact),ncol=nlevels(varX),dimnames=list(levels(fact),levels(varX)))\n")
	  cat("for (i in 1:nlevels(varX)) {\n")
	  cat("  for (j in 1:nlevels(fact)) {\n")
	  cat(paste("    proportions[j,i] <- length(na.omit(varY[varY==\"",varY.niv,"\" & varX==levels(varX)[i] & fact==levels(fact)[j]]))/length(na.omit(",
	    "varY[varX==levels(varX)[i] & fact==levels(fact)[j]]))\n",sep=""))
	  cat("  }\n}\n")
	  if (tclvalue(Env$l.var$erreur)==Env$voc[97,1]) {
	    cat("\nstd.err <- matrix(0,nrow=nrow(proportions),ncol=ncol(proportions),dimnames=list(rownames(proportions),colnames(proportions)))\n")
	    cat("for (i in 1:ncol(proportions)) {\n")
	    cat("  for (j in 1:nrow(proportions)) {\n")
	    cat("    std.err[j,i] <- sqrt((proportions[j,i]*(1-proportions[j,i]))/(length(na.omit(varY[varX==levels(varX)[i] & fact==levels(fact)[j]]))-1))\n")
	    cat("  }\n}\n\n")
	  } else if (tclvalue(Env$l.var$erreur)==Env$voc[98,1]) {
	    cat("\nci.inf <- matrix(0,nrow=nrow(proportions),ncol=ncol(proportions),dimnames=list(rownames(proportions),colnames(proportions)))\n")
	    cat("ci.sup <- matrix(0,nrow=nrow(proportions),ncol=ncol(proportions),dimnames=list(rownames(proportions),colnames(proportions)))\n")
	    cat("for (i in 1:ncol(proportions)) {\n")
	    cat(paste("  for (j in c(\"",paste(niveaux,collapse="\",\""),"\")) {\n",sep=""))
	    cat(paste("    ci.inf[j,i] <- proportions[j,i] - binom.test(length(na.omit(varY[varY==\"",varY.niv,"\" & varX==levels(varX)[i] & fact==j])),",
		"length(na.omit(varY[varX==levels(varX)[i] & fact==j])))$conf.int[1]\n",sep=""))
	    cat(paste("    ci.sup[j,i] <- binom.test(length(na.omit(varY[varY==\"",varY.niv,"\" & varX==levels(varX)[i] & fact==j])),",
		"length(na.omit(varY[varX==levels(varX)[i] & fact==j])))$conf.int[2] - proportions[j,i]\n",sep=""))
	    cat("  }\n}\n\n")
	  } else {
	    cat("\n")
	  }
	}
    }
  } else if (Env$l.var$ecran=="N") {
    varX<-tclvalue(Env$l.var$varX)
    varY<-tclvalue(Env$l.var$varY)
    facteur1<-tclvalue(Env$l.var$facteur1)
    if (tclvalue(Env$l.var$plusieurs)==0) {
	if (nchar(tclvalue(Env$l.var$facteur1))>0 & tclvalue(Env$l.var$facteur1)!=Env$voc[82,1]) {
	  niveau<-levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(tclvalue(Env$l.var$niveau))+1]
	  cat("# Preliminary data creation\n\n")
	  cat(paste("varX <- ",varX,"[",facteur1,"==\"",niveau,"\"]\n\n",sep=""))
	  cat(paste("varY <- ",varY,"[",facteur1,"==\"",niveau,"\"]\n\n",sep=""))
	}
    } else {
	cat("# Preliminary data creation\n\n")
	niveaux<-levels(Env$dataset[,tclvalue(Env$l.var$facteur1)])[as.numeric(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]])+1]
	cat(paste("varX <- ",varX,"[",facteur1," %in% c(\"",paste(niveaux,collapse="\",\""),"\")]\n\n",sep=""))
	cat(paste("varY <- ",varY,"[",facteur1," %in% c(\"",paste(niveaux,collapse="\",\""),"\")]\n\n",sep=""))
    }
  }
}


#-------------------------------------------------
# Code - Graphe (axes)
#-------------------------------------------------

code.graph.axes<-function() {
  texte<-ifelse(tclvalue(Env$l.var$graduations.orient)==Env$voc[246,1],"","par(las=1)\n")
  if (Env$l.var$ecran=="H" & tclvalue(Env$l.var$hist.type)==Env$voc[40,1]) {
    texte<-paste(texte,"axis(1, labels=middles, at=(0.5:(length(na.omit(frequencies))))",sep="")
    if (tclvalue(Env$l.var$graduations.col)!="black" & tclvalue(Env$l.var$graduations.col)!="#000000") {texte<-paste(texte,", col.axis=\"",tclvalue(Env$l.var$graduations.col),"\"",sep="")}
    if (tclvalue(Env$l.var$graduations.taille)!="1") {texte<-paste(texte,", cex.axis=",tclvalue(Env$l.var$graduations.taille),sep="")}
    texte<-paste(texte,")\n",sep="")
  } else if (Env$l.var$ecran=="M") {
    if (tclvalue(Env$l.var$box.orient)==Env$voc[68,1]) {
	if (tclvalue(Env$l.var$plusieurs)==0) {
	  texte<-paste(texte,"axis(1, labels=c(\"",paste(Env$l.var$noms1,collapse="\",\""),"\"), at=1:",length(Env$l.var$noms1),sep="")
	} else {
	  n.lab <- length(Env$l.var$noms1)
	  n.tot <- nlevels(Env$l.var$facteur.interaction)
	  n <- n.tot/n.lab
	  deb <- ((n.tot+1)-n*(n.lab-1))/2
	  texte<-paste(texte,"mtext(c(\"",paste(Env$l.var$noms1,collapse="\",\""),"\"), side=1, line=1, at=c(",paste(seq(deb,deb+n*(n.lab-1),n),collapse=","),")",sep="")
	}
    } else {
	texte<-paste(texte,"axis(1",sep="")
    }
    if (tclvalue(Env$l.var$plusieurs)==0) {
	if (tclvalue(Env$l.var$graduations.col)!="black" & tclvalue(Env$l.var$graduations.col)!="#000000") {texte<-paste(texte,", col.axis=\"",tclvalue(Env$l.var$graduations.col),"\"",sep="")}
	if (tclvalue(Env$l.var$graduations.taille)!="1") {texte<-paste(texte,", cex.axis=",tclvalue(Env$l.var$graduations.taille),sep="")}
    }
    texte<-paste(texte,")\n",sep="")
  } else if (Env$l.var$ecran=="B") {
    if (tclvalue(Env$l.var$encadre)==0) {
	ordonnee<-if (Env$l.code$y.inf>=0 & Env$l.code$y.sup>0) {
	  Env$l.code$y.inf
	} else if (Env$l.code$y.inf<0 & Env$l.code$y.sup>0){
	  0
	} else if (Env$l.code$y.inf<0 & Env$l.code$y.sup<=0) {
	  Env$l.code$y.sup
	}
	texte<-paste(texte,"abline(h=",ordonnee,")\n",sep="")
    }
    if (tclvalue(Env$l.var$nobar)==1) {
	texte<-paste(texte,"mtext(",sep="")
	if(tclvalue(Env$l.var$moyprop)=="moy"){
	  texte<-paste(texte,"c(\"",paste(Env$l.var$noms1,collapse="\",\""),"\")",sep="")
	} else {
	  texte<-paste(texte,"c(\"",paste(Env$l.var$nomsprop.fac,collapse="\",\""),"\")",sep="")
	}
	texte<-paste(texte,", side=1",sep="")
	if (nrow(t(Env$l.var$add.abscisses))==1) {
	  texte<-paste(texte,", at=graph",sep="")
	} else {
	  texte<-paste(texte,", at=c(",paste(round(apply(Env$l.var$add.abscisses,2,mean),1),collapse=","),")",sep="")
	}
	texte<-paste(texte,", line=1.1",sep="")
	texte<-paste(texte,",\n  cex=",tclvalue(Env$l.var$legendes.taille),sep="")
	texte<-paste(texte,", col=\"",tclvalue(Env$l.var$legendes.col),"\")\n",sep="")
    }
  } else {
    texte<-paste("axis(1",sep="")
    if (tclvalue(Env$l.var$graduations.col)!="black" & tclvalue(Env$l.var$graduations.col)!="#000000") {texte<-paste(texte,", col.axis=\"",tclvalue(Env$l.var$graduations.col),"\"",sep="")}
    if (tclvalue(Env$l.var$graduations.taille)!="1") {texte<-paste(texte,", cex.axis=",tclvalue(Env$l.var$graduations.taille),sep="")}
    texte<-paste(texte,")\n",sep="")
  }
  if (Env$l.var$ecran=="M") {
    if (tclvalue(Env$l.var$box.orient)==Env$voc[68,1]) {
	texte<-paste(texte,"axis(2",sep="")
	if (tclvalue(Env$l.var$graduations.col)!="black" & tclvalue(Env$l.var$graduations.col)!="#000000") {texte<-paste(texte,", col.axis=\"",tclvalue(Env$l.var$graduations.col),"\"",sep="")}
	if (tclvalue(Env$l.var$graduations.taille)!="1") {texte<-paste(texte,", cex.axis=",tclvalue(Env$l.var$graduations.taille),sep="")}
	texte<-paste(texte,")\n",sep="")
    } else {
	if (tclvalue(Env$l.var$plusieurs)==0) {
	  texte<-paste(texte,"axis(2, labels=c(\"",paste(Env$l.var$noms1,collapse="\",\""),"\"), at=1:",length(Env$l.var$noms1),sep="")
	} else {
	  n.lab <- length(Env$l.var$noms1)
	  n.tot <- nlevels(Env$l.var$facteur.interaction)
	  n <- n.tot/n.lab
	  deb <- ((n.tot+1)-n*(n.lab-1))/2
	  texte<-paste(texte,"mtext(c(\"",paste(Env$l.var$noms1,collapse="\",\""),"\"), side=2, line=1, at=c(",paste(seq(deb,deb+n*(n.lab-1),n),collapse=","),")",sep="")
	}
	if (tclvalue(Env$l.var$plusieurs)==0) {
	  if (tclvalue(Env$l.var$graduations.col)!="black" & tclvalue(Env$l.var$graduations.col)!="#000000") {texte<-paste(texte,", col.axis=\"",tclvalue(Env$l.var$graduations.col),"\"",sep="")}
	  if (tclvalue(Env$l.var$graduations.taille)!="1") {texte<-paste(texte,", cex.axis=",tclvalue(Env$l.var$graduations.taille),sep="")}
	}
	texte<-paste(texte,")\n",sep="")
    }
  } else {
    texte<-paste(texte,"axis(2",sep="")
    if (tclvalue(Env$l.var$graduations.col)!="black" & tclvalue(Env$l.var$graduations.col)!="#000000") {texte<-paste(texte,", col.axis=\"",tclvalue(Env$l.var$graduations.col),"\"",sep="")}
    if (tclvalue(Env$l.var$graduations.taille)!="1") {texte<-paste(texte,", cex.axis=",tclvalue(Env$l.var$graduations.taille),sep="")}
    texte<-paste(texte,")\n",sep="")
  }
  texte<-paste(texte,ifelse(tclvalue(Env$l.var$graduations.orient)==Env$voc[246,1],"\n","par(las=0)\n\n"),sep="")
  cat(texte)
}


#-------------------------------------------------
# Code - Graphe (titre)
#-------------------------------------------------

code.graph.titre<-function() {
  titre.x<-""
  titre.y<-""
  if (Env$l.var$ecran=="M") {
    titre.x<-if (tclvalue(Env$l.var$box.orient)==Env$voc[68,1]) {
	tclvalue(Env$l.var$titre.axenoms)
    } else {
	tclvalue(Env$l.var$titre.axevaleurs)
    }
    titre.y<-if (tclvalue(Env$l.var$box.orient)==Env$voc[68,1]) {
	tclvalue(Env$l.var$titre.axevaleurs)
    } else {
	tclvalue(Env$l.var$titre.axenoms)
    }
  } else {
    titre.x<-tclvalue(Env$l.var$titre.axehor)
    titre.y<-tclvalue(Env$l.var$titre.axever)
  }
  if (tclvalue(Env$l.var$titre)!="" | tclvalue(Env$l.var$soustitre)!="" | titre.x!="" | titre.y!="") {
    if (nchar(tclvalue(Env$l.var$soustitre))!=0) {
	if (tclvalue(Env$l.var$titre)!="") {
	  texte<-paste("title(main=\"",tclvalue(Env$l.var$titre),"\"",sep="")
	  if (tclvalue(Env$l.var$titre.col)!="black" & tclvalue(Env$l.var$titre.col)!="#000000") {texte<-paste(texte,", col.main=\"",tclvalue(Env$l.var$titre.col),"\"",sep="")}
	  texte<-paste(texte,", cex.main=",tclvalue(Env$l.var$titre.taille),sep="")
	  texte<-paste(texte,", line=2.2)\n",sep="")
	  cat(texte)
	}
	texte<-paste("title(main=\"",tclvalue(Env$l.var$soustitre),"\"",sep="")
	if (tclvalue(Env$l.var$titre.col)!="black" & tclvalue(Env$l.var$titre.col)!="#000000") {texte<-paste(texte,", col.main=\"",tclvalue(Env$l.var$titre.col),"\"",sep="")}
	texte<-paste(texte,", cex.main=",round(0.7*as.numeric(tclvalue(Env$l.var$titre.taille)),2),sep="")
	texte<-paste(texte,", line=0.9)\n",sep="")
	cat(texte)
	if (titre.x!="" | titre.y!="") {
	  texte<-"title("
	  if (titre.x!="") {
	    texte<-paste(texte,"xlab=\"",titre.x,"\"",sep="")
	  }
	  if (titre.y!="") {
	    if (titre.x!="") {
		texte<-paste(texte,", ylab=\"",titre.y,"\"",sep="")
	    } else {
		texte<-paste(texte,"ylab=\"",titre.y,"\"",sep="")
	    }
	  }
	  if (tclvalue(Env$l.var$legendes.col)!="black" & tclvalue(Env$l.var$legendes.col)!="#000000") {texte<-paste(texte,", col.lab=\"",tclvalue(Env$l.var$legendes.col),"\"",sep="")}
	  if (tclvalue(Env$l.var$legendes.taille)!="1") {texte<-paste(texte,", cex.lab=",tclvalue(Env$l.var$legendes.taille),sep="")}
	  texte<-paste(texte,")\n\n",sep="")
	  cat(texte)
	}
    } else {
	texte<-"title("
	if (tclvalue(Env$l.var$titre)!="") {
	  texte<-paste(texte,"main=\"",tclvalue(Env$l.var$titre),"\"",sep="")
	  if (tclvalue(Env$l.var$titre.col)!="black" & tclvalue(Env$l.var$titre.col)!="#000000") {texte<-paste(texte,", col.main=\"",tclvalue(Env$l.var$titre.col),"\"",sep="")}
	  if (tclvalue(Env$l.var$titre.taille)!="1.5") {texte<-paste(texte,", cex.main=",tclvalue(Env$l.var$titre.taille),sep="")}
	}
	if (titre.x!="") {
	  if (tclvalue(Env$l.var$titre)!="") {
	    texte<-paste(texte,",\n  xlab=\"",titre.x,"\"",sep="")
	  } else {
	    texte<-paste(texte,"xlab=\"",titre.x,"\"",sep="")
	  }
	}
	if (titre.y!="") {
	  if (tclvalue(Env$l.var$titre)!="" | titre.x!="") {
	    if (tclvalue(Env$l.var$titre)=="") {
		texte<-paste(texte,", ylab=\"",titre.y,"\"",sep="")
	    } else {
	      if (titre.x=="") {
		  texte<-paste(texte,",\n  ylab=\"",titre.y,"\"",sep="")
		} else {
		  texte<-paste(texte,", ylab=\"",titre.y,"\"",sep="")
		}
	    }
	  } else {
	    texte<-paste(texte,"ylab=\"",titre.y,"\"",sep="")
	  }
	}
	if (titre.x!="" | titre.y!="") {
	  if (tclvalue(Env$l.var$legendes.col)!="black" & tclvalue(Env$l.var$legendes.col)!="#000000") {texte<-paste(texte,", col.lab=\"",tclvalue(Env$l.var$legendes.col),"\"",sep="")}
	  if (tclvalue(Env$l.var$legendes.taille)!="1") {texte<-paste(texte,", cex.lab=",tclvalue(Env$l.var$legendes.taille),sep="")}
	}
	cat(paste(texte,")\n\n",sep=""))
    }
  }
}


#-------------------------------------------------
# Code - Graphe (position légende)
#-------------------------------------------------

code.graph.posleg<-function() {
  position<-if (tclvalue(Env$l.var$legende.pos)==Env$voc[102,1]) {"top"} else
    if (tclvalue(Env$l.var$legende.pos)==Env$voc[103,1]) {"topright"} else
    if (tclvalue(Env$l.var$legende.pos)==Env$voc[104,1]) {"left"} else
    if (tclvalue(Env$l.var$legende.pos)==Env$voc[105,1]) {"center"} else
    if (tclvalue(Env$l.var$legende.pos)==Env$voc[106,1]) {"right"} else
    if (tclvalue(Env$l.var$legende.pos)==Env$voc[107,1]) {"bottomleft"} else
    if (tclvalue(Env$l.var$legende.pos)==Env$voc[108,1]) {"bottom"} else
    if (tclvalue(Env$l.var$legende.pos)==Env$voc[109,1]) {"bottomright"} else
    {"topleft"}
  return(position)
}


#-------------------------------------------------
# Code - Graphe (histogramme)
#-------------------------------------------------

code.graph.hist<-function() {
  variable<-tclvalue(Env$l.var$variable)
  facteur1<-tclvalue(Env$l.var$facteur1)
  texte<-""
  if (tclvalue(Env$l.var$hist.type)==Env$voc[40,1]) {
    texte<-"barplot(frequencies, axes=FALSE, ann=FALSE, space=0"
    texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur1A),"\"",sep="")
    if (tclvalue(Env$l.var$col.borduresA)!="black" & tclvalue(Env$l.var$col.borduresA)!="#000000") {texte<-paste(texte,", border=\"",tclvalue(Env$l.var$col.borduresA),"\"",sep="")}
    texte<-paste(texte,",\n  xlim=c(",round(Env$l.code$x.inf,2),",",round(Env$l.code$x.sup,2),")",sep="")
    texte<-paste(texte,", ylim=c(0,",round(Env$l.code$y.sup,2),"))\n\n",sep="")
    cat(texte)
    code.graph.axes()
    code.graph.titre()
    if (tclvalue(Env$l.var$encadre)==1) {cat("box()\n\n")}
  } else {
    if (nchar(facteur1)>0 & facteur1!=Env$voc[82,1]) {
	texte<-paste(texte,"hist(variable, axes=FALSE, ann=FALSE",sep="")
    } else {
	texte<-paste(texte,"hist(",variable,", axes=FALSE, ann=FALSE",sep="")
    }
    texte<-paste(texte,", freq=",ifelse(tclvalue(Env$l.var$hist.type)==Env$voc[41,1],"TRUE","FALSE"),sep="")
    texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur1A),"\"",sep="")
    if (tclvalue(Env$l.var$col.borduresA)!="black" & tclvalue(Env$l.var$col.borduresA)!="#000000") {texte<-paste(texte,", border=\"",tclvalue(Env$l.var$col.borduresA),"\"",sep="")}
    if (tclvalue(Env$l.var$hist.barres)!="Auto") {texte<-paste(texte,", breaks=",as.numeric(tclvalue(Env$l.var$hist.barres))-1,sep="")}
    texte<-paste(texte,",\n  xlim=c(",round(Env$l.code$x.inf,2),",",round(Env$l.code$x.sup,2),")",sep="")
    texte<-paste(texte,", ylim=c(0,",round(Env$l.code$y.sup,2),"))\n\n",sep="")
    cat(texte)
    code.graph.axes()
    code.graph.titre()
    if (tclvalue(Env$l.var$encadre)==1) {cat("box()\n\n")}
    if (tclvalue(Env$l.var$hist.type)==Env$voc[42,1] & tclvalue(Env$l.var$hist.dens)==1) {
	texte<-""
	cat("# Distribution curve of the data\n\n")
	if (nchar(facteur1)>0 & facteur1!=Env$voc[82,1]) {
	  texte<-paste(texte,"lines(density(na.omit(variable))",sep="")
	} else {
	  texte<-paste(texte,"lines(density(na.omit(",variable,"))",sep="")
	}
	if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	if (tclvalue(Env$l.var$epaisseur1)!="1") {texte<-paste(texte,", lwd=",tclvalue(Env$l.var$epaisseur1),sep="")}
	texte<-paste(texte,", lty=",type.trait(type=tclvalue(Env$l.var$trait1)),")\n\n",sep="")
	cat(texte)
    }
  }
}


#-------------------------------------------------
# Code - Graphe (boîtes à moustaches)
#-------------------------------------------------

code.graph.moust<-function() {
  variable<-tclvalue(Env$l.var$variable)
  facteur1<-tclvalue(Env$l.var$facteur1)
  facteur2<-tclvalue(Env$l.var$facteur2)
  texte<-""
  if (nchar(facteur2)>0 & facteur2!=Env$voc[82,1]) {
    texte<-paste("boxplot(",variable," ~ interaction, axes=FALSE, ann=FALSE",sep="")
  } else {
    texte<-paste("boxplot(",variable," ~ ",facteur1,", axes=FALSE, ann=FALSE",sep="")
  }
  orient<-ifelse (tclvalue(Env$l.var$box.orient)==Env$voc[67,1],"hor","ver")
  log.axes<-if (orient=="ver" & tclvalue(Env$l.var$log.axevaleurs)==1) {
    "y"
  } else if (orient=="hor" & tclvalue(Env$l.var$log.axevaleurs)==1) {
    "x"
  } else {""}
  texte<-paste(texte,", horizontal=",ifelse(orient=="hor","TRUE","FALSE"),sep="")
  if (tclvalue(Env$l.var$plusieurs)==0) {
    texte<-paste(texte,",\n  col=\"",tclvalue(Env$l.var$couleur1A),"\"",sep="")
    texte<-paste(texte,", boxcol=\"",tclvalue(Env$l.var$col.borduresA),"\"",sep="")
    texte<-paste(texte,", medcol=\"",tclvalue(Env$l.var$col.borduresA),"\"",sep="")
    texte<-paste(texte,",\n  whiskcol=\"",tclvalue(Env$l.var$col.borduresA),"\"",sep="")
    texte<-paste(texte,", staplecol=\"",tclvalue(Env$l.var$col.borduresA),"\"",sep="")
    texte<-paste(texte,", whisklty=",type.trait(type=tclvalue(Env$l.var$trait1)),sep="")
    texte<-paste(texte,",\n  outcol=\"",tclvalue(Env$l.var$col.borduresA),"\"",sep="")
  } else {
    texte<-paste(texte,",\n  col=c(\"",paste(Env$l.var$couleur1B,collapse="\",\""),"\")",sep="")
    texte<-paste(texte,", boxcol=c(\"",paste(Env$l.var$col.borduresB,collapse="\",\""),"\")",sep="")
    texte<-paste(texte,", medcol=c(\"",paste(Env$l.var$col.borduresB,collapse="\",\""),"\")",sep="")
    texte<-paste(texte,",\n  whiskcol=c(\"",paste(Env$l.var$col.borduresB,collapse="\",\""),"\")",sep="")
    texte<-paste(texte,", staplecol=c(\"",paste(Env$l.var$col.borduresB,collapse="\",\""),"\")",sep="")
    texte<-paste(texte,", whisklty=",type.trait(type=tclvalue(Env$l.var$trait1)),sep="")
    texte<-paste(texte,",\n  outcol=c(\"",paste(Env$l.var$col.borduresB,collapse="\",\""),"\")",sep="")
  }
  texte<-paste(texte,", range=",as.numeric(tclvalue(Env$l.var$lg.moustaches)),sep="")
  texte<-paste(texte,", notch=",ifelse(tclvalue(Env$l.var$ICmediane)=="1","TRUE","FALSE"),sep="")
  if (log.axes!="") {texte<-paste(texte,", log=",log.axes,sep="")}
  texte<-paste(texte,",\n  ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),")",sep="")
  if (tclvalue(Env$l.var$varwidth)==1) {texte<-paste(texte,", varwidth=TRUE",sep="")}
  texte<-paste(texte,")\n\n",sep="")
  cat(texte)
  if (tclvalue(Env$l.var$boxmoy)==1) {
    facteur<-NULL
    if (nchar(facteur2)>0 & facteur2!=Env$voc[82,1]) {
	facteur <- factor(paste(Env$dataset[,tclvalue(Env$l.var$facteur1)],Env$dataset[,tclvalue(Env$l.var$facteur2)],sep=":"))
    } else {
	facteur <- Env$dataset[,tclvalue(Env$l.var$facteur1)]
    }
    texte<-""
    if (tclvalue(Env$l.var$box.orient)==Env$voc[67,1]) {
	texte<-paste("points(means",sep="")
	texte<-paste(texte,", 1:",nlevels(facteur),sep="")
    } else {
	texte<-paste("points(1:",nlevels(facteur),sep="")
	texte<-paste(texte,", means",sep="")
    }
    texte<-paste(texte,", cex=2",sep="")
    texte<-paste(texte,", col=c(\"",paste(Env$l.var$col.borduresB,collapse="\",\""),"\")",sep="")
    texte<-paste(texte,", pch=\"+\")\n\n",sep="")
    cat(texte)
  }
  code.graph.axes()
  code.graph.titre()
  if (tclvalue(Env$l.var$encadre)==1) {cat("box()\n\n")}
  if (tclvalue(Env$l.var$plusieurs)==1 & tclvalue(Env$l.var$legende)==1) {
    cat("# Legend\n\n")
    cat("par(xpd=TRUE)\n")
    texte<-paste("legend(\"",code.graph.posleg(),"\"",sep="")
    texte<-paste(texte,", legend=c(\"",paste(Env$l.var$noms2,collapse="\",\""),"\")",sep="")
    texte<-paste(texte,",\n  fill=c(\"",paste(Env$l.var$couleur1B,collapse="\",\""),"\")",sep="")
    if (nchar(tclvalue(Env$l.var$legende.titre))>0) {texte<-paste(texte,", title=\"",tclvalue(Env$l.var$legende.titre),"\"",sep="")}
    texte<-paste(texte,")\n",sep="")
    cat(texte)
    cat("par(xpd=FALSE)\n\n")
  }
}


#-------------------------------------------------
# Code - Graphe (barres - barres d'erreur)
#-------------------------------------------------

code.graph.barres.erreurs<-function(type) {
  if (tclvalue(Env$l.var$erreur)==Env$voc[96,1]) {
    cat("# Error bars\n\n")
    texte<-""
    texte<-paste("arrows(graph, ",type," - std.dev, graph, ",type," + std.dev, code=3, angle=90, length=0.1",sep="")
    if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}

    texte<-paste(texte,")\n\n",sep="")
    cat(texte)
  } else if (tclvalue(Env$l.var$erreur)==Env$voc[97,1]) {
    cat("# Error bars\n\n")
    texte<-""
    texte<-paste("arrows(graph, ",type," - std.err, graph, ",type," + std.err, code=3, angle=90, length=0.1",sep="")
    if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
    texte<-paste(texte,")\n\n",sep="")
    cat(texte)
  } else if (tclvalue(Env$l.var$erreur)==Env$voc[98,1]) {
    cat("# Error bars\n\n")
    texte<-""
    if (tclvalue(Env$l.var$moyprop)=="moy") {
	texte<-paste("arrows(graph, ",type," - ci, graph, ",type," + ci, code=3, angle=90, length=0.1",sep="")
	if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	texte<-paste(texte,")\n\n",sep="")
    } else {
	texte<-paste("arrows(graph, ",type," - ci.inf, graph, ",type," + ci.sup, code=3, angle=90, length=0.1",sep="")
	if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	texte<-paste(texte,")\n\n",sep="")
    }
    cat(texte)
  }
}


#-------------------------------------------------
# Code - Graphe (barres)
#-------------------------------------------------

code.graph.barres<-function() {
  texte<-""
  if (tclvalue(Env$l.var$moyprop)=="moy") {
    if (tclvalue(Env$l.var$plusieurs)==0) {
	texte<-"graph <- barplot(means"
	if (tclvalue(Env$l.var$nobar)==0) {
	  texte<-paste(texte,", axes=FALSE, ann=FALSE",sep="")
	  texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur1A),"\"",sep="")
	  if (tclvalue(Env$l.var$col.borduresA)!="black" & tclvalue(Env$l.var$col.borduresA)!="#000000") {texte<-paste(texte,", border=\"",tclvalue(Env$l.var$col.borduresA),"\"",sep="")}
	  if (graphe.log()!="") {texte<-paste(texte,", log=\"",graphe.log(),"\"",sep="")}
	  texte<-paste(texte,",\n  ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),")",sep="")
	  texte<-paste(texte,", names=c(\"",paste(Env$l.var$noms1,collapse="\",\""),"\"), xpd=FALSE)\n\n",sep="")
	  if (graphe.log()=="" & tclvalue(Env$l.var$hachuresA)!="1") {
	    hachures<-graphe.hachures(num=as.numeric(tclvalue(Env$l.var$hachuresA)))
	    texte<-paste(texte,"barplot(means, axes=FALSE, ann=FALSE",sep="")
	    texte<-paste(texte,", col=\"",tclvalue(Env$l.var$col.borduresA),"\"",sep="")
	    texte<-paste(texte,", border=\"",tclvalue(Env$l.var$col.borduresA),"\"",sep="")
	    texte<-paste(texte,", ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),")",sep="")
	    texte<-paste(texte,",\n  density=",hachures$densite,sep="")
	    texte<-paste(texte,", angle=",hachures$angle,sep="")
	    texte<-paste(texte,", names.arg=rep(\"\",",length(Env$l.var$noms1),")",sep="")
	    texte<-paste(texte,", xpd=FALSE, add=TRUE)\n\n",sep="")
	  }
	} else {
	  if (graphe.log()!="") {texte<-paste(texte,", log=\"",graphe.log(),"\"",sep="")}
	  texte<-paste(texte,", ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),")",sep="")
	  texte<-paste(texte,", plot=FALSE)\n",sep="")
	  texte<-paste(texte,"plot(graph, means, axes=FALSE, ann=FALSE, pch=16, cex=1.7",sep="")
	  texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur1A),"\"",sep="")
	  if (graphe.log()!="") {texte<-paste(texte,", log=\"",graphe.log(),"\"",sep="")}
	  texte<-paste(texte,",\n  xlim=c(",min(Env$l.var$add.abscisses)-0.5,",",max(Env$l.var$add.abscisses)+0.5,")",sep="")
	  texte<-paste(texte,", ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),"))\n\n",sep="")
	}
	cat(texte)
	code.graph.axes()
	code.graph.titre()
	if (tclvalue(Env$l.var$encadre)==1) {cat("box()\n\n")}
	if (tclvalue(Env$l.var$erreur)%in%Env$voc[96:98,1]) {
	  code.graph.barres.erreurs(type="means")
	  if (tclvalue(Env$l.var$nobar)==1) {
	    cat("# Points above error bars\n\n")
	    texte<-paste("points(graph, means, pch=16, cex=1.7, col=\"",tclvalue(Env$l.var$couleur1A),"\")\n\n",sep="")
	    cat(texte)
	  }
	}
    } else {
	texte<-"graph <- barplot(means"
	if (tclvalue(Env$l.var$nobar)==0) {
	  texte<-paste(texte,", axes=FALSE, ann=FALSE",sep="")
	  texte<-paste(texte,", col=c(\"",paste(Env$l.var$couleur1B,collapse="\",\""),"\")",sep="")
	  if (any(!Env$l.var$col.borduresB%in%c("black","#000000"))) {texte<-paste(texte,", border=c(\"",paste(Env$l.var$col.borduresB,collapse="\",\""),"\")",sep="")}
	  if (graphe.log()!="") {texte<-paste(texte,", log=\"",graphe.log(),"\"",sep="")}
	  texte<-paste(texte,",\n  ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),")",sep="")
	  texte<-paste(texte,", beside=",ifelse(tclvalue(Env$l.var$stack)==1,"FALSE","TRUE"),sep="")
	  texte<-paste(texte,", names=c(\"",paste(Env$l.var$noms1,collapse="\",\""),"\"), xpd=FALSE)\n\n",sep="")
	  if (graphe.log()=="" & any(Env$l.var$hachuresB!=1)) {
	    hachures<-graphe.hachures(num=Env$l.var$hachuresB)
	    texte<-paste(texte,"barplot(means, axes=FALSE, ann=FALSE",sep="")
	    texte<-paste(texte,", col=c(\"",paste(Env$l.var$col.borduresB,collapse="\",\""),"\")",sep="")
	    texte<-paste(texte,", border=c(\"",paste(Env$l.var$col.borduresB,collapse="\",\""),"\")",sep="")
	    texte<-paste(texte,",\n  ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),")",sep="")
	    texte<-paste(texte,", beside=",ifelse(tclvalue(Env$l.var$stack)==1,"FALSE","TRUE"),sep="")
	    texte<-paste(texte,", density=c(",paste(hachures$densite,collapse=","),")",sep="")
	    texte<-paste(texte,", angle=c(",paste(hachures$angle,collapse=","),")",sep="")
	    texte<-paste(texte,",\n  names.arg=rep(\"\",",length(Env$l.var$noms1),")",sep="")
	    texte<-paste(texte,", xpd=FALSE, add=TRUE)\n\n",sep="")
	  }
	} else {
	  if (graphe.log()!="") {texte<-paste(texte,", log=\"",graphe.log(),"\"",sep="")}
	  texte<-paste(texte,", ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),")",sep="")
	  texte<-paste(texte,", beside=",ifelse(tclvalue(Env$l.var$stack)==1,"FALSE","TRUE"),sep="")
	  texte<-paste(texte,", plot=FALSE)\n",sep="")
	  texte<-paste(texte,"plot(",sep="")
	  if(tclvalue(Env$l.var$stack)==1){
	    texte<-paste(texte,"c(",paste(rep(Env$l.var$add.abscisses,each=length(Env$l.var$add.abscisses)),collapse=","),")",sep="")
	  } else {
	    texte<-paste(texte,"graph",sep="")
	  }
	  texte<-paste(texte,", means, axes=FALSE, ann=FALSE, pch=16, cex=1.7",sep="")
	  texte<-paste(texte,", col=c(\"",paste(Env$l.var$couleur1B,collapse="\",\""),"\")",sep="")
	  if (graphe.log()!="") {texte<-paste(texte,", log=\"",graphe.log(),"\"",sep="")}
	  texte<-paste(texte,",\n  xlim=c(",min(Env$l.var$add.abscisses)-0.5,",",max(Env$l.var$add.abscisses)+0.5,")",sep="")
	  texte<-paste(texte,", ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),"))\n\n",sep="")
	}
	cat(texte)
	code.graph.axes()
	code.graph.titre()
	if (tclvalue(Env$l.var$encadre)==1) {cat("box()\n\n")}
	if (tclvalue(Env$l.var$stack)==0 & tclvalue(Env$l.var$erreur)%in%Env$voc[96:98,1]) {
	  code.graph.barres.erreurs(type="means")
	  if (tclvalue(Env$l.var$nobar)==1) {
	    cat("# Points above error bars\n\n")
	    texte<-paste("points(graph, means, pch=16, cex=1.7, col=c(\"",paste(Env$l.var$couleur1B,collapse="\",\""),"\"))\n\n",sep="")
	    cat(texte)
	  }
	}
	if (tclvalue(Env$l.var$legende)==1) {
	  cat("# Legend\n\n")
	  cat("par(xpd=TRUE)\n")
	  texte<-paste("legend(\"",code.graph.posleg(),"\"",sep="")
	  texte<-paste(texte,", legend=c(\"",paste(Env$l.var$noms2,collapse="\",\""),"\")",sep="")
	  if (tclvalue(Env$l.var$nobar)==0) {
	    texte<-paste(texte,",\n  fill=c(\"",paste(Env$l.var$couleur1B,collapse="\",\""),"\")",sep="")
	  } else {
	    texte<-paste(texte,",\n  pch=16, pt.cex=1.7",sep="")
	    texte<-paste(texte,", col=c(\"",paste(Env$l.var$couleur1B,collapse="\",\""),"\")",sep="")
	  }
	  if (nchar(tclvalue(Env$l.var$legende.titre))>0) {texte<-paste(texte,", title=\"",tclvalue(Env$l.var$legende.titre),"\"",sep="")}
	  texte<-paste(texte,")\n",sep="")
	  cat(texte)
	  cat("par(xpd=FALSE)\n\n")
	}
    }
  } else {
    if (tclvalue(Env$l.var$plusieurs)==0) {
	texte<-"graph <- barplot(proportions"
	if (tclvalue(Env$l.var$nobar)==0) {
	  texte<-paste(texte,", axes=FALSE, ann=FALSE",sep="")
	  texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur1A),"\"",sep="")
	  if (tclvalue(Env$l.var$col.borduresA)!="black" & tclvalue(Env$l.var$col.borduresA)!="#000000") {texte<-paste(texte,", border=\"",tclvalue(Env$l.var$col.borduresA),"\"",sep="")}
	  if (graphe.log()!="") {texte<-paste(texte,", log=\"",graphe.log(),"\"",sep="")}
	  texte<-paste(texte,",\n  ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),")",sep="")
	  texte<-paste(texte,", names=c(\"",paste(Env$l.var$nomsprop.fac,collapse="\",\""),"\"), xpd=FALSE)\n\n",sep="")
	  if (graphe.log()=="" & tclvalue(Env$l.var$hachuresA)!="1") {
	    hachures<-graphe.hachures(num=as.numeric(tclvalue(Env$l.var$hachuresA)))
	    texte<-paste(texte,"barplot(proportions, axes=FALSE, ann=FALSE",sep="")
	    texte<-paste(texte,", col=\"",tclvalue(Env$l.var$col.borduresA),"\"",sep="")
	    texte<-paste(texte,", border=\"",tclvalue(Env$l.var$col.borduresA),"\"",sep="")
	    texte<-paste(texte,", ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),")",sep="")
	    texte<-paste(texte,",\n  density=",hachures$densite,sep="")
	    texte<-paste(texte,", angle=",hachures$angle,sep="")
	    texte<-paste(texte,", names.arg=rep(\"\",",length(Env$l.var$nomsprop.fac),")",sep="")
	    texte<-paste(texte,", xpd=FALSE, add=TRUE)\n\n",sep="")
	  }
	} else {
	  if (graphe.log()!="") {texte<-paste(texte,", log=\"",graphe.log(),"\"",sep="")}
	  texte<-paste(texte,", ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),")",sep="")
	  texte<-paste(texte,", plot=FALSE)\n",sep="")
	  texte<-paste(texte,"plot(graph, proportions, axes=FALSE, ann=FALSE, pch=16, cex=1.7",sep="")
	  texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur1A),"\"",sep="")
	  if (graphe.log()!="") {texte<-paste(texte,", log=\"",graphe.log(),"\"",sep="")}
	  texte<-paste(texte,",\n  xlim=c(",min(Env$l.var$add.abscisses)-0.5,",",max(Env$l.var$add.abscisses)+0.5,")",sep="")
	  texte<-paste(texte,", ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),"))\n\n",sep="")
	}
	cat(texte)
	code.graph.axes()
	code.graph.titre()
	if (tclvalue(Env$l.var$encadre)==1) {cat("box()\n\n")}
	if (tclvalue(Env$l.var$erreur)%in%Env$voc[97:98,1]) {
	  code.graph.barres.erreurs(type="proportions")
	  if (tclvalue(Env$l.var$nobar)==1) {
	    cat("# Points above error bars\n\n")
	    texte<-paste("points(graph, proportions, pch=16, cex=1.7, col=\"",tclvalue(Env$l.var$couleur1A),"\")\n\n",sep="")
	    cat(texte)
	  }
	}
    } else {
	texte<-"graph <- barplot(proportions"
	if (tclvalue(Env$l.var$nobar)==0) {
	  texte<-paste(texte,", axes=FALSE, ann=FALSE",sep="")
	  texte<-paste(texte,", col=c(\"",paste(Env$l.var$couleur1B,collapse="\",\""),"\")",sep="")
	  if (any(!Env$l.var$col.borduresB%in%c("black","#000000"))) {texte<-paste(texte,", border=c(\"",paste(Env$l.var$col.borduresB,collapse="\",\""),"\")",sep="")}
	  if (graphe.log()!="") {texte<-paste(texte,", log=\"",graphe.log(),"\"",sep="")}
	  texte<-paste(texte,",\n  ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),")",sep="")
	  texte<-paste(texte,", beside=",ifelse(tclvalue(Env$l.var$stack)==1,"FALSE","TRUE"),sep="")
	  texte<-paste(texte,", names=c(\"",paste(Env$l.var$nomsprop.fac,collapse="\",\""),"\"), xpd=FALSE)\n\n",sep="")
	  if (graphe.log()=="" & any(Env$l.var$hachuresB!=1)) {
	    hachures<-graphe.hachures(num=Env$l.var$hachuresB)
	    texte<-paste(texte,"barplot(proportions, axes=FALSE, ann=FALSE",sep="")
	    texte<-paste(texte,", col=c(\"",paste(Env$l.var$col.borduresB,collapse="\",\""),"\")",sep="")
	    texte<-paste(texte,", border=c(\"",paste(Env$l.var$col.borduresB,collapse="\",\""),"\")",sep="")
	    texte<-paste(texte,",\n  ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),")",sep="")
	    texte<-paste(texte,", beside=",ifelse(tclvalue(Env$l.var$stack)==1,"FALSE","TRUE"),sep="")
	    texte<-paste(texte,", density=c(",paste(hachures$densite,collapse=","),")",sep="")
	    texte<-paste(texte,", angle=c(",paste(hachures$angle,collapse=","),")",sep="")
	    texte<-paste(texte,",\n  names.arg=rep(\"\",",length(Env$l.var$nomsprop.fac),")",sep="")
	    texte<-paste(texte,", xpd=FALSE, add=TRUE)\n\n",sep="")
	  }
	} else {
	  if (graphe.log()!="") {texte<-paste(texte,", log=\"",graphe.log(),"\"",sep="")}
	  texte<-paste(texte,", ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),")",sep="")
	  texte<-paste(texte,", beside=",ifelse(tclvalue(Env$l.var$stack)==1,"FALSE","TRUE"),sep="")
	  texte<-paste(texte,", plot=FALSE)\n",sep="")
	  texte<-paste(texte,"plot(",sep="")
	  if(tclvalue(Env$l.var$stack)==1){
	    texte<-paste(texte,"c(",paste(rep(Env$l.var$add.abscisses,each=length(Env$l.var$add.abscisses)),collapse=","),")",sep="")
	  } else {
	    texte<-paste(texte,"graph",sep="")
	  }
	  texte<-paste(texte,", proportions, axes=FALSE, ann=FALSE, pch=16, cex=1.7",sep="")
	  texte<-paste(texte,", col=c(\"",paste(Env$l.var$couleur1B,collapse="\",\""),"\")",sep="")
	  if (graphe.log()!="") {texte<-paste(texte,", log=\"",graphe.log(),"\"",sep="")}
	  texte<-paste(texte,",\n  xlim=c(",min(Env$l.var$add.abscisses)-0.5,",",max(Env$l.var$add.abscisses)+0.5,")",sep="")
	  texte<-paste(texte,", ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),"))\n\n",sep="")
	}
	cat(texte)
	code.graph.axes()
	code.graph.titre()
	if (tclvalue(Env$l.var$encadre)==1) {cat("box()\n\n")}
	if (tclvalue(Env$l.var$stack)==0 & tclvalue(Env$l.var$erreur)%in%Env$voc[97:98,1]) {
	  code.graph.barres.erreurs(type="proportions")
	  if (tclvalue(Env$l.var$nobar)==1) {
	    cat("# Points above error bars\n\n")
	    texte<-paste("points(graph, proportions, pch=16, cex=1.7, col=c(\"",paste(Env$l.var$couleur1B,collapse="\",\""),"\"))\n\n",sep="")
	    cat(texte)
	  }
	}
	if (tclvalue(Env$l.var$legende)==1) {
	  cat("# Legend\n\n")
	  cat("par(xpd=TRUE)\n")
	  texte<-paste("legend(\"",code.graph.posleg(),"\"",sep="")
	  texte<-paste(texte,", legend=c(\"",paste(Env$l.var$nomsprop,collapse="\",\""),"\")",sep="")
	  if (tclvalue(Env$l.var$nobar)==0) {
	    texte<-paste(texte,",\n  fill=c(\"",paste(Env$l.var$couleur1B,collapse="\",\""),"\")",sep="")
	  } else {
	    texte<-paste(texte,",\n  pch=16, pt.cex=1.7",sep="")
	    texte<-paste(texte,", col=c(\"",paste(Env$l.var$couleur1B,collapse="\",\""),"\")",sep="")
	  }
	  if (nchar(tclvalue(Env$l.var$legende.titre))>0) {texte<-paste(texte,", title=\"",tclvalue(Env$l.var$legende.titre),"\"",sep="")}
	  texte<-paste(texte,")\n",sep="")
	  cat(texte)
	  cat("par(xpd=FALSE)\n\n")
	}
    }
  }
}


#-------------------------------------------------
# Code - Graphe (camembert)
#-------------------------------------------------

code.graph.cam<-function() {
  texte<-"pie(variable"
  texte<-paste(texte,", col=c(\"",paste(Env$l.var$couleur1B,collapse="\",\""),"\")",sep="")
  if (tclvalue(Env$l.var$col.borduresA)!="black" & tclvalue(Env$l.var$col.borduresA)!="#000000") {texte<-paste(texte,", border=\"",tclvalue(Env$l.var$col.borduresA),"\"",sep="")}
  if (tclvalue(Env$l.var$cam.lien)==1) {
    texte<-paste(texte,",\n  labels=c(\"",paste(Env$l.var$nomsparts,collapse="\",\""),"\")",sep="")
  } else {
    texte<-paste(texte,",\n  labels=\"\"",sep="")
  }
  texte<-paste(texte,", clockwise=",ifelse(tclvalue(Env$l.var$cam.orient)==Env$voc[120,1],"TRUE","FALSE"),sep="")
  texte<-paste(texte,", init.angle=",ifelse(tclvalue(Env$l.var$cam.orient)==Env$voc[120,1],as.numeric(tclvalue(Env$l.var$cam.start))+90,
    -1*as.numeric(tclvalue(Env$l.var$cam.start))+90),")\n",sep="")
  if (any(Env$l.var$hachuresB!=1)) {
    hachures<-graphe.hachures(num=Env$l.var$hachuresB)
    texte<-paste(texte,"par(new=TRUE)\n",sep="")
    texte<-paste(texte,"pie(variable",sep="")
    if (tclvalue(Env$l.var$col.borduresA)!="black" & tclvalue(Env$l.var$col.borduresA)!="#000000") {
	texte<-paste(texte,", col=\"",tclvalue(Env$l.var$col.borduresA),"\"",sep="")
	texte<-paste(texte,", border=\"",tclvalue(Env$l.var$col.borduresA),"\"",sep="")
    }
    paste(texte,", labels=\"\"",sep="")
    texte<-paste(texte,", clockwise=",ifelse(tclvalue(Env$l.var$cam.orient)==Env$voc[120,1],"TRUE","FALSE"),sep="")
    texte<-paste(texte,", init.angle=",ifelse(tclvalue(Env$l.var$cam.orient)==Env$voc[120,1],as.numeric(tclvalue(Env$l.var$cam.start))+90,
	-1*as.numeric(tclvalue(Env$l.var$cam.start))+90),sep="")
    texte<-paste(texte,",\n  density=c(",paste(hachures$densite,collapse=","),")",sep="")
    texte<-paste(texte,", angle=c(",paste(hachures$angle,collapse=","),"))\n\n",sep="")
  }
  cat(texte)
  code.graph.titre()
  if (tclvalue(Env$l.var$cam.lien)==0 & tclvalue(Env$l.var$legende)==1) {
    cat("# Legend\n\n")
    cat("par(xpd=TRUE)\n")
    texte<-paste("legend(\"",code.graph.posleg(),"\"",sep="")
    texte<-paste(texte,", legend=c(\"",paste(Env$l.var$nomsparts,collapse="\",\""),"\")",sep="")
    texte<-paste(texte,",\n  fill=c(\"",paste(Env$l.var$couleur1B,collapse="\",\""),"\")",sep="")
    if (nchar(tclvalue(Env$l.var$legende.titre))>0) {texte<-paste(texte,", title=\"",tclvalue(Env$l.var$legende.titre),"\"",sep="")}
    paste(texte,", labels=\"\")\n",sep="")
    cat(texte)
    cat("par(xpd=FALSE)\n\n")
  }
}


#-------------------------------------------------
# Code - Graphe (courbe - barres d'erreur un)
#-------------------------------------------------

code.graph.courbe.erreurs1<-function(type) {
  if (tclvalue(Env$l.var$erreur)==Env$voc[96,1]) {
    cat("# Error bars\n\n")
    texte<-paste("width <- 0.015 * (",Env$l.code$x.sup," - ",Env$l.code$x.inf,")\n",sep="")
    texte<-paste(texte,"arrows(abscissae, ",type," - std.dev, abscissae, ",type," + std.dev, code=3, angle=90, length=width",sep="")
    if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
    texte<-paste(texte,")\n\n",sep="")
    cat(texte)
  } else if (tclvalue(Env$l.var$erreur)==Env$voc[97,1]) {
    cat("# Error bars\n\n")
    texte<-paste("width <- 0.015 * (",Env$l.code$x.sup," - ",Env$l.code$x.inf,")\n",sep="")
    texte<-paste(texte,"arrows(abscissae, ",type," - std.err, abscissae, ",type," + std.err, code=3, angle=90, length=width",sep="")
    if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
    texte<-paste(texte,")\n\n",sep="")
    cat(texte)
  } else if (tclvalue(Env$l.var$erreur)==Env$voc[98,1]) {
    cat("# Error bars\n\n")
    texte<-paste("width <- 0.015 * (",Env$l.code$x.sup," - ",Env$l.code$x.inf,")\n",sep="")
    if (tclvalue(Env$l.var$moyprop)=="moy") {
	texte<-paste(texte,"arrows(abscissae, ",type," - ci, abscissae, ",type," + ci, code=3, angle=90, length=width",sep="")
	if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	texte<-paste(texte,")\n\n",sep="")
    } else {
	texte<-paste(texte,"arrows(abscissae, ",type," - ci.inf, abscissae, ",type," + ci.sup, code=3, angle=90, length=width",sep="")
	if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	texte<-paste(texte,")\n\n",sep="")
    }
    cat(texte)
  }
}


#-------------------------------------------------
# Code - Graphe (courbe - barres d'erreur plusieurs)
#-------------------------------------------------

code.graph.courbe.erreurs2<-function(type) {
  if (tclvalue(Env$l.var$erreur)==Env$voc[96,1]) {
    cat("# Error bars\n\n")
    texte<-paste("width <- 0.015 * (",Env$l.code$x.sup," - ",Env$l.code$x.inf,")\n",sep="")
    cat(texte)
    for (i in 1:length(Env$l.var$noms1)) {
	texte<-paste("arrows(abscissae, ",type,"[",i,",] - std.dev[",i,",], abscissae, ",type,"[",i,",] + std.dev[",i,",], code=3, angle=90, length=width",sep="")
	if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	texte<-paste(texte,")\n",sep="")
	cat(texte)
    }
    cat("\n")
  } else if (tclvalue(Env$l.var$erreur)==Env$voc[97,1]) {
    cat("# Error bars\n\n")
    texte<-paste("width <- 0.015 * (",Env$l.code$x.sup," - ",Env$l.code$x.inf,")\n",sep="")
    cat(texte)
    for (i in 1:length(Env$l.var$noms1)) {
	texte<-paste("arrows(abscissae, ",type,"[",i,",] - std.err[",i,",], abscissae, ",type,"[",i,",] + std.err[",i,",], code=3, angle=90, length=width",sep="")
	if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	texte<-paste(texte,")\n",sep="")
	cat(texte)
    }
    cat("\n")
  } else if (tclvalue(Env$l.var$erreur)==Env$voc[98,1]) {
    cat("# Error bars\n\n")
    texte<-paste("width <- 0.015 * (",Env$l.code$x.sup," - ",Env$l.code$x.inf,")\n",sep="")
    cat(texte)
    if (tclvalue(Env$l.var$moyprop)=="moy") {
	for (i in 1:length(Env$l.var$noms1)) {
	  texte<-paste("arrows(abscissae, ",type,"[",i,",] - ci[",i,",], abscissae, ",type,"[",i,",] + ci[",i,",], code=3, angle=90, length=width",sep="")
	  if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	  texte<-paste(texte,")\n",sep="")
	  cat(texte)
	}
    } else {
	for (i in 1:length(Env$l.var$noms1)) {
	  texte<-paste("arrows(abscissae, ",type,"[",i,",] - ci.inf[",i,",], abscissae, ",type,"[",i,",] + ci.sup[",i,",], code=3, angle=90, length=width",sep="")
	  if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	  texte<-paste(texte,")\n",sep="")
	  cat(texte)
	}
    }
    cat("\n")
  }
}


#-------------------------------------------------
# Code - Graphe (courbe - légende)
#-------------------------------------------------

code.graph.courbe.leg<-function() {
  symboles<-graphe.symboles(num=Env$l.var$symboleB)
  lignes<-type.ligne(type=Env$l.var$type.courbeB)
  traits<-type.trait(type=Env$l.var$trait2)
  epaisseur<-Env$l.var$epaisseur2
  taille.pts<-Env$l.var$taille.ptsB
  if (any(lignes=="p")) {
    traits[which(lignes=="p")]<-NA
    epaisseur[which(lignes=="p")]<-NA
  }
  if (any(lignes%in%c("l","h"))) {
    symboles[which(lignes%in%c("l","h"))]<-NA
    taille.pts[which(lignes%in%c("l","h"))]<-NA
  }
  cat("# Legend\n\n")
  cat("par(xpd=TRUE)\n")
  texte<-paste("legend(\"",code.graph.posleg(),"\"",sep="")
  texte<-paste(texte,", legend=c(\"",paste(Env$l.var$noms1,collapse="\",\""),"\")",sep="")
  texte<-paste(texte,", col=c(\"",paste(Env$l.var$couleur2B,collapse="\",\""),"\")",sep="")
  texte<-paste(texte,",\n  pch=c(",paste(symboles,collapse=","),")",sep="")
  texte<-paste(texte,", pt.cex=c(",paste(taille.pts,collapse=","),")",sep="")
  texte<-paste(texte,", lty=c(",paste(traits,collapse=","),")",sep="")
  texte<-paste(texte,", lwd=c(",paste(epaisseur,collapse=","),")",sep="")
  if (nchar(tclvalue(Env$l.var$legende.titre))>0) {texte<-paste(texte,", title=\"",tclvalue(Env$l.var$legende.titre),"\"",sep="")}
  cat(paste(texte,")\n",sep=""))
  cat("par(xpd=FALSE)\n\n")
}


#-------------------------------------------------
# Code - Graphe (courbe)
#-------------------------------------------------

code.graph.courbe<-function() {
  texte<-""
  if (tclvalue(Env$l.var$moyprop)=="moy") {
    if (tclvalue(Env$l.var$plusieurs)==0) {
	cat("abscissae <- as.numeric(names(means))\n")
	texte<-"plot(means ~ abscissae, axes=FALSE, ann=FALSE"
	texte<-paste(texte,", xlim=c(",round(Env$l.code$x.inf,2),",",round(Env$l.code$x.sup,2),")",sep="")
	texte<-paste(texte,", ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),")",sep="")
	if (graphe.log()!="") {texte<-paste(texte,", log=\"",graphe.log(),"\"",sep="")}
	if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	texte<-paste(texte,",\n  pch=",graphe.symboles(num=as.numeric(tclvalue(Env$l.var$symboleA))),sep="")
	if (tclvalue(Env$l.var$taille.ptsA)!="1") {texte<-paste(texte,", cex=",as.numeric(tclvalue(Env$l.var$taille.ptsA)),sep="")}
	texte<-paste(texte,", type=\"",type.ligne(type=tclvalue(Env$l.var$type.courbeA)),"\"",sep="")
	if (type.ligne(type=tclvalue(Env$l.var$type.courbeA))!="p") {
	  texte<-paste(texte,", lty=",type.trait(type=tclvalue(Env$l.var$trait1)),sep="")
	  texte<-paste(texte,", lwd=",as.numeric(tclvalue(Env$l.var$epaisseur1)),sep="")
	}
	cat(paste(texte,")\n\n",sep=""))
	code.graph.axes()
	code.graph.titre()
	if (tclvalue(Env$l.var$encadre)==1) {cat("box()\n\n")}
	code.graph.courbe.erreurs1(type="means")
    } else {
	symboles<-graphe.symboles(num=Env$l.var$symboleB)
	lignes<-type.ligne(type=Env$l.var$type.courbeB)
	traits<-type.trait(type=Env$l.var$trait2)
	cat("abscissae <- as.numeric(colnames(means))\n")
	texte<-"plot(means[1,] ~ abscissae, axes=FALSE, ann=FALSE"
	texte<-paste(texte,", xlim=c(",round(Env$l.code$x.inf,2),",",round(Env$l.code$x.sup,2),")",sep="")
	texte<-paste(texte,", ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),")",sep="")
	if (graphe.log()!="") {texte<-paste(texte,", log=\"",graphe.log(),"\"",sep="")}
	texte<-paste(texte,", col=\"",Env$l.var$couleur2B[1],"\"",sep="")
	texte<-paste(texte,",\n  pch=",symboles[1],sep="")
	if (Env$l.var$taille.ptsB[1]!="1") {texte<-paste(texte,", cex=",as.numeric(Env$l.var$taille.ptsB[1]),sep="")}
	texte<-paste(texte,", type=\"",lignes[1],"\"",sep="")
	if (lignes[1]!="p") {
	  texte<-paste(texte,", lty=",traits[1],sep="")
	  texte<-paste(texte,", lwd=",as.numeric(Env$l.var$epaisseur2[1]),sep="")
	}
	cat(paste(texte,")\n",sep=""))
	for (i in 2:length(Env$l.var$noms1)) {
	  texte<-paste("lines(abscissae, means[",i,",]",sep="")
	  texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")
	  texte<-paste(texte,", pch=",symboles[i],sep="")
	  if (Env$l.var$taille.ptsB[i]!="1") {texte<-paste(texte,", cex=",as.numeric(Env$l.var$taille.ptsB[i]),sep="")}
	  texte<-paste(texte,", type=\"",lignes[i],"\"",sep="")
	  if (lignes[i]!="p") {
	    texte<-paste(texte,", lty=",traits[i],sep="")
	    texte<-paste(texte,", lwd=",as.numeric(Env$l.var$epaisseur2[i]),sep="")
	  }
	  cat(paste(texte,")\n",sep=""))
	}
	cat("\n")
	code.graph.axes()
	code.graph.titre()
	if (tclvalue(Env$l.var$encadre)==1) {cat("box()\n\n")}
	code.graph.courbe.erreurs2(type="means")
	if (tclvalue(Env$l.var$legende)==1) {code.graph.courbe.leg()}
    }
  } else {
    if (tclvalue(Env$l.var$plusieurs)==0) {
	cat("abscissae <- as.numeric(as.character(levels(varX)))\n")
	texte<-"plot(proportions ~ abscissae, axes=FALSE, ann=FALSE"
	texte<-paste(texte,", xlim=c(",round(Env$l.code$x.inf,2),",",round(Env$l.code$x.sup,2),")",sep="")
	texte<-paste(texte,", ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),")",sep="")
	if (graphe.log()!="") {texte<-paste(texte,", log=\"",graphe.log(),"\"",sep="")}
	if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	texte<-paste(texte,",\n  pch=",graphe.symboles(num=as.numeric(tclvalue(Env$l.var$symboleA))),sep="")
	if (tclvalue(Env$l.var$taille.ptsA)!="1") {texte<-paste(texte,", cex=",as.numeric(tclvalue(Env$l.var$taille.ptsA)),sep="")}
	texte<-paste(texte,", type=\"",type.ligne(type=tclvalue(Env$l.var$type.courbeA)),"\"",sep="")
	if (type.ligne(type=tclvalue(Env$l.var$type.courbeA))!="p") {
	  texte<-paste(texte,", lty=",type.trait(type=tclvalue(Env$l.var$trait1)),sep="")
	  texte<-paste(texte,", lwd=",as.numeric(tclvalue(Env$l.var$epaisseur1)),sep="")
	}
	cat(paste(texte,")\n\n",sep=""))
	code.graph.axes()
	code.graph.titre()
	if (tclvalue(Env$l.var$encadre)==1) {cat("box()\n\n")}
	code.graph.courbe.erreurs1(type="proportions")
    } else {
	symboles<-graphe.symboles(num=Env$l.var$symboleB)
	lignes<-type.ligne(type=Env$l.var$type.courbeB)
	traits<-type.trait(type=Env$l.var$trait2)
	cat("abscissae <- as.numeric(as.character(levels(varX)))\n")
	texte<-"plot(proportions[1,] ~ abscissae, axes=FALSE, ann=FALSE"
	texte<-paste(texte,", xlim=c(",round(Env$l.code$x.inf,2),",",round(Env$l.code$x.sup,2),")",sep="")
	texte<-paste(texte,", ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),")",sep="")
	if (graphe.log()!="") {texte<-paste(texte,", log=\"",graphe.log(),"\"",sep="")}
	texte<-paste(texte,", col=\"",Env$l.var$couleur2B[1],"\"",sep="")
	texte<-paste(texte,",\n  pch=",symboles[1],sep="")
	if (Env$l.var$taille.ptsB[1]!="1") {texte<-paste(texte,", cex=",as.numeric(Env$l.var$taille.ptsB[1]),sep="")}
	texte<-paste(texte,", type=\"",lignes[1],"\"",sep="")
	if (lignes[1]!="p") {
	  texte<-paste(texte,", lty=",traits[1],sep="")
	  texte<-paste(texte,", lwd=",as.numeric(Env$l.var$epaisseur2[1]),sep="")
	}
	cat(paste(texte,")\n",sep=""))
	for (i in 2:length(Env$l.var$noms1)) {
	  texte<-paste("lines(abscissae, proportions[",i,",]",sep="")
	  texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")
	  texte<-paste(texte,", pch=",symboles[i],sep="")
	  if (Env$l.var$taille.ptsB[i]!="1") {texte<-paste(texte,", cex=",as.numeric(Env$l.var$taille.ptsB[i]),sep="")}
	  texte<-paste(texte,", type=\"",lignes[i],"\"",sep="")
	  if (lignes[i]!="p") {
	    texte<-paste(texte,", lty=",traits[i],sep="")
	    texte<-paste(texte,", lwd=",as.numeric(Env$l.var$epaisseur2[i]),sep="")
	  }
	  cat(paste(texte,")\n",sep=""))
	}
	cat("\n")
	code.graph.axes()
	code.graph.titre()
	if (tclvalue(Env$l.var$encadre)==1) {cat("box()\n\n")}
	code.graph.courbe.erreurs2(type="proportions")
	if (tclvalue(Env$l.var$legende)==1) {code.graph.courbe.leg()}
    }
  }
}


#-------------------------------------------------
# Code - Graphe (nuage)
#-------------------------------------------------

code.graph.nuage<-function() {
  varX<-tclvalue(Env$l.var$varX)
  varY<-tclvalue(Env$l.var$varY)
  facteur1<-tclvalue(Env$l.var$facteur1)
  texte<-""
  if (tclvalue(Env$l.var$plusieurs)==0) {
    if (nchar(facteur1)>0 & facteur1!=Env$voc[82,1]) {
	texte<-"plot(varY ~ varX, axes=FALSE, ann=FALSE"
    } else {
	texte<-paste("plot(",varY," ~ ",varX,", axes=FALSE, ann=FALSE",sep="")
    }
    texte<-paste(texte,", xlim=c(",round(Env$l.code$x.inf,2),",",round(Env$l.code$x.sup,2),")",sep="")
    texte<-paste(texte,", ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),")",sep="")
    if (graphe.log()!="") {texte<-paste(texte,", log=\"",graphe.log(),"\"",sep="")}
    if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
    texte<-paste(texte,",\n  pch=",graphe.symboles(num=as.numeric(tclvalue(Env$l.var$symboleA))),sep="")
    if (tclvalue(Env$l.var$taille.ptsA)!="1") {texte<-paste(texte,", cex=",as.numeric(tclvalue(Env$l.var$taille.ptsA)),sep="")}
    cat(paste(texte,")\n\n",sep=""))
    if (tclvalue(Env$l.var$ptlab)==1) {
	texte<-""
	if (nchar(facteur1)>0 & facteur1!=Env$voc[82,1]) {
	  texte<-"text(varX, varY, pos=3, offset=0.4"
	} else {
	  texte<-paste("text(",varX,", ",varY,", pos=3, offset=0.4",sep="")
	}
	if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	texte<-paste(texte,", cex=",0.65*as.numeric(tclvalue(Env$l.var$taille.ptsA)),")\n\n",sep="")
	cat(texte)
    }
    code.graph.axes()
    code.graph.titre()
    if (tclvalue(Env$l.var$encadre)==1) {cat("box()\n\n")}
    if (tclvalue(Env$l.var$droiteA)==Env$voc[145,1]) {
	cat("# Regression line\n\n")
	texte<-""
	if (nchar(facteur1)>0 & facteur1!=Env$voc[82,1]) {
	  texte<-"model <- lm(varY ~ varX)\n"
	} else {
	  texte<-paste("model <- lm(",varY," ~ ",varX,")\n",sep="")
	}
	texte<-paste(texte,"abline(model$coefficients",sep="")
	if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	texte<-paste(texte,", lty=",type.trait(type=tclvalue(Env$l.var$trait1)),sep="")
	texte<-paste(texte,", lwd=",as.numeric(tclvalue(Env$l.var$epaisseur1)),sep="")
	cat(paste(texte,")\n",sep=""))
	if (tclvalue(Env$l.var$intervalA)==Env$voc[261,1]) {
	  if (nchar(facteur1)>0 & facteur1!=Env$voc[82,1]) {
	    texte<-"varX2 <- seq(min(varX, na.rm=TRUE), max(varX, na.rm=TRUE), abs(max(varX, na.rm=TRUE) - min(varX, na.rm=TRUE)) / 1000)\n"
	    texte<-paste(texte,"varY2 <- predict(model, list(varX=varX2), interval=\"confidence\")\n",sep="")
	  } else {
	    texte<-paste("varX2 <- seq(min(",varX,", na.rm=TRUE), max(",varX,", na.rm=TRUE), abs(max(",varX,", na.rm=TRUE) - min(",varX,", na.rm=TRUE)) / 1000)\n",sep="")
	    texte<-paste(texte,"varY2 <- predict(model, list(",varX,"=varX2), interval=\"confidence\")\n",sep="")
	  }
	  texte<-paste(texte,"lines(varX2, varY2[,\"lwr\"], lty=2",sep="")
	  if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	  texte<-paste(texte,")\n",sep="")
	  texte<-paste(texte,"lines(varX2, varY2[,\"upr\"], lty=2",sep="")
	  if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	  texte<-paste(texte,")\n",sep="")
	  cat(texte)
	}
	if (tclvalue(Env$l.var$intervalA)==Env$voc[262,1]) {
	  if (nchar(facteur1)>0 & facteur1!=Env$voc[82,1]) {
	    texte<-"varX2 <- seq(min(varX, na.rm=TRUE), max(varX, na.rm=TRUE), abs(max(varX, na.rm=TRUE) - min(varX, na.rm=TRUE)) / 1000)\n"
	    texte<-paste(texte,"varY2 <- predict(model, list(varX=varX2), interval=\"confidence\")\n",sep="")
	  } else {
	    texte<-paste("varX2 <- seq(min(",varX,", na.rm=TRUE), max(",varX,", na.rm=TRUE), abs(max(",varX,", na.rm=TRUE) - min(",varX,", na.rm=TRUE)) / 1000)\n",sep="")
	    texte<-paste(texte,"varY2 <- predict(model, list(",varX,"=varX2), interval=\"prediction\")\n",sep="")
	  }
	  texte<-paste(texte,"lines(varX2, varY2[,\"lwr\"], lty=3",sep="")
	  if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	  texte<-paste(texte,")\n",sep="")
	  texte<-paste(texte,"lines(varX2, varY2[,\"upr\"], lty=3",sep="")
	  if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	  texte<-paste(texte,")\n",sep="")
	  cat(texte)
	}
	if (tclvalue(Env$l.var$intervalA)==Env$voc[263,1]) {
	  if (nchar(facteur1)>0 & facteur1!=Env$voc[82,1]) {
	    texte<-"varX2 <- seq(min(varX, na.rm=TRUE), max(varX, na.rm=TRUE), abs(max(varX, na.rm=TRUE) - min(varX, na.rm=TRUE)) / 1000)\n"
	    texte<-paste(texte,"varY2.conf <- predict(model, list(varX=varX2), interval=\"confidence\")\n",sep="")
	    texte<-paste(texte,"varY2.pred <- predict(model, list(varX=varX2), interval=\"prediction\")\n",sep="")
	  } else {
	    texte<-paste("varX2 <- seq(min(",varX,", na.rm=TRUE), max(",varX,", na.rm=TRUE), abs(max(",varX,", na.rm=TRUE) - min(",varX,", na.rm=TRUE)) / 1000)\n",sep="")
	    texte<-paste(texte,"varY2.conf <- predict(model, list(",varX,"=varX2), interval=\"confidence\")\n",sep="")
	    texte<-paste(texte,"varY2.pred <- predict(model, list(",varX,"=varX2), interval=\"prediction\")\n",sep="")
	  }
	  texte<-paste(texte,"lines(varX2, varY2.conf[,\"lwr\"], lty=2",sep="")
	  if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	  texte<-paste(texte,")\n",sep="")
	  texte<-paste(texte,"lines(varX2, varY2.conf[,\"upr\"], lty=2",sep="")
	  if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	  texte<-paste(texte,")\n",sep="")
	  texte<-paste(texte,"lines(varX2, varY2.pred[,\"lwr\"], lty=3",sep="")
	  if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	  texte<-paste(texte,")\n",sep="")
	  texte<-paste(texte,"lines(varX2, varY2.pred[,\"upr\"], lty=3",sep="")
	  if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	  texte<-paste(texte,")\n",sep="")
	  cat(texte)
	}
	cat("\n")
    } else if (tclvalue(Env$l.var$droiteA)==Env$voc[146,1]) {
	cat("# Regression line\n\n")
	texte<-""
	if (nchar(facteur1)>0 & facteur1!=Env$voc[82,1]) {
	  texte<-"b <- sd(varY, na.rm=TRUE) / sd(varX, na.rm=TRUE) * sign(cov(varX, varY, use=\"complete.obs\"))\n"
	  texte<-paste(texte,"a <- mean(varY, na.rm=TRUE) - b * mean(varX, na.rm=TRUE)\n",sep="")
	} else {
	  texte<-paste("b <- sd(",varY,", na.rm=TRUE) / sd(",varX,", na.rm=TRUE) * sign(cov(",varX,", ",varY,", use=\"complete.obs\"))\n",sep="")
	  texte<-paste(texte,"a <- mean(",varY,", na.rm=TRUE) - b * mean(",varX,", na.rm=TRUE)\n",sep="")
	}
	texte<-paste(texte,"abline(a, b",sep="")
	if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	texte<-paste(texte,", lty=",type.trait(type=tclvalue(Env$l.var$trait1)),sep="")
	texte<-paste(texte,", lwd=",as.numeric(tclvalue(Env$l.var$epaisseur1)),sep="")
	cat(paste(texte,")\n\n",sep=""))
    } else if (tclvalue(Env$l.var$droiteA)==Env$voc[147,1]) {
	cat("# Regression line\n\n")
	texte<-""
	if (nchar(facteur1)>0 & facteur1!=Env$voc[82,1]) {
	  texte<-"model <- lm(varY ~ varX + I(varX ^ 2))\n"
	  texte<-paste(texte,"varX2 <- seq(min(varX, na.rm=TRUE), max(varX, na.rm=TRUE), abs(max(varX, na.rm=TRUE) - min(varX, na.rm=TRUE)) / 1000)\n",sep="")
	  texte<-paste(texte,"varY2 <- predict(model, list(varX = varX2))\n",sep="")
	} else {
	  texte<-paste("model <- lm(",varY," ~ ",varX," + I(",varX," ^ 2))\n",sep="")
	  texte<-paste(texte,"varX2 <- seq(min(",varX,", na.rm=TRUE), max(",varX,", na.rm=TRUE), abs(max(",varX,", na.rm=TRUE) - min(",varX,", na.rm=TRUE)) / 1000)\n",sep="")
	  texte<-paste(texte,"varY2 <- predict(model, list(",varX," = varX2))\n",sep="")
	}
	texte<-paste(texte,"lines(varX2, varY2",sep="")
	if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	texte<-paste(texte,", lty=",type.trait(type=tclvalue(Env$l.var$trait1)),sep="")
	texte<-paste(texte,", lwd=",as.numeric(tclvalue(Env$l.var$epaisseur1)),sep="")
	cat(paste(texte,")\n",sep=""))
	if (tclvalue(Env$l.var$intervalA)==Env$voc[261,1]) {
	  texte<-""
	  if (nchar(facteur1)>0 & facteur1!=Env$voc[82,1]) {
	    texte<-paste("varY2.conf <- predict(model, list(varX=varX2), interval=\"confidence\")\n",sep="")
	  } else {
	    texte<-paste("varY2.conf <- predict(model, list(",varX,"=varX2), interval=\"confidence\")\n",sep="")
	  }
	  texte<-paste(texte,"lines(varX2, varY2.conf[,\"lwr\"], lty=2",sep="")
	  if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	  texte<-paste(texte,")\n",sep="")
	  texte<-paste(texte,"lines(varX2, varY2.conf[,\"upr\"], lty=2",sep="")
	  if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	  texte<-paste(texte,")\n",sep="")
	  cat(texte)
	}
	if (tclvalue(Env$l.var$intervalA)==Env$voc[262,1]) {
	  texte<-""
	  if (nchar(facteur1)>0 & facteur1!=Env$voc[82,1]) {
	    texte<-paste("varY2.pred <- predict(model, list(varX=varX2), interval=\"confidence\")\n",sep="")
	  } else {
	    texte<-paste("varY2.pred <- predict(model, list(",varX,"=varX2), interval=\"prediction\")\n",sep="")
	  }
	  texte<-paste(texte,"lines(varX2, varY2.pred[,\"lwr\"], lty=3",sep="")
	  if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	  texte<-paste(texte,")\n",sep="")
	  texte<-paste(texte,"lines(varX2, varY2.pred[,\"upr\"], lty=3",sep="")
	  if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	  texte<-paste(texte,")\n",sep="")
	  cat(texte)
	}
	if (tclvalue(Env$l.var$intervalA)==Env$voc[263,1]) {
	  texte<-""
	  if (nchar(facteur1)>0 & facteur1!=Env$voc[82,1]) {
	    texte<-paste("varY2.conf <- predict(model, list(varX=varX2), interval=\"confidence\")\n",sep="")
	    texte<-paste(texte,"varY2.pred <- predict(model, list(varX=varX2), interval=\"prediction\")\n",sep="")
	  } else {
	    texte<-paste("varY2.conf <- predict(model, list(",varX,"=varX2), interval=\"confidence\")\n",sep="")
	    texte<-paste(texte,"varY2.pred <- predict(model, list(",varX,"=varX2), interval=\"prediction\")\n",sep="")
	  }
	  texte<-paste(texte,"lines(varX2, varY2.conf[,\"lwr\"], lty=2",sep="")
	  if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	  texte<-paste(texte,")\n",sep="")
	  texte<-paste(texte,"lines(varX2, varY2.conf[,\"upr\"], lty=2",sep="")
	  if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	  texte<-paste(texte,")\n",sep="")
	  texte<-paste(texte,"lines(varX2, varY2.pred[,\"lwr\"], lty=3",sep="")
	  if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	  texte<-paste(texte,")\n",sep="")
	  texte<-paste(texte,"lines(varX2, varY2.pred[,\"upr\"], lty=3",sep="")
	  if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")}
	  texte<-paste(texte,")\n",sep="")
	  cat(texte)
	}
	cat("\n")
    } else if (tclvalue(Env$l.var$droiteA)==Env$voc[148,1]) {
	cat("# Tendency curve\n\n")
	texte<-""
	if (nchar(facteur1)>0 & facteur1!=Env$voc[82,1]) {
	  texte<-"panel.smooth(varX, varY"
	} else {
	  texte<-paste("panel.smooth(",varX,", ",varY,sep="")
	}
	if (tclvalue(Env$l.var$couleur2A)!="black" & tclvalue(Env$l.var$couleur2A)!="#000000") {
	  texte<-paste(texte,", col=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")
	  texte<-paste(texte,", col.smooth=\"",tclvalue(Env$l.var$couleur2A),"\"",sep="")
	}
	texte<-paste(texte,", pch=",graphe.symboles(num=as.numeric(tclvalue(Env$l.var$symboleA))),sep="")
	if (tclvalue(Env$l.var$taille.ptsA)!="1") {texte<-paste(texte,", cex=",as.numeric(tclvalue(Env$l.var$taille.ptsA)),sep="")}
	texte<-paste(texte,", lty=",type.trait(type=tclvalue(Env$l.var$trait1)),sep="")
	texte<-paste(texte,", lwd=",as.numeric(tclvalue(Env$l.var$epaisseur1)),sep="")
	cat(paste(texte,")\n\n",sep=""))
    }
  } else {
    niveaux<-levels(Env$dataset[,facteur1])[as.numeric(strsplit(tclvalue(Env$l.var$niveau),split=" ")[[1]])+1]
    symboles<-graphe.symboles(num=Env$l.var$symboleB)
    texte<-paste("plot(varY[",facteur1,"==\"",niveaux[1],"\"] ~ varX[",facteur1,"==\"",niveaux[1],"\"], axes=FALSE, ann=FALSE",sep="")
    texte<-paste(texte,", xlim=c(",round(Env$l.code$x.inf,2),",",round(Env$l.code$x.sup,2),")",sep="")
    texte<-paste(texte,", ylim=c(",round(Env$l.code$y.inf,2),",",round(Env$l.code$y.sup,2),")",sep="")
    if (graphe.log()!="") {texte<-paste(texte,", log=\"",graphe.log(),"\"",sep="")}
    texte<-paste(texte,", col=\"",Env$l.var$couleur2B[1],"\"",sep="")
    texte<-paste(texte,",\n  pch=",symboles[1],sep="")
    if (Env$l.var$taille.ptsB[1]!="1") {texte<-paste(texte,", cex=",as.numeric(Env$l.var$taille.ptsB[1]),sep="")}
    cat(paste(texte,")\n",sep=""))
    for (i in 2:length(niveaux)) {
	texte<-paste("points(varY[",facteur1,"==\"",niveaux[i],"\"] ~ varX[",facteur1,"==\"",niveaux[i],"\"]",sep="")
	texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")
	texte<-paste(texte,", pch=",symboles[i],sep="")
	if (Env$l.var$taille.ptsB[i]!=1) {texte<-paste(texte,", cex=",Env$l.var$taille.ptsB[i],sep="")}
	cat(paste(texte,")\n",sep=""))
    }
    cat("\n")
    if (tclvalue(Env$l.var$ptlab)==1) {
	for (i in 1:length(niveaux)) {
	  texte<-paste("text(varX[",facteur1,"==\"",niveaux[i],"\"], varY[",facteur1,"==\"",niveaux[i],"\"], pos=3, offset=0.4",sep="")
	  texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")
	  texte<-paste(texte,", cex=",round(0.65*Env$l.var$taille.ptsB[i],2),")\n",sep="")
	  cat(texte)
	}
	cat("\n")
    }
    code.graph.axes()
    code.graph.titre()
    if (tclvalue(Env$l.var$encadre)==1) {cat("box()\n\n")}
    for (i in 1:length(niveaux)) {
	if (Env$l.var$droiteB[i]==Env$voc[145,1]) {
	  cat("# Regression line\n\n")
	  texte<-paste("varX.",i," <- varX[",facteur1,"==\"",niveaux[i],"\"]\n",sep="")
	  texte<-paste(texte,"varY.",i," <- varY[",facteur1,"==\"",niveaux[i],"\"]\n",sep="")
	  texte<-paste(texte,"model.",i," <- lm(varY.",i," ~ varX.",i,")\n",sep="")
	  texte<-paste(texte,"abline(model.",i,"$coefficients",sep="")
	  if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	  texte<-paste(texte,", lty=",type.trait(type=Env$l.var$trait2[i]),sep="")
	  texte<-paste(texte,", lwd=",Env$l.var$epaisseur2[i],sep="")
	  cat(paste(texte,")\n",sep=""))
	  if (Env$l.var$intervalB[i]==Env$voc[261,1]) {
	    texte<-paste("varX2.",i," <- seq(min(varX.",i,", na.rm=TRUE), max(varX.",i,", na.rm=TRUE), abs(max(varX.",i,", na.rm=TRUE) - min(varX.",i,", na.rm=TRUE)) / 1000)\n",sep="")
	    texte<-paste(texte,"varY2.",i," <- predict(model.",i,", list(varX.",i,"=varX2.",i,"), interval=\"confidence\")\n",sep="")
	    texte<-paste(texte,"lines(varX2.",i,", varY2.",i,"[,\"lwr\"], lty=2",sep="")
	    if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	    texte<-paste(texte,")\n",sep="")
	    texte<-paste(texte,"lines(varX2.",i,", varY2.",i,"[,\"upr\"], lty=2",sep="")
	    if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	    texte<-paste(texte,")\n",sep="")
	    cat(texte)
	  }
	  if (Env$l.var$intervalB[i]==Env$voc[262,1]) {
	    texte<-paste("varX2.",i," <- seq(min(varX.",i,", na.rm=TRUE), max(varX.",i,", na.rm=TRUE), abs(max(varX.",i,", na.rm=TRUE) - min(varX.",i,", na.rm=TRUE)) / 1000)\n",sep="")
	    texte<-paste(texte,"varY2.",i," <- predict(model.",i,", list(varX.",i,"=varX2.",i,"), interval=\"prediction\")\n",sep="")
	    texte<-paste(texte,"lines(varX2.",i,", varY2.",i,"[,\"lwr\"], lty=3",sep="")
	    if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	    texte<-paste(texte,")\n",sep="")
	    texte<-paste(texte,"lines(varX2.",i,", varY2.",i,"[,\"upr\"], lty=3",sep="")
	    if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	    texte<-paste(texte,")\n",sep="")
	    cat(texte)
	  }
	  if (Env$l.var$intervalB[i]==Env$voc[263,1]) {
	    texte<-paste("varX2.",i," <- seq(min(varX.",i,", na.rm=TRUE), max(varX.",i,", na.rm=TRUE), abs(max(varX.",i,", na.rm=TRUE) - min(varX.",i,", na.rm=TRUE)) / 1000)\n",sep="")
	    texte<-paste(texte,"varY2.",i,".conf <- predict(model.",i,", list(varX.",i,"=varX2.",i,"), interval=\"confidence\")\n",sep="")
	    texte<-paste(texte,"varY2.",i,".pred <- predict(model.",i,", list(varX.",i,"=varX2.",i,"), interval=\"prediction\")\n",sep="")
	    texte<-paste(texte,"lines(varX2.",i,", varY2.",i,".conf[,\"lwr\"], lty=2",sep="")
	    if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	    texte<-paste(texte,")\n",sep="")
	    texte<-paste(texte,"lines(varX2.",i,", varY2.",i,".conf[,\"upr\"], lty=2",sep="")
	    if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	    texte<-paste(texte,")\n",sep="")
	    texte<-paste(texte,"lines(varX2.",i,", varY2.",i,".pred[,\"lwr\"], lty=3",sep="")
	    if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	    texte<-paste(texte,")\n",sep="")
	    texte<-paste(texte,"lines(varX2.",i,", varY2.",i,".pred[,\"upr\"], lty=3",sep="")
	    if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	    texte<-paste(texte,")\n",sep="")
	    cat(texte)
	  }
	  cat("\n")
	} else if (Env$l.var$droiteB[i]==Env$voc[146,1]) {
	  cat("# Regression line\n\n")
	  texte<-paste("b <- sd(varY[",facteur1,"==\"",niveaux[i],"\"], na.rm=TRUE) / sd(varX[",facteur1,"==\"",niveaux[i],"\"], na.rm=TRUE) * sign(cov(varX[",facteur1,"==\"",niveaux[i],"\"],\n  varY[",facteur1,"==\"",niveaux[i],"\"], use=\"complete.obs\"))\n",sep="")
	  texte<-paste(texte,"a <- mean(varY[",facteur1,"==\"",niveaux[i],"\"], na.rm=TRUE) - b * mean(varX[",facteur1,"==\"",niveaux[i],"\"], na.rm=TRUE)\n",sep="")
	  texte<-paste(texte,"abline(a, b",sep="")
	  if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	  texte<-paste(texte,", lty=",type.trait(type=Env$l.var$trait2[i]),sep="")
	  texte<-paste(texte,", lwd=",Env$l.var$epaisseur2[i],sep="")
	  cat(paste(texte,")\n\n",sep=""))
	} else if (Env$l.var$droiteB[i]==Env$voc[147,1]) {
	  cat("# Regression line\n\n")
	  texte<-paste("varX.",i," <- varX[",facteur1,"==\"",niveaux[i],"\"]\n",sep="")
	  texte<-paste(texte,"varY.",i," <- varY[",facteur1,"==\"",niveaux[i],"\"]\n",sep="")
	  texte<-paste(texte,"model.",i," <- lm(varY.",i," ~ varX.",i," + I(varX.",i," ^ 2))\n",sep="")
	  texte<-paste(texte,"varX2.",i," <- seq(min(varX.",i,", na.rm=TRUE), max(varX.",i,", na.rm=TRUE), abs(max(varX.",i,", na.rm=TRUE) - min(varX.",i,", na.rm=TRUE)) / 1000)\n",sep="")
	  texte<-paste(texte,"varY2.",i," <- predict(model.",i,", list(varX.",i," = varX2.",i,"))\n",sep="")
	  texte<-paste(texte,"lines(varX2.",i,", varY2.",i,sep="")
	  if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	  texte<-paste(texte,", lty=",type.trait(type=Env$l.var$trait2[i]),sep="")
	  texte<-paste(texte,", lwd=",Env$l.var$epaisseur2[i],sep="")
	  cat(paste(texte,")\n",sep=""))
	  if (Env$l.var$intervalB[i]==Env$voc[261,1]) {
	    texte<-paste("varY2.",i,".conf <- predict(model.",i,", list(varX.",i,"=varX2.",i,"), interval=\"confidence\")\n",sep="")
	    texte<-paste(texte,"lines(varX2.",i,", varY2.",i,".conf[,\"lwr\"], lty=2",sep="")
	    if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	    texte<-paste(texte,")\n",sep="")
	    texte<-paste(texte,"lines(varX2.",i,", varY2.",i,".conf[,\"upr\"], lty=2",sep="")
	    if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	    texte<-paste(texte,")\n",sep="")
	    cat(texte)
	  }
	  if (Env$l.var$intervalB[i]==Env$voc[262,1]) {
	    texte<-paste("varY2.",i,".pred <- predict(model.",i,", list(varX.",i,"=varX2.",i,"), interval=\"prediction\")\n",sep="")
	    texte<-paste(texte,"lines(varX2.",i,", varY2.",i,".pred[,\"lwr\"], lty=3",sep="")
	    if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	    texte<-paste(texte,")\n",sep="")
	    texte<-paste(texte,"lines(varX2.",i,", varY2.",i,".pred[,\"upr\"], lty=3",sep="")
	    if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	    texte<-paste(texte,")\n",sep="")
	    cat(texte)
	  }
	  if (Env$l.var$intervalB[i]==Env$voc[263,1]) {
	    texte<-paste("varY2.",i,".conf <- predict(model.",i,", list(varX.",i,"=varX2.",i,"), interval=\"confidence\")\n",sep="")
	    texte<-paste(texte,"varY2.",i,".pred <- predict(model.",i,", list(varX.",i,"=varX2.",i,"), interval=\"prediction\")\n",sep="")
	    texte<-paste(texte,"lines(varX2.",i,", varY2.",i,".conf[,\"lwr\"], lty=2",sep="")
	    if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	    texte<-paste(texte,")\n",sep="")
	    texte<-paste(texte,"lines(varX2.",i,", varY2.",i,".conf[,\"upr\"], lty=2",sep="")
	    if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	    texte<-paste(texte,")\n",sep="")
	    texte<-paste(texte,"lines(varX2.",i,", varY2.",i,".pred[,\"lwr\"], lty=3",sep="")
	    if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	    texte<-paste(texte,")\n",sep="")
	    texte<-paste(texte,"lines(varX2.",i,", varY2.",i,".pred[,\"upr\"], lty=3",sep="")
	    if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")}
	    texte<-paste(texte,")\n",sep="")
	    cat(texte)
	  }
	  cat("\n")
	} else if (Env$l.var$droiteB[i]==Env$voc[148,1]) {
	  cat("# Tendency curve\n\n")
	  texte<-paste("panel.smooth(varX[",facteur1,"==\"",niveaux[i],"\"], varY[",facteur1,"==\"",niveaux[i],"\"]",sep="")
	  if (Env$l.var$couleur2B[i]!="black" & Env$l.var$couleur2B[i]!="#000000") {
	    texte<-paste(texte,", col=\"",Env$l.var$couleur2B[i],"\"",sep="")
	    texte<-paste(texte,", col.smooth=\"",Env$l.var$couleur2B[i],"\"",sep="")
	  }
	  texte<-paste(texte,", pch=",symboles[i],sep="")
	  if (Env$l.var$taille.ptsB[i]!=1) {texte<-paste(texte,", cex=",Env$l.var$taille.ptsB[i],sep="")}
	  texte<-paste(texte,", lty=",type.trait(type=Env$l.var$trait2[i]),sep="")
	  texte<-paste(texte,", lwd=",Env$l.var$epaisseur2[i],sep="")
	  cat(paste(texte,")\n\n",sep=""))
	}
    }
    if (tclvalue(Env$l.var$legende)==1) {
	cat("# Legend\n\n")
	texte<-paste("legend(\"",code.graph.posleg(),"\"",sep="")
	texte<-paste(texte,", legend=c(\"",paste(Env$l.var$noms1,collapse="\",\""),"\")",sep="")
	texte<-paste(texte,", col=c(\"",paste(Env$l.var$couleur2B,collapse="\",\""),"\")",sep="")
	texte<-paste(texte,",\n  pch=c(",paste(symboles,collapse=","),")",sep="")
	texte<-paste(texte,", pt.cex=c(",paste(Env$l.var$taille.ptsB,collapse=","),")",sep="")
	if (nchar(tclvalue(Env$l.var$legende.titre))>0) {texte<-paste(texte,", title=\"",tclvalue(Env$l.var$legende.titre),"\"",sep="")}
	cat(paste(texte,")\n\n",sep=""))
    }
  }
}


#-------------------------------------------------
# Code - Graphe
#-------------------------------------------------

code.graph<-function() {
  cat("# Graph\n\n")
  if (Env$l.var$ecran=="H") {
    code.graph.hist()
  } else if (Env$l.var$ecran=="M") {
    code.graph.moust()
  } else if (Env$l.var$ecran=="B") {
    code.graph.barres()
  } else if (Env$l.var$ecran=="Ca") {
    code.graph.cam()
  } else if (Env$l.var$ecran=="Co") {
    code.graph.courbe()
  } else if (Env$l.var$ecran=="N") {
    code.graph.nuage()
  }
}


#-------------------------------------------------
# Code - Vérification si enregistrement du code
#   avant de tracer un graphe
#-------------------------------------------------

pretracer<-function() {
  test<-tracer()
  if (test==TRUE) {
    if (Env$l.code$ask==FALSE) {
	Env$l.code$ask<-TRUE
	code.ask()
    } else {
	if (Env$l.code$save==TRUE) {
	  code.open()
	}
    }
  }
}


#-------------------------------------------------
# Fenêtre principale
#-------------------------------------------------

ouvrir.GrapheR<-function() {
  Env$Fen<-tktoplevel()
  tktitle(Env$Fen)<-"GrapheR"
  tkwm.geometry(Env$Fen, "+30+30")
  ### Barre de navigation
  Env$l.wdg$vide1<-tklabel(Env$Fen,text="  ",font=Env$police)
  Env$l.wdg$but.data<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","But_data.gif",fsep=.Platform$file.sep)),command=function() {navigation(type="data")})
  Env$l.wdg$sep1<-tklabel(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Sep.gif",fsep=.Platform$file.sep)))
  Env$l.wdg$but.his<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","But_hist.gif",fsep=.Platform$file.sep)),command=function() {navigation(type="hist")})
  Env$l.wdg$but.box<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","But_box.gif",fsep=.Platform$file.sep)),command=function() {navigation(type="moust")})
  Env$l.wdg$but.bar<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","But_bar.gif",fsep=.Platform$file.sep)),command=function() {navigation(type="barres")})
  Env$l.wdg$but.pie<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","But_cam.gif",fsep=.Platform$file.sep)),command=function() {navigation(type="cam")})
  Env$l.wdg$but.crb<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","But_courb.gif",fsep=.Platform$file.sep)),command=function() {navigation(type="courbe")})
  Env$l.wdg$but.sct<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","But_nuage.gif",fsep=.Platform$file.sep)),command=function() {navigation(type="nuage")})
  Env$l.wdg$sep2<-tklabel(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Sep.gif",fsep=.Platform$file.sep)))
  Env$l.wdg$but.win<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","But_window.gif",fsep=.Platform$file.sep)),command=new.window)
  Env$l.wdg$sep3<-tklabel(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Sep.gif",fsep=.Platform$file.sep)))
  Env$l.wdg$but.drw<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","But_draw.gif",fsep=.Platform$file.sep)),command=pretracer)
  Env$l.wdg$sep4<-tklabel(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Sep.gif",fsep=.Platform$file.sep)))
  Env$l.wdg$but.hor<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","But_hor.gif",fsep=.Platform$file.sep)),command=horizontal)
  Env$l.wdg$but.ver<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","But_ver.gif",fsep=.Platform$file.sep)),command=vertical)
  Env$l.wdg$but.drt<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","But_drt.gif",fsep=.Platform$file.sep)),command=affine)
  Env$l.wdg$but.dis<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","But_dis.gif",fsep=.Platform$file.sep)),command=distrib)
  Env$l.wdg$but.txt<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","But_text.gif",fsep=.Platform$file.sep)),command=texte)
  Env$l.wdg$but.p<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","But_p.gif",fsep=.Platform$file.sep)),command=pval)
  Env$l.wdg$sep5<-tklabel(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Sep.gif",fsep=.Platform$file.sep)))
  Env$l.wdg$but.save<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","But_save.gif",fsep=.Platform$file.sep)),command=enregistrer)
  Env$l.wdg$sep6<-tklabel(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Sep.gif",fsep=.Platform$file.sep)))
  Env$l.wdg$but.lang<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","But_lang.gif",fsep=.Platform$file.sep)),command=language.change)
  Env$l.wdg$but.help<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","But_help.gif",fsep=.Platform$file.sep)),command=aide)
  tkbind(Env$l.wdg$but.data,"<Enter>",function() {msg(text=Env$voc[221,1],type="info")})
  tkbind(Env$l.wdg$but.data,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.wdg$but.his,"<Enter>",function() {msg(text=Env$voc[222,1],type="info")})
  tkbind(Env$l.wdg$but.his,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.wdg$but.box,"<Enter>",function() {msg(text=Env$voc[223,1],type="info")})
  tkbind(Env$l.wdg$but.box,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.wdg$but.bar,"<Enter>",function() {msg(text=Env$voc[224,1],type="info")})
  tkbind(Env$l.wdg$but.bar,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.wdg$but.pie,"<Enter>",function() {msg(text=Env$voc[225,1],type="info")})
  tkbind(Env$l.wdg$but.pie,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.wdg$but.crb,"<Enter>",function() {msg(text=Env$voc[226,1],type="info")})
  tkbind(Env$l.wdg$but.crb,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.wdg$but.sct,"<Enter>",function() {msg(text=Env$voc[227,1],type="info")})
  tkbind(Env$l.wdg$but.sct,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.wdg$but.win,"<Enter>",function() {msg(text=Env$voc[228,1],type="info")})
  tkbind(Env$l.wdg$but.win,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.wdg$but.drw,"<Enter>",function() {msg(text=Env$voc[229,1],type="info")})
  tkbind(Env$l.wdg$but.drw,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.wdg$but.hor,"<Enter>",function() {msg(text=Env$voc[230,1],type="info")})
  tkbind(Env$l.wdg$but.hor,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.wdg$but.ver,"<Enter>",function() {msg(text=Env$voc[231,1],type="info")})
  tkbind(Env$l.wdg$but.ver,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.wdg$but.drt,"<Enter>",function() {msg(text=Env$voc[232,1],type="info")})
  tkbind(Env$l.wdg$but.drt,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.wdg$but.dis,"<Enter>",function() {msg(text=Env$voc[233,1],type="info")})
  tkbind(Env$l.wdg$but.dis,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.wdg$but.txt,"<Enter>",function() {msg(text=Env$voc[234,1],type="info")})
  tkbind(Env$l.wdg$but.txt,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.wdg$but.p,"<Enter>",function() {msg(text=Env$voc[235,1],type="info")})
  tkbind(Env$l.wdg$but.p,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.wdg$but.save,"<Enter>",function() {msg(text=Env$voc[236,1],type="info")})
  tkbind(Env$l.wdg$but.save,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.wdg$but.lang,"<Enter>",function() {msg(text=Env$voc[237,1],type="info")})
  tkbind(Env$l.wdg$but.lang,"<Leave>",function() {msg(text="",type="info")})
  tkbind(Env$l.wdg$but.help,"<Enter>",function() {msg(text=Env$voc[238,1],type="info")})
  tkbind(Env$l.wdg$but.help,"<Leave>",function() {msg(text="",type="info")})
  tkgrid(Env$l.wdg$vide1,Env$l.wdg$but.data,Env$l.wdg$sep1,Env$l.wdg$but.his,Env$l.wdg$but.box,Env$l.wdg$but.bar,Env$l.wdg$but.pie,Env$l.wdg$but.crb,
    Env$l.wdg$but.sct,Env$l.wdg$sep2,Env$l.wdg$but.win,Env$l.wdg$sep3,Env$l.wdg$but.drw,Env$l.wdg$sep4,Env$l.wdg$but.hor,Env$l.wdg$but.ver,
    Env$l.wdg$but.drt,Env$l.wdg$but.dis,Env$l.wdg$but.txt,Env$l.wdg$but.p,Env$l.wdg$sep5,Env$l.wdg$but.save,Env$l.wdg$sep6,Env$l.wdg$but.lang,
    Env$l.wdg$but.help)
  tkgrid(tklabel(Env$Fen,text="",font=Env$police2),row=1,column=0)
  ### Cadre de texte pour les messages
  Env$l.wdg$message.wdg<-tkentry(Env$Fen,width=100,font=Env$police5,textvariable=Env$l.var$message,state="readonly",readonlybackground="white")
  tkgrid(Env$l.wdg$message.wdg,row=2,column=0,columnspan=25)
  tkgrid(tklabel(Env$Fen,text="",font=Env$police2),row=3,column=0)
  msg(text=Env$voc[1,1],type="info")
  ### Frame gauche
  Env$l.frames$Frg<-tkframe(Env$Fen)
  Env$l.lab$lab1<-tklabel(Env$l.frames$Frg,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[1,1],fsep=.Platform$file.sep)))
  Env$l.wdg$but.lab1<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)),command=function() {
	if (Env$l.frames$Fr1.status==1) {fr1.close()} else {
	if (Env$l.var$ecran=="D") {fr1.openD()} else
	if (Env$l.var$ecran=="H") {fr1.openH()} else
	if (Env$l.var$ecran=="M") {fr1.openM()} else
	if (Env$l.var$ecran=="B") {fr1.openB()} else
	if (Env$l.var$ecran=="Ca") {fr1.openCa()} else
	if (Env$l.var$ecran=="Co") {fr1.openCo()} else
	if (Env$l.var$ecran=="N") {fr1.openN()}
    }
  })
  tkgrid(Env$l.lab$lab1,Env$l.wdg$but.lab1)
  ## Frame 1
  Env$l.frames$Fr1<-tkframe(Env$l.frames$Frg)
  Env$l.fr1$vide<-tklabel(Env$l.frames$Fr1,text="",font=Env$police2)
  tkgrid(Env$l.fr1$vide)
  fr1.openD()
  tkgrid(Env$l.frames$Fr1)
  tkgrid(tklabel(Env$l.frames$Frg,text="",font=Env$police2))
  Env$l.lab$lab2<-tklabel(Env$l.frames$Frg,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[2,1],fsep=.Platform$file.sep)))
  Env$l.wdg$but.lab2<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)),command=function() {
    if (Env$l.frames$Fr2.status==1) {fr2.close()} else {
	if (Env$l.var$ecran=="D") {fr2.openD()} else
	if (Env$l.var$ecran%in%c("H","M","B","Ca","Co","N")) {fr2.opengraphe()}
    }
  })
  tkgrid(Env$l.lab$lab2,Env$l.wdg$but.lab2)
  ## Frame 2
  Env$l.frames$Fr2<-tkframe(Env$l.frames$Frg)
  Env$l.fr2$vide<-tklabel(Env$l.frames$Fr2,text="",font=Env$police2)
  tkgrid(Env$l.fr2$vide)
  fr2.close()
  tkgrid(Env$l.frames$Fr2)
  tkgrid(tklabel(Env$l.frames$Frg,text="",font=Env$police2))
  Env$l.lab$lab3<-tklabel(Env$l.frames$Frg,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[3,1],fsep=.Platform$file.sep)))
  Env$l.wdg$but.lab3<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)),command=function() {
    if (Env$l.frames$Fr3.status==1) {fr3.close()} else {
	if (Env$l.var$ecran=="D") {fr3.openD()} else
	if (Env$l.var$ecran=="H") {fr3.openH()} else
	if (Env$l.var$ecran=="M") {fr3.openM()} else
	if (Env$l.var$ecran=="B") {fr3.openB()} else
	if (Env$l.var$ecran=="Ca") {fr3.openCa()} else
	if (Env$l.var$ecran%in%c("Co","N")) {fr3.openCoN()}
    }
  })
  tkgrid(Env$l.lab$lab3,Env$l.wdg$but.lab3)
  ## Frame 3
  Env$l.frames$Fr3<-tkframe(Env$l.frames$Frg)
  Env$l.fr3$vide<-tklabel(Env$l.frames$Fr3,text="",font=Env$police2)
  tkgrid(Env$l.fr3$vide)
  fr3.close()
  tkgrid(Env$l.frames$Fr3)
  tkgrid(tklabel(Env$l.frames$Frg,text="",font=Env$police2))
  Env$l.lab$lab4<-tklabel(Env$l.frames$Frg,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[4,1],fsep=.Platform$file.sep)))
  Env$l.wdg$but.lab4<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)),command=function() {
    if (Env$l.frames$Fr4.status==1) {fr4.close()} else {
	if (Env$l.var$ecran=="D") {fr4.openD()} else
	if (Env$l.var$ecran=="H") {fr4.openH()} else
	if (Env$l.var$ecran=="M") {fr4.openM()} else
	if (Env$l.var$ecran=="B") {fr4.openB()} else
	if (Env$l.var$ecran=="Ca") {fr4.openCa()} else
	if (Env$l.var$ecran=="Co") {fr4.openCo()} else
	if (Env$l.var$ecran=="N") {fr4.openN()}
    }
  })
  tkgrid(Env$l.lab$lab4,Env$l.wdg$but.lab4)
  ## Frame 4
  Env$l.frames$Fr4<-tkframe(Env$l.frames$Frg)
  Env$l.fr4$vide<-tklabel(Env$l.frames$Fr4,text="",font=Env$police2)
  tkgrid(Env$l.fr4$vide)
  fr4.close()
  tkgrid(Env$l.frames$Fr4)
  tkgrid(tklabel(Env$l.frames$Frg,text="",font=Env$police2))
  Env$l.lab$lab5<-tklabel(Env$l.frames$Frg,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images",Env$img[18,1],fsep=.Platform$file.sep)))
  Env$l.wdg$but.lab5<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)),command=function() {
    if (Env$l.frames$Fr5.status==1) {fr5.close()} else {
	if (Env$l.var$ecran=="D") {fr5.openD()} else
	if (Env$l.var$ecran=="H") {fr5.openH()} else
	if (Env$l.var$ecran=="M") {fr5.openM()} else
	if (Env$l.var$ecran=="B") {fr5.openB()} else
	if (Env$l.var$ecran=="Ca") {fr5.openCa()} else
	if (Env$l.var$ecran=="Co") {fr5.openCo()} else
	if (Env$l.var$ecran=="N") {fr5.openN()}
    }
  })
  tkgrid(Env$l.lab$lab5,Env$l.wdg$but.lab5)
  ## Frame 5
  Env$l.frames$Fr5<-tkframe(Env$l.frames$Frg)
  Env$l.fr5$vide<-tklabel(Env$l.frames$Fr5,text="",font=Env$police2)
  tkgrid(Env$l.fr5$vide)
  fr5.close()
  tkgrid(Env$l.frames$Fr5)
  tkgrid(tklabel(Env$l.frames$Frg,text="",font=Env$police2))
  Env$l.lab$lab6<-tklabel(Env$l.frames$Frg,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Lab0.gif",fsep=.Platform$file.sep)))
  Env$l.wdg$but.lab6<-tkbutton(Env$Fen,image=tkimage.create("photo",file=file.path(path.package("GrapheR"),"images","Fleche_haut.gif",fsep=.Platform$file.sep)),command=function() {
    if (Env$l.frames$Fr6.status==1) {fr6.close()} else {
	if (Env$l.var$ecran=="M") {fr6.openM()} else
	if (Env$l.var$ecran=="B") {fr6.openB()} else
	if (Env$l.var$ecran=="Co") {fr6.openCo()} else
	if (Env$l.var$ecran=="N") {fr6.openN()}
    }
  },state="disabled")
  tkgrid(Env$l.lab$lab6,Env$l.wdg$but.lab6)
  ## Frame 6
  Env$l.frames$Fr6<-tkframe(Env$l.frames$Frg)
  Env$l.fr6$vide<-tklabel(Env$l.frames$Fr6,text="",font=Env$police2)
  tkgrid(Env$l.fr6$vide)
  fr6.close()
  tkgrid(Env$l.frames$Fr6)
  tkgrid(tklabel(Env$l.frames$Frg,text="",font=Env$police2))
  tkgrid(Env$l.frames$Frg,row=4,column=0,columnspan=25)
  ### Frame droite
  Env$l.frames$Frd<-tkframe(Env$Fen)
  ## Frame 7
  Env$l.frames$Fr7<-tkframe(Env$l.frames$Frd,relief="ridge",borderwidth=0)
  Env$l.fr7$vide<-tklabel(Env$l.frames$Fr7,text="",font=Env$police2)
  tkgrid(Env$l.fr7$vide)
  tkgrid(Env$l.frames$Fr7)
  tkgrid(Env$l.frames$Frd,row=1,column=26,rowspan=4,sticky="n")
  ### Frame droite finale
  Env$l.frames$Frd2<-tkframe(Env$Fen)
  Env$l.wdg$vide2<-tklabel(Env$l.frames$Frd2,text="",font=Env$police)
  tkgrid(Env$l.wdg$vide2)
  tkgrid(Env$l.frames$Frd2,row=0,column=27)
}


#-------------------------------------------------
# Fermeture de l'interface
#-------------------------------------------------

fermer.GrapheR<-function() {
  tkdestroy(Env$Fen)
}
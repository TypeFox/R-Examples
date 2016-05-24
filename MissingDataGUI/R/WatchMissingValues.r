##' The Main Window of Missing Data GUI.
##'
##' This function is to open the missing data GUI. The widgets shown
##' in the GUI include: a table of all variables in the dataset, a
##' checkbox group of categorical variables to condition on, a table
##' of variables which have missing values to coloy by, a radio of
##' imputation methods, a radio of graph types, three command buttons,
##' and a graphics device. In this GUI the user can: 1)change the name
##' and class of the selected variable; 2)look at the numeric summary
##' for the missing values in the selected variables; 3)look at the
##' plot of imputed data, under one of the imputation methods and one
##' of the graph types and one color-by variable, with or without the
##' conditions; 4)export the imputed data as well as the missing
##' shadow matrix, and save them to a data file(csv).
##'
##' The missing data GUI consists of two tabs. In the summary tab,
##' there are a list of all variables, a list of variables having
##' missing values to color by, two radios for imputation methods and
##' graph types respectively, a checkbox group for the conditional
##' variables, four buttons and a graphics device. In the help tab,
##' the layout is the same as the summary tab.  But when the users
##' move their mouse on those widgets, or click any of those radios or
##' buttons, the functions of all widgets will be described at the
##' place of the graphics device. The attributes of the variables can
##' be changed. If the user double clicks on any variables in the top
##' left table of missing-data GUI, an attribute window will pop
##' up. Then the name could be edited, and the class could be changed
##' to one of the four classes: integer, numeric, factor, and
##' character. When a numeric variable is changed to a categorical
##' variable, the condtions in the bottom left checkbox group will be
##' updated. If the list of the color by variables is very long, the
##' selector allows text entry to find the variable when this widget
##' is active.
##' @param h A list with components obj referring to the button "Watch
##' Missing Values" in \code{\link{MissingDataGUI}}.
##' @param data A data frame which is shown in the missing-data
##' GUI. If it is null, then parameter gt must not be null.
##' @param gt A widget created by gtable(). It should be passed from
##' \code{\link{MissingDataGUI}}.
##' @param size.width the width of window. Default to be 1000, and the
##' minimal is 800.
##' @param size.height the height of window. Default to be 750, and the
##' minimal is 600.
##' @param ... Other parameters to be passed to this function.
##' @return NULL
##' @author Xiaoyue Cheng <\email{xycheng@@unomaha.edu}>
##' @import ggplot2 cairoDevice gWidgetsRGtk2
##' @importFrom reshape melt.data.frame
##' @importFrom GGally ggpairs
##' @importFrom GGally ggparcoord
##' @export
##' @examples
##' if(interactive()){
##' data(tao)
##' WatchMissingValues(data=tao)
##' data(brfss)
##' WatchMissingValues(data=brfss)
##' }
##'
WatchMissingValues = function(h, data=NULL, gt=NULL, size.width=1000, size.height=750, ...){
  if (is.null(data) & is.null(gt)) {
    gmessage("There's no input.", icon="error")
    return()
  }

  m = new.env()
  m$group151 = m$glay151 = list()

  if (is.null(data)) {
    if (length(svalue(gt)) == 0) {
      gtfile = gt[1,]
    } else {
      gtfile = svalue(gt)[1]
    }
    dataset = read.csv(file = gtfile, header = TRUE)
  } else {
    dataset = data
  }

  rows = nrow(dataset)
  cols = ncol(dataset)
  halfcol = as.integer(cols/2)
  dataclass = as.character(sapply(dataset, function(x) class(x)[1]))
  vNApct = sapply(dataset, function(avec){mean(is.na(avec))})

  if (sum(vNApct)==0) {
    if (cols%%2==0) {
      chr_shadowmatrix = unlist(sapply(dataset[,(halfcol+1):cols],as.character))
      if (all(dataclass[(halfcol+1):cols]=="logical") |
        all(chr_shadowmatrix %in% c("TRUE","FALSE")) |
        all(chr_shadowmatrix %in% c("0","1"))) {
        imp_dat = data.frame(dataset[,1:halfcol],row_number=1:rows)
        tmp = dataset[,1:halfcol]
        tmp[sapply(dataset[,(halfcol+1):cols],as.logical)] = NA
        dataset = tmp
        is_imputed_data = TRUE
        gmessage(paste('There are no missing values in this data set. The last', halfcol,
                       'columns are used as the shadow matrix.'), icon = "info")
      } else {
        is_imputed_data = FALSE
        gmessage('There are no missing values in this data set. If the data is combined
                  with a shadow matrix, then the number of columns of the shadow matrix
                  is not the same as original data.', icon = "error")
      }
    } else {
      is_imputed_data = FALSE
      gmessage('There are no missing values in this data set. If the data is combined
                with a shadow matrix, then the number of columns of the shadow matrix
                is not the same as original data.', icon = "error")
    }
  } else {is_imputed_data = FALSE}
  vname = as.character(names(dataset))
  dataclass = as.character(sapply(dataset, function(x) class(x)[1]))
  NAcol = which(sapply(dataset, function(avec){all(is.na(avec))}))
  vNApct = sapply(dataset, function(avec){mean(is.na(avec))})
  m$dataset = dataset

  #####------------------------------------#####
  ##  Graph and SavePlot share too much code  ##
  ##  So I made the following functions       ##
  #####------------------------------------#####
  initializ = function(){
      m$name_select=svalue(gt11, index = TRUE)
      m$n=length(m$name_select)
      m$imp_method=svalue(gr142)
      m$graphtype=svalue(gr143)
      cond=check123[svalue(check123,index=T)]
      idx_cond = if (length(cond)) (cond %in% gt11[m$name_select,2]) else FALSE
      if (any(idx_cond) && any((gt11[(gt11[,2] %in% cond[idx_cond]),3])>0) && m$graphtype != "Missingness Map") {
        gmessage("At least one of the conditional variables contains missing values and is selected to be imputed.\n\nWe will unselect the variable(s).\n\nTip: please impute the conditional variable before conditioning on it.",icon='warning')
        miss_cond = gt11[(gt11[,2] %in% cond[idx_cond]),1]
        svalue(gt11) = setdiff(m$name_select, miss_cond)
        m$name_select=svalue(gt11, index = TRUE)
        m$n=length(m$name_select)
      }
      if (length(cond)==0 || m$imp_method=='Below 10%') m$cond=NULL else {
        m$cond=m$dataset[,cond,drop=FALSE]
      }
      m$colorby=as.character(svalue(radio125))
      if (length(m$colorby)==0 | 'Missing on Selected Variables' %in% m$colorby) {
          m$colorby="Missing on Selected Variables"
      } else {
          if ("Missing Any Variables" %in% m$colorby) {
              m$colorby="Missing Any Variables"
          }
      }
      m$mi_n=svalue(text32112)
      if (is.na(m$mi_n) || m$mi_n < 1) {
        m$mi_n=3
        message("Number of imputation sets is set to 3.")
        svalue(text32112)=3
      }
      m$mi_seed=svalue(text32122)
      if (is.na(m$mi_seed)) {
        m$mi_seed=1234567
        message("Random number seed is set to 1234567.")
        svalue(text32122)=1234567
      }
      m$nn_k=svalue(text32212)
      if (is.na(m$nn_k)) {
        m$nn_k=5
        message("Number of the nearest neighbors is set to 5.")
        svalue(text32212)=5
      }
      m$parajitter=svalue(radio32311)
      if (m$imp_method == 'MI:mice') attr(m$imp_method,'method')=gt3421[m$name_select,4]
  }
  initial_plot = function(){
      if (m$n == 0) {
          if (m$graphtype!="Missingness Map"){
              gmessage("Please select at least one variable!", icon = "error")
              return(TRUE)
          } else {
              m$name_select = 1:nrow(gt11)
              m$n = length(m$name_select)
              return(FALSE)
          }
      } else if (m$graphtype=="Missingness Map") {
        if (m$n ==1) {
          gmessage("Please select at least two variables for a missingness map.", icon = "error")
          return(TRUE)
        } else {return(FALSE)}
      }

      if ( (!exists('imp_dat',envir=m))  || m$graphtype!="Below 10%" ) {
          m$dat = imputation(origdata=m$dataset[,gt11[m$name_select,2],drop=FALSE],
                      method=m$imp_method, vartype=as.character(gt11[m$name_select,3]),
                      missingpct=as.numeric(as.character(gt11[m$name_select,4])),
                      condition=m$cond,knn=m$nn_k,mi.n=m$mi_n,mi.seed=m$mi_seed)
          if (all(sapply(m$dat,nrow)==0)) return(TRUE)
          for (j in 1:length(m$dat)) {
            m$dat[[j]]=m$dat[[j]][,c(gt11[m$name_select,2],'row_number'),drop=FALSE]
          }
      } else {
          m$dat = list(Imported=data.frame(m$imp_dat[,c(gt11[m$name_select,2])],m$imp_dat[,ncol(m$imp_dat)]))
      }

#       for (j in 1:length(m$dat)){
#           for (i in 1:m$n){
#               eval(parse(text=paste("m$dat[[j]][,i]=as.",as.character(gt11[m$name_select,3])[i],"(as.character(m$dat[[j]][,i]))",sep="")))
#           }
#       }

      if (m$colorby[1]=='Missing on Selected Variables') {
          Missing <- !complete.cases(m$dataset[,gt11[m$name_select,2]])
      } else {
          if (m$colorby[1]=='Missing Any Variables') {
              Missing <- !complete.cases(m$dataset)
          } else {
              Missing <- !complete.cases(m$dataset[,m$colorby])
          }
      }
      m$Missing <- Missing[m$dat[[1]][,ncol(m$dat[[1]])]]
      return(FALSE)
  }
  graph_hist = function(j, i){
      dat = m$dat[[j]]
      if (!(is.numeric(dat[,i]) | is.character(dat[,i]) | is.factor(dat[,i]))) return()
      tmpdat = data.frame(dat,Missing=m$Missing)
      p=qplot(tmpdat[,i],data=tmpdat,geom='histogram',fill=Missing,
              position='stack', xlab=names(tmpdat)[i])
      if (is.numeric(dat[,i])) suppressMessages(print(p))
      if (is.character(dat[,i])) suppressMessages(print(p+coord_flip()))
      if (is.factor(dat[,i]) &
          as.numeric(as.character(gt11[m$name_select,4]))[i]<1) suppressMessages(print(p+coord_flip()))
  }
  graph_spine = function(j, i){
    dat = m$dat[[j]]
    if (!(is.numeric(dat[,i]) | is.character(dat[,i]) | is.factor(dat[,i]))) return()
    if (is.numeric(dat[,i])){
      itv=pretty(dat[,i],n=10)
      itvct=cut(dat[,i],itv)
    } else {
      itvct=dat[,i]
    }
    a=table(itvct,m$Missing)
    b=matrix(a,ncol=2)
    bsum=rowSums(b)
    b1=apply(b,2,`/`,bsum)
    b1sum=c(0,cumsum(bsum))
    if (is.numeric(dat[,i])){
      d1=data.frame(x1=b1sum[-length(b1sum)],x2=b1sum[-1],y1=0,y2=b1[,1],Missing='FALSE')
      d2=data.frame(x1=b1sum[-length(b1sum)],x2=b1sum[-1],y1=b1[,1],y2=1,Missing='TRUE')
      b1center=b1sum
    } else {
      b1center=b1sum[-length(b1sum)]+diff(b1sum)/2
      b1width=diff(b1sum)*9/20
      d1=data.frame(x1=b1center-b1width,x2=b1center+b1width,y1=0,y2=b1[,1],Missing='FALSE')
      d2=data.frame(x1=b1center-b1width,x2=b1center+b1width,y1=b1[,1],y2=1,Missing='TRUE')
      itv=rownames(a)
    }
    d=rbind(d1,d2)
    p=ggplot(d, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill=Missing))+geom_rect()+
      scale_x_continuous(breaks=b1center, labels=itv)+
      theme(panel.grid.minor.x = element_blank())+
      xlab(colnames(dat)[i])+ylab("Proportion")
    if (is.numeric(dat[,i])) suppressMessages(print(p)) else suppressMessages(print(p+coord_flip()))
  }
  graph_pair = function(j, legend.pos, type, contour='id'){
      if (m$n < 2) {
          gmessage("You selected less than two variables! Please re-select.", icon = "error")
          return()
      }
      if (m$n > 5) {
          gmessage("You selected more than five variables! Only the first five are displayed.", icon = "warning")
          m$n = 5
      }
      dat = m$dat[[j]]
      dat = dat[order(m$Missing),]
      if (m$n==2) {
          Missing=factor(m$Missing)[order(m$Missing)]
          rsltnx=ifelse(is.numeric(dat[,1]),resolution(dat[,1]),2.5)
          rsltny=ifelse(is.numeric(dat[,2]),resolution(dat[,2]),2.5)
          print(qplot(dat[,1],dat[,2], color=Missing, geom='jitter',alpha=I(0.7),
                      position=position_jitter(width=rsltnx/8,height=rsltny/8),
                      size=I(3),xlab=colnames(dat)[1],ylab=colnames(dat)[2]) +
                    theme(legend.position=legend.pos))
      } else {
          dat$Missings=factor(m$Missing)[order(m$Missing)]
          if (contour == 'jitter'){
            for (i in 1:m$n){
              idx = which(is.na(m$dataset[,m$name_select[i]])[order(m$Missing)])
              if (is.numeric(dat[,i]) & length(idx)) {
                rsltn = diff(range(dat[,i],na.rm=TRUE))/1000
                dat[idx,i] = dat[idx,i]+rnorm(length(idx),sd=rsltn)
              }
            }
          }
          m$p = ggpairs(dat,columns=1:m$n,
                        upper=list(continuous='density',
                                   combo='box',
                                   discrete='facetbar'),
                        mapping = ggplot2::aes(color = Missings,alpha=0.5))
          for(pmrow in 1:m$p$nrow){
            for(pmcol in 1:m$p$ncol){
              m$p[pmrow,pmcol] <- m$p[pmrow,pmcol] +
                scale_colour_discrete() +
                scale_fill_discrete()
            }
          }
          if (.Platform$OS.type!='windows' || type=='save') {print(m$p)} else {
            winprint(m$p)
          }
      }
  }
  graph_pcp = function(j,...){
      if (m$n==1) {
          gmessage('You only selected one variable. Cannot plot the parallel coordinates.',
                   icon = "error")
          return(TRUE)
      }
      library_call("reshape")
      dat = m$dat[[j]]
      if (m$parajitter == "Jitter on categorical variables") {
        rowz = nrow(dat)
        for (i in which(sapply(lapply(dat[,1:m$n],class),`[`,1) %in% c('logical','factor','ordered','character'))) {
          if (is.character(dat[,i])) dat[,i]=factor(dat[,i])
          dat[,i] = as.integer(dat[,i])
          set.seed(m$mi_seed)
          dat[,i] = dat[,i] + runif(rowz, -0.1, 0.1)
        }
      } else {
        for (i in which(sapply(lapply(dat[,1:m$n],class),`[`,1) == 'ordered')) {
          dat[,i] = as.integer(dat[,i])
        }
      }
      dat$Missing = m$Missing
      dat = dat[order(dat$Missing),]
      m$p = ggparcoord(dat,1:m$n,groupColumn='Missing',scale='uniminmax',
                       alphaLines=0.5,...)+scale_colour_discrete()
      return(FALSE)
  }
  graph_map = function(){
      Mapdat=data.frame(is.na(m$dataset[,gt11[m$name_select,2]]))
      Mapdat$observation=1:nrow(Mapdat)
      mapdat=melt.data.frame(Mapdat,"observation")
      colnames(mapdat)[3]="Missing"
      mapdat$variable=factor(mapdat$variable,levels=rev(levels(mapdat$variable)),
                             labels=rev(levels(mapdat$variable)))
      q1=qplot(observation,variable,data=mapdat,geom='tile',fill=Missing,main='Original Missingness Map')

      Mapdat2=Mapdat[order(rowSums(Mapdat[,-(m$n+1)]),decreasing=TRUE),]
      Mapdat2=Mapdat2[,order(colSums(Mapdat2),decreasing=TRUE)]
      Mapdat2$Reordered_Observation=1:nrow(Mapdat2)
      mapdat=melt.data.frame(Mapdat2[,-1],"Reordered_Observation")
      colnames(mapdat)[3]="Missing"
      mapdat$variable=factor(mapdat$variable,levels=rev(levels(mapdat$variable)),
                             labels=rev(levels(mapdat$variable)))
      q2=qplot(Reordered_Observation,variable,data=mapdat,geom='tile',fill=Missing,main='Sorted by the number of missing in variables and observations')

      dist31=dist(Mapdat[,-(m$n+1)])
      order31=hclust(dist31,"ward.D")$order
      Mapdat3=Mapdat[order31,-(m$n+1)]
      dist32=dist(t(Mapdat[,-(m$n+1)]))
      order32=hclust(dist32,"ward.D")$order
      Mapdat3=Mapdat3[,order32]
      Mapdat3$Reordered_Observation=1:nrow(Mapdat3)
      mapdat=melt.data.frame(Mapdat3[,],"Reordered_Observation")
      colnames(mapdat)[3]="Missing"
      mapdat$variable=factor(mapdat$variable,levels=rev(levels(mapdat$variable)),
                             labels=rev(levels(mapdat$variable)))
      q3=qplot(Reordered_Observation,variable,data=mapdat,geom='tile',fill=Missing,main='Sorted by hierarchical clustering of missingness')

      return(list(q1=q1,q2=q2,q3=q3))
  }

  #####------------------------------------------------------#####
  ##  VariableOptions is the handler when double clicking gt4.  ##
  ##  It gives a new window for                                 ##
  ##          editing the attributes of variables.              ##
  #####------------------------------------------------------#####
  VariableOptions = function(h, ...) {
    gt11input = gwindow("Attributes", visible = T, width = 300,
                        height = 200, parent = combo1)
    gt11input0 = ggroup(horizontal = FALSE, container = gt11input,
                        expand = TRUE)
    gt11input1 = ggroup(container = gt11input0, expand = TRUE)
    gt11input11 = glabel("Name:", container = gt11input1)
    gt11input12 = gedit(text = svalue(gt11), container = gt11input1, expand = TRUE)

    gt11input2 = ggroup(container = gt11input0, expand = TRUE)
    gt11input21 = glabel("Class:", container = gt11input2)
    gt11input22 = gcombobox(union(gt11[svalue(gt11, index = TRUE),3],
                            c("integer", "numeric", "logical", "character", "factor", "ordered")),
                            container = gt11input2, expand = TRUE)

    gt11input3 = ggroup(container = gt11input0, expand = TRUE)
    gt11input31 = gbutton("Ok", container = gt11input3, expand = TRUE,
                          handler = function(h, ...) {
                            if (svalue(gt11input12) != "") {
                              idx = svalue(gt11, index = TRUE)
                              colnames(m$dataset)[colnames(m$dataset)==as.character(gt11[idx, 2])] = svalue(gt11input12)
                              tmpcolorby = radio125[,]
                              tmpcolorby[tmpcolorby==as.character(gt11[idx, 2])]=svalue(gt11input12)
                              gt11[idx, 2] = svalue(gt11input12)
                              if (svalue(gt11input22) == 'logical') {
                                  if (all(as.character(m$dataset[,idx])
                                          %in% c('0','1','FALSE','TRUE',NA))) {
                                      m$dataset[,idx] = as.logical(as.character(m$dataset[,idx]))
                                  } else {
                                      gmessage("It cannot be converted to a logical variable.")
                                      svalue(gt11input22) = gt11[idx, 3]
                                  }
                              }
                              if (svalue(gt11input22) == 'numeric') {
                                  if (any(tryCatch(na.omit(as.numeric(as.character(m$dataset[,idx]))), warning=function(w){return('warning')})=="warning")) {
                                      gmessage("It cannot be converted to a numeric variable.")
                                      svalue(gt11input22) = gt11[idx, 3]
                                  } else {
                                      m$dataset[,idx] = as.numeric(as.character(m$dataset[,idx]))
                                  }
                              }
                              if (svalue(gt11input22) == 'integer') {
                                  if (any(tryCatch(na.omit(as.integer(as.character(m$dataset[,idx]))), warning=function(w){return('warning')})=="warning")) {
                                      gmessage("It cannot be converted to an integer variable.")
                                      svalue(gt11input22) = gt11[idx, 3]
                                  } else {
                                      m$dataset[,idx] = as.integer(as.character(m$dataset[,idx]))
                                  }
                              }
                              if (svalue(gt11input22) == 'character') {
                                  m$dataset[,idx] = as.character(m$dataset[,idx])
                              }
                              if (svalue(gt11input22) == 'factor') {
                                  m$dataset[,idx] = factor(as.character(m$dataset[,idx]))
                              }
                              if (svalue(gt11input22) == 'ordered') {
                                  m$dataset[,idx] = if (!is.factor(m$dataset[,idx])) {ordered(as.character(m$dataset[,idx]))} else {ordered(m$dataset[,idx])}
                                  lev = levels(m$dataset[,idx])
                                  gt11input31level = gwindow('Levels',visible = T, height=300, parent = gt11input)
                                  gt11input31level0 = ggroup(horizontal = FALSE, container = gt11input31level, expand = TRUE)
                                  gt11input31level1 = glabel("The current levels are:", container=gt11input31level0)
                                  gt11input31level2 = gdf(data.frame(level=lev),name="",container=gt11input31level0, expand=TRUE)
                                  gt11input31level3 = glabel("Redefine the order of the levels:", container=gt11input31level0)
                                  gt11input31level4 = ggroup(horizontal = TRUE, container = gt11input31level0)
                                  gt11input31level41 = glabel("levels =", container = gt11input31level4)
                                  gt11input31level42 = gtext(text=sprintf("c(%s)",paste(1:length(lev),collapse=',')),container=gt11input31level4, height = 20, width = 200)
                                  gt11input31level5 = ggroup(horizontal = TRUE, container = gt11input31level0)
                                  gt11input31level51 = gbutton("Execute", container = gt11input31level5, expand = TRUE, handler = function(h, ...){
                                      neworder = eval(parse(text=svalue(gt11input31level42)))
                                      if (all(sort(neworder)==1:length(lev))) {
                                          levels(m$dataset[,idx])=levels(m$dataset[,idx])[neworder]
                                      } else {
                                          gmessage('Your input cannot be recognized.')
                                      }
                                      dispose(gt11input31level)
                                      dispose(gt11input)
                                  })
                                  gt11input31level52 = gbutton("Cancel", container = gt11input31level5, expand = TRUE, handler = function(h, ...) {
                                      dispose(gt11input31level)
                                  })
                              }
                              gt11[svalue(gt11, index = TRUE), 3] = svalue(gt11input22)
                              check123[,] = gt11[gt11[,3] %in% c('factor','logical','character','ordered'),2]
                              radio125[,] = tmpcolorby
                              gt21[,] = gt11[,]
                              gt3421[,1:3] = gt11[,1:3]
                              gt3421[,4] = mice_default(gt11[,3],dataset[,gt11[,1]])
                              if (svalue(gt11input22) != 'ordered') dispose(gt11input)
                            }
                            else {
                              gmessage("Variable name could not be empty!")
                              svalue(gt11input12) = svalue(gt11)
                            }
                          })
    gt11input32 = gbutton("Cancel", container = gt11input3,
                          expand = TRUE, handler = function(h, ...) {
                            dispose(gt11input)
                          })
  }

  #####----------------------------------------------------#####
  ##  miceMethod is the handler when double clicking gt3421.  ##
  ##  It gives a new window for                               ##
  ##          editing the mice methods of variables.          ##
  #####----------------------------------------------------#####
  miceMethod = function(h, ...) {
    k = svalue(gt3421, index = TRUE)
    if (gt3421[k,4] == '') return()
    gt34input = gwindow("mice methods", visible = T, width = 300,
                        height = 200, parent = combo1)
    gt34input0 = ggroup(horizontal = FALSE, container = gt34input,
                        expand = TRUE)
    gt34input1 = ggroup(container = gt34input0, expand = TRUE, horizontal=FALSE)

    gt34input11 = glabel(paste("Name:",gt3421[k,2],"\nClass:",gt3421[k,3],"\nMethod for mice imputation:"), container = gt34input1)
    possible_methods = list(numeric=c('pmm','norm','norm.nob','norm.boot','norm.predict','mean','quadratic','cart','ri','sample'),
                            level2=c('pmm','2lonly.pmm','logreg','logreg.boot','cart','sample'),
                            factor3=c('pmm','2lonly.pmm','polyreg','lda','cart','sample'),
                            ordered3=c('pmm','2lonly.pmm','polr','cart','sample'))
    s1 = mice_default(gt3421[k,3], dataset[,gt3421[k,1],drop=FALSE])
    s2 = switch(s1,'pmm'=1,'logreg'=2,'polyreg'=3,'polr'=4)
    valid_methods = union(gt3421[k,4], possible_methods[[s2]])
    gt34input12 = gcombobox(valid_methods, container = gt34input1, expand = TRUE)

    gt34input2 = ggroup(container = gt34input0, expand = TRUE)
    gt34input21 = gbutton("Ok", container = gt34input2, expand = TRUE,
                          handler = function(h, ...) {
                            gt3421[k,4]=svalue(gt34input12)
                            dispose(gt34input)
                            })
    gt34input22 = gbutton("Cancel", container = gt34input2, expand = TRUE,
                          handler = function(h, ...) {dispose(gt34input)})
  }

  #####------------------------------#####
  ##  NumSmry is the handler of gb145.  ##
  ##  (gbutton: Numeric Summary)        ##
  #####------------------------------#####
  NumSmry = function(h,...) {
    name_select = svalue(gt11, index = TRUE)
    n = length(name_select)
    if (n == 0) {
        n = nrow(gt11)
        name_select = 1:n
    }
    tmp = m$dataset[,gt11[name_select,2],drop=FALSE]
    cond=check123[svalue(check123,index=T)]

    if (length(cond)==0) {
        l = m$dataset[1,,drop=FALSE]
    } else {
        l = unique(m$dataset[,cond,drop=FALSE])
        ml = apply(l,1,paste,collapse='--')
        tmpid = apply(m$dataset[,cond,drop=FALSE],1,function(d){
            which(ml == paste(unlist(d),collapse='--'))})
    }
    NumSumforMisVal = gwindow("Numeric Summary for Missing Values", visible = T,
                              width = min(400*nrow(l),1200), height = min(200+n*20,800),
                              parent = combo1)
    condgroup0 = ggroup(container = NumSumforMisVal, horizontal = TRUE, expand = TRUE)
    groupN=labelN11=labelN12=labelN13=labelN14=groupN15=tableN15=list()
    for (i in 1:nrow(l)){
        if (length(cond)==0) {
            tmpset = tmp
            groupN[[i]] = ggroup(container = condgroup0, horizontal = FALSE, expand = TRUE)
        } else {
            tmpset = tmp[tmpid==i,,drop=FALSE]
            groupN[[i]] = gframe(text = paste(paste(cond,unlist(l[i,]),sep='='),collapse=','),
                                 container = condgroup0, horizontal = FALSE, expand = TRUE)
        }
        numsummary = compute_missing_pct(tmpset)
        labelN11[[i]] = glabel('Missing:', container=groupN[[i]], pos=0)
        labelN12[[i]] = glabel(paste("    ",round(numsummary$totalmissingpct*100,2),
                                     "% of the numbers",sep=""), container=groupN[[i]])
        labelN13[[i]] = glabel(paste("    ",round(numsummary$varmissingpct*100,2),
                                     "% of variables",sep=""), container=groupN[[i]])
        labelN14[[i]] = glabel(paste("    ",round(numsummary$casemissingpct*100,2),
                                     "% of samples",sep=""), container=groupN[[i]])

        groupN15[[i]]= ggroup(container = groupN[[i]], use.scrollwindow = TRUE,
                              horizontal = FALSE, expand = TRUE)
        tableN15[[i]] = gtable(numsummary$missingsummary, container=groupN15[[i]], expand = TRUE)
    }
  }

  #####----------------------------#####
  ##  Graph is the handler of gb144.  ##
  ##  (gbutton: Watch the data)       ##
  #####----------------------------#####
  Graph = function(h,...) {
    # graphics.off()

    if (!exists('k', envir = m)) m$k = 0
    # l = length(glay15[1,1])
    # if (l>0) for (i in 1:l) dispose(glay15[1,1])

    initializ()
    if(initial_plot()) return()

    if (m$graphtype == "Missingness Map") {
        lab = "Map"; k=1
    } else {
        lab = names(m$dat); k = length(m$dat)
    }

    for (j in 1:k) {
        m$group151[[j+m$k]] = ggroup(container = glay15, label = lab[j], expand = TRUE, horizontal = FALSE)
        m$glay151[[j+m$k]] = glayout(container = m$group151[[j+m$k]], use.scrollwindow = TRUE)
    }

    m$size = size(combo1)*c(.1,1)-c(10,255)

    if (m$graphtype=="Histogram/Barchart") {
      for (j in 1:k) {
        size(m$glay151[[j+m$k]])=c(m$size[1],300*m$n)
          for (i in 1:m$n) {
              m$glay151[[j+m$k]][i, 1, expand = TRUE] = ggraphics(container = m$glay151[[j+m$k]])
              graph_hist(j, i)
          }
      }
    }
    if (m$graphtype=="Spinogram/Spineplot") {
        for (j in 1:k) {
          size(m$glay151[[j+m$k]])=c(m$size[1],300*m$n)
            for (i in 1:m$n) {
                m$glay151[[j+m$k]][i, 1, expand = TRUE] = ggraphics(container = m$glay151[[j+m$k]])
                graph_spine(j, i)
            }
        }
    }
    if (m$graphtype=="Pairwise Plots") {
        for (j in 1:k) {
          size(m$glay151[[j+m$k]])=m$size
            m$glay151[[j+m$k]][1, 1, expand = TRUE] = ggraphics(container = m$glay151[[j+m$k]])
          if (m$n>2 && m$imp_method %in% c('Below 10%','Simple') && lab[j]!='Random Value' ) {
            graph_pair(j, 'bottom', type='plot', contour='jitter')
          } else {
            graph_pair(j, 'bottom', type='plot')
          }
        }
    }
    if (m$graphtype=="Parallel Coordinates") {
        for (j in 1:k) {
          size(m$glay151[[j+m$k]])=c(m$size[1],900)
            if (!graph_pcp(j,title='Original variable order')) {
                m$glay151[[j+m$k]][1, 1, expand = TRUE] = ggraphics(container = m$glay151[[j+m$k]])
                print(m$p+theme(legend.position='bottom'))
            }
            if (!graph_pcp(j,order='allClass',title='Variables ordered by their variation between colors')) {
                m$glay151[[j+m$k]][2, 1, expand = TRUE] = ggraphics(container = m$glay151[[j+m$k]])
                print(m$p+theme(legend.position='bottom'))
            }
        }
    }
    if (m$graphtype=="Missingness Map"){
      size(m$glay151[[1+m$k]])=c(m$size[1],1200)
        q = graph_map()
        m$glay151[[1+m$k]][1, 1, expand = TRUE] = ggraphics(container = m$glay151[[1+m$k]])
        print(q$q1)
        m$glay151[[1+m$k]][2, 1, expand = TRUE] = ggraphics(container = m$glay151[[1+m$k]])
        print(q$q2)
        m$glay151[[1+m$k]][3, 1, expand = TRUE] = ggraphics(container = m$glay151[[1+m$k]])
        print(q$q3)
    }
    # svalue(glay15[1,1])=1
    m$k = m$k + k
  }

  #####---------------------------------#####
  ##  ExportData is the handler of gb146.  ##
  ##  (gbutton: Export the data)           ##
  #####---------------------------------#####
  ExportData = function(h,...){

      initializ()

      if (m$n == 0) {
          svalue(gt11) = 1:nrow(gt11)
          m$name_select = svalue(gt11, index = TRUE)
          m$n = length(m$name_select)
          gmessage("All variables are selected to export.")
      }

      dat = imputation(origdata=m$dataset[,gt11[m$name_select,2],drop=FALSE],
                       method=m$imp_method, vartype=as.character(gt11[m$name_select,3]),
                       missingpct=as.numeric(as.character(gt11[m$name_select,4])),
                       condition=m$cond,knn=m$nn_k,mi.n=m$mi_n,mi.seed=m$mi_seed)
      dat = lapply(dat, function(x) {
          x = data.frame(x[,gt11[m$name_select,2]],
                         is.na(m$dataset[x[,ncol(x)],gt11[m$name_select,2]]))
          colnames(x) = c(gt11[m$name_select,2],paste('Missing', gt11[m$name_select,2], sep='_'))
          return(x)
          })

      entire_dat = dat
      for (j in 1:length(dat)) {
          entire_dat[[j]] = data.frame(m$dataset,is.na(m$dataset))
          colnames(entire_dat[[j]])=c(colnames(m$dataset),paste('Missing',colnames(m$dataset),sep='_'))
          entire_dat[[j]][,colnames(dat[[j]])] = dat[[j]]
      }

    ExpData = function(opa,opb){
        opa = svalue(opa)
        opb = svalue(opb)
        if (opa=='All columns' && opb) return(entire_dat)
        if (opa=='Selected columns' && opb) return(dat)
        if (opa=='All columns' && !opb) return(lapply(entire_dat, function(x) x[,1:ncol(m$dataset)]))
        if (opa=='Selected columns' && !opb) return(lapply(dat, function(x) x[,1:m$n]))
    }

    gExport = gwindow("Export Options", visible = T, width = 300, height = 200, parent = combo1)
    ExGroup = ggroup(container = gExport, expand = TRUE, horizontal = FALSE)
    ExFrame = gframe(text = "Export", container = ExGroup)
    ExRadio = gradio(c('All columns','Selected columns'), container = ExFrame)
    ExCheck = gcheckbox(text = "Export the imputed data with a shadow matrix which consists of TRUE's and FALSE's.\nTRUE's represent for the missing values, and FALSE's represent for non-missings.", checked=TRUE, container = ExGroup)
    ExGroupB = ggroup(container = ExGroup, expand = FALSE, horizontal = TRUE)
    ExLabel2 = glabel(text = "Directory  ", container = ExGroupB)
    ExpText2 = gtext(text = getwd(), width=400, height=20, container = ExGroupB)
    ExButtoB = gbutton(text = "Browse...", container = ExGroupB, handler = function(h,...){
      ExFile = gfile(text = 'filename', type = 'selectdir')
      if (!is.na(ExFile)) svalue(ExpText2) = ExFile
    })
    ExGroupT = ggroup(container = ExGroup, expand = FALSE, horizontal = TRUE)
    ExLabel1 = glabel(text = "File name ", container = ExGroupT)
    ExpText1 = gtext(text = "", width=400, height=30, container = ExGroupT)
    ExGroupB = ggroup(container = ExGroup, expand = TRUE, horizontal = TRUE)
    addSpace(ExGroupB, 20)
    ExButto1 = gbutton(text = "Save as .csv format", container = ExGroupB, handler=function(h,...){
        fexdat = ExpData(ExRadio,ExCheck)
        filename = paste(svalue(ExpText2),'/',gsub('\\.csv$','',svalue(ExpText1)),sep='')
        fileSuffix = paste('_impute_',gsub('(MI:)|(%)|( )',"",m$imp_method),'_',gsub('[^a-zA-Z0-9]',"",names(fexdat)),sep='')
        for (j in 1:length(fexdat)) {
            write.csv(fexdat[[j]], file=paste(filename,fileSuffix[j],'.csv',sep=''), row.names=FALSE)
        }
        gmessage("The data is exported!")
        dispose(gExport)
    })
    addSpace(ExGroupB, 20)
    ExButto2 = gbutton(text = "Save as .rda format", container = ExGroupB, handler=function(h,...){
        fexdat = ExpData(ExRadio,ExCheck)
        filename = paste(svalue(ExpText2),'/',gsub('\\.rda$','',svalue(ExpText1)),sep='')
        fileSuffix = paste('_impute_',gsub('[^a-zA-Z0-9]',"",m$imp_method),sep='')
        save(fexdat, file=paste(filename,fileSuffix,'.rda',sep=''))
        gmessage("The data is exported!")
        dispose(gExport)
    })
    addSpace(ExGroupB, 20)
    ExButto3 = gbutton(text = "Save as a list of data frames", container = ExGroupB, handler=function(h,...){
        filename = svalue(ExpText1)
        fexdat = ExpData(ExRadio,ExCheck)
        eval(parse(text=paste(filename,'<<- fexdat')))
        gmessage(paste("An R data frame is saved with the name ",filename))
        dispose(gExport)
        })
  }

  #####---------------------------------#####
  ##  ExportData is the handler of gb148.  ##
  ##  (gbutton: Save the plot)             ##
  #####---------------------------------#####
  SavePlot = function(h,...){
    initializ()
    if(initial_plot()) return()
    if (m$graphtype == "Missingness Map") {
        lab = "Map"; k=1
    } else {
        lab = gsub("[^A-Za-z0-9]","",names(m$dat)); k = length(m$dat)
    }
    savename = gfile(type="save")
    if (m$graphtype=="Histogram/Barchart") {
      for (j in 1:k) {
          for (i in 1:m$n) {
              png(filename = paste(savename,'_',lab[j],'_hist_',i,'.png',sep=''),width = 7, height = 5,units = "in", res=90)
              graph_hist(j, i)
              dev.off()
          }
      }
    }
    if (m$graphtype=="Spinogram/Spineplot") {
      for (j in 1:k) {
          for (i in 1:m$n) {
              png(filename = paste(savename,'_',lab[j],'_spinogram_',i,'.png',sep=''),width = 7, height = 5,units = "in", res=90)
              graph_spine(j, i)
              dev.off()
          }
      }
    }
    if (m$graphtype=="Pairwise Plots") {
      for (j in 1:k) {
          png(filename = paste(savename,'_',lab[j],'_pairwise.png',sep=''), width = 2*m$n, height = 2*m$n, units = "in", res=90)
          if (m$n>2 && m$imp_method %in% c('Below 10%','Simple') && lab[j]!='Random Value' ) {
            graph_pair(j, 'right', type='save', contour='jitter')
          } else {
            graph_pair(j, 'right',type='save')
          }
          dev.off()
      }
    }
    if (m$graphtype=="Parallel Coordinates") {
        for (j in 1:k) {
            if (!graph_pcp(j)) {
                png(filename = paste(savename,'_',lab[j],'_pcp_orig.png',sep=''),width = (m$n+2), height = 4, units = "in", res=90)
                print(m$p+theme(legend.position='bottom'))
                dev.off()
            }
            if (!graph_pcp(j,order='allClass')) {
              png(filename = paste(savename,'_',lab[j],'_pcp_sort.png',sep=''),width = (m$n+2), height = 4, units = "in", res=90)
              print(m$p+theme(legend.position='bottom'))
              dev.off()
            }
        }
    }
    if (m$graphtype=="Missingness Map"){
        q = graph_map()
        png(filename = paste(savename,'_map_1.png',sep=''),width = 6, height = max(4,round(m$n/8)), units = "in", res=90)
        print(q$q1)
        dev.off()

        png(filename = paste(savename,'_map_2.png',sep=''),width = 6, height = max(4,round(m$n/8)), units = "in", res=90)
        print(q$q2)
        dev.off()

        png(filename = paste(savename,'_map_3.png',sep=''),width = 6, height = max(4,round(m$n/8)), units = "in", res=90)
        print(q$q3)
        dev.off()
    }
  }

  #####-------------------------------#####
  ##  New window for missing values      ##
  #####-------------------------------#####
  size.width = round(max(c(size.width, 700)))
  size.height = round(max(c(size.height, 450)))
  if (size.height/size.width > 0.8) size.height = size.width * 0.8
  if (size.width/size.height > 2) size.width = size.height * 2
  combo1 <- gwindow("Missing Values", visible = T, width = size.width, height = size.height)
  tab <- gnotebook(container = combo1)

  #####------------------------------------------------#####
  ##  In the first tab we can:                            ##
  ##  (1) Watch and change the name or type of variables. ##
  ##  (2) Numeric or graphic summary.                     ##
  ##  (3) Save the imputation and plots.                  ##
  #####------------------------------------------------#####
  group11 = ggroup(container = tab, label = "Summary", expand = TRUE, horizontal = FALSE)
  group1100 = ggroup(container = group11, expand = TRUE)
  group12 = ggroup(container = group1100, use.scrollwindow = TRUE,
                   horizontal = FALSE, expand = TRUE)
  g12width = min(max(nchar(vname))*3+200,round(size.width/2))
  size(group12) = c(g12width, size.height-25)
  nametable = data.frame(ID=1:length(vname), Variables=as.character(vname),
                         Class=as.character(dataclass), NApct=as.character(round(vNApct,3)),
                         stringsAsFactors=FALSE)
  gt11 = gtable(nametable, multiple = T, container = group12,
                expand = TRUE, chosencol = 2)
  addhandlerdoubleclick(gt11, handler = VariableOptions)

  label121 = glabel('Categorical variables to condition on',container=group12)
  check123 = gcheckboxgroup(nametable$Variables[nametable$Class %in%
    c('factor','logical','character')], container=group12, use.table=TRUE)
  size(check123) = c(g12width-5,round(size.height/5))

  group13 = ggroup(horizontal = FALSE, container = group1100, expand = TRUE)
  group14 = ggroup(horizontal = TRUE, container = group13)
  size(group14) = c(min(round(0.6*size.width),size.width-g12width-20), 160)
  tmpcolorby = data.frame(`Color by the missing of`= c('Missing Any Variables',
                                                       'Missing on Selected Variables',
                                                       nametable[vNApct>0,2]))
  tmpcolorby[,1]=as.character(tmpcolorby[,1])
  radio125 = gtable(tmpcolorby, container=group14, expand=TRUE, multiple=TRUE)
  addHandlerKeystroke(radio125, handler = function(h,...){})
  gframe142 = gframe(text = "Method", container = group14)
  gr142 = gradio(c('Below 10%','Simple','Neighbor','MI:areg',
                   'MI:norm','MI:mice','MI:mi'),
                 container = gframe142, handler = function(h,...){
                   if (svalue(gr142)=='Below 10%') {
                     svalue(check123) = FALSE
                   }
                 })
  if (is_imputed_data) {enabled(gr142) = FALSE}
  gframe143 = gframe(text = "Graph Type", container = group14)
  gr143 = gradio(c('Histogram/Barchart','Spinogram/Spineplot','Pairwise Plots',
                   'Parallel Coordinates','Missingness Map'), container = gframe143)

  group144 = ggroup(horizontal = FALSE, container = group14)
  gb145 = gbutton('Summary', container = group144,
                  handler = NumSmry)
  gb144 = gbutton("P l o t", container = group144,
                  handler = Graph)
  gb146 = gbutton('Export data', container = group144,
                  handler = ExportData)
  gb148 = gbutton('Save plot', container = group144,
                  handler = SavePlot)
  gb147 = gbutton('Q u i t', container = group144,
                  handler = function(h,...){
                    dispose(combo1)
                  })

  group15 = ggroup(horizontal = FALSE, container = group13,
                   expand = TRUE, use.scrollwindow = TRUE)
  glay15 = gnotebook(container = group15, closebuttons =TRUE, expand=TRUE)

  #####------------------------------------------------#####
  ##  In the second tab we can:                           ##
  ##  Look at the help documents.                         ##
  #####------------------------------------------------#####
  group21 = ggroup(container = tab, label = "Help", expand = TRUE, horizontal = FALSE)
  group2100 = ggroup(container = group21, expand = TRUE)
  group22 = ggroup(container = group2100, use.scrollwindow = TRUE,
                   horizontal = FALSE, expand = T)
  size(group22) = c(g12width, size.height-25)
  gt21 = gtable(nametable, multiple = T, container = group22,
                expand = TRUE, chosencol = 2)
  addHandlerMouseMotion(gt21, handler = function(h,...){
    if (exists('text25')) svalue(text25) = capture.output(cat("\n\n   This table displays all variables in the data set and reports the type and the percentage of missing values for each variable.\n\n   Users can sort the variables by NA's percent.\n\n   Doubleclicking one row can change the variable name, as well as the data type.\n\n   The table allows text entry to find a variable."))
	})

  label221 = glabel('Categorical variables to condition on',container=group22)
  check223 = gcheckboxgroup(nametable$Variables[nametable$Class %in%
                            c('factor','logical','character')],
                            container=group22, use.table=TRUE)
  size(check223) = c(g12width-5,round(size.height/5))
  addHandlerMouseMotion(check223, handler = function(h,...){
    if (exists('text25')) svalue(text25) = capture.output(cat("\n\n   This list displays all categorical variables. We can make multiple selection on them.\n\n   Once we select one or more variables, the data set will be divided into the categories. The numeric summary and imputation will be made in each category.\n\n   If the imputation method is 'Below 10%', then the selected conditioning variables are ignored.\n\n   The list allows text entry to find a variable."))
	})

  group23 = ggroup(horizontal = FALSE, container = group2100,
                   expand = TRUE)
  group24 = ggroup(horizontal = TRUE, container = group23)
  size(group24) = c(min(round(0.6*size.width),size.width-g12width-20), 160)
  help_colorlist = function(h,...){
      if (exists('text25')) svalue(text25) = capture.output(cat("\n\n   This list displays all variables which have missing values.\n\n   If the user chooses one of them, the color of the plot showing on the right will change based on whether the cases being NA on that variable or not.\n\n   The user can also choose several variables. Then the color of the plot will be based on whether the cases have missing values on any of those variable.\n\n   The first row 'Missing Any Variables' means the plot is colored by whether the observation is complete or not, regardless of the selected variables.\n\n   The second row 'Missing on Selected Variables' means whether the observations have missing values on any of the selected variable.\n\n   The widget allows text entry to find a particular variable if the list is quite long."))
  }
  radio225 = gtable(data.frame(`Color by the missing of`= c('Missing Any Variables',
                    'Missing on Selected Variables', nametable[vNApct>0,2])),
                    container = group24, expand = TRUE, handler = help_colorlist)
  addHandlerMouseMotion(radio225, handler = help_colorlist)

  gframe242 = gframe(text = "Method", container = group24)
  help_methods = function(h,...){
      if (exists('text25')) svalue(text25) = capture.output(cat("\n\n   This list displays all the imputation methods.\n\n   Users can make one selection.\n\n      (1) 'Below 10%' means NA's of one variable will be replaced by the value which equals to the minimum of the variable minus 10% of the range. For categorical variables, NA's are treated as a new category. Under this status the selected conditioning variables are ignored. If the data are already imputed, then this item will show the imputed result.\n\n      (2) 'Simple' will create three tabs: Median, Mean, and Random Value. 'Median' means NA's will be replaced by the median of this variable (omit NA's). 'Mean' means NA's will be replaced by the mean of the variable (omit NA's). The median does not apply to the nominal variable, neither does the mean to the categorical variable. In these cases the mode (omit NA's) is provided. 'Random Value' means NA's will be replaced by any values of this variable (omit NA's) which are randomly selected. \n\n      (3) 'Neighbor' contains two methods: 'Average Neighbor' and 'Random Neighbor'. 'Average Neighbor' will replace the NA's by the mean of the nearest neighbors. 'Random Neighbor' substitutes the missing for a random sample of the k nearest neighbors. The number of neighbors is default to 5, and editable in the Settings tab. The nearest neighbor method requires at lease one case to be complete, at least two variables to be selected, and no factor/character variables. The ordered factors are treated as integers. The method will return the overall mean or a global random sample value if the observation only contains NA's. \n\n     For all the multiple imputation methods below, the number of imputed sets (default to be 3) and random number seed can be set in the Settings tab.\n\n      (4) 'MI:areg' uses function 'aregImpute' from package 'Hmisc'. It requires at lease one case to be complete, and at least two variables to be selected.\n\n    (5) 'MI:norm' uses function 'imp.norm' from package 'norm'. It requires all selected variables to be numeric(at least integer), and at least two variables to be selected. Sometimes it cannot converge, then the programme will leave NA's without imputation.\n\n      (6) 'MI:mice' uses the mice package. In the Settings tab, the imputing methods for all variables are displayed. Doubleclicking the variable will allow the user to change the method.\n\n    (7) 'MI:mi' employes the mi package.\n\n "))
  }
  gr242 = gradio(c('Below 10%','Simple','Neighbor','MI:areg',
                   'MI:norm','MI:mice','MI:mi'),
                 container = gframe242, handler = help_methods)
  addHandlerMouseMotion(gr242, handler = help_methods)

  gframe243 = gframe(text = "Graph Type", container = group24)
  help_plottype = function(h,...){
      if (exists('text25')) svalue(text25) = capture.output(cat("\n\n   This frame shows all plots we can make.\n\n      (1)'Histogram/Barchart' will display histograms (numeric variables) and barcharts(categorical variables) for each variable selected.\n\n      (2)'Spinogram/Spineplot' shows the spineplot for each selected variable.\n\n      (3)'Pairwise Plots' will generate n*(n-1) plots if n variables are selected.\n\n       For two quantitative variables, scatterplots are in the lower matrix and contour density plots are in the upper matrix. For one categorical variable and one quantitative variable, histograms are in the lower matrix and boxplots are in the upper matrix. Two categorical variables will give a barchart.\n\n       When n>5, then only the first 5 variables are displayed. When n=2, a scatterplot is drawn for the two variables. When 2<n<=5, the function 'ggpairs' from package 'GGally' is used.\n\n      (4)'Parallel Coordinates' displays two parallel coordinates plots for the selected variables. One shows the axes in the original order, and the other in an optimized order, determined by whether the missings and non-missings can be separated well by the variable. (order='allClass' is used. See ?GGally::ggparcoord)\n\n       Whether jittered on the categorical variables can be set in the Settings tab.\n\n      (5)'Missingness Map' shows the positions of missing values in all the observations from the variable selected, regardless of the imputation. Three maps are given, and the variables and observations are sorted as the titles say."))
  }
  gr243 = gradio(c('Histogram/Barchart','Spinogram/Spineplot','Pairwise Plots',
                   'Parallel Coordinates','Missingness Map'), container = gframe243,
                 handler = help_plottype)
  addHandlerMouseMotion(gr243, handler = help_plottype)

  group244 = ggroup(horizontal = FALSE, container = group24)
  help_numeric_summary = function(h,...){
      if (exists('text25')) svalue(text25) = capture.output(cat("\n\n   Clicking this button will create another window which presents the numeric summaries for missing values.\n\n   In this summary window, the missing percentage of all the numbers, variables, and cases are presented.\n\n   Besides, there is a table of the missing levels. The table has n+1 rows, where n = # of selected variables. For each i in 0:n, the table gives the count of cases which have i missing values, as well as the percentage of those cases.\n\n   All the summaries and tables can be conditioned on the categorical variables."))
  }
  gb245 = gbutton('Summary', container = group244, handler = help_numeric_summary)
  addHandlerMouseMotion(gb245, handler = help_numeric_summary)

  help_plot = function(h,...){
      if (exists('text25')) svalue(text25) = capture.output(cat("\n\n   Clicking this button will draw a plot based on the data and options the user chooses."))
  }
  gb244 = gbutton("P l o t", container = group244, handler = help_plot)
  addHandlerMouseMotion(gb244, handler = help_plot)

  help_export = function(h,...){
      if (exists('text25')) svalue(text25) = capture.output(cat("\n\n   Clicking this button will export the imputed data based on the other chosen options.\n\n   A window interface will pop up to let the user change some exporting options. The imputation method will be completed in the file name automatically.\n\n   The imputed data can be saved in csv or rda formats, or to a list in the R console."))
  }
  gb246 = gbutton('Export data', container = group244, handler = help_export)
  addHandlerMouseMotion(gb246, handler = help_export)

  help_save = function(h,...){
      if (exists('text25')) svalue(text25) = capture.output(cat("\n\n   Clicking this button will save the plot(s) to png file(s) based on the options the user chooses.\n\n   A user can define the file name for the plot, and the graph type will be suffixed automatically."))
  }
  gb248 = gbutton('Save plot', container = group244, handler = help_save)
  addHandlerMouseMotion(gb248, handler = help_save)

  help_quit = function(h,...){
      if (exists('text25')) svalue(text25) = capture.output(cat("\n\n   Clicking this button will destroy the main window."))
  }
  gb247 = gbutton('Q u i t', container = group244, handler = help_quit)
  addHandlerMouseMotion(gb247, handler = help_quit)


  group25 = ggroup(horizontal = FALSE, container = group23,
                   expand = TRUE, use.scrollwindow = TRUE)
  text25 = gtext(container = group25, expand = TRUE,
                 use.scrollwindow = TRUE, font.attr=c(family="monospace"))

  #####------------------------------------------------#####
  ##  In the third tab we can:                            ##
  ##  (1) Change the default settings of the methods.     ##
  #####------------------------------------------------#####
  group31 = ggroup(container = tab, label = "Settings", expand = TRUE, horizontal=FALSE)
  group32 = glayout(container = group31, expand = TRUE, spacing = 10)

  group32[1,1] = ggroup(container = group32, expand = TRUE, horizontal=FALSE)
  frame342 = gframe(text = "MI:mice", container = group32[1,1], horizontal=FALSE)

  miceSettings = data.frame(nametable[,1:3], Method=mice_default(as.character(dataclass),dataset),
                            stringsAsFactors=FALSE)
  gt3421 = gtable(miceSettings, multiple = T, container = frame342, expand = FALSE, chosencol = 2)
  size(gt3421) = c(round(0.45 * size.width), size.height-30)
  addhandlerdoubleclick(gt3421, handler = miceMethod)

  group32[1,2] = ggroup(container = group32, expand=FALSE, horizontal=FALSE)

  frame3211 = gframe(text = "Multiple imputation",
                    container = group32[1,2], horizontal=FALSE)
  group3211 = ggroup(container = frame3211)
  label32111 = glabel(text="Number of imputed sets :  ", container=group3211)
  text32112 = gedit(text="3", container=group3211, width=2, coerce.with=as.integer)
  group3212 = ggroup(container = frame3211)
  label32121 = glabel(text="Random number seed :  ", container=group3212)
  text32122 = gedit(text="1234567", container=group3212, width=8, coerce.with=as.integer)

  label321 = glabel(container = group32[1,2])

  frame3221 = gframe(text = "Neighbor:",
                    container = group32[1,2], horizontal=FALSE)
  group3221 = ggroup(container = frame3221)
  label32211 = glabel(text="Number of neighbors :  ", container=group3221)
  text32212 = gedit(text="5", container=group3221, width=2, coerce.with=as.integer)

  label322 = glabel(container = group32[1,2])

  frame3231 = gframe(text = "Parallel coordinates plot",
                        container = group32[1,2], horizontal=FALSE)
  radio32311 = gradio(c('Jitter on categorical variables','No jitter'),container = frame3231)

  svalue(tab)=1
}

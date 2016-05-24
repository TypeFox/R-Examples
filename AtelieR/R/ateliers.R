#-----------------------------------------------------------------------------------------
#
#               Yvonnick Noel, U. of Brittany, Rennes, France, 2007-2011
#                        Statistical workshops for teaching
#
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
#                           Construction of main GUI
#-----------------------------------------------------------------------------------------

.ws = proto(

   create = function(.) {
   
     # Interface
    .$version = packageDescription("AtelieR")$Version
    .$window = gwindow(.$translate("AtelieR: Statistical workshops in R"),visible=FALSE)
    .$bigGroup = ggroup(cont = window)
    .$statusBar = gstatusbar(.$translate("Ready."),cont=.$window)
     add(.$bigGroup,.$nb <- gnotebook(tab.pos=2,closebuttons=TRUE),expand=TRUE)
     add(.$bigGroup,ggraphics())

     # Menu actions
     aClose      = gaction(label=.$translate("Quit"),icon="quit", handler=function(h,...)  dispose(.$window))
     aNormal     = gaction(label=.$translate("Construction of the gaussian distribution"), handler=.ws1$create)
     aScale      = gaction(label=.$translate("Change of origin and scale"),                handler=.ws2$create)
     aMean       = gaction(label=.$translate("Distribution of a sample mean"),             handler=.ws3$create)
     aVar        = gaction(label=.$translate("Distribution of a sample variance"),         handler=.ws5$create)
     aCalc       = gaction(label=.$translate("Probability calculator"),                    handler=.ws4$create)
     aProp       = gaction(label=.$translate("On a proportion"),                           handler=.ws6$create)
     aKprop      = gaction(label=.$translate("On several proportions"),                    handler=.ws9$create)
     aCatdist    = gaction(label=.$translate("On a categorical distribution"),             handler=.ws16$create)
     aTable      = gaction(label=.$translate("On a contingency table"),                    handler=.ws10$create)
     aBayesvar   = gaction(label=.$translate("On a sample variance"),                      handler=.ws11$create)
     aBayesmean  = gaction(label=.$translate("On a sample mean"),                          handler=.ws8$create)
     aBayesboth  = gaction(label=.$translate("On both a mean and a variance"),             handler=.ws12$create)
     # aBayesdelta = gaction(label=.$translate("On a standardized effect"),                  handler=.ws14$create)
     # a2means     = gaction(label=.$translate("On two means"),                              handler=.ws15$create)
     aKmeans     = gaction(label=.$translate("On several means"),                          handler=.ws13$create)
     aAbout      = gaction(label=.$translate("About..."),                                  handler=.$aboutAtelieR)
     
     # Menu tree
     tmp = list(Session = list(Quit=aClose),
                Modules = list(Understand = list(normal=aNormal,scale=aScale,mean=aMean,var=aVar),
                               Compute    = list(calc=aCalc),
                               "Bayesian inference" = list(prop    = aProp,
                                                           kprop   = aKprop,
                                                           catdist = aCatdist,
                                                           tab     = aTable,
                                                           bvar    = aBayesvar,
                                                           bmean   = aBayesmean,
                                                           bboth   = aBayesboth,
                                                           # bdelta  = aBayesdelta,
                                                           # m2      = a2means,
                                                           km      = aKmeans)),
                Help    = list(About=aAbout))
                
     names(tmp$Session) = .$translate(names(tmp$Session))
     names(tmp$Modules) = .$translate(names(tmp$Modules))
     names(tmp) = .$translate(names(tmp))
    .$menu = gmenu(tmp,cont=.$window)
   
   },
   
   show = function(.) {
     svalue(.$nb) = 1
     visible(.$window) = TRUE
   },
   
   setStatus = function(.,text) {
     svalue(.$statusBar)
     svalue(.$statusBar) = text
   },

   ### Gettext utility for translating messages
   translate = function(.,...) {
     gettext(..., domain="R-AtelieR")
   },

  aboutAtelieR = function(.,h,...) {
    aboutMessage = gbasicdialog(title="About...",do.buttons=FALSE)
    messageFrame = gframe(cont=aboutMessage,horizontal=FALSE)
    add(messageFrame,glabel("<span foreground='blue' size='x-large' weight='ultrabold'>AtelieR</span>",markup=TRUE))
    add(messageFrame,glabel(paste("version",.$version)))
    add(messageFrame,glabel("\n   A GTK GUI for elementary Bayesian statistics   \n"))
    add(messageFrame,glabel("<b>Yvonnick Noel</b>",markup=TRUE))
    add(messageFrame,glabel("University of Brittany at Rennes, France"))
    add(messageFrame,glabel("<i>yvonnick.noel@uhb.fr</i>\n",markup=TRUE))
    visible(aboutMessage,set=TRUE)
  },

  #---------------------------------------------------------------------------------------
  #  SLOT       INITIAL VALUE                CONTENT
  #---------------------------------------------------------------------------------------
  window      = NULL,                   # Main window
  bigGroup    = NULL,                   # Main group
  nb          = NULL,                   # Main notebook
  menu        = NULL,                   # Menu
  version     = NULL,                   # package version
  statusBar   = NULL                    # Status bar
 )

AtelieR = function() {

  .ws$create()
  .ws$show()

}

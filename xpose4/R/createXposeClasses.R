# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.

"createXposeClasses" <-
  function(nm7=F) {

    setClassUnion("character or NULL",c("character","NULL"),where=.GlobalEnv
)
    setClassUnion("character or numeric",c("character","numeric"),where=.GlobalEnv)
    setClassUnion("numeric or NULL",c("numeric","NULL"),where=.GlobalEnv)
    setClassUnion("data.frame or NULL",c("data.frame","NULL"),where=.GlobalEnv)
    setClassUnion("list or NULL",c("list","NULL"),where=.GlobalEnv)
    setClassUnion("lang or numeric",c("vector","numeric","list"),where=.GlobalEnv)
    setClassUnion("logical or numeric",c("logical","numeric"),where=.GlobalEnv)

    if(nm7) ipred.def <- "IPRED" else ipred.def <- "IPRE"
    if(nm7) iwres.def <- "IWRES" else iwres.def <- "IWRE"

    if(nm7) {
      labels.list  = list(
        OCC  = "Occasion",
        TIME = "Time",
        PRED = "Population predictions",
        IPRED = "Individual predictions",
        WRES = "Weighted residuals",
        CWRES = "Conditional weighted residuals",
        IWRES = "Individual weighted residuals",
        DV   = "Observations",
        RES  = "Residuals",
        CL = "Clearance",
        V  = "Volume",
        TAD  = "Time after dose"
        )
    } else {
      labels.list  = list(
        OCC  = "Occasion",
        TIME = "Time",
        PRED = "Population predictions",
        IPRE = "Individual predictions",
        WRES = "Weighted residuals",
        CWRES = "Conditional weighted residuals",
        IWRE = "Individual weighted residuals",
        DV   = "Observations",
        RES  = "Residuals",
        CL = "Clearance",
        V  = "Volume",
        TAD  = "Time after dose"
        )
    }


    setClass("xpose.prefs",where=.GlobalEnv,
             representation(Xvardef       = "list",
                            Labels        = "list",
                            Graph.prefs   = "list",
                            Miss          = "numeric",
                            Cat.levels    = "numeric",
                            DV.Cat.levels = "numeric",
                            Subset        = "character or NULL",
                            Gam.prefs     = "list",
                            Bootgam.prefs = "list"
                            ),
             prototype(
                       Xvardef = list(
                         id      = "ID",
                         idlab   = "ID",
                         idv     = "TIME",
                         occ     = "OCC",
                         dv      = "DV",
                         pred    = "PRED",
                         ipred   = ipred.def,
                         iwres   = iwres.def,
                         wres    = "WRES",
                         cwres   = "CWRES",
                         res     = "RES",
                         parms   = c("CL","V","V1","V2","V3","Q","Q1","Q2","Q3","KA",
                           "ETA1","ETA2","ETA3","ETA4","ETA5","ETA6","ETA7",
                           "ETA8","ETA9","ET10","ET11","ET12","ET13","ET14",
                           "ET15","ET16","ET17","ET18","ET19","ET20"),
                         covariates = c("GENO","SEX","RACE","DOSE","FLAG","DAY","PAT",
                           "GEND","AGE","WT","HT","CRCL","CLCR"),
                         ranpar  = c("ETA1","ETA2","ETA3","ETA4","ETA5","ETA6","ETA7",
                           "ETA8","ETA9","ET10","ET11","ET12","ET13","ET14","ET15",
                           "ET16","ET17","ET18","ET19","ET20"),
                         tvparms = c("TVCL","TVV","TVV1","TVV2","TVV3","TVQ","TVQ1",
                           "TVQ2","TVQ3","TVKA")
                         ),

                       Labels  = labels.list,

                       Graph.prefs = list(
                         type   = "b" ,
                         pch    = 1   ,
                         cex    = 0.8 ,
                         lty    = 1   ,
                         lwd    = 1   ,
                         col    = 4   ,
                         fill   = "lightblue",
                         grid   = FALSE ,
                         aspect = "fill"   ,

                         ## By arguments
                         condvar   = NULL,
                         byordfun  = "median" ,
                         ordby     = NULL     ,
                         shingnum  = 6        ,
                         shingol   = 0.5      ,

                         ## Abline settings
                         abline = NULL ,
                         abllwd = 1    ,
                         ablcol = 1    ,
                         abllty = 1,

                         ## Smooth settings
                         smooth = NULL ,
                         smlwd  = 2    ,
                         smcol  = "red" ,
                         smlty  = 1    ,
                         smspan = 2/3  ,
                         smdegr = 1,

                         ## Lm settings
                         lmline = NULL,
                         lmlwd  = 2,
                         lmcol  = 2,
                         lmlty  = 1,

                         ## Superpose line settings
                         suline = NULL,
                         sulwd  = 2,
                         sucol  = 3,
                         sulty  = 1,
                         suspan = 2/3,
                         sudegr = 1,

                         ## Text label settings,
                         ids    = FALSE,
                         idsmode= NULL,
                         idsext = 0.05, ## In each end
                         idscex = 0.7,
                         idsdir = "both",

                         ## Dilution stuff
                         dilfrac = 0.7,
                         diltype = NULL,
                         dilci   = 0.95,

                         ## Prediction interval stuff
                         PIuplty = 2,
                         PIdolty = 2,
                         PImelty = 1,
                         PIuptyp = "l",
                         PIdotyp = "l",
                         PImetyp = "l",
                         PIupcol = "black",
                         PIdocol = "black",
                         PImecol = "black",
                         PIuplwd = 2,
                         PIdolwd = 2,
                         PImelwd = 2,
                         PIupltyR = 1,
                         PIdoltyR = 1,
                         PImeltyR = 2,
                         PIuptypR = "l",
                         PIdotypR = "l",
                         PImetypR = "l",
                         PIupcolR = "blue",
                         PIdocolR = "blue",
                         PImecolR = "blue",
                         PIuplwdR = 2,
                         PIdolwdR = 2,
                         PImelwdR = 2,
                         PIupltyM = 1,
                         PIdoltyM = 1,
                         PImeltyM = 2,
                         PIuptypM = "l",
                         PIdotypM = "l",
                         PImetypM = "l",
                         PIupcolM = "darkgreen",
                         PIdocolM = "darkgreen",
                         PImecolM = "darkgreen",
                         PIuplwdM = 0.5,
                         PIdolwdM = 0.5,
                         PImelwdM = 0.5,
                         PIarcol = "lightgreen",
                         PIlimits=c(0.025,0.975),

                         ## Categorical x-variable
                         bwhoriz = FALSE,
                         bwratio = 1.5,
                         bwvarwid = FALSE,
                         bwdotpch = 16,
                         bwdotcol = "black",
                         bwdotcex = 1,
                         bwreccol = "blue",
                         bwrecfill= "transparent",
                         bwreclty = 1,
                         bwreclwd = 1,
                         bwumbcol = "blue",
                         bwumblty = 1,
                         bwumblwd = 1,
                         bwoutcol  ="blue" ,
                         bwoutcex  = 0.8,
                         bwoutpch  = 1,

                         ##Histogram settings
                         hicol     = 5,#"blue",
                         hiborder  = "black",
                         hilty     = 1,
                         hilwd     = 1,
                         hidlty    = 2,
                         hidlwd    = 2,
                         hidcol    = 1
                         ),

                       Miss       = -99,
                       Cat.levels = 4,
                       DV.Cat.levels = 7,
                       Subset     = NULL,

                       Gam.prefs  = list(
                         onlyfirst=TRUE,
                         wts=FALSE,
                         start.mod=NULL,
                         steppit=TRUE,
                         disp = NULL,
                         nmods=3,
                         smoother1=0,
                         smoother2=1,
                         smoother3="ns",
                         smoother4="ns",
                         arg1=NULL,
                         arg2=NULL,
                         arg3="df=2",
                         arg4="df=3",
                         excl1=NULL,
                         excl2=NULL,
                         excl3=NULL,
                         excl4=NULL,
                         extra=NULL,
                         plot.ids=TRUE,
                         medianNorm=TRUE
                         ),
                       Bootgam.prefs = list(n = 100,
                         algo = "fluct.ratio",
                         conv.value = as.numeric(1.04),
                         check.interval = as.numeric(20),
                         start.check = as.numeric(50),
                         liif = as.numeric(0.2),
                         ljif.conv = as.numeric(25),
                         seed = NULL,
                         start.mod = NULL,
                         excluded.ids = NULL
                         )                  
                       
                       )
             )


    setClass("xpose.data",where=.GlobalEnv,
             representation(Data      = "data.frame or NULL",
                            SData     = "data.frame or NULL",
                            Data.firstonly = "data.frame or NULL",
                            SData.firstonly = "data.frame or NULL",
                            Runno     = "character or numeric",
                            Nsim      = "numeric or NULL",
                            Doc       = "character or NULL",
                            Prefs     = "xpose.prefs"
                            ),
             prototype(Data    = NULL,
                       SData   = NULL,
                       Data.firstonly    = NULL,
                       SData.firstonly   = NULL,
                       Nsim    = NULL,
                       Runno   = NULL,
                       Doc     = NULL),
             validity = test.xpose.data
             )




    invisible()
  }





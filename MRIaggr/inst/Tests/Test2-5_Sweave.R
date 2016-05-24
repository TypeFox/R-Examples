#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%
#%%%%% Fichier test 2.5 : Test de la fonction constSweave
#%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
require(MRIaggr)

options(error=function() traceback(2)) 
options(max.print=10000)

####  constSweave ####
# ?constSweave
# constSweave
# function (dir, identifier = NULL, param = NULL, table = NULL, 
#           extra_text = NULL, subsection = NULL, index_subsection = NULL, 
#           subsubsection = NULL, index_subsubsection = NULL, legend = NULL, 
#           trace = FALSE, width = list(0.9, 0.9, 0.9), trim = list(c(0,0, 0, 0), c(0, 0, 0, 0), c(0, 160, 0, 0)), width.legend = 0.35, 
#           trim.legend = c(0, 0, 0, 0), title = "", date = "", author = "") 

#### test

#### example 

## load an MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## directories for storage
if(("Display" %in% list.files())  == FALSE){dir.create("Display")}
if(("DWI" %in% list.files("Display"))  == FALSE){dir.create("Display/DWI")}
if(("DWI_lesion" %in% list.files("Display"))  == FALSE){dir.create("Display/DWI_lesion")}
if(("T2" %in% list.files("Display"))  == FALSE){dir.create("Display/T2")}

## plot generation
multiplot(MRIaggr.Pat1_red,param="DWI_t0",
             window="png",path="Display/DWI/")
multiplot(MRIaggr.Pat1_red,param="DWI_t0",
             index1=list(coords="MASK_T2_FLAIR_t2",outline=TRUE),as.logical=TRUE,
             window="png",path="Display/DWI_lesion/")
multiplot(MRIaggr.Pat1_red,param="T2_FLAIR_t2",
             window="png",path="Display/T2/")

MRIaggr.Pat1_red@identifier <- "Pat2"

multiplot(MRIaggr.Pat1_red,param="DWI_t0",
             window="png",path="Display/DWI/")
multiplot(MRIaggr.Pat1_red,param="DWI_t0",
             index1=list(coords="MASK_T2_FLAIR_t2",outline=TRUE),as.logical=TRUE,
             window="png",path="Display/DWI_lesion/")
multiplot(MRIaggr.Pat1_red,param="T2_FLAIR_t2",
             window="png",path="Display/T2/")

## Sweave generation
tablePat <- list(cbind(Id0=MRIaggr.Pat1_red@identifier,
                       selectClinic(MRIaggr.Pat1_red,param=c("Age","Gender"))),
                 cbind(Id0=MRIaggr.Pat1_red@identifier,
                       selectClinic(MRIaggr.Pat1_red,param=c("FinalStroke_volume","AcuteStroke_volume")))
)

res <- constSweave(dir="Display",
                   table=tablePat)

cat(res$ls.text[[1]],sep="")

#### Sweave doc ########################
# <<label=test,results=tex,echo=FALSE>>=
#  require(MRIaggr)
#
#### Sweave generation
# data("MRIaggr.Pat1_red", package="MRIaggr")
#
# 
# tablePat <- list(cbind(Id0=MRIaggr.Pat1_red@identifier,
#                     selectClinic(MRIaggr.Pat1_red,param=c("Age","Gender"))),
#     cbind(Id0=MRIaggr.Pat1_red@identifier,
#           selectClinic(MRIaggr.Pat1_red,param=c("FinalStroke_volume","AcuteStroke_volume")))
# )
#
# extra_text <- c("blablabla \\n \\n",
#                "\\\bigskip \\n \\n",
#                "blablabla \\n \\n",
#                "\\\[Y = \\\beta X + \\\varepsilon \\\] \\n \\n"
# )
# legend <- list("DWI alone","DWI with clinician segmentation","final contrast")
# 
# res <- constSweave(dir="Display",
#                    extra_text=list(extra_text,extra_text),
#                    table=tablePat)
#
#
#### preamble
# cat(res$text.preamble,sep="")
#
#### document 
# cat(res$text.begin,sep="") 
#
#
#
# for(iter_list in 1:length(res$ls.text)){
#   cat(res$ls.text[[iter_list]],sep="")
# }
#
#  cat(res$text.end,sep="")
# @
########################################

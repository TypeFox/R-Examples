require("NNTbiomarker")
# options(error=function(){
#   cat("Recovering!\n"); recover()
#   })

#stepsTableInitial = read.csv("SteppingStoneTable.csv", stringsAsFactors=F)[-1,]  ## TODO
#data(stepsTableInitial)
source("conveniences.R", local=F)
source("stepsTableInitial.R", local=F)
rownames(stepsTableInitial) = NULL
stepsTableInitial[[1]] = as.character(stepsTableInitial[[1]])

disableActionButton <- function(id,session) {
  catn("Disabling button " %&% id)
  jsCode = list(code=
                  paste("$('#", id,
                        "').prop('disabled',true);"
                        , sep=""))
  catn("jsCode ", jsCode[[1]])
  #$("#reportButton").prop('disabled',true)
  session$sendCustomMessage(type="jsCode", jsCode)
}
enableActionButton <- function(id,session) {
  catn("Enabling button " %&% id)
  session$sendCustomMessage(type="jsCode",
                            list(code= paste("$('#",id,"').prop('disabled',false)"
                                             ,sep="")))
}



completedToggle = function(number) {
  toggleLabelStyle =
    "color:#0099FF;vertical-align:top;text-align:left;font-size:150%;font-style:italic"
  div(
    hr(),
    fluidRow(
      column(2, ""),
      #HTML(str_dup("&nbsp;", 15)),
      column(8, div(style=toggleLabelStyle,
                    HTML("Is this step done?")
      ),
      HTML(stringr::str_dup("&nbsp;", 1)),
      radioButtons(  inline=TRUE,
                     ### inline requires in shinyUI a tags$head adding script for Shiny.addCustomMessageHandler("jsCode")
                     "stepStatus" %&% number,
                     label=NULL, choices=c("Not yet", "Done"))
      )
    )
  )
}

sectionHeader = function(number, content) {
  div(#class = "well container-fluid",
    style="border: 2px solid;
        padding: 10px;
        background: #dddddd;
        border-top-left-radius: 2em",
    list(hr(),
         h2(paste0('(', number, ') ',
                   stepsTableInitial[number, "Stepping stone"])),
         h3("Question " %&% number, ": ",
            stepsTableInitial[number, "Question"]),
         div(style="margin-bottom:0.5cm",
             content),
         completedToggle(number),
         hr()
    )
  )
}


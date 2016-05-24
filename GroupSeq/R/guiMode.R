"guiMode" <-
function()
{
  #requireNamespace(tcltk)
  backupScipen <- options(scipen=10)[[1]]

  taskWindow<-tktoplevel()
  tkwm.title(taskWindow,"Choose a Task")
  listBoxTasks<-tklistbox(taskWindow,height=4,width=40,selectmode="single",background="white")
  tkgrid(tklabel(taskWindow,text="Select a Task!"))
  tkgrid(listBoxTasks)
  tasks <- c("-1- Compute Bounds.","-2- Compute Drift given Power and Bounds.",
    "-3- Compute Probabilities given Bounds and Drift.","-4- Compute Confidence Interval.")
  for (i in (1:4))
  {
    tkinsert(listBoxTasks,"end",tasks[i])
  }
  tkselection.set(listBoxTasks,0)  # Default task is Task -1-.  Indexing starts at zero.

  OnOKtaskWindow <- function()
  {
    taskChoice <- as.numeric(tkcurselection(listBoxTasks))+1
    if(length(taskChoice)<1)
    {
     tkmessageBox(message="You must select a task!",icon="error",type="ok")
    }
    else
    {
     #call according function
     switch(taskChoice, guiInputTask1(taskWindow), guiInputTask2(taskWindow),
                       guiInputTask3(taskWindow), guiInputTask4(taskWindow) )
    }
  }

  quitGroupSeq <- function()
  {
   tkdestroy(taskWindow)

   # restore scipen value
   options(scipen=backupScipen)
   cat("GroupSeq closed by user.\n\n")
   return()
  }

  OK.button <-tkbutton(taskWindow,text="   Perform Selected Task   ",command=OnOKtaskWindow)
  Quit.buttton <- tkbutton(taskWindow,text="  QUIT GroupSeq  ",command=quitGroupSeq)

  # place buttons
  tkgrid(OK.button)
  tkgrid(tklabel(taskWindow,text="")) # Blank line
  tkgrid(Quit.buttton)
  tkgrid(tklabel(taskWindow,text="")) # Blank line

  tkfocus(taskWindow)


}

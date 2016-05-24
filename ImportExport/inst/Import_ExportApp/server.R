shinyServer(function(input, output,session) {

################ Format detection ################################# 
format<-reactive({ validate(need(!is.null(input$file),"File not selected"))
  if(input$correctFormat){format_detector(input$file$name)}
  else{ input$newformat}})
  
  
  
output$format<-renderText({format() })


output$listtables<- renderUI({
  if(format()=="MICROSOFT_OFFICE_ACCES_.mdb"){
  mycon<-RODBC::odbcConnectAccess(input$file$datapath,uid=input$uid_access_import,pwd=input$pwd_access_import)
  tables<- sqlTables(mycon)$TABLE_NAME 
  close(mycon)
  selectInput("tables","Tables in the file",choices=tables,multiple=F,selectize=F) 
  }
  else{NULL}
})


outputOptions(output, "format", suspendWhenHidden=FALSE)

#################### Reading data #########################################  

data_inic<-eventReactive(input$goButton,{
  
  
  #if for each extension
  
if(format()=="SPSS_.sav"){
  if(input$arguments_spss == ""){spss_import(input$file$datapath,use.value.labels=input$button_spss_import,out.format=input$date_spss_import)}
  else {eval(parse(text=paste("spss_import(",input$arguments_spss,")",sep="")))}
  }

else if(format()=="TXT" | format()=="CSV"){
  if(input$arguments_table == ""){table_import(input$file$datapath,sep=input$sep_txt_import,dec=input$dec_txt_import,header=input$button_txt_import)}
  else {eval(parse(text=paste("table_import(",input$arguments_table,")",sep="")))}
  }

  
else if(format()=="XLS" | format()=="XLSX"){
  if(input$arguments_excel == ""){read.xlsx(input$file$datapath,sheetIndex=input$sheetIndex_excel_import,sheetName=if(input$sheetName_excel_import != ""){input$sheetName_excel_import},header=input$button_excel_import )}
  else {eval(parse(text=paste("read.xlsx(",input$arguments_excel,")",sep="")))}
  } 
  
 

  
else if(format()=="MICROSOFT_OFFICE_ACCES_.mdb"){
  if(input$arguments_acces == ""){
    access_import(input$file$datapath,uid=input$uid_access_import,pwd=input$pwd_access_import,table_names=
                    ifelse(input$tableNames_access_import == "",input$tables,input$tableNames_access_import),
                  SQL_query=input$tableNames_access_import!="",out.format=input$date_access_import)}
  else {eval(parse(text=paste("access_import(",input$arguments_acces,")",sep="")))}
  }  
  
else if(format()=="STATA"){
  if(input$arguments_stata == ""){read_dta(input$file$datapath)}
  else {eval(parse(text=paste("read_dta(",input$arguments_stata,")",sep="")))}
  }   
  
else if(format()=="SAS_.sas7bdat"){if(input$arguments_SAS == ""){
  read_sas(input$file$datapath)}
  else {eval(parse(text=paste("read_sas(",input$arguments_SAS,")",sep="")))}
  }   
  
else if(format()=="R_.Rda"){get(load(as.character(input$file$datapath)))}
  
else if(format()=="SAS_.xport"){
  if(input$arguments_SASxport == ""){read.xport(input$file$datapath)}
  else {eval(parse(text=paste("read.xport(",input$arguments_SASxport,")",sep="")))}
  } 

    })

data<-reactive({
  if(input$format_corrector){format_corrector(data_inic())}
  else{data_inic()}
})

observeEvent(input$goButton,updateCollapse(session,id="collapse",open="Step 2) Database checkup") )



######## Database modification #####################################

output$format_corrector2<- renderText({ input$format_corrector})
output$variables <- renderText({ paste( input$vars,collapse=", ")})
output$unselected_variables <- renderText({ paste( setdiff(names(data()),input$vars),collapse=", ")})



output$var_view<-renderTable({ var_view(data()[input$vars]) })

output$summary <- renderPrint({
createTable(compareGroups(~.,data()[setdiff(names(data()[input$vars]),as.character(subset(var_view(data()[input$vars]),formats=="dates times" | formats=="date")$varname))],method=NA,alpha=as.numeric(input$summary_radiobuttons)))
  })
output$datesummary<- renderPrint({ if(length(as.character(subset(var_view(data()[input$vars]),formats=="dates times" | formats=="date")$varname) != 0)){
                                      summary(data()[as.character(subset(var_view(data()[input$vars]),formats=="dates times" | formats=="date")$varname)])
                                       } 
                                   else{ "No date variables in the data frame"}
  })


## datatable##
output$datatable<-renderDataTable(data()[input$vars] , options = list(pageLength = 10,scrollCollapse=T,scrollX=T) )

output$listvars<- renderUI({
  vars<-(if(length(names(data())) <100){ names(data())}
         else{names(data())[1:100]})
  
selectInput("vars","Included variables",vars,multiple=T,selected=if(input$allvar){names(data())},selectize=F) 
  
})




########### Data writting ######################################

file_name<-reactive({input$fileName})
outputFormat<-reactive({input$outputFormat})



### Access #
observeEvent(input$access_export_button,{
  if(input$arguments_access_export == ""){
    
  if(input$varview_access){access_export(file=file.choose(),list(data()[input$vars],var_view(data()[input$vars])), tablename=c(input$tableNames_access_export,"var_view")) }
  else{access_export(file=file.choose(),date_to_char(as.data.frame(data()[input$vars])),tablename=input$tableNames_access_export)}
  } 
  else{eval(parse(text=paste("access_export(",input$arguments_access_export,")",sep="")))}
})


#spss export##
observeEvent(input$SPSS_export_button,{
  dir<-choose.dir();
  setwd(dir);
  if(input$arguments_SPSS_export == ""){spss_export(as.data.frame(data()[input$vars]),file.save=input$final_spss_export,var.keep=input$var.keep_spss,file.runsyntax=input$filerunsyntax,file.data=input$data_spss_export,run.spss=input$run_spss) }
  else{eval(parse(text=paste("spss_export(",input$arguments_access_export,")",sep="")))}
})
  

output$file <- downloadHandler( 
  
  filename=function() { 
    if(outputFormat()=="XLSX"){ paste(file_name(), '.xlsx', sep='') }
    else if(outputFormat()=="XLS"){ paste(file_name(), '.xls', sep='') }
    else if(outputFormat()=="TXT"){ paste(file_name(), '.txt', sep='') }
    else if(outputFormat()=="CSV"){ paste(file_name(), '.csv', sep='') }
    else if(outputFormat()=="SPSS_.sav"){ paste(file_name(), '.sav', sep='') }
    else if(outputFormat()=="MICROSOFT_OFFICE_ACCES_.mdb"){ input$accessfile_export$datapath }
    else if(outputFormat()=="STATA"){ paste(file_name(), '.dta', sep='') }
    else if(outputFormat()=="R_.Rda"){ paste(file_name(), '.Rda', sep='') }
    else if(outputFormat()=="SAS_.xport"){ paste(file_name(), '.xport', sep='') }

  
  }
  , 

  content=function(file) {
    if(outputFormat()=="XLSX" | outputFormat()=="XLS"){
      if(input$arguments_excel_export == ""){ 
        if(input$varview_excel){excel_export(list(data()[input$vars],var_view(data()[input$vars])), file,table_names=c(input$tablename_excel,"var_view")) }
        else{excel_export(date_to_char(data()[input$vars]),file,table_names=input$tablename_excel)}
    } else{eval(parse(text=paste("excel_export(",input$arguments_excel_export,")",sep="")))}} 
    

    else if(outputFormat()=="TXT"){
      if(input$arguments_TXT_export == ""){write.table(data()[input$vars],file,sep=input$sep_txt_export,dec=input$dec_txt_export,na=input$NA_txt_export)}
      else {eval(parse(text=paste("write.table(",input$arguments_TXT_export,")",sep="")))}}
    
    else if(outputFormat()=="CSV"){
      if(input$arguments_csv_export == ""){write.csv(data()[input$vars],file,dec=input$dec_csv_export,na=input$NA_csv_export)}
      else {eval(parse(text=paste("write.csv(",input$arguments_csv_export,")",sep="")))}}
  
    else if(outputFormat()=="STATA"){write_dta(date_to_char(as.data.frame(data()[input$vars])),file)}
    
    else if(outputFormat()=="R_.Rda"){x<-as.data.frame(data()[input$vars]);save(x,file=file)}
    
    else if(outputFormat()=="SAS_.xport"){
      if(input$arguments_SAS_export == ""){write.foreign(data()[input$vars],file,codefile=input$codefile,package="SAS")}
      else{ eval(parse(text=paste("write.foreign(",input$arguments_SAS_export,")",sep="")))}}
    
  }
  
)






})






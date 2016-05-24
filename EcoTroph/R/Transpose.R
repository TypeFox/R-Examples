Transpose <-
function(tab_smooth,tab_input,column)
{
if(missing(tab_smooth))
cat("tab_smooth is missing\n")
if(missing(tab_input))
cat("tab_input is missing\n")
if(missing(column))
cat("column is missing\n")

tab_Trans <- array(dim=c(length(rownames(tab_smooth)),length(tab_input$group_name)),dimnames=list(rownames(tab_smooth),tab_input$group_name))  
pas <- round(0.00+(as.double(rownames(tab_smooth)[length(rownames(tab_smooth))])-as.double(rownames(tab_smooth)[length(rownames(tab_smooth))-1])),3)       ## recalculation of the 'pas', argument of the create.smooth function

for(groupe in tab_input$group_name)
{
tmp_tl=tab_input[tab_input$group_name==groupe,]$TL
tmp_tl <- seq(0,7,pas)[as.numeric(cut(tmp_tl+0.0000000000001,seq(0,7,pas),right=FALSE))]  ## assignment of a trophic class to the trophic groups
tab_Trans[,groupe]=tab_input[tab_input$group_name==groupe,column]*tab_smooth[,as.character(tmp_tl)]
}

class(tab_Trans)<-"Transpose"
return (tab_Trans)
}


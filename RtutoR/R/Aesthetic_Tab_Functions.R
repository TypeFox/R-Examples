mapAes2Var <- function(plot_type,var_list,layer_count) {

  Geom <- Aes <- Type <- NULL

  applicable_aes <- aesthetics %>% filter(Geom == plot_type & Type == "Aesthetic") %>% select(Aes)
  id_val = paste("map_aes_controls",layer_count,sep="")
  div(id = id_val,
      h4(strong(paste("LAYER ",layer_count, " - MAP AESTHETIC TO VARIABLE")),style="color:darkblue"),
      selectInput(inputId = paste("select_aes",layer_count,sep=""),"Select an Aesthetic",choices = c("",as.character(unlist(applicable_aes)))),
      selectInput(inputId = paste("aes_mapping_var",layer_count,sep=""),"Select Mapping Variable",choices = c("",var_list))
  )

}

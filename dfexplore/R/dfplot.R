# generic dfplot method 
setGeneric(name="dfplot",
           def=function(dfdescription,title = NULL){standardGeneric("dfplot")})

# dfplot for data.frame class. Create a data.frame.description object and then plot it
setMethod(f="dfplot",
          signature="data.frame",
          def = function(dfdescription, title = NULL){
            
            if(is.null(title)) title = as.character(as.list(match.call())$dfdescription)
            # Create a Data.frame.description object
            object_df_description <- data.frame.description(dfdescription,title )
            
            # Plot it
            return(dfplot(object_df_description))
          }
)

# dfplot method for Data.frame.description class. Main dfplot method.
setMethod(f="dfplot",
          signature="Data.frame.description",
          def = function(dfdescription,title=dfdescription@title){
            
            ### Choose good colors ###
            modalities <- c("NA", "integer",   "numeric",   "Date",         "character", "factor", "ordered",  "other")
            couleurs <- c("red","palegreen1","palegreen3","rosybrown1","paleturquoise1","paleturquoise3","paleturquoise4","snow2")
            
            used_classes <- modalities %in% unique(dfdescription@df_pour_plot$cell_content)
            
            ### Create a ggplot object ###
            # init vectors used in ggplot to make R CMD check happy and avoid
            # the note " no visible binding for global variable for ‘variable’
            variable <- NULL
            observation <- NULL
            cell_content <- NULL
            
            gg_tableplot <- 
              ggplot(data=dfdescription@df_pour_plot, aes( x=variable, y=observation) ) +
              geom_tile(aes(fill=cell_content)) +
              theme(axis.text.x = element_text(angle = 90)) +
              theme(panel.background = element_blank()) +
              scale_fill_manual(name="Class of data", values=couleurs[used_classes]) +
              scale_y_reverse() + # draw obs from first up to last down
              ggtitle(paste("Classes of data in '",dfdescription@title,"'",sep="")) +
              coord_cartesian(ylim=c(nrow(dfdescription@df_source),1))
            
            # Return a ggplot object. Directly draw if not assigned
            return(gg_tableplot)
          }
)


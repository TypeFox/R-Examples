readhelpers.formatmeta <-
function(metaD_unformat, metaD_n_cols, brLens, sTreeFlg) {
  # Descr:  formatting metadata
  # Deps:   -
  # I/p:    metaD_unformat
  #         metaD_n_cols
  #         brLens
  #         sTreeFlg

  debugBool = get("P2C2M_flg_dbgBool", envir=P2C2M_globalVars)
  if (debugBool) {
    cat("\n", xtermStyle::style("DEBUG> readhelpers.formatmeta", fg="red"), 
        sep="")
  }
  
  # DEBUGLINES:
  #cat("\nmetaD_unformat\n"); print(metaD_unformat)
  #cat("\nmetaD_n_cols\n"); print(metaD_n_cols)

################################
# 1. Setting up empty matrices #
################################
  if(metaD_n_cols==1) {
    metaD_formatted = array(dim=c(length(metaD_unformat), metaD_n_cols+2))
  }
  if(sTreeFlg && metaD_n_cols==3) {
    metaD_formatted = array(dim=c(length(metaD_unformat), metaD_n_cols+3))
  }

###################################################################
# 2. Save node name, branch length and demographic info to matrix #
###################################################################
  for(i in 1:length(metaD_unformat)) {
    brInfo = unlist(strsplit(brLens[i], ":"))
    metaInfo = unlist(strsplit(metaD_unformat[i], ","))

    if(metaD_n_cols==1) {
      metaD_formatted[i,] = c(brInfo, metaInfo)
    }
    if(sTreeFlg && metaD_n_cols==3) {
      dmvValue = mean(as.numeric(c(metaInfo[2], metaInfo[3])))          # Calculation of dmvValue as mean between dmv_start and dmv_end
      metaD_formatted[i,] = c(brInfo, metaInfo, dmvValue)               # Fill the first two columns of the matrix "metaD_formatted" with node name and branch length info,and fill the remaining three columns with the metadata info
    }
  }

  # DEBUGLINES:
  #cat("\nmetaD_formatted\n"); print(metaD_formatted)

##########################
# 3. Formatting metadata #
##########################
  metaD_formatted = metaD_formatted[which(metaD_formatted[,1]!=""),]    # Remove all those metaD_formatted entries that do not contain a branch name (i.e. that represent nodes)

  if(sTreeFlg && metaD_n_cols==1) {
      colnames(metaD_formatted) = c("br", "length", "dmv")              # Name the columns of matrix "metaD_formatted" depending on the number of metadata infos
  }
  if(!sTreeFlg && metaD_n_cols==1) {
      colnames(metaD_formatted) = c("br", "length", "rate")
  }
  if(sTreeFlg && metaD_n_cols==3) {
      colnames(metaD_formatted) = c("br", "length", "dmt",
                                    "dmv_start", "dmv_end", "dmv")
  }

  rownames(metaD_formatted) = metaD_formatted[,1]                       # name the rows of of matrix "metaD_formatted" by the node name

  return(metaD_formatted)
}

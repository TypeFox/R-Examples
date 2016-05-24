
#' @title VMS File Class
#'
#' @description  VMS File reference Class
#'\itemize{
#'\item{\code{path} character path of the vms file}
#'\item{\code{data} data.frame data from the vms file}
#'}
#'
#' @details \code{vms_File} class is a reference class for the VMS file.
#'
#' @docType class
#' @aliases vms_File-class
#'
#' @name vms_File

vms_File <- setRefClass("vms_File", 
                        fields = list(path = "character",
                                      data = "data.frame"))


#' @title VMS DataBase class 
#'
#' @description VMS DataBase reference Class
#'
#'\itemize{
#'\item{\code{dir} character path of the VMS DataBase}
#'\item{\code{db} character name of the VMS DataBase}
#'\item{\code{tab} character tables in the VMS DataBase}
#'
#'}
#' @details \code{vms_DB} class is a reference class for the VMS DataBase.
#'
#' @docType class
#' @aliases vms_DB-class
#' 
#' @name vms_DB

vms_DB <- setRefClass("vms_DB", 
                      fields = list(dir = "character",
                                    db = "character",
                                    tab = "character"))


#' @title VMS DataBase Query class 
#'
#' @description VMS DataBase Query reference Class
#'
#'\itemize{
#'\item{\code{dir} character path of the VMS DataBase}
#'\item{\code{db} character name of the VMS DataBase}
#'\item{\code{que} character query for the VMS DataBase}
#'
#'}
#' @details \code{que_vms_DB} class is a reference class for the VMS DataBase.
#'
#' @docType class
#' @aliases que_vms_DB-class
#' 
#' @name que_vms_DB

que_vms_DB <- setRefClass("que_vms_DB", 
                          fields = list(dir = "character",
                                        db = "character",
                                        que = "character"))


#' @title LogBook File Class
#'
#' @description LogBook File reference Class
#' 
#'\itemize{
#'\item{\code{path} character path of the logbook file}
#'\item{\code{data} data.frame data from the logbook file}
#'
#'}
#' @details \code{log_File} class is a reference class for the LogBook file.
#'
#' @docType class
#' @aliases log_File-class
#'
#'
#' @name log_File

log_File <- setRefClass("log_File", 
                        fields = list(path = "character", 
                                      data = "data.frame"))


#' @title LogBook DataBase class 
#'
#' @description LogBook DataBase reference Class
#' 
#'\itemize{
#'\item{\code{dir} character path of the LogBook DataBase}
#'\item{\code{db} character name of the LogBook DataBase}
#'\item{\code{tab} character tables in the LogBook DataBase}
#'
#'}
#' @details \code{log_DB} class is a reference class for the Logbook DataBase.
#'
#' @docType class
#' @aliases log_DB-class
#'
#' 
#' @name log_DB
#' 

log_DB <- setRefClass("log_DB", 
                      fields = list(dir = "character", 
                                    db = "character",
                                    tab = "character"))


#' @title LogBook Clustering class 
#'
#' @description LogBook Clustering rData file reference Class
#' 
#'\itemize{
#'\item{\code{data} character The Clustering result data}
#'\item{\code{options} character Metier manual annotation}
#'
#'}
#' @details \code{log_Cla} class is a reference class for the Logbook Clustering rData file.
#'
#' @docType class
#' @aliases log_Cla-class
#'
#' 
#' @name log_Cla
#' 

log_Cla <- setRefClass("log_Cla", fields = c("data", "options"))


#' @title Harbours Coordinates Shape File class 
#'
#' @description Harbours coordinates Shape File reference Class
#' 
#'\itemize{
#'\item{\code{path} character The Path to the harbours coordinates shape file}
#'\item{\code{data} SpatialPointsDataFrame The loaded harbours coordinates data}
#'
#'}
#' @details \code{harbCoo} class is a reference class for the Harbours Coordinates shape file.
#'
#' @docType class
#' @aliases harbCoo-class
#' 
#' @name harbCoo
#' 

harbCoo <- setRefClass("harbCoo",
                       fields = list(path = "character",
                                     data = "SpatialPointsDataFrame"))


#' @title Land Map Shape File class 
#'
#' @description Land Map Shape File reference Class
#' 
#'\itemize{
#'\item{\code{path} character The Path to the Land Map shape file}
#'\item{\code{data} SpatialPolygonsDataFrame The loaded Land Map data}
#'
#'}
#' @details \code{polymap} class is a reference class for the Land Map shape file.
#'
#' @docType class
#' @aliases polymap-class
#'
#' @param path character The Path to the Land Map shape file
#' @param data SpatialPolygonsDataFrame The loaded Land Map data
#' 
#' @name polymap
#' 

polymap <- setRefClass("polymap",
                       fields = list(path = "character",
                                     data = "SpatialPolygonsDataFrame"))


#' @title Bathymetry rData File class 
#'
#' @description Bathymetry rData File reference Class
#' 
#'\itemize{
#'\item{\code{path} character The Path to the Bathymetry rData file}
#'\item{\code{data} XYZ-matrix The loaded Bathymetry data}
#'
#'}
#' @details \code{bathymetry} class is a reference class for the Bathymetry rData file.
#'
#' @docType class
#' @aliases bathymetry-class
#'
#' 
#' @name bathymetry
#' 

bathymetry <- setRefClass("bathymetry", fields = c("path", "data"))

miso_list <- R6Class("miso_list",
                     portable = FALSE,
                     class = TRUE,
                     public = list(
                       th_lst = NULL,
                       
                       add_sou = function(new_dataset){
                         
                         self$th_lst <- c(th_lst, new_dataset)
                         
                         upd_sou_num()
                       },
                       
                       rem_sou = function(){
                         
                         self$th_lst[length(th_lst)] <- NULL
                         
                         upd_sou_num()
                         
                       },
                       
                       upd_sou_num = function(){
                         len_lis <- length(self$th_lst)
                         if(len_lis != 0){
                           for(j in 1:len_lis){
                             th_lst[[j]]$set_objNum(j)
                           }
                         }
                       },
                       
                       upd_sou_idt = function(num, value){
                         th_lst[[num]]$set_ids_typ(value)
                       }
                       
                     ))


mix_sou <- R6Class("mixsou",
                   portable = FALSE,
                   class = TRUE,
                   public = list(
                     dbPath = NULL,
                     objNum = NULL,
                     
                     sourceName = NULL,
                     vessIds = NULL,
                     ids_typ = NULL,
                     
                     vessDate = NULL,
                     n_ping = NULL,
                     n_vess = NULL,
                     n_day = NULL,
                     n_pida = NULL,
                     flee_ids = NULL,
                     ids_tab = NULL,
                     res_mat = NULL,
                     widFrame = NULL,
                     #                      widGroup = NULL,
                     widEdit = NULL,
                     
                     set_ids_typ = function(value){
                       self$ids_typ <- value
                     },
                     set_objNum = function(value){
                       self$objNum <- value
                     },
                     
                     initialize = function(db_path, main_cont = "", sel_wid, cur_sou){
                       obj_lst <- cur_sou$th_lst
                       
                       self$dbPath <- db_path
                       self$objNum <- length(obj_lst)+1
                       self$vessIds <- as.character(read.csv.sql(dbPath, sql = "select distinct I_NCEE from file", sep = ";", eol = "\n")[,1])
                       self$vessDate <- read.csv.sql(dbPath, sql = "select DATE from file", sep = ";", eol = "\n")[,1]
                       self$n_ping <- length(vessDate)
                       self$n_vess <- length(vessIds)
                       self$n_day <- length(unique(round(vessDate)))
                       self$n_pida <- as.numeric(table(floor(vessDate)))
                       
                       if(class(main_cont) == "gGroup"){
                         self$widFrame <- gframe(container = main_cont, horizontal = TRUE)
                         
                         addSpring(widFrame)
                         self$widEdit <- gedit(text = paste("Source ", length(obj_lst)+1, sep = ""), width = 10, container = widFrame, handler = function(h,...){
                           self$sourceName <- svalue(widEdit)
                           plot_comp(svalue(sel_wid), cur_sou$th_lst)
                           
                         })
                         
                         self$sourceName <- svalue(widEdit)
                         addSpring(widFrame)
                         gimage(system.file("ico/view-refresh-5.ico", package="vmsbase"), container = widFrame, handler = function(h,...){
                           self$sourceName <- svalue(widEdit)
                           plot_comp(svalue(sel_wid), cur_sou$th_lst)
                           
                         })
                         
                         gimage(system.file("ico/list-remove-4.ico", package="vmsbase"), container = widFrame, handler = function(h,...){
                           if(gconfirm("Remove this dataset from the list?", title="Confirm", icon = c("warning"), parent=NULL)){
                             delete(main_cont, widFrame)
                             
                             cur_sou$rem_sou()
                             
                             if(length(cur_sou$th_lst) != 0){
                               plot_comp(svalue(sel_wid), cur_sou$th_lst)
                             }else{
                               plot.new()
                             }
                             
                           }
                         })
                         
                         addSpring(widFrame)
                         
                         if(objNum == 1){
                           self$flee_ids <- vessIds
                           
                         }else{
                           
                           items <- cbind(obj_lst[[objNum-1]]$flee_ids[1:5], vessIds[1:5])
                           #                            colnames(items) <- c(obj_lst[[objNum-1]]$sourceName, sourceName)
                           temp_dia <- gbasicdialog(title="Please, select an option", do.buttons = FALSE)
                           size(temp_dia) <- c(550, 350)
                           up_g <- ggroup(horizontal = FALSE, container = temp_dia)
                           
                           glabel("\nThis is a sample of the vessel IDs found in your dataset:", container = up_g)
                           addSpace(up_g, 15, horizontal = FALSE)
                           tab_g <- ggroup(horizontal = TRUE, container = up_g)
                           size(tab_g) <- c(450, 150)
                           
                           addSpring(tab_g)
                           id_ty_1 <- gradio(c("CFR", "MMSI", "IMO", "Other"), selected = 1, container = tab_g)
                           svalue(id_ty_1) <- "CFR"
                           t_tab_1 <- gtable(items[,1], container = tab_g, expand = TRUE, fill = TRUE)
                           
                           names(t_tab_1) <- obj_lst[[objNum-1]]$sourceName
                           
                           addSpring(tab_g)
                           t_tab_2 <- gtable(items[,2], container = tab_g, expand = TRUE, fill = TRUE)
                           id_ty_2 <- gradio(c("CFR", "MMSI", "IMO", "Other"), selected = 1, container = tab_g)
                           
                           names(t_tab_2) <- sourceName
                           
                           svalue(id_ty_2) <- "CFR"
                           addSpring(tab_g)
                           
                           addSpace(up_g, 15, horizontal = FALSE)
                           dia_g <- ggroup(horizontal = TRUE, container = up_g)
                           addSpring(dia_g)
                           bu1_g <- ggroup(horizontal = FALSE, container = dia_g)
                           gimage(system.file("ico/document-save-2.ico", package="vmsbase"), container = bu1_g)
                           gbutton("\t    Merge\nvessels ID 'as it is'", container = bu1_g, handler = function(h,...){
                             
                             self$flee_ids <- unique(c(cur_sou$th_lst[[objNum-1]]$flee_ids, vessIds))
                             cur_sou$upd_sou_idt(objNum-1, svalue(id_ty_1))
                             set_ids_typ(svalue(id_ty_2))
                             
                             dispose(temp_dia)
                           })
                           addSpring(dia_g)
                           bu2_g <- ggroup(horizontal = FALSE, container = dia_g)
                           gimage(system.file("ico/address-book-new.ico", package="vmsbase"), container = bu2_g)
                           gbutton("\t      Load\nAssociation Table", container = bu2_g, handler = function(h,...){
                             
                             tmp_tab <- gsub('\\\\', '/', gfile(text = "Select Pairs",type = "open",
                                                                filter = list("ID Pairs data" = list(patterns = c("*")),
                                                                              "All files" = list(patterns = c("*")))))
                             self$ids_tab <- read.csv(tmp_tab, header = FALSE)
                             
                             tab_matc <- ids_tab[ids_tab[,1] %in% vessIds,1]
                             tab_lst <- as.list(ids_tab[ids_tab[,1] %in% vessIds,2])
                             names(tab_lst) <- tab_matc
                             
                             for(i in names(tab_lst)){
                               self$vessIds[which(vessIds == i)] <-  as.character(tab_lst[i])
                             }
                             
                             self$flee_ids <- unique(c(cur_sou$th_lst[[objNum-1]]$flee_ids, vessIds))
                             
                             cur_sou$upd_sou_idt(objNum-1, svalue(id_ty_1))
                             set_ids_typ(svalue(id_ty_2))
                             
                             dispose(temp_dia)
                           })
                           addSpring(dia_g)
                           bu3_g <- ggroup(horizontal = FALSE, container = dia_g)
                           gimage(system.file("ico/insert-link-2.ico", package="vmsbase"), container = bu3_g)
                           gbutton("\t  Run\nPairing Code", container = bu3_g, handler = function(h,...){
                             
                             cur_sou$upd_sou_idt(objNum-1, svalue(id_ty_1))
                             set_ids_typ(svalue(id_ty_2))
                             
                             merge_dia <- gbasicdialog(title="Set the parameters", do.buttons = FALSE, parent = temp_dia)
                             glabel("\nWarning!\n\n", container = merge_dia)
                             
                             bi_g <- ggroup(horizontal = FALSE, container = merge_dia)
                             
                             bi2_g <- gframe("Parameters", horizontal = TRUE, container = bi_g)
                             addSpring(bi2_g)
                             par1_g <- ggroup(horizontal = FALSE, container = bi2_g)
                             par1_l <- glabel("Rounding Level", container = par1_g)
                             par1_s <- gslider(from = 2, to = 4, by = 1, container = par1_g)
                             addSpring(bi2_g)
                             par2_g <- ggroup(horizontal = FALSE, container = bi2_g)
                             par2_l <- glabel("Overlap Threshold", container = par2_g)
                             par2_s <- gslider(from = 1, to = 30, by = 1, value = 10, container = par2_g)
                             addSpring(bi2_g)
                             
                             start_b <- gbutton("   Start\nMerging", container = bi2_g, handler = function(h,...){
                               
                               data_2 <- read.csv.sql(obj_lst[[1]]$dbPath, "select I_NCEE, LON, LAT, DATE from file", sep = ";")
                               data_1 <- read.csv.sql(dbPath, "select I_NCEE, LON, LAT, DATE from file", sep = ";")
                               if(nrow(data_1) != 0 & nrow(data_2) != 0){
                                 
                                 self$res_mat <- merge_source(data_source1 = data_1, data_source2 = data_2, rnd_level = svalue(par1_s), minover = svalue(par2_s))
                                 Obs_med <- as.numeric(self$res_mat)
                                 Obs_med <- Obs_med[which(Obs_med < 20)]
                                 hist(Obs_med, 100)
                                 abline(v = 10, col = 2, lwd = 2, lty = 2)
                                 
                               }
                               enabled(bigra_g) <- TRUE
                             })
                             
                             bigra_g <- ggroup(horizontal = FALSE, container = bi_g)
                             
                             par3_g <- ggroup(horizontal = FALSE, container = bigra_g)
                             par3_l <- glabel("Distance Threshold", container = par3_g)
                             par3_s <- gslider(from = 0.1, to = 20, by = 0.1, value = 5, container = par3_g, handler = function(h,...){
                               
                               Obs_med <- as.numeric(res_mat)
                               Obs_med <- Obs_med[which(Obs_med < 20)]
                               hist(Obs_med, 31)
                               abline(v = svalue(par3_s), col = 2, lwd = 2, lty = 2)
                               pos_pairs_D <- which(res_mat < svalue(par3_s), arr.ind=TRUE)
                               cat("\n", nrow(unique(cbind(rownames(res_mat)[pos_pairs_D[,1]],colnames(res_mat)[pos_pairs_D[,2]]))), " pairs with a ", svalue(par3_s), "-kilometer threshold", sep = "")
                               
                             })
                             
                             ggraphics(width = 500, height = 300, container = bigra_g)
                             
                             save_b <- gbutton("   Save\nMerging", container = bigra_g, handler = function(h,...){
                               
                               pos_pairs_D <- which(res_mat < svalue(par3_s), arr.ind=TRUE)
                               self$ids_tab <- unique(cbind(rownames(res_mat)[pos_pairs_D[,1]],colnames(res_mat)[pos_pairs_D[,2]]))
                               
                               tab_matc <- ids_tab[ids_tab[,1] %in% vessIds,1]
                               tab_lst <- as.list(ids_tab[ids_tab[,1] %in% vessIds,2])
                               names(tab_lst) <- tab_matc
                               for(i in names(tab_lst)){
                                 self$vessIds[which(vessIds == i)] <-  as.character(tab_lst[i])
                               }
                               self$flee_ids <- unique(c(cur_sou$th_lst[[objNum-1]]$flee_ids, vessIds))
                               
                               dispose(merge_dia)
                             })
                             addSpring(bi2_g)
                             enabled(bigra_g) <- FALSE
                             visible(merge_dia, set=TRUE)
                             dispose(temp_dia)
                             
                           })
                           addSpring(dia_g)
                           visible(temp_dia, set=TRUE)
                         }
                       }
                     })
)
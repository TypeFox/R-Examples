


setMethod ('show' , 'EotMode', 
           function(object) {
             cat('class                :', class(object), '\n')
             cat('name                 :', names(object), '\n')
             cat('base point (x, y)    :', object@coords_bp, '\n')
             cat('cum. expl. variance  :', object@cum_exp_var, '\n')
             cat('dimensions           : ', 
                 raster::nrow(object@r_predictor), ', ', 
                 raster::ncol(object@r_predictor), ', ', 
                 raster::ncell(object@r_predictor),
                 '  (nrow, ncol, ncell)\n', sep="") 
             cat('resolution           : ', 
                 raster::xres(object@r_predictor), ', ', 
                 raster::yres(object@r_predictor), '  (x, y)\n', sep="")
             cat('extent               : ', object@r_predictor@extent@xmin, 
                 ', ', object@r_predictor@extent@xmax, ', ', 
                 object@r_predictor@extent@ymin, ', ', 
                 object@r_predictor@extent@ymax, 
                 '  (xmin, xmax, ymin, ymax)\n', sep="")
             cat('coord. ref.          :', 
                 raster::projection(object@r_predictor, TRUE), '\n')
           }
)	

setMethod ('show' , 'EotStack', 
           function(object) {
             cat('class                :', class(object), '\n')
             cat('cum. expl. variance  :', 
                 object[[nmodes(object)]]@cum_exp_var, '\n')
             cat('names                :', names(object), '\n')
             cat('dimensions           : ', 
                 raster::nrow(object[[1]]@r_predictor), ', ', 
                 raster::ncol(object[[1]]@r_predictor), ', ', 
                 raster::ncell(object[[1]]@r_predictor), ', ', 
                 nmodes(object), '  (nrow, ncol, ncell, nmodes)\n', sep="") 
             cat('resolution           : ', 
                 raster::xres(object[[1]]@r_predictor), ', ', 
                 raster::yres(object[[1]]@r_predictor), 
                 '  (x, y)\n', sep="")
             cat('extent               : ', 
                 object[[1]]@r_predictor@extent@xmin, ', ', 
                 object[[1]]@r_predictor@extent@xmax, ', ', 
                 object[[1]]@r_predictor@extent@ymin, ', ', 
                 object[[1]]@r_predictor@extent@ymax, 
                 '  (xmin, xmax, ymin, ymax)\n', sep="")
             cat('coord. ref.          :', 
                 raster::projection(object[[1]]@r_predictor, TRUE), '\n')
           }
)
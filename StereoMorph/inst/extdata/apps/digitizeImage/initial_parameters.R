init_params <- list()
init_params[['app_dir']] <- c('/Applications/XAMPP/xamppfiles/htdocs/data_analysis/r_package_development/StereoMorph/inst/extdata/apps/digitizeImage')

init_params[['prev_wd']] <- c('/Users/aaron/Documents/Research/Software Manuals/StereoMorph/Digitizing App Tutorial/Example project')

init_params[['img_name']] <- c('mug_001.jpg')

init_params[['image_id']] <- c('mug_001')

init_params[['img_size']] <- c('110136')

init_params[['img_file']] <- c('Mugs/mug_001.jpg')

init_params[['scaling_units']] <- c('cm')

init_params[['scaling']] <- c('NA')

init_params[['ruler_interval']] <- c('1')

init_params[['checkerboard_nx']] <- c('NA')

init_params[['checkerboard_ny']] <- c('NA')

init_params[['shapes_file']] <- c('Shapes/mug_001.txt')

init_params[['img_width']] <- c('2144')

init_params[['img_height']] <- c('1424')

init_params[['landmark_color_blur']] <- c('blue')

init_params[['landmark_color_focus']] <- c('green')

init_params[['curve_color_blur']] <- c('purple')

init_params[['control_point_color_blur']] <- c('purple')

init_params[['control_point_color_focus']] <- c('red')

init_params[['landmark_radius']] <- c('4')

init_params[['control_point_radius']] <- c('4')

init_params[['marker_stroke_width']] <- c('2')

init_params[['landmarks_ref']] <- c('Handle_Cup_Out_Top','Handle_Cup_In_Top','Handle_Cup_Out_Bottom','Handle_Cup_In_Bottom','Cup_Rim_Handle','Cup_Rim_Opposite','Cup_Base_Handle','Cup_Base_Opposite')

init_params[['curves_ref']][[1]] <- c('Handle_Cup_Out','Handle_Cup_Out_Top','Handle_Cup_Out_Bottom')
init_params[['curves_ref']][[2]] <- c('Handle_Cup_In','Handle_Cup_In_Top','Handle_Cup_In_Bottom')
init_params[['curves_ref']][[3]] <- c('Side_Handle_In','Handle_Cup_In_Top','Handle_Cup_In_Bottom')

init_params[['unsaved_landmarks']] <- c('FALSE')

init_params[['unsaved_curves']] <- c('FALSE')

init_params[['prev_img']] <- c('FALSE')

init_params[['next_img']] <- c('TRUE')

init_params[['landmarks']][[1]] <- c('Handle_Cup_Out_Top','1041','324')
init_params[['landmarks']][[2]] <- c('Handle_Cup_In_Top','1033','531')
init_params[['landmarks']][[3]] <- c('Handle_Cup_Out_Bottom','1027','1088')
init_params[['landmarks']][[4]] <- c('Handle_Cup_In_Bottom','1027','884')
init_params[['landmarks']][[5]] <- c('Cup_Rim_Handle','1053','182')
init_params[['landmarks']][[6]] <- c('Cup_Rim_Opposite','1970','196')
init_params[['landmarks']][[7]] <- c('Cup_Base_Handle','1027','1209')
init_params[['landmarks']][[8]] <- c('Cup_Base_Opposite','1943','1248')

init_params[['control_points']][[1]] <- c('Handle_Cup_Out','1041','324','1039','335','995','321','870','296','784','342','649','421','611','629','583','885','700','1017','817','1139','1027','1088')
init_params[['control_points']][[2]] <- c('Handle_Cup_In','1033','531','1028','523','1018','522','1003','516','996','482','960','396','893','407','762','418','710','595','683','696','703','804','726','947','821','991','930','1040','977','943','993','904','1008','897','1022','897','1027','884')
init_params[['control_points']][[3]] <- c('Side_Handle_In','1033','531','1030','689','1027','884')

init_params[['ruler_pixel']] <- c('118.404222755233')

init_params[['ruler_points']][[1]] <- c('Ruler point 1','260','1138')
init_params[['ruler_points']][[2]] <- c('Ruler point 2','259','1022')
init_params[['ruler_points']][[3]] <- c('Ruler point 3','260','906')
init_params[['ruler_points']][[4]] <- c('Ruler point 4','259','789')
init_params[['ruler_points']][[5]] <- c('Ruler point 5','257','671')
init_params[['ruler_points']][[6]] <- c('Ruler point 6','256','553')
init_params[['ruler_points']][[7]] <- c('Ruler point 7','254','433')
init_params[['ruler_points']][[8]] <- c('Ruler point 8','252','314')
init_params[['ruler_points']][[9]] <- c('Ruler point 9','252','193')
init_params[['ruler_points']][[10]] <- c('Ruler point 10','251','71')


json_string <- '{"app_dir":["/Applications/XAMPP/xamppfiles/htdocs/data_analysis/r_package_development/StereoMorph/inst/extdata/apps/digitizeImage"], "prev_wd":["/Users/aaron/Documents/Research/Software Manuals/StereoMorph/Digitizing App Tutorial/Example project"], "img_name":["mug_001.jpg"], "image_id":["mug_001"], "img_size":["110136"], "img_file":["Mugs/mug_001.jpg"], "scaling_units":["cm"], "scaling":["NA"], "ruler_interval":["1"], "checkerboard_nx":["NA"], "checkerboard_ny":["NA"], "shapes_file":["Shapes/mug_001.txt"], "img_width":["2144"], "img_height":["1424"], "landmark_color_blur":["blue"], "landmark_color_focus":["green"], "curve_color_blur":["purple"], "control_point_color_blur":["purple"], "control_point_color_focus":["red"], "landmark_radius":["4"], "control_point_radius":["4"], "marker_stroke_width":["2"], "landmarks_ref":["Handle_Cup_Out_Top", "Handle_Cup_In_Top", "Handle_Cup_Out_Bottom", "Handle_Cup_In_Bottom", "Cup_Rim_Handle", "Cup_Rim_Opposite", "Cup_Base_Handle", "Cup_Base_Opposite"], "curves_ref":[["Handle_Cup_Out", "Handle_Cup_Out_Top", "Handle_Cup_Out_Bottom"], ["Handle_Cup_In", "Handle_Cup_In_Top", "Handle_Cup_In_Bottom"], ["Side_Handle_In", "Handle_Cup_In_Top", "Handle_Cup_In_Bottom"]],"unsaved_landmarks":["FALSE"], "unsaved_curves":["FALSE"], "prev_img":["FALSE"], "next_img":["TRUE"], "landmarks":[["Handle_Cup_Out_Top", "1041", "324"], ["Handle_Cup_In_Top", "1033", "531"], ["Handle_Cup_Out_Bottom", "1027", "1088"], ["Handle_Cup_In_Bottom", "1027", "884"], ["Cup_Rim_Handle", "1053", "182"], ["Cup_Rim_Opposite", "1970", "196"], ["Cup_Base_Handle", "1027", "1209"], ["Cup_Base_Opposite", "1943", "1248"]],"control_points":[["Handle_Cup_Out", "1041", "324", "1039", "335", "995", "321", "870", "296", "784", "342", "649", "421", "611", "629", "583", "885", "700", "1017", "817", "1139", "1027", "1088"], ["Handle_Cup_In", "1033", "531", "1028", "523", "1018", "522", "1003", "516", "996", "482", "960", "396", "893", "407", "762", "418", "710", "595", "683", "696", "703", "804", "726", "947", "821", "991", "930", "1040", "977", "943", "993", "904", "1008", "897", "1022", "897", "1027", "884"], ["Side_Handle_In", "1033", "531", "1030", "689", "1027", "884"]],"ruler_pixel":["118.404222755233"], "ruler_points":[["Ruler point 1", "260", "1138"], ["Ruler point 2", "259", "1022"], ["Ruler point 3", "260", "906"], ["Ruler point 4", "259", "789"], ["Ruler point 5", "257", "671"], ["Ruler point 6", "256", "553"], ["Ruler point 7", "254", "433"], ["Ruler point 8", "252", "314"], ["Ruler point 9", "252", "193"], ["Ruler point 10", "251", "71"]]}'

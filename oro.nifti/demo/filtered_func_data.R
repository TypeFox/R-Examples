ffd <- readNIfTI(file.path(system.file("nifti", package="oro.nifti"),
                           "filtered_func_data"))
image(ffd)
orthographic(ffd)

print.clDeviceID <- function(x, ...) {
  i <- .Call("ocl_get_device_info", x)
  cat(" OpenCL device '", i$name, "'\n", sep='')
  x
}

print.clPlatformID <- function(x, ...) {
  i <- .Call("ocl_get_platform_info", x)
  cat(" OpenCL platform '", i$name, "'\n", sep='')
  x
}

print.clContext <- function(x, ...) {
  cat(" OpenCL context ")
  print.default(x, ...)
  x
}

names.clKernel <- function(x) names(attributes(x))
`$.clKernel` <- function(x, name) attr(x, name)
`$<-.clKernel` <- function(x, name, value) stop("clKernel properies are read-only")
print.clKernel <- function(x, ...) {
  cat(" OpenCL kernel '", attr(x, "name"),"'\n", sep='')
  a <- attributes(x)
  a$class <- NULL
  a$name <- NULL
  print(a)
  x
}
  
oclPlatforms <- function() .Call("ocl_platforms")
oclDevices <- function(platform = oclPlatforms()[[1]], type="default") .Call("ocl_devices", platform, type)
oclSimpleKernel <- function(device, name, code, precision = c("single", "double", "best")) {
  precision <- match.arg(precision)
  if (precision == "best") { # detect supported precision from the device
    precision <- if (any(grepl("cl_khr_fp64", oclInfo(device)$exts))) {
      code <- c("#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n", gsub("\\bfloat\\b", "double", code))
      "double"
    } else "single"
  }
  .Call("ocl_ez_kernel", device, name, code, precision)
}
oclRun <- function(kernel, size, ..., native.result=FALSE, wait=TRUE, dim=size) .External("ocl_call", kernel, size, native.result, wait, dim, ...)
oclResult <- function(context, wait = TRUE) .Call("ocl_collect_call", context, wait)

oclInfo <- function(item) UseMethod("oclInfo")
oclInfo.clDeviceID <- function(item) .Call("ocl_get_device_info", item)
oclInfo.clPlatformID <- function(item) .Call("ocl_get_platform_info", item)
oclInfo.list <- function(item) lapply(item, oclInfo)

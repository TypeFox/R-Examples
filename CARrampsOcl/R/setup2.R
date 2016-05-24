setup2 <- function(dev)
{

# kronVect 3Q

code2 <- c(
"#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n",
"__kernel void krony3(\n",
"__global double *result,\n",
"const unsigned int count,\n",
"__global double *a,\n",
"__global double *b,\n",
"__global double *c,\n",
"__global double *y,\n",
"const unsigned int acols,\n",
"const unsigned int arows,\n",
"const unsigned int bcols,\n",
"const unsigned int brows,\n",
"const unsigned int ccols,\n",
"const unsigned int crows)\n",
"{\n",
"// Vector element index\n",
"int nIndex = get_global_id(0);\n",
"result[nIndex] = 0;\n",
"double Csub = 0.0;\n",
"int arow = (int)(nIndex / (brows*crows) );\n",
"int brow = (int)  (nIndex - arow * brows * crows)  /  crows;\n",
"int crow =  nIndex % crows;\n",
"int acol = 0;\n",
"int bcol = 0;\n",
"int ccol = 0;\n",
"for (int i = 0; i < acols * bcols * ccols; i++) {\n",
"acol = (int)(i / (bcols*ccols));\n",
"bcol = (int) (i - acol * bcols * ccols)  /  ccols;\n",
"ccol = i % ccols;\n",
"Csub += ((a[arow + acol * arows] * b[brow + bcol * brows] * c[crow + ccol*crows]) * y[i]);\n",
"}\n",
"result[nIndex] = Csub ;\n",
"};\n"
)

oclSimpleKernel(dev[[1]], "krony3", code2, "double")
}

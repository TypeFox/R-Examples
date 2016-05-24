setup1 <- function(dev)
{
# kronVect 2Q
code <- c(
"#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n",
"__kernel void krony(\n",
"__global double *result,\n",
"const unsigned int count,\n",
"__global double *a,\n",
"__global double *b,\n",
"__global double *y,\n",
"const unsigned int acols,\n",
"const unsigned int arows,\n",
"const unsigned int bcols,\n",
"const unsigned int brows)\n",
"{\n",
"// Vector element index\n",
"int nIndex = get_global_id(0);\n",
"result[nIndex] = 0;\n",
"double Csub = 0.0;\n",
"int arow = (int)(nIndex / brows);\n",
"int brow = nIndex % brows;\n",
"int acol = 0;\n",
"int bcol = 0;\n",
"for (int i = 0; i < acols * bcols; i++) {\n",
"acol = (int)(i / bcols);\n",
"bcol = i % bcols;\n",
"Csub += ((a[arow + acol * arows] * b[brow + bcol * brows]) * y[i]);\n",
"}\n",
"result[nIndex] = Csub ;\n",
"};\n"
)

oclSimpleKernel(dev[[1]], "krony", code, "double")
}

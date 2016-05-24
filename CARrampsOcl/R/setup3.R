setup3 <- function(dev)
{
# sampling


 code4 = c(
      "#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n",
      "__kernel void sampling(\n",
      "  __global double* output,\n",
      " const unsigned int count,\n",
      "  __global double* smat,\n",
      "  __global double* D,\n",
      "  __global double* By,\n",
      "  __global double* alpha,\n",
      "  __global double* beta,\n",
      " const int nsamp,\n",
      " const int N,\n",
      " const int Fm1,\n",
      " const int k1)\n",
      "{\n",
      "  int idxtot = get_global_id(0), j, k;\n",
      "  double s0, whole = 0.0, logpostdensnumer = 0.0, logpostdensdenom ;\n",
      "  double neweigennumer, neweigen, newbeta ;\n",
      "if( idxtot < nsamp ) \n",
      "  { \n",
      " /* next line corrected to add alpha[0] 08/30/12 MKC */ \n" ,
      " double newalpha = (double) (N - k1) / 2.0 + alpha[0] ; \n" ,
      "   s0 = 1.0 ;\n",
      "   for( j = 0; j < Fm1 ; j++) \n",
      "      { \n",
      "   /* next line added 08/30/12 MKC */ \n" ,
      "        newalpha += alpha[j+1] ; \n",
      "        s0 -= smat[idxtot * Fm1 + j] ; \n",
      "      } \n",
      "   for( k = 0; k < N; k++) \n",
      "      { \n",
      "        neweigennumer = 0.0 ; \n",
      "        for( j = 0; j < Fm1; j++) \n",
      "            neweigennumer += smat[idxtot * Fm1 + j] * D[k*Fm1 + j] ; \n",
      "        neweigen = s0 * neweigennumer / (s0 + neweigennumer) ; \n",
      "        logpostdensnumer += neweigen > 0.0 ? log(neweigen) / 2.0 : 0.0 ;\n", 
      "        whole += neweigen * By[k] * By[k] ; \n",
      "      } \n",
      "   newbeta = whole / 2.0 + s0 * beta[0] ; \n",
      "   logpostdensnumer += log(s0) * (alpha[0] - 1.0) ; \n",
      "   for( j = 0; j < Fm1; j++) \n",
      "      { \n",
      "       logpostdensnumer += log(smat[idxtot * Fm1 + j]) * (alpha[j+1] - 1.0) ; \n",
      "       newbeta += smat[idxtot * Fm1 + j] * beta[j+1] ; \n",
      "      } \n",
      "   logpostdensdenom = log(newbeta) * newalpha ; \n",
      "   output[idxtot] = logpostdensnumer - logpostdensdenom ; \n",
      "   output[idxtot + nsamp] = newbeta ;\n",

      "  } \n",  

      "};")

oclSimpleKernel(dev[[1]], "sampling", code4, "double")


}

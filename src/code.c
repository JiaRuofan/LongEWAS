//NR method
//Current
#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>
#include <gsl/gsl_spblas.h>


double**** make4Darray(int a1, int a2, int a3, int a4){
  double ****tmp;
  tmp = (double ****)malloc(a1*sizeof(double ***));
  for(int i=0; i<a1; i++){
    tmp[i] = (double ***) malloc(a2*sizeof(double **));
    for(int j=0; j<a2; j++){
      tmp[i][j] = (double **) malloc(a3*sizeof(double*));
      for(int k=0; k<a3; k++){
        tmp[i][j][k] = (double *)malloc(a4*sizeof(double));
      }
    }
  }
  return tmp;
}

void delet4Darray(double ****tmp, int a1, int a2, int a3, int a4){
  for(int i=0; i<a1; i++){
    for(int j=0; j<a2; j++){
      for(int k=0; k<a3; k++){
        free(tmp[i][j][k]);
      }
      free(tmp[i][j]);
    }
    free(tmp[i]);
  }
  free(tmp);
}


double** make2Darray(int a1, int a2){
  double **tmp;
  tmp = (double **)malloc(a1*sizeof(double *));
  for(int i=0; i<a1; i++){
    tmp[i] = (double *) malloc(a2*sizeof(double));
  }
  return tmp;
}

void delet2Darray(double **tmp, int a1, int a2){
  for(int i=0; i<a1; i++){
    free(tmp[i]);
  }
  free(tmp);
}

double*** make3Darray(int a1, int a2, int a3){
  double ***tmp;
  tmp = (double ***)malloc(a1*sizeof(double **));
  for(int i=0; i<a1; i++){
    tmp[i] = (double **) malloc(a2*sizeof(double *));
    for(int j=0; j<a2; j++){
      tmp[i][j] = (double *) malloc(a3*sizeof(double));
    }
  }
  return tmp;
}

void delet3Darray(double ***tmp, int a1, int a2, int a3){
  for(int i=0; i<a1; i++){
    for(int j=0; j<a2; j++){
      free(tmp[i][j]);
    }
    free(tmp[i]);
  }
  free(tmp);
}

void gsl_matrix_inv(gsl_matrix *a)
{
  size_t n=a->size1;
  size_t m=a->size2;
  
  gsl_matrix *temp1=gsl_matrix_calloc(n,n);
  gsl_matrix_memcpy(temp1,a);
  
  gsl_permutation *p=gsl_permutation_calloc(n);
  int sign=0;
  gsl_linalg_LU_decomp(temp1,p,&sign);
  gsl_matrix *inverse=gsl_matrix_calloc(n,n);
  
  gsl_linalg_LU_invert(temp1,p,inverse);
  gsl_matrix_memcpy(a,inverse);
  
  gsl_permutation_free(p);
  gsl_matrix_free(temp1);
  gsl_matrix_free(inverse);
  
}


void MultiplyMatrix(double **a,double **b,double **c,int ROW, int COL,int RC){
  int i,j,k;
  
  for (i=0; i<ROW; i++){
    for(j=0; j<COL; j++){
      c[i][j]= 0 ;
      for(k=0; k<RC; k++){
        c[i][j]= c[i][j] + a[i][k]*b[k][j] ;
        //c[i][j] += a[i][k]*b[k][j] ;
      }
    }
  }
}

double Determinant(double **a,int n){
  int i,j,j1,j2;
  double det = 0;
  double **m=NULL;
  
  if (n < 1) { /* Error */
    
  } else if (n == 1) { /* Shouldn't get used */
    det = a[0][0];
  } else if (n == 2) {
    det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
  } else {
    det = 0;
    for (j1=0;j1<n;j1++) {
      m = (double **) malloc((n-1)*sizeof(double *));
      for (i=0;i<n-1;i++)
        m[i] = (double *) malloc((n-1)*sizeof(double));
      for (i=1;i<n;i++) {
        j2 = 0;
        for (j=0;j<n;j++) {
          if (j == j1)
            continue;
          m[i-1][j2] = a[i][j];
          j2++;
        }
      }
      det += pow(-1.0,j1+2.0) * a[0][j1] * Determinant(m,n-1);
      for (i=0;i<n-1;i++)
        free(m[i]);
      free(m);
    }
  }
  return(det);
}

/*Find the cofactor matrix of a square matrix*/

void CoFactor(double **a, int n, double **b){
  int i,j,ii,jj,i1,j1;
  double det;
  double **c;
  
  c = (double **) malloc((n-1)*sizeof(double *));
  for (i=0;i<n-1;i++)
    c[i] = (double *) malloc((n-1)*sizeof(double));
  
  for (j=0;j<n;j++) {
    for (i=0;i<n;i++) {
      
      /* Form the adjoint a_ij */
      i1 = 0;
      for (ii=0;ii<n;ii++) {
        if (ii == i)
          continue;
        j1 = 0;
        for (jj=0;jj<n;jj++) {
          if (jj == j)
            continue;
          c[i1][j1] = a[ii][jj];
          j1++;
        }
        i1++;
      }
      
      /* Calculate the determinate */
      det = Determinant(c,n-1);
      
      /* Fill in the elements of the cofactor */
      b[i][j] = pow(-1.0,i+j+2.0) * det;
    }
  }
  for (i=0;i<n-1;i++)
    free(c[i]);
  free(c);
}

/*Transpose of a square matrix, do it in place*/
void Transpose(double **a,int n){
  int i,j;
  double tmp;
  
  for (i=1;i<n;i++) {
    for (j=0;j<i;j++) {
      tmp = a[i][j];
      a[i][j] = a[j][i];
      a[j][i] = tmp;
    }
  }
}

void transpose(double **x, double **y,int p,int q)//Transpose
{
  
  for (int i = 0; i < p; i++)
  {
    for (int j = 0; j < q; j++)
    {
      y[j][i] = x[i][j];
    }
  }
}

/*calculate the inverse*/
void inverse(double **a, int n, double **a_inv){
  double det;
  double  **cofac_a;
  cofac_a = (double **)malloc(n*sizeof(double*));
  for(int i =0; i <n;i++){
    cofac_a[i] = (double *)malloc(n*sizeof(double));
  }
  CoFactor(a, n, cofac_a);
  Transpose(cofac_a, n); //turn the cofacotor matrix into the adjoint matrix
  det = Determinant(a, n); 
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      a_inv[i][j] = cofac_a[i][j] / det;
    }
  }
  
  for(int i=0; i<n; i++){
    free(cofac_a[i]);
  }
  free(cofac_a);
}

void gsl_hat(gsl_matrix *a, int row, int col, gsl_matrix *b){
  
}

void hat(double **B, int row, int col, double **Theta){
  double** BT=make2Darray(col,row);
  transpose(B,BT,row,col);
  double** BTB=make2Darray(col,col);
  MultiplyMatrix(BT,B,BTB,col,col,row);
  double** BTBI=make2Darray(col,col);
  inverse(BTB,col,BTBI);
  MultiplyMatrix(BTBI,BT,Theta,col,row,col);
  
  delet2Darray(BT,col,row);
  delet2Darray(BTB,col,col);
  delet2Darray(BTBI,col,col);
}

void hat1(double **B, int row, int col, double **Theta){
  gsl_matrix *BTB = gsl_matrix_alloc(col,col); 
  double** BT=make2Darray(col,row);
  transpose(B,BT,row,col);
  double value=0;
  double** BTBI=make2Darray(col,col);
  for (int i=0;i<col;i++){
    for (int j=0;j<(i+1);j++){
      value=0;
      for (int l=0;l<row;l++){
        value=value+B[l][i]*B[l][j];
      }
      gsl_matrix_set(BTB, i, j, value);
      gsl_matrix_set(BTB, j, i, value);
    }
  }
  gsl_matrix_inv(BTB);
  for (int i=0;i<col;i++){
    for (int j=0;j<col;j++){
      BTBI[i][j]=gsl_matrix_get(BTB,i,j);
    }
  }
  MultiplyMatrix(BTBI,BT,Theta,col,row,col);
  
  delet2Darray(BT,col,row);
  gsl_matrix_free(BTB);
  delet2Darray(BTBI,col,col);
}

void MatrixVector(double **X,double *y,double *b,int row, int col){
  
  for (int i=0; i<row; i++){
    b[i]=0;
    for(int j=0; j<col; j++){
      b[i]=b[i]+X[i][j]*y[j];
    }
  }
}

//========================================================================================================================
//quadratic programming
//========================================================================================================================
void quadprog(double **Dmat, double *dvec, double *xvec, int K){ //K is the dimension of xvec
  if(K==1){
    xvec[0] = 1;
  }else{
    double **Dmat_inv, *x_star, s1, s2;
    double Dmat_inv_sum=0, x_star_sum=0, *Dmat_inv_rowsum;
    int num_negatives=0;
    double *ind_negative;
    Dmat_inv = (double **)malloc(K*sizeof(double*));
    x_star = (double *)malloc(K*sizeof(double));
    Dmat_inv_rowsum = (double *)malloc(K*sizeof(double));
    ind_negative = (double *)malloc(K*sizeof(double));
    for(int k=0; k<K; k++){
      Dmat_inv[k] = (double *)malloc(K*sizeof(double));
    }
    inverse(Dmat, K, Dmat_inv);
    for(int k=0; k<K; k++){
      s1 = 0;
      s2 = 0;
      for(int k1=0; k1<K; k1++){
        s1 += Dmat_inv[k][k1]*dvec[k1];
        Dmat_inv_sum += Dmat_inv[k][k1];
        s2 += Dmat_inv[k][k1];
      }
      x_star[k] = s1;
      Dmat_inv_rowsum[k] = s2;
      x_star_sum += s1;
    }
    for(int k=0; k<K; k++){
      xvec[k] = x_star[k] + (1-x_star_sum)/Dmat_inv_sum*Dmat_inv_rowsum[k];
      if(xvec[k]<0){
        num_negatives++;
        ind_negative[k] = 1; 
      }else{
        ind_negative[k] = 0;
      }
    }
    free(x_star);
    free(Dmat_inv_rowsum);
    for(int k=0; k<K; k++){
      free(Dmat_inv[k]);
    }
    free(Dmat_inv);
    
    if(num_negatives == 0){
      free(ind_negative);
    }else{
      int Knew = K-num_negatives, i, j;
      double ** Dmat_new, *dvec_new, *xvec_sub;
      Dmat_new = (double **)malloc((Knew)*sizeof(double*));
      for(int k=0; k<Knew; k++){
        Dmat_new[k] = (double *) malloc((Knew)*sizeof(double));
      }
      dvec_new = (double *) malloc((Knew)*sizeof(double));
      xvec_sub = (double *) malloc((Knew)*sizeof(double));
      i = 0;
      
      for(int k1=0; k1<K; k1++){
        if(ind_negative[k1]==0){
          dvec_new[i] = dvec[k1];
          j = 0;
          for(int k2=0; k2<K; k2++){
            if(ind_negative[k2]==0){
              Dmat_new[i][j] = Dmat[k1][k2];
              j++;
            }
          }
          i++;
        }else{
          xvec[k1] = 0;
        }
      }
      
      quadprog(Dmat_new, dvec_new, xvec_sub, Knew);
      i=0;
      for(int k=0; k<K; k++){
        if(ind_negative[k]==0){
          xvec[k] = xvec_sub[i];
          i++;
        }
      }
      free(dvec_new);
      free(xvec_sub);
      for(int k=0; k<Knew;k++){
        free(Dmat_new[k]);
      }
      free(Dmat_new);
    }
    
  }
  
} 

double gsl_matrix_trace(gsl_matrix *a)
{
  double tr=0;
  for (size_t i=0;i<a->size1;i++)
  {
    tr=tr+gsl_matrix_get(a,i,i);
  }
  return(tr);
}


void gsl_matrix_mul(gsl_matrix *a,gsl_matrix *b,gsl_matrix *c)
{
  for (size_t i=0;i<a->size1;i++)
  {
    for (size_t j=0;j<b->size2;j++)
    {
      double sum=0.0;
      for (size_t k=0;k<b->size1;k++)
      {				   sum+=gsl_matrix_get(a,i,k)*gsl_matrix_get(b,k,j);
      }
      gsl_matrix_set(c,i,j,sum);
    }
  }
}

void gsl_matrix_mul_pa(gsl_matrix *a,gsl_matrix *b,gsl_matrix *c)
{
  
  double** ab=make2Darray(a->size1,b->size2);
#pragma omp parallel for
  for (size_t i=0;i<a->size1;i++)
  {
    for (size_t j=0;j<b->size2;j++)
    {
      ab[i][j]=0;
      for (size_t k=0;k<b->size1;k++)
      {				   ab[i][j]+=gsl_matrix_get(a,i,k)*gsl_matrix_get(b,k,j);
      }
      gsl_matrix_set(c,i,j,ab[i][j]);
    }
  }
  delet2Darray(ab,a->size1,b->size2);
  
}

void gsl_matrixvector_mul(gsl_matrix *a,gsl_vector *b,gsl_vector *c,int row,int col)
{
  double sum=0.0;
  for (int i=0;i<row;i++)
  {
    sum=0.0;
    for (int j=0;j<col;j++)
    {
      sum+=gsl_matrix_get(a,i,j)*gsl_vector_get(b,j);
      
    }
    gsl_vector_set(c,i,sum);
  }
}





double get_det(gsl_matrix * A)
{
  double det=0.0; 
  int n = A->size1;
  gsl_permutation *p = gsl_permutation_calloc(n);
  gsl_matrix *tmpA = gsl_matrix_calloc(n, n);
  int signum;
  gsl_matrix_memcpy(tmpA, A);
  gsl_linalg_LU_decomp(tmpA, p, &signum);
  det = gsl_linalg_LU_det(tmpA, signum);
  gsl_permutation_free(p);
  gsl_matrix_free(tmpA);
  return det;
}

double VMV(double *t1,double *t2,double **Dmat,int a, int b)
{
  double v=0;
  for (int i=0;i<a;i++){
    for (int j=0;j<b;j++){
      v=v+t1[i]*Dmat[i][j]*t2[j];
    }
  }
  return(v);
}


double gsl_VMV(gsl_vector *b,gsl_matrix *a,gsl_vector *c,int row,int col)
{
  double v=0;
  for (int i=0;i<row;i++){
    for (int j=0;j<col;j++){
      v=v+gsl_vector_get(b,i)*gsl_matrix_get(a,i,j)*gsl_vector_get(c,j);
    }
  }
  return(v);
}


void gsl_matrix_trans(gsl_matrix *a,gsl_matrix *b)
{
  for (size_t i=0;i<a->size1;i++)
  {
    for (size_t j=0;j<a->size2;j++)
    {
      gsl_matrix_set(b,j,i,gsl_matrix_get(a,i,j));
    }
  }
}

void IterativeSolver(gsl_spmatrix *A,gsl_vector *f,gsl_vector *u,int n){
  /* convert to compressed column format */
  gsl_spmatrix *C;
  C = gsl_spmatrix_ccs(A);
  
  /* now solve the system with the GMRES iterative solver */
  {
    const double tol = 1.0e-6;  /* solution relative tolerance */
  const size_t max_iter = 10; /* maximum iterations */
  const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
  gsl_splinalg_itersolve *work =
    gsl_splinalg_itersolve_alloc(T, n, 0);
  size_t iter = 0;
  double residual;
  int status;
  
  /* initial guess u = 0 */
  gsl_vector_set_zero(u);
  
  /* solve the system A u = f */
  do
  {
    status = gsl_splinalg_itersolve_iterate(C, f, tol, u, work);
    
    /* print out residual norm ||A*u - f|| */
    residual = gsl_splinalg_itersolve_normr(work);
    //fprintf(stderr, "iter %zu residual = %.12e\n", iter, residual);
    
    //if (status == GSL_SUCCESS)
    // fprintf(stderr, "Converged\n");
  }
  while (status == GSL_CONTINUE && ++iter < max_iter);
  
  /* output solution */
  //for (i = 0; i < n; ++i)
  //  {
  //    double xi = (i + 1) * h;
  //    double u_exact = sin(M_PI * xi);
  //   double u_gsl = gsl_vector_get(u, i);
  
  //   printf("%f %.12e %.12e\n", xi, u_gsl, u_exact);
  //  }
  
  gsl_splinalg_itersolve_free(work);
  }
}


int main(void){
  
  //unknown parameter
  int n=300;//sample size
  int type=6; //number of cell type
  int T=3; //period
  int S=1;  //number of covariates for y
  int G=20000; //number of CpG
  int Q=1; //number of covariates for x
  int nT=900;
  int m=3;
  int *t_num=(int *) malloc(n*sizeof(int));
  //double** tc=make2Darray(n,T);
  int **tc;
  tc = (int **)malloc(n*sizeof(int *));
  for(int i=0; i<n; i++){
    tc[i] = (int *) malloc(T*sizeof(int));
  }
  double*** O=make3Darray(n,G,T);
  double*** x=make3Darray(n,T,Q);
  double*** y=make3Darray(n,T,1+S);
  double*** z=make3Darray(n,T,m);
  double*** p=make3Darray(n,type,T);
  double*** p_old=make3Darray(n,type,T);
  double**** beta=make4Darray(G,type,T,Q);
  double*** gamma=make3Darray(G,type,1+S);
  double**** beta_var=make4Darray(G,type,T,Q);
  double*** gamma_var=make3Darray(G,type,1+S);
  double**** sigma_eta=make4Darray(G,type,m,m);
  double** sigma_delta=make2Darray(G,type);
  double *sigma_eps=(double *) malloc(G*sizeof(double));
  
  double**** mu=make4Darray(n,G,type,T);
  double*** Emu=make3Darray(n,G,(T+m)*type);
  double**** Emu2=make4Darray(n,G,(T+m)*type,(T+m)*type);
  int maxit=2000;
  int it=0;
  double value;
  int thre_num;
  double sigma_delta_old;
  
  
  
  gsl_vector *y_tilde_jk = gsl_vector_alloc(nT);
  gsl_matrix *X_tilde = gsl_matrix_alloc(nT,1+S+T*Q);
  gsl_matrix *X_tilde_T = gsl_matrix_alloc(1+S+T*Q,nT);
  gsl_matrix *XX_tilde = gsl_matrix_alloc(1+S+T*Q,1+S+T*Q);
  gsl_matrix *X_tilde_hat = gsl_matrix_alloc(1+S+T*Q,nT);
  gsl_vector *beta_hat = gsl_vector_alloc(1+S+T*Q);
  
  gsl_matrix_set_zero(X_tilde);
  gsl_matrix_set_zero(X_tilde_T);
  gsl_matrix_set_zero(XX_tilde);
  gsl_matrix_set_zero(X_tilde_hat);
  
  
  
  double LL;
  double LL_old;
  double *LL_vec=(double *) malloc(n*sizeof(double));
  
  
  double*** Dmat=make3Darray(n,type,type);
  double** d=make2Darray(n,type);
  
  double** pvec=make2Darray(n,type);
  
  double *value_vec=(double *) malloc(n*sizeof(double));
  double *value_vec_G=(double *) malloc(G*sizeof(double));
  
  double *value_temp1=(double *) malloc(G*sizeof(double));
  double *value_temp2=(double *) malloc(G*sizeof(double));
  int *nt_vec=(int *) malloc(G*sizeof(int));
  
  gsl_matrix *XP_tilde = gsl_matrix_alloc(nT,(1+S+T*Q)*type);
  gsl_matrix *XP_tilde_T = gsl_matrix_alloc((1+S+T*Q)*type,nT);
  gsl_matrix *W = gsl_matrix_alloc(nT,nT);
  gsl_matrix *XW = gsl_matrix_alloc((1+S+T*Q)*type,nT);
  gsl_matrix *XWX = gsl_matrix_alloc((1+S+T*Q)*type,(1+S+T*Q)*type);
  
  
  
  
  
  double**** beta_z=make4Darray(G,type,T,Q);
  double*** gamma_z=make3Darray(G,type,1+S);
  
  
  
  
  FILE* fp11=fopen("t_num.txt","r"); //Ã¦â€°â€œÃ¥Â¼â‚¬Ã¦â€“â€¡Ã¤Â»Â¶
  if(fp11==NULL)
  {
    printf("no file");
    return -1;
  }
  for(int i=0;i<n;i++)
  {
    
    fscanf(fp11,"%d",&t_num[i]);/*Ã¦Â¯ÂÃ¦Â¬Â¡Ã¨Â¯Â»Ã¥Ââ€“Ã¤Â¸â‚¬Ã¤Â¸ÂªÃ¦â€¢Â°Ã¯Â¼Å’fscanfÃ¥â€¡Â½Ã¦â€¢Â°Ã©Ââ€¡Ã¥Ë†Â°Ã§Â©ÂºÃ¦Â Â¼Ã¦Ë†â€“Ã¨â‚¬â€¦Ã¦ÂÂ¢Ã¨Â¡Å’Ã§Â»â€œÃ¦ÂÅ¸*/
  
  }
  fclose(fp11);
  
  
  
  FILE* fp2=fopen("tc.txt","r"); //Ã¦â€°â€œÃ¥Â¼â‚¬Ã¦â€“â€¡Ã¤Â»Â¶
  if(fp2==NULL)
  {
    printf("no file");
    return -1;
  }
  for(int i=0;i<n;i++)
  {
    for(int t=0;t<t_num[i];t++)
    {
      fscanf(fp2,"%d",&tc[i][t]);/*Ã¦Â¯ÂÃ¦Â¬Â¡Ã¨Â¯Â»Ã¥Ââ€“Ã¤Â¸â‚¬Ã¤Â¸ÂªÃ¦â€¢Â°Ã¯Â¼Å’fscanfÃ¥â€¡Â½Ã¦â€¢Â°Ã©Ââ€¡Ã¥Ë†Â°Ã§Â©ÂºÃ¦Â Â¼Ã¦Ë†â€“Ã¨â‚¬â€¦Ã¦ÂÂ¢Ã¨Â¡Å’Ã§Â»â€œÃ¦ÂÅ¸*/
    }
  }
  fclose(fp2);
  
  
  
  FILE* fp4=fopen("X.txt","r"); //Ã¦â€°â€œÃ¥Â¼â‚¬Ã¦â€“â€¡Ã¤Â»Â¶
  if(fp4==NULL)
  {
    printf("no file");
    return -1;
  }
  for(int i=0;i<n;i++)
  {
    for(int t=0;t<t_num[i];t++)
    {
      for(int j=0;j<Q;j++)
      {
        fscanf(fp4,"%lf",&x[i][tc[i][t]][j]);/*Ã¦Â¯ÂÃ¦Â¬Â¡Ã¨Â¯Â»Ã¥Ââ€“Ã¤Â¸â‚¬Ã¤Â¸ÂªÃ¦â€¢Â°Ã¯Â¼Å’fscanfÃ¥â€¡Â½Ã¦â€¢Â°Ã©Ââ€¡Ã¥Ë†Â°Ã§Â©ÂºÃ¦Â Â¼Ã¦Ë†â€“Ã¨â‚¬â€¦Ã¦ÂÂ¢Ã¨Â¡Å’Ã§Â»â€œÃ¦ÂÅ¸*/
      }
    }
  }
  fclose(fp4);
  
  
  for(int i=0;i<n;i++)
  {
    for(int t=0;t<t_num[i];t++)
    {
      
        z[i][tc[i][t]][2]=x[i][tc[i][t]][0];
      }}
  
  
  FILE* fp1=fopen("Y.txt","r"); //Ã¦â€°â€œÃ¥Â¼â‚¬Ã¦â€“â€¡Ã¤Â»Â¶
  if(fp1==NULL)
  {
    printf("no file");
    return -1;
  }
  for(int i=0;i<n;i++)
  {
    for(int t=0;t<t_num[i];t++)
    {
      for(int j=0;j<(1+S);j++)
      {
        fscanf(fp1,"%lf",&y[i][tc[i][t]][j]);/*Ã¦Â¯ÂÃ¦Â¬Â¡Ã¨Â¯Â»Ã¥Ââ€“Ã¤Â¸â‚¬Ã¤Â¸ÂªÃ¦â€¢Â°Ã¯Â¼Å’fscanfÃ¥â€¡Â½Ã¦â€¢Â°Ã©Ââ€¡Ã¥Ë†Â°Ã§Â©ÂºÃ¦Â Â¼Ã¦Ë†â€“Ã¨â‚¬â€¦Ã¦ÂÂ¢Ã¨Â¡Å’Ã§Â»â€œÃ¦ÂÅ¸*/
      }
    }
  }
  fclose(fp1);
  
  
  for(int i=0;i<n;i++)
  {
    for(int t=0;t<t_num[i];t++)
    {
      for(int j=0;j<(1+S);j++)
      {
        z[i][tc[i][t]][j]=y[i][tc[i][t]][j];
      }}}
  
  //   FILE* fp3=fopen("O.txt","r"); //Ã¦â€°â€œÃ¥Â¼â‚¬Ã¦â€“â€¡Ã¤Â»Â¶
  //   if(fp3==NULL)
  //   {
  //     printf("no file");
  //    return -1;
  //  }
  //  for(int i=0;i<n;i++)
  //  {
  
  //     for(int j=0;j<G;j++)
  //    {
  //      for(int t=0;t<T;t++)
  //      {
  //      fscanf(fp3,"%lf",&O[i][j][t]);/*Ã¦Â¯ÂÃ¦Â¬Â¡Ã¨Â¯Â»Ã¥Ââ€“Ã¤Â¸â‚¬Ã¤Â¸ÂªÃ¦â€¢Â°Ã¯Â¼Å’fscanfÃ¥â€¡Â½Ã¦â€¢Â°Ã©Ââ€¡Ã¥Ë†Â°Ã§Â©ÂºÃ¦Â Â¼Ã¦Ë†â€“Ã¨â‚¬â€¦Ã¦ÂÂ¢Ã¨Â¡Å’Ã§Â»â€œÃ¦ÂÅ¸*/
  //    }
  //  }
  //  }
  //   fclose(fp3);
  
  
  
  char Ofilename[10];
  FILE *fpSet[n];
  for (int i=0;i<(n);i++){
    sprintf(Ofilename,"O%d.txt",(i+1));
    fpSet[i]=fopen(Ofilename,"r");  
    if(!fpSet[i])
    {
      printf("no file");
      return -1;
    }
    
    
    
    for(int j=0;j<G;j++)
    {
      for(int t=0;t<t_num[i];t++)
      {
        fscanf(fpSet[i],"%lf",&O[i][j][tc[i][t]]);/*Ã¦Â¯ÂÃ¦Â¬Â¡Ã¨Â¯Â»Ã¥Ââ€“Ã¤Â¸â‚¬Ã¤Â¸ÂªÃ¦â€¢Â°Ã¯Â¼Å’fscanfÃ¥â€¡Â½Ã¦â€¢Â°Ã©Ââ€¡Ã¥Ë†Â°Ã§Â©ÂºÃ¦Â Â¼Ã¦Ë†â€“Ã¨â‚¬â€¦Ã¦ÂÂ¢Ã¨Â¡Å’Ã§Â»â€œÃ¦ÂÅ¸*/
      }
    }
    fclose(fpSet[i]);
  }
  
  
  //printf("%d\n",1);
  
  FILE* fp6=fopen("pinit.txt","r"); //Ã¦â€°â€œÃ¥Â¼â‚¬Ã¦â€“â€¡Ã¤Â»Â¶
  if(fp6==NULL)
  {
    printf("no file");
    return -1;
  }
  for(int i=0;i<n;i++)
  {
    for(int t=0;t<t_num[i];t++)
    {
      for(int j=0;j<type;j++)
      {
        
        fscanf(fp6,"%lf",&p[i][j][tc[i][t]]);/*Ã¦Â¯ÂÃ¦Â¬Â¡Ã¨Â¯Â»Ã¥Ââ€“Ã¤Â¸â‚¬Ã¤Â¸ÂªÃ¦â€¢Â°Ã¯Â¼Å’fscanfÃ¥â€¡Â½Ã¦â€¢Â°Ã©Ââ€¡Ã¥Ë†Â°Ã§Â©ÂºÃ¦Â Â¼Ã¦Ë†â€“Ã¨â‚¬â€¦Ã¦ÂÂ¢Ã¨Â¡Å’Ã§Â»â€œÃ¦ÂÅ¸*/
        
      }
    }
  }
  fclose(fp6);
  
  
   for (int j=0;j<G;j++){
    
    sigma_eps[j]=0.001;
    
  }
  
  
  
  
  for (int j=0;j<G;j++){
    for (int ty=0;ty<type;ty++){
      sigma_delta[j][ty]=0.01;
    }
  }
  
  
  for (int j=0;j<G;j++){
    for (int ty=0;ty<type;ty++){
      for (int ii=0;ii<m;ii++){
        sigma_eta[j][ty][ii][ii]=0.01;
      }
    }
  }
  
  
  for (int j=0;j<G;j++){
    for (int t=0;t<T;t++){
      gamma[j][t][0]=0.5;
    }
  }
  
  
  
  int nt=0;
  
  for (int i=0;i<n;i++){
    for (int t=0;t<t_num[i];t++){
      
      
      for (int s=0;s<(1+S);s++){
        gsl_matrix_set(X_tilde, nt, s,y[i][tc[i][t]][s]);
      }
      for (int l=0;l<Q;l++){
        gsl_matrix_set(X_tilde, nt, 1+S+tc[i][t]*Q+l,x[i][tc[i][t]][l]);
      }
      
      
      nt=nt+1;
    }
  }
  
  
  
  
  gsl_matrix_trans(X_tilde,X_tilde_T);
  gsl_matrix_mul(X_tilde_T,X_tilde,XX_tilde);
  gsl_matrix_inv(XX_tilde);
  gsl_matrix_mul(XX_tilde,X_tilde_T,X_tilde_hat);
  
  
  
  while(1){
    
    it=it+1;
    if (it>maxit){break;}
    
    
    ////E step
    
    for (int i=0;i<n;i++){
      for (int j=0;j<G;j++){
        
        
        for (int k=0;k<type;k++){
          for (int t=0;t<t_num[i];t++){
            mu[i][j][k][tc[i][t]]=0;
            for (int s=0;s<(S+1);s++){
              mu[i][j][k][tc[i][t]]=mu[i][j][k][tc[i][t]]+gamma[j][k][s]*y[i][tc[i][t]][s];
            }
            
            for (int l=0;l<Q;l++){
              mu[i][j][k][tc[i][t]]=mu[i][j][k][tc[i][t]]+beta[j][k][tc[i][t]][l]*x[i][tc[i][t]][l];
            }
            
          }  
        }
      }}
    
    // printf("%lf\n",mu[6][6][3][2]);
#pragma omp parallel for    
    for (int i=0;i<n;i++){
      for (int j=0;j<G;j++){
        
        
        gsl_vector *L_ij = gsl_vector_alloc(type*(t_num[i]+m));
        gsl_vector *Emu_ij = gsl_vector_alloc(type*(t_num[i]+m));
        gsl_matrix *Psi_ij = gsl_matrix_alloc(type*(t_num[i]+m),type*(t_num[i]+m));//v,w
        double** v_inverse=make2Darray(m,m);
        
        gsl_matrix_set_zero(Psi_ij);
        gsl_vector_set_zero(L_ij);
        gsl_vector_set_zero(Emu_ij);
        /////////////////L_ij
        for (int t=0;t<t_num[i];t++){
          
          for (int k=0;k<type;k++){
            
            value_vec[i]=((O[i][j][tc[i][t]]*p[i][k][tc[i][t]]/sigma_eps[j])+(mu[i][j][k][tc[i][t]]/sigma_delta[j][k]));
            //printf("%lf\n",p[i][k][tc[i][t]]);
            gsl_vector_set(L_ij, (t*type+k), value_vec[i]);
            
          }
          
          
        }
        
        
        //printf("%lf\n",gsl_vector_get(L_ij,1));
        
        
        for (int s=0;s<m;s++){
          
          for (int k=0;k<type;k++){
            value_vec[i]=0;
            for (int t=0;t<t_num[i];t++){
              value_vec[i]=value_vec[i]+mu[i][j][k][tc[i][t]]*z[i][tc[i][t]][s];
            }
            value_vec[i]=-(value_vec[i]/sigma_delta[j][k]);
            gsl_vector_set(L_ij, (type*t_num[i]+k*m+s), value_vec[i]);
          }
          
        }
        
        
        /////////////////Psi_ij
        
        //diag
        for (int t=0;t<t_num[i];t++){
          
          for (int k=0;k<type;k++){
            value_vec[i]=(p[i][k][tc[i][t]]*p[i][k][tc[i][t]]/sigma_eps[j]+1/sigma_delta[j][k]);
            gsl_matrix_set(Psi_ij, (t*type+k),(t*type+k),value_vec[i]);
            
          }
          
        }
        
        
        for (int t=0;t<t_num[i];t++){
          
          for (int k1=1;k1<type;k1++){
            
            for (int k2=0;k2<k1;k2++){
              value_vec[i]=(p[i][k1][tc[i][t]]*p[i][k2][tc[i][t]]/sigma_eps[j]);
              gsl_matrix_set(Psi_ij, t*type+k1,t*type+k2,value_vec[i]);
              gsl_matrix_set(Psi_ij, t*type+k2,t*type+k1,value_vec[i]);
            }
          }
          
        }
        
        
        
        
        for (int s=0;s<m;s++){
          for (int k=0;k<type;k++){
            
            for (int t=0;t<t_num[i];t++){
              value_vec[i]=-z[i][tc[i][t]][s]/sigma_delta[j][k];
              gsl_matrix_set(Psi_ij, t*type+k,t_num[i]*type+k*m+s,value_vec[i]);
              gsl_matrix_set(Psi_ij, t_num[i]*type+k*m+s, t*type+k,value_vec[i]);
            }
          }
        }  
        
        
        
        for (int k=0;k<type;k++){
          
          inverse(sigma_eta[j][k],m,v_inverse);
          
          //gsl_matrix_set(Psi_ij, t_num[i]*type+k*(1+S+Q), t_num[i]*type+k*(1+S+Q),v_inverse[0][0]+t_num[i]/sigma_delta[j][k]);
          
          for (int s=0;s<m;s++){
            
            for (int s1=0;s1<m;s1++){
              
              value_vec[i]=0;
              for (int t=0;t<t_num[i];t++){
                value_vec[i]=value_vec[i]+z[i][tc[i][t]][s]*z[i][tc[i][t]][s1]/sigma_delta[j][k];
              }
              
              gsl_matrix_set(Psi_ij, t_num[i]*type+k*m+s1, t_num[i]*type+k*m+s,v_inverse[s][s1]+value_vec[i]);
              
              
            }
            
            
          }
          
          
        }
        
        
        //
        
        
        
        
        
        
        
        gsl_matrix_inv(Psi_ij);
        gsl_matrixvector_mul(Psi_ij,L_ij,Emu_ij,(t_num[i]+m)*type,(t_num[i]+m)*type);
        
        
        
        for (int t=0;t<t_num[i];t++){
          for (int k=0;k<type;k++){
            Emu[i][j][tc[i][t]*type+k]=gsl_vector_get(Emu_ij,t*type+k);
          }
        }
        
        for (int k=0;k<type;k++){
          for (int ii=0;ii<m;ii++){
            Emu[i][j][T*type+k*m+ii]=gsl_vector_get(Emu_ij,t_num[i]*type+k*m+ii);
          }
        }
        
        for (int t1=0;t1<t_num[i];t1++){
          for (int k1=0;k1<type;k1++){
            for (int t2=0;t2<t_num[i];t2++){
              for (int k2=0;k2<type;k2++){
                Emu2[i][j][tc[i][t1]*type+k1][tc[i][t2]*type+k2]=Emu[i][j][tc[i][t1]*type+k1]*Emu[i][j][tc[i][t2]*type+k2]+gsl_matrix_get(Psi_ij,t1*type+k1,t2*type+k2);
              }
            }
          }
        }
        
        
        for (int k1=0;k1<type;k1++){
          for (int k2=0;k2<type;k2++){
            
            for (int s1=0;s1<m;s1++){
              for (int s2=0;s2<m;s2++){
                Emu2[i][j][T*type+k1*m+s1][T*type+k2*m+s2]=Emu[i][j][T*type+k1*m+s1]*Emu[i][j][T*type+k2*m+s2]+gsl_matrix_get(Psi_ij,t_num[i]*type+k1*m+s1,t_num[i]*type+k2*m+s2);
              }
            }
            
            
          }
        }
        
        for (int k1=0;k1<type;k1++){
          for (int t2=0;t2<t_num[i];t2++){
            for (int k2=0;k2<type;k2++){
              for (int s=0;s<m;s++){
                Emu2[i][j][T*type+k1*m+s][tc[i][t2]*type+k2]=Emu[i][j][T*type+k1*m+s]*Emu[i][j][tc[i][t2]*type+k2]+gsl_matrix_get(Psi_ij,t_num[i]*type+k1*m+s,t2*type+k2);
                Emu2[i][j][tc[i][t2]*type+k2][T*type+k1*m+s]=Emu2[i][j][T*type+k1*m+s][tc[i][t2]*type+k2];
              }
              
            }
          }
        }
        
        gsl_matrix_free(Psi_ij);
        gsl_vector_free(Emu_ij);
        gsl_vector_free(L_ij);
        delet2Darray(v_inverse,m,m);
      }}
    
    
#pragma omp parallel for    
    for (int j=0;j<G;j++){
      
      gsl_vector *y_tilde_jk = gsl_vector_alloc(nT);
      gsl_vector *beta_hat = gsl_vector_alloc(1+S+T*Q);
      
      for (int k=0;k<type;k++){
        nt_vec[j]=0;
        for (int i=0;i<n;i++){
          for (int t=0;t<t_num[i];t++){
            value_vec_G[j]=0;
            for (int s=0;s<m;s++){
              value_vec_G[j]=value_vec_G[j]-Emu[i][j][T*type+k*m+s]*z[i][tc[i][t]][s];
            }
            gsl_vector_set(y_tilde_jk, nt_vec[j], Emu[i][j][tc[i][t]*type+k]+value_vec_G[j]);
            nt_vec[j]=nt_vec[j]+1;
          }
        }
        
        
        gsl_matrixvector_mul(X_tilde_hat,y_tilde_jk,beta_hat,1+S+T*Q,nT);
        
        
        //beta=make4Darray(G,type,T,S);
        
        for (int t=0;t<T;t++){
          for (int q=0;q<Q;q++){
            beta[j][k][t][q]=gsl_vector_get(beta_hat, 1+S+t*Q+q);
          }
        }
        
        
        for (int s=0;s<(S+1);s++){
          gamma[j][k][s]=gsl_vector_get(beta_hat, s);
        }
        
        
      }
      
      gsl_vector_free(y_tilde_jk);
      gsl_vector_free(beta_hat);
    }
    
    
    ///p
    
    
    for (int i=0;i<n;i++){
      for (int k=0;k<type;k++){
        for (int t=0;t<t_num[i];t++){
          p_old[i][k][tc[i][t]]=p[i][k][tc[i][t]];
        }
      }
    }
    
    
#pragma omp parallel for     
    for (int i=0;i<n;i++){
      
      for (int t=0;t<t_num[i];t++){
        
        
        for (int k1=0;k1<type;k1++){
          for (int k2=0;k2<type;k2++){
            
            Dmat[i][k1][k2]=0;
            for (int j=0;j<G;j++){
              Dmat[i][k1][k2]=Dmat[i][k1][k2]+Emu2[i][j][tc[i][t]*type+k1][tc[i][t]*type+k2]/sigma_eps[j];
            }
          }
        }
        
        
        for (int k2=0;k2<type;k2++){
          
          d[i][k2]=0;
          for (int j=0;j<G;j++){
            d[i][k2]=d[i][k2]+O[i][j][tc[i][t]]*Emu[i][j][tc[i][t]*type+k2]/sigma_eps[j];
          }
        }
        
        
        
        quadprog(Dmat[i],d[i],pvec[i],type);
        
        
        for (int k=0;k<type;k++){
          p[i][k][tc[i][t]]=pvec[i][k];
        }
        
      }
    }
    
    
    
    
    
    ///sigma_delta
#pragma omp parallel for   
    for (int i=0;i<n;i++){
      for (int j=0;j<G;j++){
        
        
        for (int k=0;k<type;k++){
          for (int t=0;t<t_num[i];t++){
            mu[i][j][k][tc[i][t]]=0;
            for (int s=0;s<(S+1);s++){
              mu[i][j][k][tc[i][t]]=mu[i][j][k][tc[i][t]]+gamma[j][k][s]*y[i][tc[i][t]][s];
            }
            
            for (int l=0;l<Q;l++){
              mu[i][j][k][tc[i][t]]=mu[i][j][k][tc[i][t]]+beta[j][k][tc[i][t]][l]*x[i][tc[i][t]][l];
            }
            
          }  
        }
      }}
    
    
    
    // thre_num=0;
#pragma omp parallel for   
    for (int j=0;j<G;j++){
      
      for (int k=0;k<type;k++){
        
        sigma_delta_old=sigma_delta[j][k];
        
        value_vec_G[j]=0;
        for (int i=0;i<n;i++){
          for (int t=0;t<t_num[i];t++){
            
            value_vec_G[j]=value_vec_G[j]+mu[i][j][k][tc[i][t]]*mu[i][j][k][tc[i][t]]+Emu2[i][j][tc[i][t]*type+k][tc[i][t]*type+k]-2*Emu[i][j][tc[i][t]*type+k]*mu[i][j][k][tc[i][t]];
            
            for (int s=0;s<m;s++){
              value_vec_G[j]=value_vec_G[j]-2*z[i][tc[i][t]][s]*Emu2[i][j][tc[i][t]*type+k][T*type+k*m+s]+2*z[i][tc[i][t]][s]*Emu[i][j][T*type+k*m+s]*mu[i][j][k][tc[i][t]];
            }
            
            
            for (int s1=0;s1<m;s1++){
              for (int s2=0;s2<m;s2++){
                value_vec_G[j]=value_vec_G[j]+z[i][tc[i][t]][s1]*z[i][tc[i][t]][s2]*Emu2[i][j][T*type+k*m+s1][T*type+k*m+s2];
              }
            }
            
            
            
          }
          
        }
        
        sigma_delta[j][k]=value_vec_G[j]/(nT);
        //printf("%lf\n",value/(n*T));
        // if (fabs(sigma_delta_old-sigma_delta[j][k])>0.00001){thre_num=thre_num+1;}
      }
    }
    //printf("%lf\n",sigma_delta[67][4]);
    // if (thre_num<1){break;}
    
    ///sigma_eta
#pragma omp parallel for     
    for (int j=0;j<G;j++){
      
      for (int k=0;k<type;k++){
        
        
        for (int ii1=0;ii1<m;ii1++){
          for (int ii2=0;ii2<m;ii2++){
            
            value_vec_G[j]=0;
            for (int i=0;i<n;i++){
              
              value_vec_G[j]=value_vec_G[j]+Emu2[i][j][T*type+k*m+ii1][T*type+k*m+ii2];
            }
            
            sigma_eta[j][k][ii1][ii2]=value_vec_G[j]/n;
            
          }
        }
        
        
      }
    }
    
    
#pragma omp parallel for   
    for (int i=0;i<n;i++){
      LL_vec[i]=0;
      for (int j=0;j<G;j++){
        
        gsl_matrix *L_mat = gsl_matrix_alloc(t_num[i],t_num[i]);
        gsl_vector *L_vec = gsl_vector_alloc(t_num[i]);
        gsl_matrix_set_zero(L_mat);
        gsl_vector_set_zero(L_vec);
        double *zit1= (double *)malloc(m*sizeof(double));
        double *zit2= (double *)malloc(m*sizeof(double));
        
        for (int t=0;t<t_num[i];t++){
          
          value_vec[i]=O[i][j][tc[i][t]];
          for (int k=0;k<type;k++){
            
            value_vec[i]=value_vec[i]-p[i][k][tc[i][t]]*mu[i][j][k][tc[i][t]];
            
          }
          
          gsl_vector_set(L_vec, t, value_vec[i]);
          
        }
        
        for (int t1=1;t1<t_num[i];t1++){
          for (int t2=0;t2<t1;t2++){
            
            for (int ii=0;ii<m;ii++){
              zit1[ii]=z[i][tc[i][t1]][ii];
              zit2[ii]=z[i][tc[i][t2]][ii];
            }
            
            
            value_vec[i]=0;
            for (int k=0;k<type;k++){
              
              value_vec[i]=value_vec[i]+p[i][k][tc[i][t1]]*p[i][k][tc[i][t2]]*VMV(zit1,zit2,sigma_eta[j][k],m,m);
              
            }
            gsl_matrix_set(L_mat, t1, t2, value_vec[i]);
            gsl_matrix_set(L_mat, t2, t1, value_vec[i]);
          }
        }
        
        
        for (int t=0;t<t_num[i];t++){
          
          
          for (int ii=0;ii<m;ii++){
            zit1[ii]=z[i][tc[i][t]][ii];
          }
          
          
          value_vec[i]=0;
          for (int k=0;k<type;k++){
            
            value_vec[i]=value_vec[i]+p[i][k][tc[i][t]]*p[i][k][tc[i][t]]*VMV(zit1,zit1,sigma_eta[j][k],m,m)+p[i][k][tc[i][t]]*p[i][k][tc[i][t]]*sigma_delta[j][k];
            
          }
          value_vec[i]=value_vec[i]+sigma_eps[j];
          gsl_matrix_set(L_mat, t, t, value_vec[i]);
        }
        gsl_matrix_inv(L_mat);
        LL_vec[i]=LL_vec[i]+0.5*get_det(L_mat)-gsl_VMV(L_vec,L_mat,L_vec,t_num[i],t_num[i]);
        
        gsl_matrix_free(L_mat);
        gsl_vector_free(L_vec);
        free(zit1);
        free(zit2);
      }
      
      
    }
    
    LL=0;
    
    for (int i=0;i<n;i++){
      
      
      LL=LL_vec[i]+LL;
      
    }
    
    
    
#pragma omp parallel for   
    for (int i=0;i<n;i++){
      LL_vec[i]=0;
      for (int j=0;j<G;j++){
        
        gsl_matrix *L_mat = gsl_matrix_alloc(t_num[i],t_num[i]);
        gsl_vector *L_vec = gsl_vector_alloc(t_num[i]);
        gsl_matrix_set_zero(L_mat);
        gsl_vector_set_zero(L_vec);
        double *zit1= (double *)malloc(m*sizeof(double));
        double *zit2= (double *)malloc(m*sizeof(double));
        
        for (int t=0;t<t_num[i];t++){
          
          value_vec[i]=O[i][j][tc[i][t]];
          for (int k=0;k<type;k++){
            
            value_vec[i]=value_vec[i]-p_old[i][k][tc[i][t]]*mu[i][j][k][tc[i][t]];
            
          }
          
          gsl_vector_set(L_vec, t, value_vec[i]);
          
        }
        
        for (int t1=1;t1<t_num[i];t1++){
          for (int t2=0;t2<t1;t2++){
            
            for (int ii=0;ii<m;ii++){
              zit1[ii]=z[i][tc[i][t1]][ii];
              zit2[ii]=z[i][tc[i][t2]][ii];
            }
            
            value_vec[i]=0;
            for (int k=0;k<type;k++){
              
              value_vec[i]=value_vec[i]+p_old[i][k][tc[i][t1]]*p_old[i][k][tc[i][t2]]*VMV(zit1,zit2,sigma_eta[j][k],m,m);
              
            }
            gsl_matrix_set(L_mat, t1, t2, value_vec[i]);
            gsl_matrix_set(L_mat, t2, t1, value_vec[i]);
          }
        }
        
        
        for (int t=0;t<t_num[i];t++){
          
          
          for (int ii=0;ii<m;ii++){
            zit1[ii]=z[i][tc[i][t]][ii];
          }
          
          
          
          value_vec[i]=0;
          for (int k=0;k<type;k++){
            
            value_vec[i]=value_vec[i]+p_old[i][k][tc[i][t]]*p_old[i][k][tc[i][t]]*VMV(zit1,zit1,sigma_eta[j][k],m,m)+p_old[i][k][tc[i][t]]*p_old[i][k][tc[i][t]]*sigma_delta[j][k];
            
          }
          value_vec[i]=value_vec[i]+sigma_eps[j];
          gsl_matrix_set(L_mat, t, t, value_vec[i]);
        }
        gsl_matrix_inv(L_mat);
        LL_vec[i]=LL_vec[i]+0.5*get_det(L_mat)-gsl_VMV(L_vec,L_mat,L_vec,t_num[i],t_num[i]);
        
        gsl_matrix_free(L_mat);
        gsl_vector_free(L_vec);
        free(zit1);
        free(zit2);
      }
      
      
    }
    
    
    
    LL_old=0;
    
    for (int i=0;i<n;i++){
      
      
      LL_old=LL_vec[i]+LL_old;
      
    }
    
    
    if (LL_old>LL){
      for (int i=0;i<n;i++){
        for (int k=0;k<type;k++){
          for (int t=0;t<t_num[i];t++){
            p[i][k][tc[i][t]]=p_old[i][k][tc[i][t]];
          }
        }
      }
      
      LL=LL_old;
    }
    
    printf("%lf\n",LL);
    
    ///sigma_eps
#pragma omp parallel for    
    for (int j=0;j<G;j++){
      value_vec_G[j]=0;
      for (int i=0;i<n;i++){
        for (int t=0;t<t_num[i];t++){
          
          for (int k=0;k<type;k++){
            value_vec_G[j]=value_vec_G[j]+Emu2[i][j][tc[i][t]*type+k][tc[i][t]*type+k]*p[i][k][tc[i][t]]*p[i][k][tc[i][t]];
          }
          
          
          for (int k=0;k<type;k++){
            value_vec_G[j]=value_vec_G[j]-2*O[i][j][tc[i][t]]*Emu[i][j][tc[i][t]*type+k]*p[i][k][tc[i][t]];
          }
          
          
          for (int k1=1;k1<type;k1++){
            for (int k2=0;k2<k1;k2++){
              value_vec_G[j]=value_vec_G[j]+2*Emu2[i][j][tc[i][t]*type+k1][tc[i][t]*type+k2]*p[i][k1][tc[i][t]]*p[i][k2][tc[i][t]];
            }
          }
          
          value_vec_G[j]=value_vec_G[j]+O[i][j][tc[i][t]]*O[i][j][tc[i][t]];
        }
      }
      sigma_eps[j]=value_vec_G[j]/(nT);
      //printf("%lf\n",sigma_eps[j]);
    }
    //printf("%lf\n",sigma_eps[4]);
    printf("%d\n",it);
    
    
    
    if (it%500<1){
      
      FILE * pFile019;
      const char *pFileName019="gamma_temp1.txt";
      pFile019 = fopen(pFileName019, "w");//Ã¨Â¿â„¢Ã¤Â¸ÂªÃ§â€Â¨Ã¢â‚¬Å“wÃ¢â‚¬ÂÃ¦ËœÂ¯Ã¥â€ â„¢Ã¦â€“â€¡Ã¤Â»Â¶Ã¯Â¼Å’Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Å½Å¸Ã¥â€ â€¦Ã¥Â®Â¹Ã¯Â¼Å’Ã¨â€¹Â¥Ã¤Â¸ÂÃ¦Æ’Â³Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Ë†â„¢Ã§â€Â¨Ã¢â‚¬Å“aÃ¢â‚¬Â
      
      for (int j = 0; j < G;j++){
        for (int k = 0; k < type; k++){
          for (int s = 0; s < (1+S); s++){
            fprintf(pFile019,"%.3f\n", gamma[j][k][s]);//Ã¨Â¿â„¢Ã©â€¡Å’Ã¥Â¾ÂªÃ§Å½Â¯Ã¥â€ â„¢Ã¥â€¦Â¥Ã¦â€“â€¡Ã¤Â»Â¶ 3Ã¤Â¸Âª 
          }
        }
      }
      
      fclose(pFile019);
      
      
      
      FILE * pFile014;
      const char *pFileName014="beta_temp1.txt";
      pFile014 = fopen(pFileName014, "w");//Ã¨Â¿â„¢Ã¤Â¸ÂªÃ§â€Â¨Ã¢â‚¬Å“wÃ¢â‚¬ÂÃ¦ËœÂ¯Ã¥â€ â„¢Ã¦â€“â€¡Ã¤Â»Â¶Ã¯Â¼Å’Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Å½Å¸Ã¥â€ â€¦Ã¥Â®Â¹Ã¯Â¼Å’Ã¨â€¹Â¥Ã¤Â¸ÂÃ¦Æ’Â³Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Ë†â„¢Ã§â€Â¨Ã¢â‚¬Å“aÃ¢â‚¬Â
      for (int j = 0; j < G; j++){
        
        for (int ty = 0; ty < type; ty++){
          for (int t = 0; t < T;t++){
            for (int q = 0; q < Q; q++){
              fprintf(pFile014,"%.3f\n", beta[j][ty][t][q]);//Ã¨Â¿â„¢Ã©â€¡Å’Ã¥Â¾ÂªÃ§Å½Â¯Ã¥â€ â„¢Ã¥â€¦Â¥Ã¦â€“â€¡Ã¤Â»Â¶ 3Ã¤Â¸Âª 
            }
          }
        }
      }
      fclose(pFile014);
      
      
      FILE * pFile20;
      const char *pFileName20="Pest1.txt";
      pFile20 = fopen(pFileName20, "w");//Ã¨Â¿â„¢Ã¤Â¸ÂªÃ§â€Â¨Ã¢â‚¬Å“wÃ¢â‚¬ÂÃ¦ËœÂ¯Ã¥â€ â„¢Ã¦â€“â€¡Ã¤Â»Â¶Ã¯Â¼Å’Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Å½Å¸Ã¥â€ â€¦Ã¥Â®Â¹Ã¯Â¼Å’Ã¨â€¹Â¥Ã¤Â¸ÂÃ¦Æ’Â³Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Ë†â„¢Ã§â€Â¨Ã¢â‚¬Å“aÃ¢â‚¬Â
      if (NULL == pFile20)
      {
        printf("error");
        return 0;
      }
      
      
      for(int i=0;i<n;i++)
      {
        for(int t=0;t<t_num[i];t++)
        {
          for(int j=0;j<type;j++)
          {
            
            fprintf(pFile20,"%lf\n", p[i][j][tc[i][t]]);
            
          }
        }
      }
      
      fclose(pFile20);
      
      
      
      FILE * pFile21;
      const char *pFileName21="sigma_eps1.txt";
      pFile21 = fopen(pFileName21, "w");//Ã¨Â¿â„¢Ã¤Â¸ÂªÃ§â€Â¨Ã¢â‚¬Å“wÃ¢â‚¬ÂÃ¦ËœÂ¯Ã¥â€ â„¢Ã¦â€“â€¡Ã¤Â»Â¶Ã¯Â¼Å’Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Å½Å¸Ã¥â€ â€¦Ã¥Â®Â¹Ã¯Â¼Å’Ã¨â€¹Â¥Ã¤Â¸ÂÃ¦Æ’Â³Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Ë†â„¢Ã§â€Â¨Ã¢â‚¬Å“aÃ¢â‚¬Â
      if (NULL == pFile21)
      {
        printf("error");
        return 0;
      }
      
      
      for(int j=0;j<G;j++)
      {
        
        fprintf(pFile21,"%lf\n", sigma_eps[j]);
        
      }
      
      fclose(pFile21);
      
      FILE * pFile22;
      const char *pFileName22="sigma_eta1.txt";
      pFile22 = fopen(pFileName22, "w");//Ã¨Â¿â„¢Ã¤Â¸ÂªÃ§â€Â¨Ã¢â‚¬Å“wÃ¢â‚¬ÂÃ¦ËœÂ¯Ã¥â€ â„¢Ã¦â€“â€¡Ã¤Â»Â¶Ã¯Â¼Å’Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Å½Å¸Ã¥â€ â€¦Ã¥Â®Â¹Ã¯Â¼Å’Ã¨â€¹Â¥Ã¤Â¸ÂÃ¦Æ’Â³Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Ë†â„¢Ã§â€Â¨Ã¢â‚¬Å“aÃ¢â‚¬Â
      if (NULL == pFile22)
      {
        printf("error");
        return 0;
      }
      
      
      for(int j=0;j<G;j++)
      {
        
        for(int k=0;k<type;k++)
        {
          for(int ii1=0;ii1<m;ii1++)
          {
            for(int ii2=0;ii2<m;ii2++)
            {
              
              fprintf(pFile22,"%lf\n", sigma_eta[j][k][ii1][ii2]);
              
            }
            
          }
        }}
      fclose(pFile22);
      
      
      
      FILE * pFile23;
      const char *pFileName23="sigma_delta1.txt";
      pFile23 = fopen(pFileName23, "w");//Ã¨Â¿â„¢Ã¤Â¸ÂªÃ§â€Â¨Ã¢â‚¬Å“wÃ¢â‚¬ÂÃ¦ËœÂ¯Ã¥â€ â„¢Ã¦â€“â€¡Ã¤Â»Â¶Ã¯Â¼Å’Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Å½Å¸Ã¥â€ â€¦Ã¥Â®Â¹Ã¯Â¼Å’Ã¨â€¹Â¥Ã¤Â¸ÂÃ¦Æ’Â³Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Ë†â„¢Ã§â€Â¨Ã¢â‚¬Å“aÃ¢â‚¬Â
      if (NULL == pFile23)
      {
        printf("error");
        return 0;
      }
      
      
      for(int j=0;j<G;j++)
      {
        
        for(int k=0;k<type;k++)
        {
          
          fprintf(pFile23,"%lf\n", sigma_delta[j][k]);
          
        }
        
      }
      
      fclose(pFile23);
      
      
      
      ////E step
      
      for (int i=0;i<n;i++){
        for (int j=0;j<G;j++){
          
          
          for (int k=0;k<type;k++){
            for (int t=0;t<t_num[i];t++){
              mu[i][j][k][tc[i][t]]=0;
              for (int s=0;s<(S+1);s++){
                mu[i][j][k][tc[i][t]]=mu[i][j][k][tc[i][t]]+gamma[j][k][s]*y[i][tc[i][t]][s];
              }
              
              for (int l=0;l<Q;l++){
                mu[i][j][k][tc[i][t]]=mu[i][j][k][tc[i][t]]+beta[j][k][tc[i][t]][l]*x[i][tc[i][t]][l];
              }
              
            }  
          }
        }}
      
      
#pragma omp parallel for    
      for (int i=0;i<n;i++){
        for (int j=0;j<G;j++){
          
          
          gsl_vector *L_ij = gsl_vector_alloc(type*(t_num[i]+m));
          gsl_vector *Emu_ij = gsl_vector_alloc(type*(t_num[i]+m));
          gsl_matrix *Psi_ij = gsl_matrix_alloc(type*(t_num[i]+m),type*(t_num[i]+m));//v,w
          double** v_inverse=make2Darray(m,m);
          
          gsl_matrix_set_zero(Psi_ij);
          gsl_vector_set_zero(L_ij);
          gsl_vector_set_zero(Emu_ij);
          /////////////////L_ij
          for (int t=0;t<t_num[i];t++){
            
            for (int k=0;k<type;k++){
              
              value_vec[i]=((O[i][j][tc[i][t]]*p[i][k][tc[i][t]]/sigma_eps[j])+(mu[i][j][k][tc[i][t]]/sigma_delta[j][k]));
              //printf("%lf\n",p[i][k][tc[i][t]]);
              gsl_vector_set(L_ij, (t*type+k), value_vec[i]);
              
            }
            
            
          }
          
          
          //printf("%lf\n",gsl_vector_get(L_ij,1));
          
          
          for (int s=0;s<m;s++){
            
            for (int k=0;k<type;k++){
              value_vec[i]=0;
              for (int t=0;t<t_num[i];t++){
                value_vec[i]=value_vec[i]+mu[i][j][k][tc[i][t]]*z[i][tc[i][t]][s];
              }
              value_vec[i]=-(value_vec[i]/sigma_delta[j][k]);
              gsl_vector_set(L_ij, (type*t_num[i]+k*m+s), value_vec[i]);
            }
            
          }
          
          
          
          
          /////////////////Psi_ij
          
          //diag
          for (int t=0;t<t_num[i];t++){
            
            for (int k=0;k<type;k++){
              value_vec[i]=(p[i][k][tc[i][t]]*p[i][k][tc[i][t]]/sigma_eps[j]+1/sigma_delta[j][k]);
              gsl_matrix_set(Psi_ij, (t*type+k),(t*type+k),value_vec[i]);
              
            }
            
          }
          
          
          for (int t=0;t<t_num[i];t++){
            
            for (int k1=1;k1<type;k1++){
              
              for (int k2=0;k2<k1;k2++){
                value_vec[i]=(p[i][k1][tc[i][t]]*p[i][k2][tc[i][t]]/sigma_eps[j]);
                gsl_matrix_set(Psi_ij, t*type+k1,t*type+k2,value_vec[i]);
                gsl_matrix_set(Psi_ij, t*type+k2,t*type+k1,value_vec[i]);
              }
            }
            
          }
          
          
          
          
          for (int s=0;s<m;s++){
            for (int k=0;k<type;k++){
              
              for (int t=0;t<t_num[i];t++){
                value_vec[i]=-z[i][tc[i][t]][s]/sigma_delta[j][k];
                gsl_matrix_set(Psi_ij, t*type+k,t_num[i]*type+k*m+s,value_vec[i]);
                gsl_matrix_set(Psi_ij, t_num[i]*type+k*m+s, t*type+k,value_vec[i]);
              }
            }
          }  
          
          
          
          
          for (int k=0;k<type;k++){
            
            inverse(sigma_eta[j][k],m,v_inverse);
            
            //gsl_matrix_set(Psi_ij, t_num[i]*type+k*(1+S+Q), t_num[i]*type+k*(1+S+Q),v_inverse[0][0]+t_num[i]/sigma_delta[j][k]);
            
            for (int s=0;s<m;s++){
              
              for (int s1=0;s1<m;s1++){
                
                value_vec[i]=0;
                for (int t=0;t<t_num[i];t++){
                  value_vec[i]=value_vec[i]+z[i][tc[i][t]][s]*z[i][tc[i][t]][s1]/sigma_delta[j][k];
                }
                
                gsl_matrix_set(Psi_ij, t_num[i]*type+k*m+s1, t_num[i]*type+k*m+s,v_inverse[s][s1]+value_vec[i]);
                
                
              }
              
              
            }
            
            
            
            
            
          }
          
          
          //
          
          
          
          
          
          
          
          gsl_matrix_inv(Psi_ij);
          gsl_matrixvector_mul(Psi_ij,L_ij,Emu_ij,(t_num[i]+m)*type,(t_num[i]+m)*type);
          
          
          
          for (int t=0;t<t_num[i];t++){
            for (int k=0;k<type;k++){
              Emu[i][j][tc[i][t]*type+k]=gsl_vector_get(Emu_ij,t*type+k);
            }
          }
          
          for (int k=0;k<type;k++){
            for (int ii=0;ii<m;ii++){
              Emu[i][j][T*type+k*m+ii]=gsl_vector_get(Emu_ij,t_num[i]*type+k*m+ii);
            }
          }
          
          for (int t1=0;t1<t_num[i];t1++){
            for (int k1=0;k1<type;k1++){
              for (int t2=0;t2<t_num[i];t2++){
                for (int k2=0;k2<type;k2++){
                  Emu2[i][j][tc[i][t1]*type+k1][tc[i][t2]*type+k2]=Emu[i][j][tc[i][t1]*type+k1]*Emu[i][j][tc[i][t2]*type+k2]+gsl_matrix_get(Psi_ij,t1*type+k1,t2*type+k2);
                }
              }
            }
          }
          
          
          for (int k1=0;k1<type;k1++){
            for (int k2=0;k2<type;k2++){
              
              for (int s1=0;s1<m;s1++){
                for (int s2=0;s2<m;s2++){
                  Emu2[i][j][T*type+k1*m+s1][T*type+k2*m+s2]=Emu[i][j][T*type+k1*m+s1]*Emu[i][j][T*type+k2*m+s2]+gsl_matrix_get(Psi_ij,t_num[i]*type+k1*m+s1,t_num[i]*type+k2*m+s2);
                }
              }
              
              
            }
          }
          
          for (int k1=0;k1<type;k1++){
            for (int t2=0;t2<t_num[i];t2++){
              for (int k2=0;k2<type;k2++){
                for (int s=0;s<m;s++){
                  Emu2[i][j][T*type+k1*m+s][tc[i][t2]*type+k2]=Emu[i][j][T*type+k1*m+s]*Emu[i][j][tc[i][t2]*type+k2]+gsl_matrix_get(Psi_ij,t_num[i]*type+k1*m+s,t2*type+k2);
                  Emu2[i][j][tc[i][t2]*type+k2][T*type+k1*m+s]=Emu2[i][j][T*type+k1*m+s][tc[i][t2]*type+k2];
                }
                
              }
            }
          }
          
          gsl_matrix_free(Psi_ij);
          gsl_vector_free(Emu_ij);
          gsl_vector_free(L_ij);
          delet2Darray(v_inverse,m,m);
        }}      
      
      
      
      
#pragma omp parallel for   
      for (int j=0;j<G;j++){
        
        
        gsl_matrix *SOC = gsl_matrix_alloc((1+S+T*Q),(1+S+T*Q));
        gsl_matrix *FOC2 = gsl_matrix_alloc((1+S+T*Q),(1+S+T*Q));
        gsl_vector *FOC = gsl_vector_alloc(1+S+T*Q);
        gsl_matrix *cov_mat = gsl_matrix_alloc((1+S+T*Q),(1+S+T*Q));
        
        for (int k=0;k<type;k++){
          
          gsl_matrix_set_zero(SOC);
          gsl_vector_set_zero(FOC);
          gsl_matrix_set_zero(FOC2);
          gsl_matrix_set_zero(cov_mat);
          
          //FOC
          for (int s=0;s<(S+1);s++){
            
            value_vec_G[j]=0;
            
            for (int i=0;i<n;i++){
              for (int t=0;t<T;t++){
                value_vec_G[j]=value_vec_G[j]+(Emu[i][j][t*type+k]-mu[i][j][k][t])*y[i][t][s]/sigma_delta[j][k];
                
                for (int ss=0;ss<m;ss++){
                  
                  value_vec_G[j]=value_vec_G[j]-z[i][t][ss]*Emu[i][j][T*type+k*m+ss]*y[i][t][s]/sigma_delta[j][k];
                  
                }
                
                
              }
              
            }
            
            
            gsl_vector_set(FOC, s, value_vec_G[j]);
            
            
          }
          
          
          for (int q=0;q<Q;q++){
            
            for (int t=0;t<T;t++){
              value_vec_G[j]=0;
              
              for (int i=0;i<n;i++){
                
                value_vec_G[j]=value_vec_G[j]+(Emu[i][j][t*type+k]-mu[i][j][k][t])*x[i][t][q]/sigma_delta[j][k];
                
                
                for (int ss=0;ss<m;ss++){
                  
                  value_vec_G[j]=value_vec_G[j]-z[i][t][ss]*Emu[i][j][T*type+k*m+ss]*x[i][t][q]/sigma_delta[j][k];
                  
                }
                
                
                
              }
              
              gsl_vector_set(FOC, 1+S+t*Q+q, value_vec_G[j]);
            }
          } 
          
          
          
          
          //SOC
          
          //diag
          
          for (int s=0;s<(S+1);s++){
            
            value_vec_G[j]=0;
            
            for (int i=0;i<n;i++){
              for (int t=0;t<T;t++){
                value_vec_G[j]=value_vec_G[j]-y[i][t][s]*y[i][t][s]/sigma_delta[j][k];
                
              }
              
            }
            
            
            gsl_matrix_set(SOC, s, s,value_vec_G[j]);
            
            
          }
          
          
          for (int q=0;q<Q;q++){
            
            for (int t=0;t<T;t++){
              value_vec_G[j]=0;
              
              for (int i=0;i<n;i++){
                
                value_vec_G[j]=value_vec_G[j]-x[i][t][q]*x[i][t][q]/sigma_delta[j][k];
                
                
              }
              
              gsl_matrix_set(SOC, 1+S+t*Q+q, 1+S+t*Q+q, value_vec_G[j]);
            }
          } 
          
          
          
          
          //offdiag
          
          
          value_vec_G[j]=0;
          
          for (int s1=0;s1<(S+1);s1++){
            for (int s2=0;s2<(S+1);s2++){
              value_vec_G[j]=0;
              
              for (int i=0;i<n;i++){
                for (int t=0;t<T;t++){
                  value_vec_G[j]=value_vec_G[j]-y[i][t][s1]*y[i][t][s2]/sigma_delta[j][k];
                  
                }
                
              }
              
              
              gsl_matrix_set(SOC, s1, s2,value_vec_G[j]);
              gsl_matrix_set(SOC, s2, s1,value_vec_G[j]);
            }
          }
          
          
          for (int s1=0;s1<(S+1);s1++){
            for (int s2=0;s2<Q;s2++){
              for (int t=0;t<T;t++){
                value_vec_G[j]=0;
                
                for (int i=0;i<n;i++){
                  
                  value_vec_G[j]=value_vec_G[j]-y[i][t][s1]*x[i][t][s2]/sigma_delta[j][k];
                  
                }
                
                gsl_matrix_set(SOC, s1, 1+S+t*Q+s2,value_vec_G[j]);
                gsl_matrix_set(SOC, 1+S+t*Q+s2, s1,value_vec_G[j]);
              }
              
            }
          }
          
          
          
          
          
          
          for (int q1=0;q1<(Q);q1++){
            for (int q2=0;q2<Q;q2++){
              
              for (int t=0;t<T;t++){
                
                value_vec_G[j]=0;
                
                for (int i=0;i<n;i++){
                  
                  value_vec_G[j]=value_vec_G[j]-x[i][t][q1]*x[i][t][q2]/sigma_delta[j][k];
                  
                }
                
                
                gsl_matrix_set(SOC, 1+S+t*Q+q1, 1+S+t*Q+q2,value_vec_G[j]);
                gsl_matrix_set(SOC, 1+S+t*Q+q1, 1+S+t*Q+q2,value_vec_G[j]);
              }
            }
            
          }
          
          
          
          
          for (int s1=0;s1<(S+1);s1++){
            for (int s2=0;s2<(S+1);s2++){
              value_vec_G[j]=0;
              
              
              for (int i=0;i<n;i++){
                for (int t1=0;t1<T;t1++){
                  for (int t2=0;t2<T;t2++){
                    value_vec_G[j]=value_vec_G[j]+y[i][t1][s1]*y[i][t2][s2]*(Emu2[i][j][t1*type+k][t2*type+k]-mu[i][j][k][t2]*Emu[i][j][t1*type+k]-mu[i][j][k][t1]*Emu[i][j][t2*type+k]+mu[i][j][k][t1]*mu[i][j][k][t2]);
                    
                    for (int ss=0;ss<m;ss++){
                      value_vec_G[j]=value_vec_G[j]-y[i][t1][s1]*y[i][t2][s2]*(z[i][t1][ss]*Emu2[i][j][t2*type+k][T*type+k*m+ss]+z[i][t2][ss]*Emu2[i][j][t1*type+k][T*type+k*m+ss]-z[i][t1][ss]*mu[i][j][k][t2]*Emu[i][j][T*type+k*m+ss]-z[i][t2][ss]*mu[i][j][k][t1]*Emu[i][j][T*type+k*m+ss]);
                    }
                    
                    
                    
                    for (int ss1=0;ss1<m;ss1++){
                      for (int ss2=0;ss2<m;ss2++){
                        value_vec_G[j]=value_vec_G[j]+y[i][t1][s1]*y[i][t2][s2]*z[i][t1][ss1]*z[i][t2][ss2]*Emu2[i][j][T*type+k*m+ss1][T*type+k*m+ss2];
                      }
                    }
                    
                    
                    //-Emu2[i][j][t2*type+k][T*type+k]-Emu2[i][j][t1*type+k][T*type+k]+mu[i][j][k][t2]*Emu[i][j][T*type+k]+mu[i][j][k][t1]*Emu[i][j][T*type+k]+Emu2[i][j][T*type+k][T*type+k]
                    //printf("%lf\n",Emu2[i][j][t1*type+k][t2*type+k]-mu[i][j][k][t2]*Emu[i][j][t1*type+k]-mu[i][j][k][t1]*Emu[i][j][t2*type+k]+mu[i][j][k][t1]*mu[i][j][k][t2]-Emu2[i][j][t2*type+k][T*type+k]-Emu2[i][j][t1*type+k][T*type+k]+mu[i][j][k][t2]*Emu[i][j][T*type+k]+mu[i][j][k][t1]*Emu[i][j][T*type+k]+Emu2[i][j][T*type+k][T*type+k]);
                  }
                }
                
              }
              
              
              for (int i1=1;i1<n;i1++){
                for (int i2=0;i2<i1;i2++){
                  
                  
                  for (int t1=0;t1<T;t1++){
                    for (int t2=0;t2<T;t2++){
                      
                      value_temp1[j]=Emu[i1][j][t1*type+k]-mu[i1][j][k][t1];
                      
                      for (int ss=0;ss<m;ss++){
                        value_temp1[j]=value_temp1[j]-Emu[i1][j][T*type+k*m+ss]*z[i1][t1][ss];
                      }
                      
                      
                      
                      value_temp2[j]=Emu[i2][j][t2*type+k]-mu[i2][j][k][t2];
                      
                      for (int ss=0;ss<m;ss++){
                        value_temp2[j]=value_temp2[j]-Emu[i2][j][T*type+k*m+ss]*z[i2][t2][ss];
                      }
                      
                      
                      value_vec_G[j]=value_vec_G[j]+y[i1][t1][s1]*y[i2][t2][s2]*value_temp1[j]*value_temp2[j];
                      
                      
                    }}
                  
                }
              }
              
              
              for (int i2=1;i2<n;i2++){
                for (int i1=0;i1<i2;i1++){
                  
                  
                  for (int t1=0;t1<T;t1++){
                    for (int t2=0;t2<T;t2++){
                      
                      value_temp1[j]=Emu[i1][j][t1*type+k]-mu[i1][j][k][t1];
                      
                      for (int ss=0;ss<m;ss++){
                        value_temp1[j]=value_temp1[j]-Emu[i1][j][T*type+k*m+ss]*z[i1][t1][ss];
                      }
                      
                      
                      
                      value_temp2[j]=Emu[i2][j][t2*type+k]-mu[i2][j][k][t2];
                      
                      for (int ss=0;ss<m;ss++){
                        value_temp2[j]=value_temp2[j]-Emu[i2][j][T*type+k*m+ss]*z[i2][t2][ss];
                      }
                      
                      
                      value_vec_G[j]=value_vec_G[j]+y[i1][t1][s1]*y[i2][t2][s2]*value_temp1[j]*value_temp2[j];
                      
                      
                    }}
                  
                }
              }
              
              
              
              
              value_vec_G[j]=value_vec_G[j]/sigma_delta[j][k];
              value_vec_G[j]=value_vec_G[j]/sigma_delta[j][k];
              gsl_matrix_set(FOC2, s1, s2,value_vec_G[j]);
              
            }
          }
          
          
          
          
          for (int s1=0;s1<(S+1);s1++){
            for (int s2=0;s2<Q;s2++){
              for (int t=0;t<T;t++){
                value_vec_G[j]=0;
                
                for (int i=0;i<n;i++){
                  for (int t1=0;t1<T;t1++){
                    value_vec_G[j]=value_vec_G[j]+y[i][t1][s1]*x[i][t][s2]*(Emu2[i][j][t1*type+k][t*type+k]-mu[i][j][k][t]*Emu[i][j][t1*type+k]-mu[i][j][k][t1]*Emu[i][j][t*type+k]+mu[i][j][k][t1]*mu[i][j][k][t]);
                    
                    for (int ss=0;ss<m;ss++){
                      value_vec_G[j]=value_vec_G[j]-y[i][t1][s1]*x[i][t][s2]*(z[i][t1][ss]*Emu2[i][j][t*type+k][T*type+k*m+ss]+z[i][t][ss]*Emu2[i][j][t1*type+k][T*type+k*m+ss]-z[i][t1][ss]*mu[i][j][k][t]*Emu[i][j][T*type+k*m+ss]-z[i][t][ss]*mu[i][j][k][t1]*Emu[i][j][T*type+k*m+ss]);
                    }
                    
                    
                    
                    for (int ss1=0;ss1<m;ss1++){
                      for (int ss2=0;ss2<m;ss2++){
                        value_vec_G[j]=value_vec_G[j]+y[i][t1][s1]*x[i][t][s2]*z[i][t1][ss1]*z[i][t][ss2]*Emu2[i][j][T*type+k*m+ss1][T*type+k*m+ss2];
                      }
                    }
                    
                    
                    
                    
                  }
                }
                
                
                for (int i1=0;i1<(n-1);i1++){
                  for (int i2=(i1+1);i2<n;i2++){
                    
                    
                    for (int t1=0;t1<T;t1++){
                      
                      
                      value_temp1[j]=Emu[i1][j][t1*type+k]-mu[i1][j][k][t1];
                      
                      for (int ss=0;ss<m;ss++){
                        value_temp1[j]=value_temp1[j]-Emu[i1][j][T*type+k*m+ss]*z[i1][t1][ss];
                      }
                      
                      
                      value_temp2[j]=Emu[i2][j][t*type+k]-mu[i2][j][k][t];
                      
                      for (int ss=0;ss<m;ss++){
                        value_temp2[j]=value_temp2[j]-Emu[i2][j][T*type+k*m+ss]*z[i2][t][ss];
                      }
                      
                      
                      value_vec_G[j]=value_vec_G[j]+y[i1][t1][s1]*x[i2][t][s2]*value_temp1[j]*value_temp2[j];
                      
                      
                    }
                    
                  }
                }
                
                for (int i2=0;i2<(n-1);i2++){
                  for (int i1=(i2+1);i1<n;i1++){
                    
                    
                    for (int t1=0;t1<T;t1++){
                      
                      
                      value_temp1[j]=Emu[i1][j][t1*type+k]-mu[i1][j][k][t1];
                      
                      for (int ss=0;ss<m;ss++){
                        value_temp1[j]=value_temp1[j]-Emu[i1][j][T*type+k*m+ss]*z[i1][t1][ss];
                      }
                      
                      
                      value_temp2[j]=Emu[i2][j][t*type+k]-mu[i2][j][k][t];
                      
                      for (int ss=0;ss<m;ss++){
                        value_temp2[j]=value_temp2[j]-Emu[i2][j][T*type+k*m+ss]*z[i2][t][ss];
                      }
                      
                      
                      value_vec_G[j]=value_vec_G[j]+y[i1][t1][s1]*x[i2][t][s2]*value_temp1[j]*value_temp2[j];
                      
                      
                    }
                    
                  }
                }
                
                value_vec_G[j]=value_vec_G[j]/sigma_delta[j][k];
                value_vec_G[j]=value_vec_G[j]/sigma_delta[j][k];
                
                gsl_matrix_set(FOC2, s1, 1+S+t*Q+s2,value_vec_G[j]);
                gsl_matrix_set(FOC2, 1+S+t*Q+s2, s1,value_vec_G[j]);
              }
              
            }
          }
          
          
          for (int q1=0;q1<Q;q1++){
            for (int q2=0;q2<Q;q2++){
              for (int t1=0;t1<T;t1++){
                for (int t2=0;t2<T;t2++){
                  value_vec_G[j]=0;
                  
                  for (int i=0;i<n;i++){
                    
                    value_vec_G[j]=value_vec_G[j]+x[i][t1][q1]*x[i][t2][q2]*(Emu2[i][j][t1*type+k][t2*type+k]-mu[i][j][k][t2]*Emu[i][j][t1*type+k]-mu[i][j][k][t1]*Emu[i][j][t2*type+k]+mu[i][j][k][t1]*mu[i][j][k][t2]);
                    
                    
                    for (int ss=0;ss<m;ss++){
                      value_vec_G[j]=value_vec_G[j]-x[i][t1][q1]*x[i][t2][q2]*(z[i][t1][ss]*Emu2[i][j][t2*type+k][T*type+k*m+ss]+z[i][t2][ss]*Emu2[i][j][t1*type+k][T*type+k*m+ss]-z[i][t1][ss]*mu[i][j][k][t2]*Emu[i][j][T*type+k*m+ss]-z[i][t2][ss]*mu[i][j][k][t1]*Emu[i][j][T*type+k*m+ss]);
                    }
                    
                    
                    
                    
                    for (int ss1=0;ss1<m;ss1++){
                      for (int ss2=0;ss2<m;ss2++){
                        value_vec_G[j]=value_vec_G[j]+x[i][t1][q1]*x[i][t2][q2]*z[i][t1][ss1]*z[i][t2][ss2]*Emu2[i][j][T*type+k*m+ss1][T*type+k*m+ss2];
                      }
                    }
                    
                    
                    
                    
                    //-Emu2[i][j][t2*type+k][T*type+k]-Emu2[i][j][t1*type+k][T*type+k]+mu[i][j][k][t2]*Emu[i][j][T*type+k]+mu[i][j][k][t1]*Emu[i][j][T*type+k]+Emu2[i][j][T*type+k][T*type+k]
                    
                  }
                  
                  
                  for (int i1=0;i1<(n-1);i1++){
                    for (int i2=(i1+1);i2<n;i2++){
                      
                      
                      value_temp1[j]=Emu[i1][j][t1*type+k]-mu[i1][j][k][t1];
                      
                      for (int ss=0;ss<m;ss++){
                        value_temp1[j]=value_temp1[j]-Emu[i1][j][T*type+k*m+ss]*z[i1][t1][ss];
                      }
                      
                      
                      
                      value_temp2[j]=Emu[i2][j][t2*type+k]-mu[i2][j][k][t2];
                      
                      for (int ss=0;ss<m;ss++){
                        value_temp2[j]=value_temp2[j]-Emu[i2][j][T*type+k*m+ss]*z[i2][t2][ss];
                      }
                      
                      
                      
                      value_vec_G[j]=value_vec_G[j]+x[i1][t1][q1]*x[i2][t2][q2]*value_temp1[j]*value_temp2[j];
                      
                      
                      
                      
                    }
                  }
                  
                  
                  
                  for (int i2=0;i2<(n-1);i2++){
                    for (int i1=(i2+1);i1<n;i1++){
                      
                      
                      
                      
                      
                      value_temp1[j]=Emu[i1][j][t1*type+k]-mu[i1][j][k][t1];
                      
                      for (int ss=0;ss<m;ss++){
                        value_temp1[j]=value_temp1[j]-Emu[i1][j][T*type+k*m+ss]*z[i1][t1][ss];
                      }
                      
                      
                      value_temp2[j]=Emu[i2][j][t2*type+k]-mu[i2][j][k][t2];
                      
                      for (int ss=0;ss<m;ss++){
                        value_temp2[j]=value_temp2[j]-Emu[i2][j][T*type+k*m+ss]*z[i2][t2][ss];
                      }
                      
                      
                      
                      value_vec_G[j]=value_vec_G[j]+x[i1][t1][q1]*x[i2][t2][q2]*value_temp1[j]*value_temp2[j];
                      
                      
                      
                      
                    }
                  }
                  
                  value_vec_G[j]=value_vec_G[j]/sigma_delta[j][k];
                  value_vec_G[j]=value_vec_G[j]/sigma_delta[j][k];
                  
                  gsl_matrix_set(FOC2, 1+S+t1*Q+q1, 1+S+t2*Q+q2, value_vec_G[j]);
                }
              } 
            } }
          
          
          
          
          
          for (int i1=0;i1<(1+S+T*Q);i1++){
            for (int i2=0;i2<(1+S+T*Q);i2++){
              
              gsl_matrix_set(cov_mat, i1, i2,(-gsl_matrix_get(SOC, i1, i2)-gsl_matrix_get(FOC2, i1, i2)+gsl_vector_get(FOC, i1)*gsl_vector_get(FOC, i2)));
              
            }
          }
          
          
          for (int i2=0;i2<(1+S+T*Q);i2++){
            
            gsl_matrix_set(cov_mat, i2, i2,(-gsl_matrix_get(SOC, i2, i2)-gsl_matrix_get(FOC2, i2, i2)+gsl_vector_get(FOC, i2)*gsl_vector_get(FOC, i2)));
            
          }
          
          
          
          //   for (int i1=0;i1<(1+S);i1++){
          
          
          //   printf("%lf\n",gsl_matrix_get(cov_mat,i1,i1));
          
          //  }
          
          
          gsl_matrix_inv(cov_mat);
          
          
          
          for (int s=0;s<(1+S);s++){
            
            gamma_var[j][k][s]=gsl_matrix_get(cov_mat,s,s);
            gamma_z[j][k][s]=gamma[j][k][s]/sqrt(gamma_var[j][k][s]);
            
          }
          
          
          for (int q=0;q<Q;q++){
            
            for (int t=0;t<T;t++){
              beta_var[j][k][t][q]=gsl_matrix_get(cov_mat,1+S+t*Q+q,1+S+t*Q+q);
              beta_z[j][k][t][q]=beta[j][k][t][q]/sqrt(beta_var[j][k][t][q]);
            }
          }
        }
        
        gsl_matrix_free(cov_mat); 
        gsl_matrix_free(SOC); 
        gsl_matrix_free(FOC2); 
        gsl_vector_free(FOC); 
        
      }
      
      FILE * pFile14;
      const char *pFileName14="beta1.txt";
      pFile14 = fopen(pFileName14, "w");//Ã¨Â¿â„¢Ã¤Â¸ÂªÃ§â€Â¨Ã¢â‚¬Å“wÃ¢â‚¬ÂÃ¦ËœÂ¯Ã¥â€ â„¢Ã¦â€“â€¡Ã¤Â»Â¶Ã¯Â¼Å’Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Å½Å¸Ã¥â€ â€¦Ã¥Â®Â¹Ã¯Â¼Å’Ã¨â€¹Â¥Ã¤Â¸ÂÃ¦Æ’Â³Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Ë†â„¢Ã§â€Â¨Ã¢â‚¬Å“aÃ¢â‚¬Â
      if (NULL == pFile14)
      {
        printf("error");
        return 0;
      }
      for (int j = 0; j < G; j++){
        
        for (int ty = 0; ty < type; ty++){
          for (int t = 0; t < T;t++){
            for (int q = 0; q < Q; q++){
              fprintf(pFile14,"%.3f\n", beta[j][ty][t][q]);//Ã¨Â¿â„¢Ã©â€¡Å’Ã¥Â¾ÂªÃ§Å½Â¯Ã¥â€ â„¢Ã¥â€¦Â¥Ã¦â€“â€¡Ã¤Â»Â¶ 3Ã¤Â¸Âª 
            }
          }
        }
      }
      fclose(pFile14);
      
      
      FILE * pFile15;
      const char *pFileName15="gamma1.txt";
      pFile15 = fopen(pFileName15, "w");//Ã¨Â¿â„¢Ã¤Â¸ÂªÃ§â€Â¨Ã¢â‚¬Å“wÃ¢â‚¬ÂÃ¦ËœÂ¯Ã¥â€ â„¢Ã¦â€“â€¡Ã¤Â»Â¶Ã¯Â¼Å’Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Å½Å¸Ã¥â€ â€¦Ã¥Â®Â¹Ã¯Â¼Å’Ã¨â€¹Â¥Ã¤Â¸ÂÃ¦Æ’Â³Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Ë†â„¢Ã§â€Â¨Ã¢â‚¬Å“aÃ¢â‚¬Â
      if (NULL == pFile15)
      {
        printf("error");
        return 0;
      }
      
      for (int g = 0; g < G; g++){
        for (int ty = 0; ty < type;ty++){
          for (int s = 0; s < (1+S); s++){
            fprintf(pFile15,"%.3f\n", gamma[g][ty][s]);//Ã¨Â¿â„¢Ã©â€¡Å’Ã¥Â¾ÂªÃ§Å½Â¯Ã¥â€ â„¢Ã¥â€¦Â¥Ã¦â€“â€¡Ã¤Â»Â¶ 3Ã¤Â¸Âª 
          }
        }
      }
      
      fclose(pFile15);
      
      
      
      FILE * pFile16;
      const char *pFileName16="beta_var1.txt";
      pFile16 = fopen(pFileName16, "w");//Ã¨Â¿â„¢Ã¤Â¸ÂªÃ§â€Â¨Ã¢â‚¬Å“wÃ¢â‚¬ÂÃ¦ËœÂ¯Ã¥â€ â„¢Ã¦â€“â€¡Ã¤Â»Â¶Ã¯Â¼Å’Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Å½Å¸Ã¥â€ â€¦Ã¥Â®Â¹Ã¯Â¼Å’Ã¨â€¹Â¥Ã¤Â¸ÂÃ¦Æ’Â³Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Ë†â„¢Ã§â€Â¨Ã¢â‚¬Å“aÃ¢â‚¬Â
      if (NULL == pFile16)
      {
        printf("error");
        return 0;
      }
      
      
      for (int j = 0; j < G;j++){
        for (int k = 0; k < type; k++){
          for (int t = 0; t < T; t++){
            for (int q = 0; q < Q; q++){
              fprintf(pFile16,"%.3f\n", beta_var[j][k][t][q]);//Ã¨Â¿â„¢Ã©â€¡Å’Ã¥Â¾ÂªÃ§Å½Â¯Ã¥â€ â„¢Ã¥â€¦Â¥Ã¦â€“â€¡Ã¤Â»Â¶ 3Ã¤Â¸Âª 
            }
          }
        }
      }
      fclose(pFile16);
      
      
      FILE * pFile18;
      const char *pFileName18="beta_z2.txt";
      pFile18 = fopen(pFileName18, "w");//Ã¨Â¿â„¢Ã¤Â¸ÂªÃ§â€Â¨Ã¢â‚¬Å“wÃ¢â‚¬ÂÃ¦ËœÂ¯Ã¥â€ â„¢Ã¦â€“â€¡Ã¤Â»Â¶Ã¯Â¼Å’Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Å½Å¸Ã¥â€ â€¦Ã¥Â®Â¹Ã¯Â¼Å’Ã¨â€¹Â¥Ã¤Â¸ÂÃ¦Æ’Â³Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Ë†â„¢Ã§â€Â¨Ã¢â‚¬Å“aÃ¢â‚¬Â
      if (NULL == pFile18)
      {
        printf("error");
        return 0;
      }
      
      
      for (int j = 0; j < G;j++){
        for (int k = 0; k < type; k++){
          for (int t = 0; t < T; t++){
            for (int q = 0; q < Q; q++){
              fprintf(pFile18,"%.3f\n", beta_z[j][k][t][q]);//Ã¨Â¿â„¢Ã©â€¡Å’Ã¥Â¾ÂªÃ§Å½Â¯Ã¥â€ â„¢Ã¥â€¦Â¥Ã¦â€“â€¡Ã¤Â»Â¶ 3Ã¤Â¸Âª 
            }
          }
        }
      }
      fclose(pFile18);
      
      
      
      FILE * pFile17;
      const char *pFileName17="gamma_var1.txt";
      pFile17 = fopen(pFileName17, "w");//Ã¨Â¿â„¢Ã¤Â¸ÂªÃ§â€Â¨Ã¢â‚¬Å“wÃ¢â‚¬ÂÃ¦ËœÂ¯Ã¥â€ â„¢Ã¦â€“â€¡Ã¤Â»Â¶Ã¯Â¼Å’Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Å½Å¸Ã¥â€ â€¦Ã¥Â®Â¹Ã¯Â¼Å’Ã¨â€¹Â¥Ã¤Â¸ÂÃ¦Æ’Â³Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Ë†â„¢Ã§â€Â¨Ã¢â‚¬Å“aÃ¢â‚¬Â
      if (NULL == pFile17)
      {
        printf("error");
        return 0;
      }
      
      for (int j = 0; j < G;j++){
        for (int k = 0; k < type; k++){
          for (int s = 0; s < (1+S); s++){
            fprintf(pFile17,"%.3f\n", gamma_var[j][k][s]);//Ã¨Â¿â„¢Ã©â€¡Å’Ã¥Â¾ÂªÃ§Å½Â¯Ã¥â€ â„¢Ã¥â€¦Â¥Ã¦â€“â€¡Ã¤Â»Â¶ 3Ã¤Â¸Âª 
          }
        }
      }
      
      fclose(pFile17);
      
      
      FILE * pFile19;
      const char *pFileName19="gamma_z2.txt";
      pFile19 = fopen(pFileName19, "w");//Ã¨Â¿â„¢Ã¤Â¸ÂªÃ§â€Â¨Ã¢â‚¬Å“wÃ¢â‚¬ÂÃ¦ËœÂ¯Ã¥â€ â„¢Ã¦â€“â€¡Ã¤Â»Â¶Ã¯Â¼Å’Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Å½Å¸Ã¥â€ â€¦Ã¥Â®Â¹Ã¯Â¼Å’Ã¨â€¹Â¥Ã¤Â¸ÂÃ¦Æ’Â³Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Ë†â„¢Ã§â€Â¨Ã¢â‚¬Å“aÃ¢â‚¬Â
      if (NULL == pFile19)
      {
        printf("error");
        return 0;
      }
      
      for (int j = 0; j < G;j++){
        for (int k = 0; k < type; k++){
          for (int s = 0; s < (1+S); s++){
            fprintf(pFile19,"%.3f\n", gamma_z[j][k][s]);//Ã¨Â¿â„¢Ã©â€¡Å’Ã¥Â¾ÂªÃ§Å½Â¯Ã¥â€ â„¢Ã¥â€¦Â¥Ã¦â€“â€¡Ã¤Â»Â¶ 3Ã¤Â¸Âª 
          }
        }
      }
      
      fclose(pFile19);
      
      
      FILE * pFile190;
      const char *pFileName190="profile.txt";
      pFile190 = fopen(pFileName190, "w");//Ã¨Â¿â„¢Ã¤Â¸ÂªÃ§â€Â¨Ã¢â‚¬Å“wÃ¢â‚¬ÂÃ¦ËœÂ¯Ã¥â€ â„¢Ã¦â€“â€¡Ã¤Â»Â¶Ã¯Â¼Å’Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Å½Å¸Ã¥â€ â€¦Ã¥Â®Â¹Ã¯Â¼Å’Ã¨â€¹Â¥Ã¤Â¸ÂÃ¦Æ’Â³Ã¨Â¦â€ Ã§â€ºâ€“Ã¥Ë†â„¢Ã§â€Â¨Ã¢â‚¬Å“aÃ¢â‚¬Â
      if (NULL == pFile190)
      {
        printf("error");
        return 0;
      }
      for (int i = 0; i < n;i++){
        for (int j = 0; j < G;j++){
          for (int k = 0; k < (T)*type; k++){
            fprintf(pFile190,"%.3f\n", Emu[i][j][k]);//Ã¨Â¿â„¢Ã©â€¡Å’Ã¥Â¾ÂªÃ§Å½Â¯Ã¥â€ â„¢Ã¥â€¦Â¥Ã¦â€“â€¡Ã¤Â»Â¶ 3Ã¤Â¸Âª 
          }
        }
      }
      
      fclose(pFile190);
      
    }
    
  }
  
}

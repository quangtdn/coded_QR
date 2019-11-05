#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <assert.h>
#include <math.h>
//#include "mkl.h"

void generate_matrix_1(double* matrix, int dim1, int dim2, int rank){
	int i, j;
	for (i = 0; i < dim1; i++){
		for(j = 0; j < dim2; j++){
			matrix[i*dim2 + j] = i+j+rank+2 ;
		}
	}
}

void transpose_matrix(double* original_mat, double* trans_mat, int dim1, int dim2){
  int i, j;
	for (i = 0; i < dim1; i++){
		for(j = 0; j < dim2; j++){
			trans_mat[i*dim2 + j] = original_mat[j*dim2 + i];
		}
	}
}

void matrix_multiply(double* matA, double* matB, double* res, int m, int n, int p){
  //matA is m x n, mat B is n x p
  for(int i = 0; i < m; i++){
    for(int j = 0; j < p; j++){
      res[i*p + j] = 0;
      for(int k = 0; k < n; k++){
        res[i*p + j] += matA[i*n + k] * matB[k*p + j];
      }
    }
  }

}


void print_matrix_d(double* matrix, int dim1, int dim2, int rank, int xrank, int yrank){
	int i, j;
  printf("From proc %d, or proc[%d][%d] \n", rank, xrank, yrank);

	for (i = 0; i < dim1; i++){
		for(j =0; j<dim2; j++){
			printf("%.3f ", matrix[i*dim2 + j]);
		}
		printf("\n");
	}
	printf("\n");
}

void print_matrix(double* matrix, int dim1, int dim2){
	int i, j;

	for (i = 0; i < dim1; i++){
		for(j =0; j<dim2; j++){
			printf("%.3f ", matrix[i*dim2 + j]);
		}
		printf("\n");
	}
	printf("\n");
}

void assign_matrix(double* mat_for_replica, double* empty_mat, int m, int n){
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      empty_mat[i * n + j] = mat_for_replica[i * n + j];
    }
  }
}
void identity_matrix(double* matrix, int dim){
  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      if(i == j){
        matrix[i * dim + j] = 1;
      }else{
        matrix[i * dim + j] = 0;
      }
    }
  }
}



void testing_QR_result(double* matA, double* tau, int m, int n){
  //the input matA here has size m x n and is the output
  //from the QR algorithm
  double* matR = (double *) malloc(m*n*sizeof(double));
  double* matQ = (double *) malloc(m*m*sizeof(double));

  int i, j;
  //wriing R
  for(i = 0; i < m; i++){
    for(j = 0; j < n; j++){
      if(i < j || i == j){
        matR[i * n + j] = matA[i * n + j];
      } else{
        matR[i * n + j] = 0;
      }
    }
  }


  //Computing Q matrix from the Y vector
  double* temp_matrix = (double *) malloc(m*m*sizeof(double));
  identity_matrix(matQ, m);
  identity_matrix(temp_matrix, m);
  for(int iter = 0; iter < n; iter++){
    double* y_col = (double *) malloc(m*sizeof(double));
    double* y_col_transpose = (double *) malloc(m*sizeof(double));
    for(i = 0; i < m; i++){
      //writing y_col and its transpose- y_col_transpose
      if(i < iter){
        y_col[i] = 0;
        y_col_transpose[i] = 0;
      }else if(i == iter){
        y_col[i] = 1;
        y_col_transpose[i] = -tau[iter];
      }else{
        y_col[i] = matA[i * n + iter];
        y_col_transpose[i] = -matA[i * n + iter]*tau[iter];
      }
    }
    /*
    printf("Heii this is tau[%d]=%.3f \n", iter, tau[iter]);
    printf("And the y vector is: \n");
    print_matrix(y_col, m, 1);
    printf("And the -tau[i]*y^T vector is: \n");
    print_matrix(y_col_transpose, 1, m);
    */

    double* vec_product = (double *) malloc(m*m*sizeof(double));
    //vec_product = -tau[i] * y_col * y_col^T
    matrix_multiply(y_col, y_col_transpose, vec_product, m, 1, m);

    //Now we do vec_product = I -tau[i] * y_col * y_col^T
    for(int dim1 = 0; dim1 < m; dim1++){
      for(int dim2 = 0; dim2 < m; dim2++){
        if(dim1 == dim2){
          vec_product[dim1 * m + dim2] += 1;
        }
      }
    }


    //We are now ready to update Q
    matrix_multiply(temp_matrix, vec_product, matQ, m, m, m);
    assign_matrix(matQ, temp_matrix, m, m);

  }



  printf("The Q matrix is:\n");
  print_matrix(matQ, m, m);
  printf("The R matrix is:\n");
  print_matrix(matR, m, n);

  double* Q_transpose = (double *) malloc(m*m*sizeof(double));
  transpose_matrix(matQ, Q_transpose, m, m);
  double* QQT = (double *) malloc(m*m*sizeof(double));
  matrix_multiply(matQ, Q_transpose, QQT, m, m, m);
  printf("Checking I = QQ^T = Q^T Q :\n QQ^T = \n");
  print_matrix(QQT, m, m);


  double* reconstruct_A = (double *) malloc(m*n*sizeof(double));
  matrix_multiply(matQ, matR, reconstruct_A, m, m, n);
  printf("The Q*R matrix is:\n");
  print_matrix(reconstruct_A, m, n);



}





int main(int argc, char** argv) {

  // implementations might need the arguments.
  MPI_Init(&argc, NULL);
/*
  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Get the name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  // Print off a hello world message
  printf("Hello world from processor %s, rank %d out of %d processors\n",
   processor_name, world_rank, world_size);
   */

  //There are r*c processors in total
  //int r = 1, c = 3;
  int debug_mode = 1;
  //each block has size b*b
  //int b = 1;
	int r, c, b;
  //default value f_r = f_c = 1.
  //matrix height = f_r * r, matrix width = f_c * c
  int f_r =1 , f_c = 1;

  int master_node = 0;

  /* For iterators */
  int i, j, k, t;




  /* Get the rank and size in the original communicator */
  int world_rank, world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);


	if(argc < 4){
		if(world_rank == master_node){
			printf("Input: P_r P_c block_size\n Note that P_r*P_c = no. of processors");
		}
	} else {
		r = atoi(argv[1]);
		c = atoi(argv[2]);
		b = atoi(argv[3]);
	}


  int color_x = world_rank / c;

//	printf("2/1 = %d and color_x = %d, r=%d, c = %d\n", 2/1, color_x, r, c);
  int color_y = world_rank % c + 999 ;// + r*(world_rank/(r*c));
  // int color_z = world_rank % (n*n);

  MPI_Comm x_comm;
  MPI_Comm y_comm;
  // MPI_Comm z_comm;

  MPI_Comm_split(MPI_COMM_WORLD, color_x, world_rank, &x_comm);
  MPI_Comm_split(MPI_COMM_WORLD, color_y, world_rank, &y_comm);
  // MPI_Comm_split(MPI_COMM_WORLD, color_z, world_rank, &z_comm);

  //x_comm: for vertical communication, i.e. comm on a column
  int x_rank, x_size;
  MPI_Comm_rank(x_comm, &x_rank);
  MPI_Comm_size(x_comm, &x_size);


  //y_comm: for horizontal communication, i.e. comm on a row
  int y_rank, y_size;
  MPI_Comm_rank(y_comm, &y_rank);
  MPI_Comm_size(y_comm, &y_size);

	int px = y_rank, py = x_rank;
	//processor has rank (px, py) in 2D coordinate

  // int z_rank, z_size;
  // MPI_Comm_rank(z_comm, &z_rank);
  // MPI_Comm_size(z_comm, &z_size);


//Testing purposes (on vertical vs hozirontal communication)
/*
  int broad = 100;
  if(world_rank ==0) {
    broad= 229;
  }

  MPI_Barrier(x_comm);

  MPI_Bcast(&broad, 1 , MPI_BYTE, 0, x_comm);
  MPI_Barrier(x_comm);

  printf("broadcast: %d Hello world from processor wolrd_rank %d / %d, px = %d / %d, py = %d / %d \n " ,broad, world_rank, world_size, px, x_size, py, y_size);

*/


	printf("world_rank = %d, r= %d, c=%d, color_x = %d, color_y = %d, px= %d, py = %d \n", world_rank,r, c, color_x, color_y, px, py);

  /*
  The input matrix A_global is distributed among the processors, each of which
  stores A_local comprising parts of A_global.
  */

  int localRow = b * f_r;
  int localColumn = b * f_c;

  double *A_local = (double *) malloc(b*b*f_r*f_c*sizeof(double));
  generate_matrix_1(A_local, b*f_r, b*f_c, world_rank);

  //double A_local[9] = {12, -51, 4,  6, 167, -68,  -4, 24, -41};
  //double A_local[9] = {1, 0, 0,   0, 1, 0,   0, 0, 1};
  print_matrix_d(A_local, b*f_r, b*f_c, world_rank, px, py);


  double *tau = (double *) malloc(b*sizeof(double));

  //This is the case of block-row distribution
  if(c == 1 && f_c == 1){
    for(i = 0; i < b; i++){
      //For this case, A_global[i, i] is always held by Processor 0
      double alpha, beta, sum_square = 0;
      if(px == 0 && py == 0){
        alpha = A_local[i * localColumn + i];
        for(j = i; j < localRow; j++){
          sum_square += A_local[j * localColumn + i] * A_local[j * localColumn + i];
          //printf("Hei sum_square on proc 0 = %.3f", sum_square);
        }
      } else {
        alpha = 0;
        for(j = 0; j < localRow; j++){
          sum_square += A_local[j * localColumn + i] * A_local[j * localColumn + i];
        }
      }

      double local_square[2] = {alpha, sum_square};
      double *global_square = (double *) malloc(2*sizeof(double));
      //MPI_Allreduce(local_square, global_square, 2,MPI_DOUBLE, MPI_SUM, x_comm);
			MPI_Allreduce(local_square, global_square, 2,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      alpha = global_square[0];
      //printf("heii alpha is: %.3f\n", alpha);
      beta = sqrt(global_square[1]);

      if(alpha > 0){
        beta = -beta;
      }
      //check carefully formula for tau --doubt that 2(beta-alpha)^2/beta
      tau[i] = (beta - alpha) / (beta);

      //Using alpha and beta, we now do panel update in parallel
      double *tempVector = (double *) malloc((localColumn - i)*sizeof(double));
      if(px == 0 && py == 0){
        A_local[i * localColumn + i] = beta;

        for(j = i+1; j < localRow; j++){
          A_local[j * localColumn + i] = A_local[j * localColumn + i]/(1.0*(alpha-beta));
        }

        //now prepare data for allreduce to compute v[i:m,i]^T A_global[i:m, i+1:b]
        // this step is part of updating trailing matrix
        //double *tempVector = (double *) malloc((localColumn - i)*sizeof(double));
        //check this malloc
        for(k = i+1; k < localColumn; k++){
          tempVector[k-i-1] = 0;
          for(j = i; j < localRow; j++){
            //v[j,i]^T * A_global[j, k]
            if(j == i) {
              //note that first index of v is normalized to 1
              tempVector[k-i-1] += A_local[j * localColumn + k];
            }else{
              tempVector[k-i-1] += A_local[j * localColumn + i] * A_local[j * localColumn + k];
            }

          }
        }


      } else {
        //the case for other processors != proc[0, 0]

        for(j = 0; j < localRow; j++){
          A_local[j * localColumn + i] = A_local[j * localColumn + i]/(1.0*(alpha-beta));
        }


        //now prepare data for allreduce to compute v[i:m,i]^T *A_global[i:m, i+1:b]
        // this step is part of updating trailing matrix
        //check this malloc
        for(k = i+1; k < localColumn; k++){
          tempVector[k-i-1] = 0;
          for(j = 0; j < localRow; j++){
            //v[j,i]^T * A_global[j, k]
            tempVector[k-i-1] += A_local[j * localColumn + i] * A_local[j * localColumn + k];
          }
        }

      }

      /*//debugging purposes
      printf("hei, tau[%d] is: %.3f, here alpha= %.3f,  beta=  %.3f  \n", i, tau[i], alpha, beta);
      printf("And the y vector is: \n");

      print_matrix(A_local, localRow, localColumn);
      */


      //Now perform all-reduce to compute z=tau[i] *v[i:m,i]^T *A_global[i:m, i+1:b]
      double *tempVector_received = (double *) malloc((localColumn - i)*sizeof(double));
      //MPI_Allreduce(tempVector, tempVector_received, localColumn - i,MPI_DOUBLE, MPI_SUM, x_comm);
			MPI_Allreduce(tempVector, tempVector_received, localColumn - i,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      //note that now the variable z = tau[i] *tempVector_received

      //The final step is to update the trailing print_matrix
      //by doing A_global[i:m, i+1:b] = A_global[i:m, i+1:b] - v[i:m, i] * z
      if(px == 0 && py == 0){
        for(j = i; j < localRow; j++){
          for(k = i+1; k < localColumn; k++){
            if(j == i){
              //Again, first index of v[] is normalized to 1
              A_local[j * localColumn + k] -= tau[i] * tempVector_received[k-i-1];
            } else {
              A_local[j * localColumn + k] -= A_local[j * localColumn + i] * tau[i] * tempVector_received[k-i-1];
            }
          }
        }
      } else {
        for(j = 0; j < localRow; j++){
          for(k = i+1; k < localColumn; k++){

            A_local[j * localColumn + k] -= A_local[j * localColumn + i] * tau[i] * tempVector_received[k-i-1];

          }
        }
      }



    }
    //print_matrix_d(A_local, localRow, localColumn, world_rank, px, py);
  }
//printf("r= %d, c=%d, color_x = %d, color_y = %d, px= %d, py = %d \n", r, c, color_x, color_y, px, py);



	if(debug_mode){
		if(world_rank == master_node){
			double *A_global = (double *) malloc(b*b*f_r*f_c*r*c*sizeof(double));

			int *xyrank = (int *) malloc(2*sizeof(int));
			double *temp_matrix = (double *) malloc(b*b*f_r*f_c*sizeof(double));
			for(int iter = 0; iter < world_size; iter++){
				xyrank[0] = 0;
				xyrank[1] = 0;

				if(iter != master_node){
					MPI_Recv(xyrank, 2, MPI_INT, iter, iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Recv(temp_matrix, b*b*f_r*f_c,MPI_DOUBLE, iter, iter+9999, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					//printf("Hei I, proc master, receive from proc(%d,  %d).\n", xyrank[0], xyrank[1]);
				} else {
					for (i = 0; i < localRow; i++){
						for(j =0; j<localColumn; j++){
							temp_matrix[i*localColumn + j] = A_local[i*localColumn + j];
						}
					}
				}
				int x_cor = xyrank[0];
				int y_cor = xyrank[1];

				//Now we fill-in A_global
				//We use t1, t2 to denote the current block (t1, t2) in A_local of
				//the sending processor. We have to allocate this (t1, t2) block into
				//its correct place in A_global
				for(int t1 = 0; t1 < f_r; t1++){
					for(int t2 = 0; t2 < f_c; t2++){
						int globalRow = b*f_r*r;
						int globalColumn = b*f_c*c;
						//now write the block (t1, t2) of size b x b into A_global
						for(i = 0; i < b; i++){
							for(j = 0; j < b; j++){
								A_global[((x_cor+t1*r)*b+i)*globalColumn + (y_cor+t2*c)*b+j] = temp_matrix[(i+t1*b)*localColumn+(j+t2*b)];
								//printf("hei %.3f \n",   A_global[((x_cor+t1*r)*b+i)*globalColumn + (y_cor+t2*c)*b+j] );
							}
						}

					}
				}


			}
			int globalRow = b*f_r*r;
			int globalColumn = b*f_c*c;
			printf("Proc master print output of QR algorithm here: \n");
			print_matrix(A_global, globalRow, globalColumn);
			printf("--------------------------------------------------------\n");

			testing_QR_result(A_global, tau, globalRow, globalColumn);

		} else {
			int xyrank[2] = {px, py};
			MPI_Send(xyrank, 2, MPI_INT, 0, world_rank, MPI_COMM_WORLD);
			MPI_Send(A_local, b*b*f_r*f_c,MPI_DOUBLE, 0, world_rank+9999, MPI_COMM_WORLD);
			//printf("Hei I, proc(%d,  %d), have just sent to proc master.\n", xyrank[0], xyrank[1]);
		}

	}




/*
  if(world_rank == 0){
    //dumb testing
    double matA[6] = {1 , 2, 3 ,4, 5, 6};
    double matB[6] = {1, 0, 1, 0, 0, 1};
    double *res = (double *) malloc(4*sizeof(double));
    matrix_multiply(matA, matB, res, 2, 3, 2);
    print_matrix_d(matA, 2, 3, world_rank, px, py);
    print_matrix_d(matB, 3, 2, world_rank, px, py);
    printf("heiii finally");
    print_matrix_d(res, 2, 2, world_rank, px, py);
    //double matA[4] = res;
    //print_matrix_d(matA, 2, 2, world_rank, px, py);
  }*/

  //printf("Processor distribution = r * c = %d * %d \n", r, c);
  MPI_Finalize();
}

/*
 * master_result_writing.h
 *
 *  Created on: 12-Feb-2017
 *      Author: poonkodi
 */

#ifndef MASTER_RESULT_WRITING_H_
#define MASTER_RESULT_WRITING_H_

int result_file_write(int n, int start_iter, int no_of_iter, int prs,int ns
		MPI_Datatype trans_primitive, int position1, int d_ele,
		char filename[100], int ele, int node, MPI_Datatype trans_derv,
		FILE* fp4, int* i, double* max_id, MPI_Status* status, double* mxm,
		double* maximum, double* delt, int* count, int* maxsize, char* buffer,
		int* j, int** glelem, struct node5* primitive, FILE* fp2, FILE* fp5,
		struct node1* xycoord, struct node2* eleminfo, struct node6* derv,
		FILE* fp3)
{
	fp4 = fopen("error.txt", "w+");
	for (n = start_iter; n <= no_of_iter; n++)
	{
		for (*i = 1; *i < prs; *i++)
		{
			MPI_Recv(&max_id[*i], 1, MPI_DOUBLE, *i, *i, MPI_COMM_WORLD,
					&*status);
		}
		*mxm = max_id[1];
		for (*i = 1; *i < prs; *i++)
		{
			if (max_id[*i] > *mxm)
			{
				*mxm = max_id[*i];
			}
		}
		*maximum = *mxm;
		//print the maximum error in a file
		fprintf(fp4, "max error is %lf for %d iteration %lf delt\n", *maximum,
				n, *delt);
		printf("max error is %lf for %d iteration %lf delt\n", *maximum, n,
				*delt);
		//results files r wriiten
		if (n % ns == 0)
		{
			for (*i = 1; *i < prs; *i++)
			{
				MPI_Pack_size(count[*i], trans_primitive, MPI_COMM_WORLD,
						&*maxsize);
				buffer = malloc(*maxsize);
				if (buffer == NULL)
					printf("buffer memory failed\n");

				MPI_Recv(buffer, *maxsize, MPI_PACKED, *i, *i, MPI_COMM_WORLD,
						&*status);
				position1 = 0;
				for (*j = 1; *j <= count[*i]; *j++)
				{
					d_ele = glelem[*i][*j];
					MPI_Unpack(buffer, *maxsize, &position1, &primitive[d_ele],
							1, trans_primitive, MPI_COMM_WORLD);
				}
				free(buffer);
				buffer = NULL;
			}
			if (n % ns == 0)
			{
				sprintf(filename, "result_%d.txt", n);
				fp2 = fopen(filename, "a");
				for (*j = 1; *j <= ele; *j++)
				{
					fprintf(fp2, "%lf %lf %lf %lf %lf\n", primitive[*j].u1,
							primitive[*j].u2, primitive[*j].u3,
							primitive[*j].u4, primitive[*j].u5);
				}
				fclose(fp2);
			}
		}
		if (n % ns == 0)
		{
			sprintf(filename, "resultplt_%d.plt", n);
			fp5 = fopen(filename, "a");
			fprintf(fp5,
					"TITLE=\"3D\"\nVARIABLES=\t\"X\",\t\"Y\",\t\"Z\",\t\"dens\",\t\"u\",\t\"v\",\t\"w\",\t\"p\"\nZone N =%d, E=%d, DATAPACKING=BLOCK, VARLOCATION=(4=CELLCENTERED,5=CELLCENTERED,6=CELLCENTERED,7=CELLCENTERED,8=CELLCENTERED), ZONETYPE=FETETRAHEDRON\n",
					node, ele);
			for (*i = 1; *i <= node; *i++)
			{
				fprintf(fp5, "%lf\n", xycoord[*i].xstart);
			}
			for (*i = 1; *i <= node; *i++)
			{
				fprintf(fp5, "%lf\n", xycoord[*i].ystart);
			}
			for (*i = 1; *i <= node; *i++)
			{
				fprintf(fp5, "%lf\n", xycoord[*i].zstart);
			}
			for (*i = 1; *i <= ele; *i++)
			{
				fprintf(fp5, "%lf\n", primitive[*i].u1);
			}
			for (*i = 1; *i <= ele; *i++)
			{
				fprintf(fp5, "%lf\n", primitive[*i].u2 / primitive[*i].u1);
			}
			for (*i = 1; *i <= ele; *i++)
			{
				fprintf(fp5, "%lf\n", primitive[*i].u3 / primitive[*i].u1);
			}
			for (*i = 1; *i <= ele; *i++)
			{
				fprintf(fp5, "%lf\n", primitive[*i].u4 / primitive[*i].u1);
			}
			for (*i = 1; *i <= ele; *i++)
			{
				fprintf(fp5, "%lf\n",
						primitive[*i].u5
								- (0.5
										* (primitive[*i].u2 * primitive[*i].u2
												+ primitive[*i].u3
														* primitive[*i].u3
												+ primitive[*i].u4
														* primitive[*i].u4)
										/ primitive[*i].u1));
			}
			for (*i = 1; *i <= ele; *i++)
			{
				fprintf(fp5, "%d %d %d %d\n", eleminfo[*i].gele_nodes[1],
						eleminfo[*i].gele_nodes[2], eleminfo[*i].gele_nodes[3],
						eleminfo[*i].gele_nodes[4]);
			}
			fclose(fp5);
		}
		//bac up files r written
		if (n % 10000 == 0)
		{
			for (*i = 1; *i < prs; *i++)
			{
				MPI_Pack_size(count[*i], trans_derv, MPI_COMM_WORLD, &*maxsize);
				buffer = malloc(*maxsize);
				if (buffer == NULL)
					printf("buffer memory failed\n");

				MPI_Recv(buffer, *maxsize, MPI_PACKED, *i, *i, MPI_COMM_WORLD,
						&*status);
				position1 = 0;
				for (*j = 1; *j <= count[*i]; *j++)
				{
					d_ele = glelem[*i][*j];
					MPI_Unpack(buffer, *maxsize, &position1, &derv[d_ele], 1,
							trans_derv, MPI_COMM_WORLD);
				}
				free(buffer);
				buffer = NULL;
			}
			fp3 = fopen("initial_cond.txt", "w");
			for (*j = 1; *j <= ele; *j++)
			{
				fprintf(fp3,
						"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
						primitive[*j].u1, primitive[*j].u2, primitive[*j].u3,
						primitive[*j].u4, primitive[*j].u5, derv[*j].ux1,
						derv[*j].ux2, derv[*j].ux3, derv[*j].ux4, derv[*j].ux5,
						derv[*j].uy1, derv[*j].uy2, derv[*j].uy3, derv[*j].uy4,
						derv[*j].uy5, derv[*j].uz1, derv[*j].uz2, derv[*j].uz3,
						derv[*j].uz4, derv[*j].uz5);
			}
			fclose(fp3);
			fp3 = fopen("recovery.txt", "a");
			fprintf(fp3, "Back up file is written for %d iteration", n);
			fclose(fp3);
		}
	}
	fclose(fp4);
	return n;
}

#endif /* MASTER_RESULT_WRITING_H_ */

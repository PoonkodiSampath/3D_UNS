/*
 * master_preproc.h
 *
 *  Created on: 12-Feb-2017
 *      Author: poonkodi
 */

#ifndef MASTER_PREPROC_H_
#define MASTER_PREPROC_H_

//neigh_iden => function for neighbour identifiation

int neigh_iden(int iele, int ele, int jj, int node, int dstore, int dstore1,
		int point, int jele, int d_node[4], int ii, int prs, char filename[100],
		int* s1, int* s2, int* s3, struct node2* eleminfo, int* i,
		int* loc_ele_surr_pts, int* iface, int** elem_neighbour, FILE* fp)
{
	//loc_ele_surr_pts[ele_no+2]
	int ipoint, jpoint, kpoint, lpoint, istore, jstore, kstore, lstore;
	//node allocation has been done according to gambit conventions
	for (iele = 1; iele <= ele; iele++)
	{
		for (jj = 1; jj <= 4; jj++)
		{
			if (jj == 1)
			{
				*s1 = 2;
				*s2 = 1;
				*s3 = 3;
			}
			else if (jj == 2)
			{
				*s1 = 1;
				*s2 = 2;
				*s3 = 4;
			}
			else if (jj == 3)
			{
				*s1 = 2;
				*s2 = 3;
				*s3 = 4;
			}
			else if (jj == 4)
			{
				*s1 = 3;
				*s2 = 1;
				*s3 = 4;
			}

			eleminfo[iele].node1_localface[jj] = eleminfo[iele].gele_nodes[*s1];
			eleminfo[iele].node2_localface[jj] = eleminfo[iele].gele_nodes[*s2];
			eleminfo[iele].node3_localface[jj] = eleminfo[iele].gele_nodes[*s3];
		}
	}
	//Identification of elements surrounding points
	for (*i = 1; *i <= node + 1; *i++)
	{
		loc_ele_surr_pts[*i] = 0;
	}
	for (*i = 1; *i <= ele; *i++)
	{
		ipoint = eleminfo[*i].gele_nodes[1] + 1;
		jpoint = eleminfo[*i].gele_nodes[2] + 1;
		kpoint = eleminfo[*i].gele_nodes[3] + 1;
		lpoint = eleminfo[*i].gele_nodes[4] + 1;
		loc_ele_surr_pts[ipoint] = loc_ele_surr_pts[ipoint] + 1;
		loc_ele_surr_pts[jpoint] = loc_ele_surr_pts[jpoint] + 1;
		loc_ele_surr_pts[kpoint] = loc_ele_surr_pts[kpoint] + 1;
		loc_ele_surr_pts[lpoint] = loc_ele_surr_pts[lpoint] + 1;
	}
	for (*i = 2; *i <= node + 1; *i++)
	{
		loc_ele_surr_pts[*i] = loc_ele_surr_pts[*i] + loc_ele_surr_pts[*i - 1];
	}
	int npoint;
	npoint = loc_ele_surr_pts[node + 1] + 2;
	printf("node is %d \n", node);
	int* ele_surr_pts;
	ele_surr_pts = (int*) malloc((npoint + 1) * sizeof(int));
	if (ele_surr_pts == NULL)
		printf("ele_surr_pts memory failed4\n");

	for (*i = 1; *i <= ele; *i++)
	{
		ipoint = eleminfo[*i].gele_nodes[1];
		jpoint = eleminfo[*i].gele_nodes[2];
		kpoint = eleminfo[*i].gele_nodes[3];
		lpoint = eleminfo[*i].gele_nodes[4];
		istore = loc_ele_surr_pts[ipoint] + 1;
		jstore = loc_ele_surr_pts[jpoint] + 1;
		kstore = loc_ele_surr_pts[kpoint] + 1;
		lstore = loc_ele_surr_pts[lpoint] + 1;
		loc_ele_surr_pts[ipoint] = istore;
		loc_ele_surr_pts[jpoint] = jstore;
		loc_ele_surr_pts[kpoint] = kstore;
		loc_ele_surr_pts[lpoint] = lstore;
		ele_surr_pts[istore] = *i;
		ele_surr_pts[jstore] = *i;
		ele_surr_pts[kstore] = *i;
		ele_surr_pts[lstore] = *i;
	}
	loc_ele_surr_pts[0] = 0.0;

	//identification of element neighbours
	for (iele = 1; iele <= ele; iele++)
	{
		//j is the number of faces in the element
		for (*iface = 1; *iface <= 4; *iface++)
		{
			elem_neighbour[*iface][iele] = 0;
		}
	}
	for (iele = 1; iele <= ele; iele++)
	{
		//j is the number of faces in the element
		for (*iface = 1; *iface <= 4; *iface++)
		{
			ipoint = eleminfo[iele].node1_localface[*iface];
			dstore = loc_ele_surr_pts[ipoint - 1];
			dstore1 = loc_ele_surr_pts[ipoint];

			for (*i = dstore + 1; *i <= dstore1; *i++)
			{
				point = 0;
				if (iele != ele_surr_pts[*i])
				{
					jele = ele_surr_pts[*i];
					d_node[1] = eleminfo[iele].node1_localface[*iface];
					d_node[2] = eleminfo[iele].node2_localface[*iface];
					d_node[3] = eleminfo[iele].node3_localface[*iface];
					for (ii = 1; ii <= 3; ii++)
					{
						for (jj = 1; jj <= 4; jj++)
						{
							if (d_node[ii] == eleminfo[jele].gele_nodes[jj])
							{
								point = point + 1;
								break;
							}
						}
					}
				}
				if (point == 3)
				{
					elem_neighbour[*iface][iele] = jele;
					break;
				}
			}
		}
	}
	fp = fopen("elem_neigh.txt", "w");
	for (iele = 1; iele <= ele; iele++)
	{
		fprintf(fp, "%d %d %d %d %d\n", iele, elem_neighbour[1][iele],
				elem_neighbour[2][iele], elem_neighbour[3][iele],
				elem_neighbour[4][iele]);
	}
	fclose(fp);
	for (*i = 1; *i < prs; *i++)
	{
		sprintf(filename, "neigh_elem_%d.dat", *i);
		fp = fopen(filename, "wb");
		for (iele = 1; iele <= ele; iele++)
		{
			fwrite(&elem_neighbour[1][iele], sizeof(int), 1, fp);
			fwrite(&elem_neighbour[2][iele], sizeof(int), 1, fp);
			fwrite(&elem_neighbour[3][iele], sizeof(int), 1, fp);
			fwrite(&elem_neighbour[4][iele], sizeof(int), 1, fp);
		}
		fclose(fp);
	}
	fp = fopen("elem_conn.txt", "w");
	for (iele = 1; iele <= ele; iele++)
	{
		fprintf(fp, "%d %d %d %d\n", eleminfo[iele].gele_nodes[1],
				eleminfo[iele].gele_nodes[2], eleminfo[iele].gele_nodes[3],
				eleminfo[iele].gele_nodes[4]);
	}
	fclose(fp);
	for (*i = 1; *i < prs; *i++)
	{
		sprintf(filename, "conn_elem_%d.dat", *i);
		fp = fopen(filename, "wb");
		for (iele = 1; iele <= ele; iele++)
		{
			fwrite(&eleminfo[iele].gele_nodes[1], sizeof(int), 1, fp);
			fwrite(&eleminfo[iele].gele_nodes[2], sizeof(int), 1, fp);
			fwrite(&eleminfo[iele].gele_nodes[3], sizeof(int), 1, fp);
			fwrite(&eleminfo[iele].gele_nodes[4], sizeof(int), 1, fp);
		}
		fclose(fp);
	}
	int* elem2proc;
	elem2proc = (int*) malloc((ele + 2) * sizeof(int));
	if (elem2proc == NULL)
		printf("elem2proc memory failed\n");

	fp = fopen("elem2proc.txt", "r");
	for (iele = 1; iele <= ele; iele++)
	{
		fscanf(fp, "%d\n", &elem2proc[iele]);
	}
	fclose(fp);
	for (*i = 1; *i < prs; *i++)
	{
		sprintf(filename, "proc2elem_%d.txt", *i);
		fp = fopen(filename, "w");
		for (iele = 1; iele <= ele; iele++)
		{
			fprintf(fp, "%d\n", elem2proc[iele]);
		}
		fclose(fp);
	}

	free(elem2proc);
	return iele;
}

//initial_cond_read_write => function for writing initial condition in a file if simulations are to begin from the iteration 1 or
// read from the previous iteration file for restarting

int initial_cond_read_write(int fg, int iele, int ele, FILE* fp, double* dummy,
		double* u1initial, double* u2initial, double* u3initial,
		double* u4initial, double* u5initial, struct node5* primitive,
		struct node6* derv)
{
	if (fg == 1)
	{
		fp = fopen("initial_cond.txt", "w");
		*dummy = 0.0;
		for (iele = 1; iele <= ele; iele++)
		{
			if (iele == ele)
			{
				fprintf(fp,
						"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
						*u1initial, *u2initial, *u3initial, *u4initial,
						*u5initial, *dummy, *dummy, *dummy, *dummy, *dummy,
						*dummy, *dummy, *dummy, *dummy, *dummy, *dummy, *dummy,
						*dummy, *dummy, *dummy);
			}
			else
			{
				fprintf(fp,
						"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
						*u1initial, *u2initial, *u3initial, *u4initial,
						*u5initial, *dummy, *dummy, *dummy, *dummy, *dummy,
						*dummy, *dummy, *dummy, *dummy, *dummy, *dummy, *dummy,
						*dummy, *dummy, *dummy);
			}
		}
		fclose(fp);
	}
	fp = fopen("initial_cond.txt", "r");
	for (iele = 1; iele <= ele; iele++)
	{
		fscanf(fp,
				"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
				&primitive[iele].u1, &primitive[iele].u2, &primitive[iele].u3,
				&primitive[iele].u4, &primitive[iele].u5, &derv[iele].ux1,
				&derv[iele].ux2, &derv[iele].ux3, &derv[iele].ux4,
				&derv[iele].ux5, &derv[iele].uy1, &derv[iele].uy2,
				&derv[iele].uy3, &derv[iele].uy4, &derv[iele].uy5,
				&derv[iele].uz1, &derv[iele].uz2, &derv[iele].uz3,
				&derv[iele].uz4, &derv[iele].uz5);
	}
	fclose(fp);
	return iele;
}

FILE* pre_meshfile_read(char oneword[1000], int temp, int result, int prs,
		char filename[100], FILE* fp, int* i, int* no_of_nodes, int* int_ele,
		int* bc_set, int* ele, int* node, struct node1* xycoord,
		struct node2* eleminfo, FILE* fp1, int* d1, int* periodic_cnt, int* j,
		int* iele, int* jface, int* k)
{
	//Reading gambit file
	fp = fopen("mesh_file.txt", "r");
	if (access("mesh_file.txt", F_OK) == -1)
	{
		printf("unable to open mesh file\n");
	}
	for (*i = 1; *i <= 6; *i++)
	{
		fgets(oneword, 500, fp);
	}
	fscanf(fp, "%d %d %d %d %d %d\n", &*no_of_nodes, &*int_ele, &temp, &*bc_set,
			&temp, &temp);
	*ele = *int_ele;
	*node = *no_of_nodes;
	*no_of_nodes = *no_of_nodes + 1;
	*int_ele = *int_ele + 1;
	xycoord = (struct node1*) malloc(*no_of_nodes * sizeof(struct node1));
	if (xycoord == NULL)
		printf("memory failed4\n");

	eleminfo = (struct node2*) malloc(*int_ele * sizeof(struct node2));
	if (eleminfo == NULL)
		printf("memory failed4\n");

	//reading of mesh file contd
	for (*i = 1; *i <= 2; *i++)
	{
		fgets(oneword, 500, fp);
	}
	for (*i = 1; *i < *no_of_nodes; *i++)
	{
		fscanf(fp, "%d ", &temp);
		fscanf(fp, "%lf %lf %lf\n", &xycoord[temp].xstart,
				&xycoord[temp].ystart, &xycoord[temp].zstart);
	}
	for (*i = 1; *i <= 2; *i++)
	{
		fgets(oneword, 500, fp);
	}
	for (*i = 1; *i < *int_ele; *i++)
	{
		fscanf(fp, "%d %d %d %d %d %d %d\n", &temp, &temp, &temp,
				&eleminfo[*i].gele_nodes[1], &eleminfo[*i].gele_nodes[2],
				&eleminfo[*i].gele_nodes[3], &eleminfo[*i].gele_nodes[4]);
	}
	fgets(oneword, 500, fp);
	do
	{
		fgets(oneword, 500, fp);
		result = strncmp(oneword, "ENDOFSECTION", 11);
	} while (result != 0);
	fgets(oneword, 500, fp);
	fp1 = fopen("boundary_data.txt", "w");
	for (*i = 1; *i <= *bc_set; *i++)
	{
		fscanf(fp, "%s %d %d %d %d\n", oneword, &temp, &*d1, &temp, &temp);
		fprintf(fp1, "%s %d\n", oneword, *d1);
		if (strncmp(oneword, "periodic", 8) == 0)
		{
			*periodic_cnt = *periodic_cnt + *d1;
		}
		for (*j = 1; *j <= *d1; *j++)
		{
			fscanf(fp, "%d %d %d\n", &*iele, &temp, &*jface);
			fprintf(fp1, "%d %d\n", *iele, *jface);
		}
		for (*j = 1; *j <= 2; *j++)
		{
			fgets(oneword, 500, fp);
		}
	}
	fclose(fp);
	fclose(fp1);
	for (*k = 1; *k < prs; *k++)
	{
		fp1 = fopen("boundary_data.txt", "r");
		sprintf(filename, "boundary_data_%d.txt", *k);
		fp = fopen(filename, "w");
		for (*i = 1; *i <= *bc_set; *i++)
		{
			fscanf(fp1, "%s %d\n", oneword, &*d1);
			fprintf(fp, "%s %d\n", oneword, *d1);
			for (*j = 1; *j <= *d1; *j++)
			{
				fscanf(fp1, "%d %d\n", &*iele, &*jface);
				fprintf(fp, "%d %d\n", *iele, *jface);
			}
		}
		fclose(fp);
		fclose(fp1);
	}
	return fp;
}

FILE* pre_meshdata_periodic(int periodic_cnt, char oneword[1000], char c1,
		int t1, int b1, int b2, int t2, int period, int temp1, int temp2,
		int t3, int t4, int check, int prs, char filename[100], int ele,
		FILE* fp, int* zonetype, int* periodic_pair1, int* periodic_pair2,
		int** periodic_face, int* i, int* j, int* periodic_elem_pair1,
		int* periodic_elem_pair2, int* periodic_mapping, int* iele)
{
	if ((periodic_cnt % 2) != 0)
	{
		printf("periodic cnt in nuetral file has not been alloted properly\n");
		exit(0);
	}
	//reading data from msh file
	fp = fopen("msh_gambitfile.txt", "r");
	if (access("msh_gambitfile.txt", F_OK) == -1)
	{
		printf("unable to open mesh file\n");
	}
	while (!(feof(fp)))
	{
		oneword[0] = '\0';
		fgets(oneword, 500, fp);
		if (strlen(oneword) != 0)
		{
			if (strncmp(oneword, "(18", 3) == 0)
			{
				sscanf(oneword, "%c %d %c %x %x %x %x", &c1, &t1, &c1, &t1, &t1,
						&zonetype[b1], &zonetype[b2]);
				while (1)
				{
					fgets(oneword, 500, fp);
					if (strncmp(oneword, "))", 2) == 0)
					{
						break;
					}
					sscanf(oneword, "%x %x", &t1, &t2);
					period = period + 1;
					periodic_pair1[period] = t1;
					periodic_pair2[period] = t2;
					*periodic_face[t1] = period;
					*periodic_face[t2] = period;

				}
				b1 = b1 + 2;
				b2 = b2 + 2;
			}
			else
			{
				if (strncmp(oneword + strlen(oneword) - 2, "(", 1) == 0)
				{
					sscanf(oneword, "%c %d %c %x %x %x", &c1, &t1, &c1, &t2,
							&temp1, &temp2);

					for (*i = temp1; *i <= temp2 + 1; *i++)
					{
						fgets(oneword, 500, fp);
					}
					if (strncmp(oneword, "))", 2) != 0)
					{
						printf("Problem in reading msh file\n");
						exit(0);
					}
				}
				else
				{
					if (strncmp(oneword, "(13", 3) == 0)
					{

						sscanf(oneword, "%c %d %c %x %x %x", &c1, &t1, &c1, &t2,
								&temp1, &temp2);

						if (t2 == 0)
						{

							if (temp1 == 1)
							{
								*periodic_face = (int*) malloc(
										(temp2 + 1) * sizeof(int));
							}
							else
							{
								printf(
										"error in memory allocation of variable periodic_face\n");
								exit(0);
							}
						}
					}
				}
			}
		}
	}
	//while loop end
	if ((periodic_cnt / 2.0) != period)
	{
		printf(
				"periodic_cnt and period values r %d %d and they should be equal. Data from nuertal/msh file has not been read properly\n",
				periodic_cnt, period);
		exit(0);
	}
	fclose(fp);
	fp = fopen("msh_gambitfile.txt", "r");
	while (!(feof(fp)))
	{
		oneword[0] = '\0';
		fgets(oneword, 500, fp);
		if (strlen(oneword) != 0)
		{
			if (strncmp(oneword, "(13", 3) == 0)
			{
				sscanf(oneword, "%c %d %c %x %x %x\n", &c1, &t1, &c1, &t2,
						&temp1, &temp2);
				if (t2 != 0)
				{
					for (*i = 1; *i <= b2; *i++)
					{
						if (t2 == zonetype[*i])
						{
							for (*j = temp1; *j <= temp2; *j++)
							{
								fscanf(fp, "%d %x %x %x %x %x\n", &t1, &t1, &t1,
										&t1, &t3, &t4);
								check = 0;
								if ((t4 == 0) || (t3 == 0))
								{
									//atleast one of the neghbouring cell should be zero indicating that it is a boundary cell
									check = 1;
								}
								if ((t4 == 0) && (t3 == 0))
								{
									//both the neighbouring cells should not be zero
									check = 0;
								}
								if (check == 0)
								{
									printf(
											"Problem in reading face data of periodic boundary condition\n");
									exit(0);
								}
								if (t3 != 0)
								{
									if (periodic_elem_pair1[*periodic_face[*j]]
											== 0)
									{
										periodic_elem_pair1[*periodic_face[*j]] =
												t3;
									}
									else
									{
										periodic_elem_pair2[*periodic_face[*j]] =
												t3;
									}
								}
								if (t4 != 0)
								{
									if (periodic_elem_pair1[*periodic_face[*j]]
											== 0)
									{
										periodic_elem_pair1[*periodic_face[*j]] =
												t4;
									}
									else
									{
										periodic_elem_pair2[*periodic_face[*j]] =
												t4;
									}
								}
							}
							//j loop
							fgets(oneword, 500, fp);
							if (strncmp(oneword, "))", 2) != 0)
							{
								printf(
										"Periodic data has not been read properly\n");
								exit(0);
							}
							break;
						}
					}
				}
			} //if strncmp
			else
			{
				if (strncmp(oneword + strlen(oneword) - 2, "(", 1) == 0)
				{
					sscanf(oneword, "%c %d %x %d %x %x", &c1, &t1, &c1, &t2,
							&temp1, &temp2);
					for (*i = temp1; *i <= temp2 + 1; *i++)
					{
						fgets(oneword, 500, fp);
					}
					if (strncmp(oneword, "))", 2) != 0)
					{
						printf("Problem in reading msh file\n");
						exit(0);
					}
				}
			}
		}
	}
	fclose(fp);
	for (*j = 1; *j < prs; *j++)
	{
		sprintf(filename, "periodic_elements_%d.txt", *j);
		fp = fopen(filename, "w");
		for (*i = 1; *i <= period; *i++)
		{
			periodic_mapping[periodic_elem_pair1[*i]] = periodic_elem_pair2[*i];
			periodic_mapping[periodic_elem_pair2[*i]] = periodic_elem_pair1[*i];
			fprintf(fp, "%d\n", periodic_elem_pair1[*i]);
			fprintf(fp, "%d\n", periodic_elem_pair2[*i]);
		}
		fclose(fp);
	}
	for (*j = 1; *j < prs; *j++)
	{
		sprintf(filename, "periodic_ele_pairing_%d.dat", *j);
		fp = fopen(filename, "wb");
		for (*iele = 1; *iele <= ele; *iele++)
		{
			fwrite(&periodic_mapping[*iele], sizeof(int), 1, fp);
		}
		fclose(fp);
	}
	return fp;
}

#endif /* MASTER_PREPROC_H_ */

/*
 * slave_preproc.h
 *
 *  Created on: 12-Feb-2017
 *      Author: poonkodi
 */

#ifndef SLAVE_PREPROC_H_
#define SLAVE_PREPROC_H_

//vol_tetra:volume of the tetrahedron

void vol_tetra(double *xstart_d1, double *ystart_d1, double *zstart_d1,
		double *xstart_d2, double *ystart_d2, double *zstart_d2,
		double *xstart_d3, double *ystart_d3, double *zstart_d3,
		double *xstart_d4, double *ystart_d4, double *zstart_d4, double *vol)
{
	a1 = *xstart_d1 - *xstart_d2;
	a2 = *xstart_d2 - *xstart_d3;
	a3 = *xstart_d3 - *xstart_d4;
	a4 = *ystart_d1 - *ystart_d2;
	a5 = *ystart_d2 - *ystart_d3;
	a6 = *ystart_d3 - *ystart_d4;
	a7 = *zstart_d1 - *zstart_d2;
	a8 = *zstart_d2 - *zstart_d3;
	a9 = *zstart_d3 - *zstart_d4;

	*vol = ((a1 * (a5 * a9 - a8 * a6)) - (a2 * (a4 * a9 - a7 * a6))
			+ (a3 * (a4 * a8 - a7 * a5))) / 6.0;
	if (*vol < 0.0)
	{
		*vol = -(*vol);
	}

}

//outward normal for any face

void triarea_normal(double *xstart_d1, double *ystart_d1, double *zstart_d1,
		double *xstart_d2, double *ystart_d2, double *zstart_d2,
		double *xstart_d3, double *ystart_d3, double *zstart_d3, double *xc,
		double *yc, double *zc, double *nx, double *ny, double *nz,
		double *area, double *deno, double *z1, double *z2, double *z3,
		double *sum)
{
	*nx = ((*ystart_d2 - *ystart_d1) * (*zstart_d3 - *zstart_d1))
			- ((*ystart_d3 - *ystart_d1) * (*zstart_d2 - *zstart_d1));

	*ny = -(((*zstart_d3 - *zstart_d1) * (*xstart_d2 - *xstart_d1))
			- ((*xstart_d3 - *xstart_d1) * (*zstart_d2 - *zstart_d1)));

	*nz = ((*xstart_d2 - *xstart_d1) * (*ystart_d3 - *ystart_d1))
			- ((*xstart_d3 - *xstart_d1) * (*ystart_d2 - *ystart_d1));

	*deno = sqrt((*nx * *nx) + (*ny * *ny) + (*nz * *nz));

	*z1 = *xc - *xstart_d1;
	*z2 = *yc - *ystart_d1;
	*z3 = *zc - *zstart_d1;

	*sum = (*nx * *z1) + (*ny * *z2) + (*nz * *z3);

	if (*sum > 0.0)
	{
		*nx = -*nx / (*deno);
		*ny = -*ny / (*deno);
		*nz = -*nz / (*deno);
	}
	else
	{
		*nx = *nx / (*deno);
		*ny = *ny / (*deno);
		*nz = *nz / (*deno);
	}

	*area = *deno / 2.0;

}

//sol_pt_cal : function for calculating solution point in the slave

int sol_pt_cal(int iele, int lele_cnt, int iface, int d1, int d2, int d3,
		int d4, struct node2* eleminfo, int** lelem_neighbour,
		double** xc_neighface, struct node1* xycoord, struct node3* centroid,
		double** yc_neighface, double** zc_neighface, double** vol_neighface,
		double** xc_face, double** yc_face, double** zc_face, double** vol_face,
		double** vol_hexahedra_face, double** xc_hexahedra_face,
		double** yc_hexahedra_face, double** zc_hexahedra_face,
		double** xc_lateral_side1, double** yc_lateral_side1,
		double** zc_lateral_side1, double** xc_lateral_side2,
		double** yc_lateral_side2, double** zc_lateral_side2,
		double** xc_lateral_side3, double** yc_lateral_side3,
		double** zc_lateral_side3, double** nx_lateral_side1,
		double** ny_lateral_side1, double** nz_lateral_side1,
		double** area_lateral_side1, double* deno, double* z1, double* z2,
		double* z3, double* sum, double** nx_lateral_side2,
		double** ny_lateral_side2, double** nz_lateral_side2,
		double** area_lateral_side2, double** nx_lateral_side3,
		double** ny_lateral_side3, double** nz_lateral_side3,
		double** area_lateral_side3, double* vol_polyhedra, struct node4* sol)
{
	//calculation of solution point
	for (iele = 1; iele <= lele_cnt; iele++)
	{
		for (iface = 1; iface <= 4; iface++)
		{
			d1 = eleminfo[iele].node1_localface[iface];
			d2 = eleminfo[iele].node2_localface[iface];
			d3 = eleminfo[iele].node3_localface[iface];
			d4 = lelem_neighbour[iface][iele];
			xc_neighface[iface][iele] = (xycoord[d1].xstart + xycoord[d2].xstart
					+ xycoord[d3].xstart + centroid[d4].xc) / 4;
			yc_neighface[iface][iele] = (xycoord[d1].ystart + xycoord[d2].ystart
					+ xycoord[d3].ystart + centroid[d4].yc) / 4;
			zc_neighface[iface][iele] = (xycoord[d1].zstart + xycoord[d2].zstart
					+ xycoord[d3].zstart + centroid[d4].zc) / 4;
			vol_tetra(&xycoord[d1].xstart, &xycoord[d1].ystart,
					&xycoord[d1].zstart, &xycoord[d2].xstart,
					&xycoord[d2].ystart, &xycoord[d2].zstart,
					&xycoord[d3].xstart, &xycoord[d3].ystart,
					&xycoord[d3].zstart, &centroid[d4].xc, &centroid[d4].yc,
					&centroid[d4].zc, &vol_neighface[iface][iele]);
			xc_face[iface][iele] = (xycoord[d1].xstart + xycoord[d2].xstart
					+ xycoord[d3].xstart + centroid[iele].xc) / 4;
			yc_face[iface][iele] = (xycoord[d1].ystart + xycoord[d2].ystart
					+ xycoord[d3].ystart + centroid[iele].yc) / 4;
			zc_face[iface][iele] = (xycoord[d1].zstart + xycoord[d2].zstart
					+ xycoord[d3].zstart + centroid[iele].zc) / 4;
			vol_tetra(&xycoord[d1].xstart, &xycoord[d1].ystart,
					&xycoord[d1].zstart, &xycoord[d2].xstart,
					&xycoord[d2].ystart, &xycoord[d2].zstart,
					&xycoord[d3].xstart, &xycoord[d3].ystart,
					&xycoord[d3].zstart, &centroid[iele].xc, &centroid[iele].yc,
					&centroid[iele].zc, &vol_face[iface][iele]);
			vol_hexahedra_face[iface][iele] = vol_neighface[iface][iele]
					+ vol_face[iface][iele];
			xc_hexahedra_face[iface][iele] = (vol_neighface[iface][iele]
					* xc_neighface[iface][iele]
					+ vol_face[iface][iele] * xc_face[iface][iele])
					/ vol_hexahedra_face[iface][iele];
			yc_hexahedra_face[iface][iele] = (vol_neighface[iface][iele]
					* yc_neighface[iface][iele]
					+ vol_face[iface][iele] * yc_face[iface][iele])
					/ vol_hexahedra_face[iface][iele];
			zc_hexahedra_face[iface][iele] = (vol_neighface[iface][iele]
					* zc_neighface[iface][iele]
					+ vol_face[iface][iele] * zc_face[iface][iele])
					/ vol_hexahedra_face[iface][iele];

			//lateral volumes (faces calculation)
			xc_lateral_side1[iface][iele] = (xycoord[d1].xstart
					+ xycoord[d2].xstart + centroid[d4].xc) / 3;
			yc_lateral_side1[iface][iele] = (xycoord[d1].ystart
					+ xycoord[d2].ystart + centroid[d4].yc) / 3;
			zc_lateral_side1[iface][iele] = (xycoord[d1].zstart
					+ xycoord[d2].zstart + centroid[d4].zc) / 3;
			xc_lateral_side2[iface][iele] = (xycoord[d2].xstart
					+ xycoord[d3].xstart + centroid[d4].xc) / 3;
			yc_lateral_side2[iface][iele] = (xycoord[d2].ystart
					+ xycoord[d3].ystart + centroid[d4].yc) / 3;
			zc_lateral_side2[iface][iele] = (xycoord[d2].zstart
					+ xycoord[d3].zstart + centroid[d4].zc) / 3;
			xc_lateral_side3[iface][iele] = (xycoord[d3].xstart
					+ xycoord[d1].xstart + centroid[d4].xc) / 3;
			yc_lateral_side3[iface][iele] = (xycoord[d3].ystart
					+ xycoord[d1].ystart + centroid[d4].yc) / 3;
			zc_lateral_side3[iface][iele] = (xycoord[d3].zstart
					+ xycoord[d1].zstart + centroid[d4].zc) / 3;
			triarea_normal(&xycoord[d1].xstart, &xycoord[d1].ystart,
					&xycoord[d1].zstart, &xycoord[d2].xstart,
					&xycoord[d2].ystart, &xycoord[d2].zstart, &centroid[d4].xc,
					&centroid[d4].yc, &centroid[d4].zc,
					&xc_hexahedra_face[iface][iele],
					&yc_hexahedra_face[iface][iele],
					&zc_hexahedra_face[iface][iele],
					&nx_lateral_side1[iface][iele],
					&ny_lateral_side1[iface][iele],
					&nz_lateral_side1[iface][iele],
					&area_lateral_side1[iface][iele], deno, z1, z2, z3, sum);
			triarea_normal(&xycoord[d2].xstart, &xycoord[d2].ystart,
					&xycoord[d2].zstart, &xycoord[d3].xstart,
					&xycoord[d3].ystart, &xycoord[d3].zstart, &centroid[d4].xc,
					&centroid[d4].yc, &centroid[d4].zc,
					&xc_hexahedra_face[iface][iele],
					&yc_hexahedra_face[iface][iele],
					&zc_hexahedra_face[iface][iele],
					&nx_lateral_side2[iface][iele],
					&ny_lateral_side2[iface][iele],
					&nz_lateral_side2[iface][iele],
					&area_lateral_side2[iface][iele], deno, z1, z2, z3, sum);
			triarea_normal(&xycoord[d3].xstart, &xycoord[d3].ystart,
					&xycoord[d3].zstart, &xycoord[d1].xstart,
					&xycoord[d1].ystart, &xycoord[d1].zstart, &centroid[d4].xc,
					&centroid[d4].yc, &centroid[d4].zc,
					&xc_hexahedra_face[iface][iele],
					&yc_hexahedra_face[iface][iele],
					&zc_hexahedra_face[iface][iele],
					&nx_lateral_side3[iface][iele],
					&ny_lateral_side3[iface][iele],
					&nz_lateral_side3[iface][iele],
					&area_lateral_side3[iface][iele], deno, z1, z2, z3, sum);
		}

		sol[iele].spx = (vol_hexahedra_face[1][iele]
				* xc_hexahedra_face[1][iele]
				+ vol_hexahedra_face[2][iele] * xc_hexahedra_face[2][iele]
				+ vol_hexahedra_face[3][iele] * xc_hexahedra_face[3][iele]
				+ vol_hexahedra_face[4][iele] * xc_hexahedra_face[4][iele])
				/ vol_polyhedra[iele];
		sol[iele].spy = (vol_hexahedra_face[1][iele]
				* yc_hexahedra_face[1][iele]
				+ vol_hexahedra_face[2][iele] * yc_hexahedra_face[2][iele]
				+ vol_hexahedra_face[3][iele] * yc_hexahedra_face[3][iele]
				+ vol_hexahedra_face[4][iele] * yc_hexahedra_face[4][iele])
				/ vol_polyhedra[iele];
		sol[iele].spz = (vol_hexahedra_face[1][iele]
				* zc_hexahedra_face[1][iele]
				+ vol_hexahedra_face[2][iele] * zc_hexahedra_face[2][iele]
				+ vol_hexahedra_face[3][iele] * zc_hexahedra_face[3][iele]
				+ vol_hexahedra_face[4][iele] * zc_hexahedra_face[4][iele])
				/ vol_polyhedra[iele];
	}
	return iele;
}

//function centroid_cal : calculation of centroids in slaves

int centroid_cal(int iele, int lele_cnt, int* d1, struct node2* eleminfo,
		int* d2, int* d3, int* d4, struct node3* centroid,
		struct node1* xycoord)
{
	//centroid of tetrahedrons
	for (iele = 1; iele <= lele_cnt; iele++)
	{
		*d1 = eleminfo[iele].gele_nodes[1];
		*d2 = eleminfo[iele].gele_nodes[2];
		*d3 = eleminfo[iele].gele_nodes[3];
		*d4 = eleminfo[iele].gele_nodes[4];
		centroid[iele].xc = (xycoord[*d1].xstart + xycoord[*d2].xstart
				+ xycoord[*d3].xstart + xycoord[*d4].xstart) / 4;
		centroid[iele].yc = (xycoord[*d1].ystart + xycoord[*d2].ystart
				+ xycoord[*d3].ystart + xycoord[*d4].ystart) / 4;
		centroid[iele].zc = (xycoord[*d1].zstart + xycoord[*d2].zstart
				+ xycoord[*d3].zstart + xycoord[*d4].zstart) / 4;

	}
	return iele;
}

int face_normal_calc(int iele, int lele_cnt, int d2, int d3, int* jface,
		int* s1, int* s2, int* s3, int* d1, int** lelem_neighbour, double* aa1,
		struct node4* sol, double* bb1, double* cc1, double* dd1, double* deno,
		double** num_dd, double* normal_x, double* normal_y, double* normal_z,
		double* t, double** xr, double** yr, double** zr)
{
	for (iele = 1; iele <= lele_cnt; iele++)
	{
		for (*jface = 1; *jface <= 4; *jface++)
		{
			if (*jface == 1)
			{
				*s1 = 2;
				*s2 = 3;
				*s3 = 4;
			}
			else if (*jface == 2)
			{
				*s1 = 1;
				*s2 = 3;
				*s3 = 4;
			}
			else if (*jface == 3)
			{
				*s1 = 1;
				*s2 = 2;
				*s3 = 4;
			}
			else if (*jface == 4)
			{
				*s1 = 1;
				*s2 = 2;
				*s3 = 3;
			}

			*d1 = lelem_neighbour[*s1][iele];
			d2 = lelem_neighbour[*s2][iele];
			d3 = lelem_neighbour[*s3][iele];
			*aa1 = ((sol[d2].spy - sol[*d1].spy) * (sol[d3].spz - sol[*d1].spz))
					- ((sol[d3].spy - sol[*d1].spy)
							* (sol[d2].spz - sol[*d1].spz));
			*bb1 = -(((sol[d2].spx - sol[*d1].spx)
					* (sol[d3].spz - sol[*d1].spz))
					- ((sol[d3].spx - sol[*d1].spx)
							* (sol[d2].spz - sol[*d1].spz)));
			*cc1 = ((sol[d2].spx - sol[*d1].spx) * (sol[d3].spy - sol[*d1].spy))
					- ((sol[d3].spx - sol[*d1].spx)
							* (sol[d2].spy - sol[*d1].spy));
			*dd1 = -*aa1 * sol[*d1].spx - *bb1 * sol[*d1].spy
					- *cc1 * sol[*d1].spz;
			//num_dd=> numerical domain of dependence
			*deno = pow(((*aa1 * *aa1) + (*bb1 * *bb1) + (*cc1 * *cc1)), 0.5);
			num_dd[*jface][iele] = (*aa1 * sol[iele].spx + *bb1 * sol[iele].spy
					+ *cc1 * sol[iele].spz + *dd1) / (*deno);
			num_dd[*jface][iele] = fabs(num_dd[*jface][iele]);
			*normal_x = *aa1;
			*normal_y = *bb1;
			*normal_z = *cc1;
			*normal_x = -*normal_x / (*deno);
			*normal_y = -*normal_y / (*deno);
			*normal_z = -*normal_z / (*deno);
			*t = (-*dd1 - *aa1 * sol[iele].spx - *bb1 * sol[iele].spy
					- *cc1 * sol[iele].spz)
					/ (*aa1 * *normal_x + *bb1 * *normal_y + *cc1 * *normal_z);
			xr[*jface][iele] = sol[iele].spx + *normal_x * *t;
			yr[*jface][iele] = sol[iele].spy + *normal_y * *t;
			zr[*jface][iele] = sol[iele].spz + *normal_z * *t;
		}
	}
	return iele;
}

int slave_global2local_numbering(int lele_cnt, char filename[100], int id,
		int ele, int n4, int d_store, int check, FILE* fp1, int* i,
		int* elem2proc, int* lelem, int* gb, int* n1, int* n2, int* n3,
		int** lelem_neighbour, double* dummy, int* sd, int* kk, int* id_proc,
		int* k2id, int* sd_cnt, FILE* fp, int* j)
{
	lele_cnt = 0;
	//variables for calculations
	sprintf(filename, "neigh_elem_%d.dat", id);
	fp1 = fopen(filename, "r");
	for (*i = 1; *i <= ele; *i++)
	{
		if (elem2proc[*i] == id - 1)
		{
			lele_cnt = lele_cnt + 1;
			lelem[lele_cnt] = *i;
			//gb -inverse mapping from global to local mapping
			gb[*i] = lele_cnt;
			//finding the neighbours and boundaries of subdomain
			fseek(fp1, (*i - 1) * sizeof(int) * 4, SEEK_SET);
			fread(&*n1, sizeof(int), 1, fp1);
			fread(&*n2, sizeof(int), 1, fp1);
			fread(&*n3, sizeof(int), 1, fp1);
			fread(&n4, sizeof(int), 1, fp1);
			lelem_neighbour[1][lele_cnt] = *n1;
			lelem_neighbour[2][lele_cnt] = *n2;
			lelem_neighbour[3][lele_cnt] = *n3;
			lelem_neighbour[4][lele_cnt] = n4;

			if (*n1 == 0)
			{
				//gb-global boundary
				dummy = dummy + 1;
			}
			else if (elem2proc[*n1] != (id - 1))
			{
				//sd-sub domain
				*sd = *sd + 1;
				if (*sd == 1)
				{
					*kk = *kk + 1;
					id_proc[*kk] = elem2proc[*n1];
					d_store = id_proc[*kk];
					k2id[d_store] = *kk;
					sd_cnt[*kk] = 1;
					sprintf(filename, "sd_boun_%d.txt", id);
					fp = fopen(filename, "a");
					fprintf(fp, " %d %d %d\n", *i, *n1, elem2proc[*n1]);
					fclose(fp);
				} //if sd==1 loop
				else
				{
					check = 0;
					for (*j = 1; *j <= *kk; *j++)
					{
						if (elem2proc[*n1] == id_proc[*j])
						{
							check = 1;
							sd_cnt[*j] = sd_cnt[*j] + 1;
							sprintf(filename, "sd_boun_%d.txt", id);
							fp = fopen(filename, "a");
							fprintf(fp, " %d %d %d\n", *i, *n1, elem2proc[*n1]);
							fclose(fp);
							break;
						}
					}
					if (check == 0)
					{
						*kk = *kk + 1;
						id_proc[*kk] = elem2proc[*n1];
						sd_cnt[*kk] = 1;
						d_store = id_proc[*kk];
						k2id[d_store] = *kk;
						sprintf(filename, "sd_boun_%d.txt", id);
						fp = fopen(filename, "a");
						fprintf(fp, " %d %d %d\n", *i, *n1, elem2proc[*n1]);
						fclose(fp);
					}
				}
			}

			if (*n2 == 0)
			{
				//gb-global boundary
				dummy = dummy + 1;
			}
			else if (elem2proc[*n2] != (id - 1))
			{
				//sd-sub domain
				*sd = *sd + 1;
				if (*sd == 1)
				{
					*kk = *kk + 1;
					id_proc[*kk] = elem2proc[*n2];
					d_store = id_proc[*kk];
					k2id[d_store] = *kk;
					sd_cnt[*kk] = 1;
					sprintf(filename, "sd_boun_%d.txt", id);
					fp = fopen(filename, "a");
					fprintf(fp, " %d %d %d\n", *i, *n2, elem2proc[*n2]);
					fclose(fp);
				} //if sd==1 loop
				else
				{
					check = 0;
					for (*j = 1; *j <= *kk; *j++)
					{
						if (elem2proc[*n2] == id_proc[*j])
						{
							check = 1;
							sd_cnt[*j] = sd_cnt[*j] + 1;
							sprintf(filename, "sd_boun_%d.txt", id);
							fp = fopen(filename, "a");
							fprintf(fp, " %d %d %d\n", *i, *n2, elem2proc[*n2]);
							fclose(fp);
							break;
						}
					}
					if (check == 0)
					{
						*kk = *kk + 1;
						id_proc[*kk] = elem2proc[*n2];
						sd_cnt[*kk] = 1;
						d_store = id_proc[*kk];
						k2id[d_store] = *kk;
						sprintf(filename, "sd_boun_%d.txt", id);
						fp = fopen(filename, "a");
						fprintf(fp, " %d %d %d\n", *i, *n2, elem2proc[*n2]);
						fclose(fp);
					}
				}
			}

			if (*n3 == 0)
			{
				//gb-global boundary
				dummy = dummy + 1;
			}
			else if (elem2proc[*n3] != (id - 1))
			{
				//sd-sub domain
				*sd = *sd + 1;
				if (*sd == 1)
				{
					*kk = *kk + 1;
					id_proc[*kk] = elem2proc[*n3];
					d_store = id_proc[*kk];
					k2id[d_store] = *kk;
					sd_cnt[*kk] = 1;
					sprintf(filename, "sd_boun_%d.txt", id);
					fp = fopen(filename, "a");
					fprintf(fp, " %d %d %d\n", *i, *n3, elem2proc[*n3]);
					fclose(fp);
				} //if sd==1 loop
				else
				{
					check = 0;
					for (*j = 1; *j <= *kk; *j++)
					{
						if (elem2proc[*n3] == id_proc[*j])
						{
							check = 1;
							sd_cnt[*j] = sd_cnt[*j] + 1;
							sprintf(filename, "sd_boun_%d.txt", id);
							fp = fopen(filename, "a");
							fprintf(fp, " %d %d %d\n", *i, *n3, elem2proc[*n3]);
							fclose(fp);
							break;
						}
					}
					if (check == 0)
					{
						*kk = *kk + 1;
						id_proc[*kk] = elem2proc[*n3];
						sd_cnt[*kk] = 1;
						d_store = id_proc[*kk];
						k2id[d_store] = *kk;
						sprintf(filename, "sd_boun_%d.txt", id);
						fp = fopen(filename, "a");
						fprintf(fp, " %d %d %d\n", *i, *n3, elem2proc[*n3]);
						fclose(fp);
					}
				}
			}

			if (n4 == 0)
			{
				//gb-global boundary
				dummy = dummy + 1;
			}
			else if (elem2proc[n4] != (id - 1))
			{
				//sd-sub domain
				*sd = *sd + 1;
				if (*sd == 1)
				{
					*kk = *kk + 1;
					id_proc[*kk] = elem2proc[n4];
					d_store = id_proc[*kk];
					k2id[d_store] = *kk;
					sd_cnt[*kk] = 1;
					sprintf(filename, "sd_boun_%d.txt", id);
					fp = fopen(filename, "a");
					fprintf(fp, " %d %d %d\n", *i, n4, elem2proc[n4]);
					fclose(fp);
				} //if sd==1 loop
				else
				{
					check = 0;
					for (*j = 1; *j <= *kk; *j++)
					{
						if (elem2proc[n4] == id_proc[*j])
						{
							check = 1;
							sd_cnt[*j] = sd_cnt[*j] + 1;
							sprintf(filename, "sd_boun_%d.txt", id);
							fp = fopen(filename, "a");
							fprintf(fp, " %d %d %d\n", *i, n4, elem2proc[n4]);
							fclose(fp);
							break;
						}
					}
					if (check == 0)
					{
						*kk = *kk + 1;
						id_proc[*kk] = elem2proc[n4];
						sd_cnt[*kk] = 1;
						d_store = id_proc[*kk];
						k2id[d_store] = *kk;
						sprintf(filename, "sd_boun_%d.txt", id);
						fp = fopen(filename, "a");
						fprintf(fp, " %d %d %d\n", *i, n4, elem2proc[n4]);
						fclose(fp);
					}
				}
			}
		}
	}
	fclose(fp1);
	return lele_cnt;
}

int* slave_subdomain_search(int kk, int sd, char filename[100], int id,
		int check, int d_ele, int lele_cnt, int d_store, int*** lelem_send,
		int* i, int* sd_cnt, int* send_cnt, int* sd_elem, FILE* fp, int* n1,
		int* n2, int* n3, int* id_proc, int* j, int* gb, int* temp,
		int* subdomain_cnt, int* lelem)
{

	*lelem_send = (int**) malloc((kk + 1) * sizeof(int*));
	if (*lelem_send == NULL)
		printf("memory failed\n");

	for (*i = 1; *i <= kk; *i++)
	{
		*lelem_send[*i] = (int*) malloc((sd_cnt[*i] + 1) * sizeof(int));
	}
	for (*i = 1; *i <= kk; *i++)
	{
		send_cnt[*i] = 0;
	}
	sd_elem = (int*) malloc((sd + 1) * sizeof(int));
	sd = 0;
	sprintf(filename, "sd_boun_%d.txt", id);
	fp = fopen(filename, "r");
	while (!feof(fp))
	{
		fscanf(fp, "%d %d %d\n", &*n1, &*n2, &*n3);
		for (*i = 1; *i <= kk; *i++)
		{
			if (id_proc[*i] == *n3)
			{
				if (send_cnt[*i] == 0)
				{
					send_cnt[*i] = 1;
					*lelem_send[*i][send_cnt[*i]] = *n1;
				}
				else
				{
					check = 0;
					for (*j = 1; *j <= send_cnt[*i]; *j++)
					{
						if (*n1 == *lelem_send[*i][*j])
						{
							check = 1;
							break;
						}
					}
					if (check == 0)
					{
						send_cnt[*i] = send_cnt[*i] + 1;
						*lelem_send[*i][send_cnt[*i]] = *n1;
					}
				}
			}
		}
	}
	fclose(fp);
	//converting global numbers of lelem_send to local numbers
	for (*i = 1; *i <= kk; *i++)
	{
		for (*j = 1; *j <= send_cnt[*i]; *j++)
		{
			d_ele = *lelem_send[*i][*j];
			*lelem_send[*i][*j] = gb[d_ele];
		}
	}
	//to find the actual number of subdomains
	sprintf(filename, "sd_boun_%d.txt", id);
	fp = fopen(filename, "r");
	while (!feof(fp))
	{
		fscanf(fp, "%d %d %d\n", &*n1, &*n2, &*n3);
		if (sd == 0)
		{
			sd = 1;
			sd_elem[sd] = *n2;
		}
		else
		{
			check = 0;
			for (*j = 1; *j <= sd; *j++)
			{
				if (*n2 == sd_elem[*j])
				{
					check = 1;
					break;
				}
			}
			if (check == 0)
			{
				sd = sd + 1;
				sd_elem[sd] = *n2;
			}
		}
	}
	fclose(fp);
	int* group;
	group = (int*) malloc(sizeof(int) * (kk + 1));
	for (*i = 1; *i <= kk; *i++)
	{
		group[*i] = id_proc[*i];
	}
	//caution id_proc is 1 less than neigh id since we compare as id-1
	//arranging processor rank in ascending order
	for (*i = 0; *i < (kk - 1); ++*i)
	{
		for (*j = 1; *j < kk - *i; ++*j)
		{
			if (group[*j] > group[*j + 1])
			{
				*temp = group[*j + 1];
				group[*j + 1] = group[*j];
				group[*j] = *temp;
			}
		}
	}
	//subdomain boundary
	for (*i = 1; *i <= sd; *i++)
	{
		*subdomain_cnt = lele_cnt + *i;
		lelem[*subdomain_cnt] = sd_elem[*i];
		d_store = sd_elem[*i];
		gb[d_store] = *subdomain_cnt;
	}
	return group;
}

int sd_boun_data_find1(int i, int bc_set, int belem_cnt[bc_set + 1],
		int bcnt_iden[bc_set + 1], char filename[100], int id,
		char oneword[1000], int check, FILE* fp, int* d1, int* j, int* d2,
		int* temp, int* elem2proc, int* b_cnt)
{
	for (i = 1; i <= bc_set; i++)
	{
		belem_cnt[i] = 0;
		bcnt_iden[i] = 0;
	}
	sprintf(filename, "boundary_data_%d.txt", id);
	fp = fopen(filename, "r");
	if (access(filename, F_OK) == -1)
	{
		printf("unable to open boundary file\n");
	}
	for (i = 1; i <= bc_set; i++)
	{
		fscanf(fp, "%s %d\n", oneword, &*d1);
		check = 0;
		for (*j = 1; *j <= *d1; *j++)
		{
			fscanf(fp, "%d %d\n", &*d2, &*temp);
			if (elem2proc[*d2] == (id - 1))
			{
				if (check == 0)
				{
					if (belem_cnt[*b_cnt + 1] == 0)
					{
						*b_cnt = *b_cnt + 1;
						bcnt_iden[i] = 1;
						check = 1;
					}
				}
				belem_cnt[*b_cnt] = belem_cnt[*b_cnt] + 1;
			}
		}
	}
	fclose(fp);
	return i;
}

int sd_boun_data_find2(char filename[100], int id, int bc_set,
		int bcnt_iden[bc_set + 1], char oneword[1000], int temp,
		int belem_cnt[bc_set + 1], int subdomain_cnt, int pg, FILE* fp, int* i,
		int* d1, int* j, int* k, int* b_cnt, int* cnt_boun_type,
		int** boun_elem, int** boun_face, int* boun_type, int* elem2proc,
		int* gb, int* tot_ele, int* n1)
{
	sprintf(filename, "boundary_data_%d.txt", id);
	fp = fopen(filename, "r");
	if (access(filename, F_OK) == -1)
	{
		printf("unable to open boundary file2\n");
	}
	for (*i = 1; *i <= bc_set; *i++)
	{
		if (bcnt_iden[*i] == 0)
		{
			fscanf(fp, "%s %d\n", oneword, &*d1);
			for (*j = 1; *j <= *d1; *j++)
			{
				fscanf(fp, "%d %d\n", &temp, &temp);
			}
		}
		else
		{
			*k = 0;
			*b_cnt = *b_cnt + 1;
			fscanf(fp, "%s %d\n", oneword, &cnt_boun_type[*b_cnt]);
			boun_elem[*b_cnt] = (int*) malloc(
					(belem_cnt[*b_cnt] + 1) * sizeof(int));
			boun_face[*b_cnt] = (int*) malloc(
					(belem_cnt[*b_cnt] + 1) * sizeof(int));
			if (strncmp(oneword, "curved_wall", 11) == 0)
			{
				boun_type[*b_cnt] = 1;
			}
			else if (strncmp(oneword, "curved_inflow", 13) == 0)
			{
				boun_type[*b_cnt] = 2;
			}
			else if (strncmp(oneword, "curved_outflow", 14) == 0)
			{
				boun_type[*b_cnt] = 3;
			}
			else if (strncmp(oneword, "curved_upper", 12) == 0)
			{
				boun_type[*b_cnt] = 4;
			}
			else if (strncmp(oneword, "flat_wall_nx", 12) == 0)
			{
				boun_type[*b_cnt] = 5;
			}
			else if (strncmp(oneword, "flat_wall_ny", 12) == 0)
			{
				boun_type[*b_cnt] = 6;
			}
			else if (strncmp(oneword, "flat_wall_nz", 12) == 0)
			{
				boun_type[*b_cnt] = 7;
			}
			else if (strncmp(oneword, "flat_inflow_nx", 14) == 0)
			{
				boun_type[*b_cnt] = 8;
			}
			else if (strncmp(oneword, "flat_inflow_ny", 14) == 0)
			{
				boun_type[*b_cnt] = 9;
			}
			else if (strncmp(oneword, "flat_inflow_nz", 14) == 0)
			{
				boun_type[*b_cnt] = 10;
			}
			else if (strncmp(oneword, "flat_outflow_nx", 15) == 0)
			{
				boun_type[*b_cnt] = 11;
			}
			else if (strncmp(oneword, "flat_outflow_ny", 15) == 0)
			{
				boun_type[*b_cnt] = 12;
			}
			else if (strncmp(oneword, "flat_outflow_nz", 15) == 0)
			{
				boun_type[*b_cnt] = 13;
			}
			else if (strncmp(oneword, "flat_upper_nx", 13) == 0)
			{
				boun_type[*b_cnt] = 14;
			}
			else if (strncmp(oneword, "flat_upper_ny", 13) == 0)
			{
				boun_type[*b_cnt] = 15;
			}
			else if (strncmp(oneword, "flat_upper_nz", 13) == 0)
			{
				boun_type[*b_cnt] = 16;
			}
			else if (strncmp(oneword, "periodic", 8) == 0)
			{
				boun_type[*b_cnt] = 17;
			}

			for (*j = 1; *j <= cnt_boun_type[*b_cnt]; *j++)
			{
				fscanf(fp, "%d %d\n", &temp, &*d1);
				if (elem2proc[temp] == (id - 1))
				{
					*k = *k + 1;
					boun_face[*b_cnt][*k] = *d1;
					boun_elem[*b_cnt][*k] = gb[temp];
				}
			}
			if (*k == belem_cnt[*b_cnt])
			{
				cnt_boun_type[*b_cnt] = *k;
			}
			else
			{
				printf(
						"error in calculating local boundary elements in processor no. %d\n",
						id);
			}
		}
	}
	fclose(fp);
	temp = 0;
	for (*i = 1; *i <= *b_cnt; *i++)
	{
		temp = temp + cnt_boun_type[*i];
	}
	//tot_ele gives the total count of elements in a particular node
	*tot_ele = subdomain_cnt + temp + 1;
	int l_periodic_cnt;
	l_periodic_cnt = 0;
	if (pg == 1)
	{
		sprintf(filename, "periodic_elements_%d.txt", id);
		fp = fopen(filename, "r");
		if (access(filename, F_OK) == -1)
		{
			printf("unable to open periodic1 file\n");
		}
		while (!feof(fp))
		{
			fscanf(fp, "%d \n", &*n1);
			if (elem2proc[*n1] == (id - 1))
			{
				l_periodic_cnt = l_periodic_cnt + 1;
			}
		}
		fclose(fp);
	}
	return l_periodic_cnt;
}

int slave_node_local_numbering(char filename[100], int id, int i, int lele_cnt,
		int r[5], int check, FILE* fp, int* j, int* lelem, int* lnode_cnt,
		int* lnode, int* gb_node, int* k, int* n1)
{
	sprintf(filename, "conn_elem_%d.dat", id);
	fp = fopen(filename, "r");
	//fp=fopen("elem_conn.dat","r");
	for (i = 1; i <= lele_cnt; i++)
	{
		*j = lelem[i];
		fseek(fp, (*j - 1) * sizeof(int) * 4, SEEK_SET);
		fread(&r[1], sizeof(int), 1, fp);
		fread(&r[2], sizeof(int), 1, fp);
		fread(&r[3], sizeof(int), 1, fp);
		fread(&r[4], sizeof(int), 1, fp);
		if (*lnode_cnt == 0)
		{
			*lnode_cnt = 1;
			lnode[*lnode_cnt] = r[1];
			gb_node[r[1]] = *lnode_cnt;
			*lnode_cnt = *lnode_cnt + 1;
			lnode[*lnode_cnt] = r[2];
			gb_node[r[2]] = *lnode_cnt;
			*lnode_cnt = *lnode_cnt + 1;
			lnode[*lnode_cnt] = r[3];
			gb_node[r[3]] = *lnode_cnt;
			*lnode_cnt = *lnode_cnt + 1;
			lnode[*lnode_cnt] = r[4];
			gb_node[r[4]] = *lnode_cnt;
		}
		else
		{
			for (*k = 1; *k <= 4; *k++)
			{
				check = 0;
				for (*n1 = 1; *n1 <= *lnode_cnt; *n1++)
				{
					if (r[*k] == lnode[*n1])
					{
						check = 1;
						break;
					}
				}
				if (check == 0)
				{
					*lnode_cnt = *lnode_cnt + 1;
					lnode[*lnode_cnt] = r[*k];
					gb_node[r[*k]] = *lnode_cnt;
				}
			}
		}
	}
	fclose(fp);
	return i;
}

int slave_global_2_local_no(int i, int lele_cnt, int d_store, int* j,
		int** lelem_neighbour, int* gb)
{
	//converting global numbers into local numbering in lelem_neighbour
	for (i = 1; i <= lele_cnt; i++)
	{

		for (*j = 1; *j <= 4; *j++)
		{
			d_store = lelem_neighbour[*j][i];
			if (d_store == 0)
			{
				lelem_neighbour[*j][i] = d_store;
			}
			else
			{
				lelem_neighbour[*j][i] = gb[d_store];
			}
		}
	}
	return i;
}

int periodic_proc_rank_arrangement(int periodic_id, int i, int pp, int temp,
		int* p_group, int* id_periodic_proc, int* j)
{
	if (periodic_id == 1)
	{
		for (i = 1; i <= pp; i++)
		{
			p_group[i] = id_periodic_proc[i];
		}
		//caution id_proc is 1 less than neigh id since we compare as id-1
		//arranging processor rank in ascending order
		for (i = 0; i < (pp - 1); ++i)
		{
			for (*j = 1; *j < pp - i; ++*j)
			{
				if (p_group[*j] > p_group[*j + 1])
				{
					temp = p_group[*j + 1];
					p_group[*j + 1] = p_group[*j];
					p_group[*j] = temp;
				}
			}
		}
	}
	return i;
}

int slave_periodic_boun_iden(char filename[100], int id, int i,
		int l_periodic_cnt, int subdomain_cnt, int periodic_proc, int d_store,
		int check, FILE* fp, int* n1, int* elem2proc, int* g_period,
		int* l_period, int* gb, FILE* fp1, int* d1, int* l_periodic_mapping,
		int* j, int** lelem_neighbour, int* nele, int* lelem, int* pp,
		int* id_periodic_proc, int* pp2id, int* pp_cnt)
{
	sprintf(filename, "periodic_elements_%d.txt", id);
	fp = fopen(filename, "r");
	if (access(filename, F_OK) == -1)
	{
		printf("unable to open periodic elements file %d\n");
	}
	i = 0;
	while (!feof(fp))
	{
		fscanf(fp, "%d \n", &*n1);
		if (elem2proc[*n1] == (id - 1))
		{
			i = i + 1;
			//g_period:gives global number of periodic element. l_period:gives the local number of periodic element
			g_period[i] = *n1;
			l_period[i] = gb[*n1];
		}
	}
	fclose(fp);
	int check1;
	sprintf(filename, "periodic_ele_pairing_%d.dat", id);
	fp1 = fopen(filename, "r");
	if (access(filename, F_OK) == -1)
	{
		printf("unable to open periodic element pairing file\n");
	}
	for (i = 1; i <= l_periodic_cnt; i++)
	{
		check1 = 0;
		*d1 = g_period[i];
		fseek(fp1, (*d1 - 1) * sizeof(int), SEEK_SET);
		fread(&*n1, sizeof(int), 1, fp1);
		if (elem2proc[*n1] == (id - 1))
		{
			l_periodic_mapping[gb[*d1]] = gb[*n1];
			check1 = 1;
		}
		if (check1 == 0)
		{
			if (gb[*n1] != 0)
			{
				for (*j = 1; *j <= subdomain_cnt; *j++)
				{
					if ((gb[*n1] == lelem_neighbour[1][*j])
							|| (gb[*n1] == lelem_neighbour[2][*j])
							|| (gb[*n1] == lelem_neighbour[3][*j])
							|| (gb[*n1] == lelem_neighbour[4][*j]))
					{
						l_periodic_mapping[gb[*d1]] = gb[*n1];
						check1 = 1;
						break;
					}
				}
			}
		}
		if (check1 == 0)
		{
			*nele = *nele + 1;
			lelem[*nele] = *n1;
			gb[*n1] = *nele;
			l_periodic_mapping[gb[*d1]] = *nele;
			periodic_proc = periodic_proc + 1;
			if (periodic_proc == 1)
			{
				//pp is similar to kk. pp is specified, id_periodic_proc gives its respective processor
				*pp = *pp + 1;
				id_periodic_proc[*pp] = elem2proc[*n1];
				d_store = id_periodic_proc[*pp];
				//processor id is given, it gives pp
				pp2id[d_store] = *pp;
				//pp_cnt: No of elements in a particular "pp"
				pp_cnt[*pp] = 1;
				sprintf(filename, "periodic_boun_%d.txt", id);
				fp = fopen(filename, "a");
				fprintf(fp, "%d %d %d\n", *d1, *n1, elem2proc[*n1]);
				fclose(fp);
			}
			else
			{
				check = 0;
				for (*j = 1; *j <= *pp; *j++)
				{
					if (elem2proc[*n1] == id_periodic_proc[*j])
					{
						check = 1;
						pp_cnt[*j] = pp_cnt[*j] + 1;
						sprintf(filename, "periodic_boun_%d.txt", id);
						fp = fopen(filename, "a");
						fprintf(fp, "%d %d %d\n", *d1, *n1, elem2proc[*n1]);
						fclose(fp);
						break;
					}
				}
				if (check == 0)
				{
					*pp = *pp + 1;
					id_periodic_proc[*pp] = elem2proc[*n1];
					pp_cnt[*pp] = 1;
					d_store = id_periodic_proc[*pp];
					pp2id[d_store] = *pp;
					sprintf(filename, "periodic_boun_%d.txt", id);
					fp = fopen(filename, "a");
					fprintf(fp, "%d %d %d\n", *d1, *n1, elem2proc[*n1]);
					fclose(fp);
				}
			}
		}
	}
	fclose(fp1);
	return i;
}

#endif /* SLAVE_PREPROC_H_ */

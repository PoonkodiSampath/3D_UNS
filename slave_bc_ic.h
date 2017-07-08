/*
 * slave_bc_ic.h
 *
 *  Created on: 12-Feb-2017
 *      Author: poonkodi
 */

#ifndef SLAVE_BC_IC_H_
#define SLAVE_BC_IC_H_

//to calculate normals

void normal_eqn(double *xstart_d1, double *ystart_d1, double *zstart_d1,
		double *xstart_d2, double *ystart_d2, double *zstart_d2,
		double *xstart_d3, double *ystart_d3, double *zstart_d3, double *xc,
		double *yc, double *zc, double *nx, double *ny, double *nz,
		double *deno, double *z1, double *z2, double *z3, double *sum,
		double *a, double *b, double *c, double *d, double *dx, double *dy,
		double *dz, double *mx, double *my, double *mz)
{
	*nx = ((*ystart_d2 - *ystart_d1) * (*zstart_d3 - *zstart_d1))
			- ((*ystart_d3 - *ystart_d1) * (*zstart_d2 - *zstart_d1));

	*ny = -(((*zstart_d3 - *zstart_d1) * (*xstart_d2 - *xstart_d1))
			- ((*xstart_d3 - *xstart_d1) * (*zstart_d2 - *zstart_d1)));

	*nz = ((*xstart_d2 - *xstart_d1) * (*ystart_d3 - *ystart_d1))
			- ((*xstart_d3 - *xstart_d1) * (*ystart_d2 - *ystart_d1));

	*a = *nx;

	*b = *ny;

	*c = *nz;

	*d = -*a * *xstart_d1 - *b * *ystart_d1 - *c * *zstart_d1;

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

	*deno = sqrt(
			pow((*xstart_d2 - *xstart_d1), 2)
					+ pow((*ystart_d2 - *ystart_d1), 2)
					+ pow((*zstart_d2 - *zstart_d1), 2));

	*dx = (*xstart_d2 - *xstart_d1) / (*deno);
	*dy = (*ystart_d2 - *ystart_d1) / (*deno);
	*dz = (*zstart_d2 - *zstart_d1) / (*deno);

	*mx = *ny * *dz - *nz * *dy;
	*my = *nz * *dx - *nx * *dz;
	*mz = *nx * *dy - *ny * *dx;

}

//Transformation matrix

void trans_matrix(double *xstart_d1, double *ystart_d1, double *zstart_d1,
		double *nx, double *ny, double *nz, double *dx, double *dy, double *dz,
		double *mx, double *my, double *mz, double *deno, double *z1,
		double *z2, double *z3, double *xc, double *yc, double *zc, double *xp,
		double *yp, double *zp)
{
	*deno = sqrt(
			pow((*xp - *xc), 2) + pow((*yp - *yc), 2) + pow((*zp - *zc), 2));

//z1,z2,z3 r tramsformed coordinates

	*z1 = *dx * (*xp - *xstart_d1) + *dy * (*yp - *ystart_d1)
			+ *dz * (*zp - *zstart_d1);
	*z2 = *mx * (*xp - *xstart_d1) + *my * (*yp - *ystart_d1)
			+ *mz * (*zp - *zstart_d1);
	*z3 = *nx * (*xp - *xstart_d1) + *ny * (*yp - *ystart_d1)
			+ *nz * (*zp - *zstart_d1);

	*z3 = *z3 + *deno;

	*xp = *dx * *z1 + *mx * *z2 + *nx * *z3 + *xstart_d1;
	*yp = *dy * *z1 + *my * *z2 + *ny * *z3 + *ystart_d1;
	*zp = *dz * *z1 + *mz * *z2 + *nz * *z3 + *zstart_d1;

}

int boun_ele_calc(int i, int b_cnt, double** ghost_elem, int* cnt_boun_type,
		int* boun_type, int* j, int* iele, int** boun_elem, int* iface,
		int** boun_face, int* nele, int** lelem_neighbour, int* d1,
		struct node2* eleminfo, int* d2, int* d3, struct node1* xycoord,
		struct node3* centroid, double** nx_boun_face, double** ny_boun_face,
		double** nz_boun_face, double* deno, double* z1, double* z2, double* z3,
		double* sum, double** aa, double** bb, double** cc, double** dd,
		double** trans_coord_dx, double** trans_coord_dy,
		double** trans_coord_dz, double** trans_coord_mx,
		double** trans_coord_my, double** trans_coord_mz, double* t, double* xp,
		double* yp, double* zp, int** periodic_elem)
{
	//boundary element calculation
	for (i = 1; i <= b_cnt; i++)
	{
		ghost_elem[i] = (double*) malloc(
				(cnt_boun_type[i] + 1) * sizeof(double));
		if (ghost_elem[i] == NULL)
			printf("ghost_elem[i]_slave memory failed4\n");

		switch (boun_type[i])
		{
		//1: curved wall
		case 1:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*iface = boun_face[i][*j];
				*nele = *nele + 1;
				ghost_elem[i][*j] = *nele;
				lelem_neighbour[*iface][*iele] = *nele;
				*d1 = eleminfo[*iele].node1_localface[*iface];
				*d2 = eleminfo[*iele].node2_localface[*iface];
				*d3 = eleminfo[*iele].node3_localface[*iface];
				normal_eqn(&xycoord[*d1].xstart, &xycoord[*d1].ystart,
						&xycoord[*d1].zstart, &xycoord[*d2].xstart,
						&xycoord[*d2].ystart, &xycoord[*d2].zstart,
						&xycoord[*d3].xstart, &xycoord[*d3].ystart,
						&xycoord[*d3].zstart, &centroid[*iele].xc,
						&centroid[*iele].yc, &centroid[*iele].zc,
						&nx_boun_face[i][*j], &ny_boun_face[i][*j],
						&nz_boun_face[i][*j], deno, z1, z2, z3, sum, &aa[i][*j],
						&bb[i][*j], &cc[i][*j], &dd[i][*j],
						&trans_coord_dx[i][*j], &trans_coord_dy[i][*j],
						&trans_coord_dz[i][*j], &trans_coord_mx[i][*j],
						&trans_coord_my[i][*j], &trans_coord_mz[i][*j]);
				*t = (-dd[i][*j] - aa[i][*j] * centroid[*iele].xc
						- bb[i][*j] * centroid[*iele].yc
						- cc[i][*j] * centroid[*iele].zc)
						/ (aa[i][*j] * nx_boun_face[i][*j]
								+ bb[i][*j] * ny_boun_face[i][*j]
								+ cc[i][*j] * nz_boun_face[i][*j]);
				centroid[*nele].xc = centroid[*iele].xc
						+ (nx_boun_face[i][*j] * *t);
				centroid[*nele].yc = centroid[*iele].yc
						+ (ny_boun_face[i][*j] * *t);
				centroid[*nele].zc = centroid[*iele].zc
						+ (nz_boun_face[i][*j] * *t);
			}
			break;
			//2: curved_inflow
		case 2:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*iface = boun_face[i][*j];
				*nele = *nele + 1;
				ghost_elem[i][*j] = *nele;
				lelem_neighbour[*iface][*iele] = *nele;
				*d1 = eleminfo[*iele].node1_localface[*iface];
				*d2 = eleminfo[*iele].node2_localface[*iface];
				*d3 = eleminfo[*iele].node3_localface[*iface];
				normal_eqn(&xycoord[*d1].xstart, &xycoord[*d1].ystart,
						&xycoord[*d1].zstart, &xycoord[*d2].xstart,
						&xycoord[*d2].ystart, &xycoord[*d2].zstart,
						&xycoord[*d3].xstart, &xycoord[*d3].ystart,
						&xycoord[*d3].zstart, &centroid[*iele].xc,
						&centroid[*iele].yc, &centroid[*iele].zc,
						&nx_boun_face[i][*j], &ny_boun_face[i][*j],
						&nz_boun_face[i][*j], deno, z1, z2, z3, sum, &aa[i][*j],
						&bb[i][*j], &cc[i][*j], &dd[i][*j],
						&trans_coord_dx[i][*j], &trans_coord_dy[i][*j],
						&trans_coord_dz[i][*j], &trans_coord_mx[i][*j],
						&trans_coord_my[i][*j], &trans_coord_mz[i][*j]);
				*t = (-dd[i][*j] - aa[i][*j] * centroid[*iele].xc
						- bb[i][*j] * centroid[*iele].yc
						- cc[i][*j] * centroid[*iele].zc)
						/ (aa[i][*j] * nx_boun_face[i][*j]
								+ bb[i][*j] * ny_boun_face[i][*j]
								+ cc[i][*j] * nz_boun_face[i][*j]);
				*xp = centroid[*iele].xc + (nx_boun_face[i][*j] * *t);
				*yp = centroid[*iele].yc + (ny_boun_face[i][*j] * *t);
				*zp = centroid[*iele].zc + (nz_boun_face[i][*j] * *t);
				trans_matrix(&xycoord[*d1].xstart, &xycoord[*d1].ystart,
						&xycoord[*d1].zstart, &nx_boun_face[i][*j],
						&ny_boun_face[i][*j], &nz_boun_face[i][*j],
						&trans_coord_dx[i][*j], &trans_coord_dy[i][*j],
						&trans_coord_dz[i][*j], &trans_coord_mx[i][*j],
						&trans_coord_my[i][*j], &trans_coord_mz[i][*j], deno,
						z1, z2, z3, &centroid[*iele].xc, &centroid[*iele].yc,
						&centroid[*iele].zc, xp, yp, zp);
				centroid[*nele].xc = *xp;
				centroid[*nele].yc = *yp;
				centroid[*nele].zc = *zp;
			}
			break;
			//3: curved_outflow
		case 3:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*iface = boun_face[i][*j];
				*nele = *nele + 1;
				ghost_elem[i][*j] = *nele;
				lelem_neighbour[*iface][*iele] = *nele;
				*d1 = eleminfo[*iele].node1_localface[*iface];
				*d2 = eleminfo[*iele].node2_localface[*iface];
				*d3 = eleminfo[*iele].node3_localface[*iface];
				normal_eqn(&xycoord[*d1].xstart, &xycoord[*d1].ystart,
						&xycoord[*d1].zstart, &xycoord[*d2].xstart,
						&xycoord[*d2].ystart, &xycoord[*d2].zstart,
						&xycoord[*d3].xstart, &xycoord[*d3].ystart,
						&xycoord[*d3].zstart, &centroid[*iele].xc,
						&centroid[*iele].yc, &centroid[*iele].zc,
						&nx_boun_face[i][*j], &ny_boun_face[i][*j],
						&nz_boun_face[i][*j], deno, z1, z2, z3, sum, &aa[i][*j],
						&bb[i][*j], &cc[i][*j], &dd[i][*j],
						&trans_coord_dx[i][*j], &trans_coord_dy[i][*j],
						&trans_coord_dz[i][*j], &trans_coord_mx[i][*j],
						&trans_coord_my[i][*j], &trans_coord_mz[i][*j]);
				*t = (-dd[i][*j] - aa[i][*j] * centroid[*iele].xc
						- bb[i][*j] * centroid[*iele].yc
						- cc[i][*j] * centroid[*iele].zc)
						/ (aa[i][*j] * nx_boun_face[i][*j]
								+ bb[i][*j] * ny_boun_face[i][*j]
								+ cc[i][*j] * nz_boun_face[i][*j]);
				*xp = centroid[*iele].xc + (nx_boun_face[i][*j] * *t);
				*yp = centroid[*iele].yc + (ny_boun_face[i][*j] * *t);
				*zp = centroid[*iele].zc + (nz_boun_face[i][*j] * *t);
				trans_matrix(&xycoord[*d1].xstart, &xycoord[*d1].ystart,
						&xycoord[*d1].zstart, &nx_boun_face[i][*j],
						&ny_boun_face[i][*j], &nz_boun_face[i][*j],
						&trans_coord_dx[i][*j], &trans_coord_dy[i][*j],
						&trans_coord_dz[i][*j], &trans_coord_mx[i][*j],
						&trans_coord_my[i][*j], &trans_coord_mz[i][*j], deno,
						z1, z2, z3, &centroid[*iele].xc, &centroid[*iele].yc,
						&centroid[*iele].zc, xp, yp, zp);
				centroid[*nele].xc = *xp;
				centroid[*nele].yc = *yp;
				centroid[*nele].zc = *zp;
			}
			break;
			//4: curved_upper
		case 4:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*iface = boun_face[i][*j];
				*nele = *nele + 1;
				ghost_elem[i][*j] = *nele;
				lelem_neighbour[*iface][*iele] = *nele;
				*d1 = eleminfo[*iele].node1_localface[*iface];
				*d2 = eleminfo[*iele].node2_localface[*iface];
				*d3 = eleminfo[*iele].node3_localface[*iface];
				normal_eqn(&xycoord[*d1].xstart, &xycoord[*d1].ystart,
						&xycoord[*d1].zstart, &xycoord[*d2].xstart,
						&xycoord[*d2].ystart, &xycoord[*d2].zstart,
						&xycoord[*d3].xstart, &xycoord[*d3].ystart,
						&xycoord[*d3].zstart, &centroid[*iele].xc,
						&centroid[*iele].yc, &centroid[*iele].zc,
						&nx_boun_face[i][*j], &ny_boun_face[i][*j],
						&nz_boun_face[i][*j], deno, z1, z2, z3, sum, &aa[i][*j],
						&bb[i][*j], &cc[i][*j], &dd[i][*j],
						&trans_coord_dx[i][*j], &trans_coord_dy[i][*j],
						&trans_coord_dz[i][*j], &trans_coord_mx[i][*j],
						&trans_coord_my[i][*j], &trans_coord_mz[i][*j]);
				*t = (-dd[i][*j] - aa[i][*j] * centroid[*iele].xc
						- bb[i][*j] * centroid[*iele].yc
						- cc[i][*j] * centroid[*iele].zc)
						/ (aa[i][*j] * nx_boun_face[i][*j]
								+ bb[i][*j] * ny_boun_face[i][*j]
								+ cc[i][*j] * nz_boun_face[i][*j]);
				*xp = centroid[*iele].xc + (nx_boun_face[i][*j] * *t);
				*yp = centroid[*iele].yc + (ny_boun_face[i][*j] * *t);
				*zp = centroid[*iele].zc + (nz_boun_face[i][*j] * *t);
				trans_matrix(&xycoord[*d1].xstart, &xycoord[*d1].ystart,
						&xycoord[*d1].zstart, &nx_boun_face[i][*j],
						&ny_boun_face[i][*j], &nz_boun_face[i][*j],
						&trans_coord_dx[i][*j], &trans_coord_dy[i][*j],
						&trans_coord_dz[i][*j], &trans_coord_mx[i][*j],
						&trans_coord_my[i][*j], &trans_coord_mz[i][*j], deno,
						z1, z2, z3, &centroid[*iele].xc, &centroid[*iele].yc,
						&centroid[*iele].zc, xp, yp, zp);
				centroid[*nele].xc = *xp;
				centroid[*nele].yc = *yp;
				centroid[*nele].zc = *zp;
			}
			break;
			//5:FLAT_WALL perpendicular to x direction
		case 5:
			*iele = boun_elem[i][1];
			*iface = boun_face[i][1];
			*d1 = eleminfo[*iele].node1_localface[*iface];
			*xp = xycoord[*d1].xstart;
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*iface = boun_face[i][*j];
				*nele = *nele + 1;
				ghost_elem[i][*j] = *nele;
				lelem_neighbour[*iface][*iele] = *nele;
				centroid[*nele].xc = *xp;
				centroid[*nele].yc = centroid[*iele].yc;
				centroid[*nele].zc = centroid[*iele].zc;
			}
			break;
			//6:FLAT_WALL perpendicular to y direction
		case 6:
			*iele = boun_elem[i][1];
			*iface = boun_face[i][1];
			*d1 = eleminfo[*iele].node1_localface[*iface];
			*yp = xycoord[*d1].ystart;
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*iface = boun_face[i][*j];
				*nele = *nele + 1;
				ghost_elem[i][*j] = *nele;
				lelem_neighbour[*iface][*iele] = *nele;
				centroid[*nele].xc = centroid[*iele].xc;
				centroid[*nele].yc = *yp;
				centroid[*nele].zc = centroid[*iele].zc;
			}
			break;
			//7:FLAT_WALL perpendicular to z direction
		case 7:
			*iele = boun_elem[i][1];
			*iface = boun_face[i][1];
			*d1 = eleminfo[*iele].node1_localface[*iface];
			*zp = xycoord[*d1].zstart;
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*iface = boun_face[i][*j];
				*nele = *nele + 1;
				ghost_elem[i][*j] = *nele;
				lelem_neighbour[*iface][*iele] = *nele;
				centroid[*nele].xc = centroid[*iele].xc;
				centroid[*nele].yc = centroid[*iele].yc;
				centroid[*nele].zc = *zp;
			}
			break;
			//8:FLAT_INFLOW perpendicular to x direction
		case 8:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*iface = boun_face[i][*j];
				*d1 = eleminfo[*iele].node1_localface[*iface];
				*xp = 2.0 * xycoord[*d1].xstart - centroid[*iele].xc;
				*nele = *nele + 1;
				ghost_elem[i][*j] = *nele;
				lelem_neighbour[*iface][*iele] = *nele;
				centroid[*nele].xc = *xp;
				centroid[*nele].yc = centroid[*iele].yc;
				centroid[*nele].zc = centroid[*iele].zc;
			}
			break;
			//9:FLAT_INFLOW perpendicular to y direction
		case 9:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*iface = boun_face[i][*j];
				*d1 = eleminfo[*iele].node1_localface[*iface];
				*yp = 2.0 * xycoord[*d1].ystart - centroid[*iele].yc;
				*nele = *nele + 1;
				ghost_elem[i][*j] = *nele;
				lelem_neighbour[*iface][*iele] = *nele;
				centroid[*nele].xc = centroid[*iele].xc;
				centroid[*nele].yc = *yp;
				centroid[*nele].zc = centroid[*iele].zc;
			}
			break;
			//10:FLAT_INFLOW perpendicular to z direction
		case 10:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*iface = boun_face[i][*j];
				*d1 = eleminfo[*iele].node1_localface[*iface];
				*zp = 2.0 * xycoord[*d1].zstart - centroid[*iele].zc;
				*nele = *nele + 1;
				ghost_elem[i][*j] = *nele;
				lelem_neighbour[*iface][*iele] = *nele;
				centroid[*nele].xc = centroid[*iele].xc;
				centroid[*nele].yc = centroid[*iele].yc;
				centroid[*nele].zc = *zp;
			}
			break;
			//11:FLAT_OUTFLOW perpendicular to x direction
		case 11:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*iface = boun_face[i][*j];
				*d1 = eleminfo[*iele].node1_localface[*iface];
				*xp = 2.0 * xycoord[*d1].xstart - centroid[*iele].xc;
				*nele = *nele + 1;
				ghost_elem[i][*j] = *nele;
				lelem_neighbour[*iface][*iele] = *nele;
				centroid[*nele].xc = *xp;
				centroid[*nele].yc = centroid[*iele].yc;
				centroid[*nele].zc = centroid[*iele].zc;
			}
			break;
			//12:FLAT_OUTFLOW perpendicular to y direction
		case 12:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*iface = boun_face[i][*j];
				*d1 = eleminfo[*iele].node1_localface[*iface];
				*yp = 2.0 * xycoord[*d1].ystart - centroid[*iele].yc;
				*nele = *nele + 1;
				ghost_elem[i][*j] = *nele;
				lelem_neighbour[*iface][*iele] = *nele;
				centroid[*nele].xc = centroid[*iele].xc;
				centroid[*nele].yc = *yp;
				centroid[*nele].zc = centroid[*iele].zc;
			}
			break;
			//13:FLAT_OUTFLOW perpendicular to z direction
		case 13:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*iface = boun_face[i][*j];
				*d1 = eleminfo[*iele].node1_localface[*iface];
				*zp = 2.0 * xycoord[*d1].zstart - centroid[*iele].zc;
				*nele = *nele + 1;
				ghost_elem[i][*j] = *nele;
				lelem_neighbour[*iface][*iele] = *nele;
				centroid[*nele].xc = centroid[*iele].xc;
				centroid[*nele].yc = centroid[*iele].yc;
				centroid[*nele].zc = *zp;
			}
			break;
			//14:FLAT_UPPER perpendicular to x direction
		case 14:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*iface = boun_face[i][*j];
				*d1 = eleminfo[*iele].node1_localface[*iface];
				*xp = 2.0 * xycoord[*d1].xstart - centroid[*iele].xc;
				*nele = *nele + 1;
				ghost_elem[i][*j] = *nele;
				lelem_neighbour[*iface][*iele] = *nele;
				centroid[*nele].xc = *xp;
				centroid[*nele].yc = centroid[*iele].yc;
				centroid[*nele].zc = centroid[*iele].zc;
			}
			break;
			//15:FLAT_UPPER perpendicular to y direction
		case 15:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*iface = boun_face[i][*j];
				*d1 = eleminfo[*iele].node1_localface[*iface];
				*yp = 2.0 * xycoord[*d1].ystart - centroid[*iele].yc;
				*nele = *nele + 1;
				ghost_elem[i][*j] = *nele;
				lelem_neighbour[*iface][*iele] = *nele;
				centroid[*nele].xc = centroid[*iele].xc;
				centroid[*nele].yc = *yp;
				centroid[*nele].zc = centroid[*iele].zc;
			}
			break;
			//16:FLAT_UPPER perpendicular to z direction
		case 16:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*iface = boun_face[i][*j];
				*d1 = eleminfo[*iele].node1_localface[*iface];
				*zp = 2.0 * xycoord[*d1].zstart - centroid[*iele].zc;
				*nele = *nele + 1;
				ghost_elem[i][*j] = *nele;
				lelem_neighbour[*iface][*iele] = *nele;
				centroid[*nele].xc = centroid[*iele].xc;
				centroid[*nele].yc = centroid[*iele].yc;
				centroid[*nele].zc = *zp;
			}
			break;
			//periodic
		case 17:
			periodic_elem[i] = (int*) malloc(
					(cnt_boun_type[i] + 1) * sizeof(int));
			//printf("cnt of periodic element is %d\n",cnt_boun_type[i]);
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*iface = boun_face[i][*j];
				*nele = *nele + 1;
				ghost_elem[i][*j] = *nele;
				//printf("periodic is %d\n",periodic_mapping[iele]);
				lelem_neighbour[*iface][*iele] = *nele;
				*d1 = eleminfo[*iele].node1_localface[*iface];
				*d2 = eleminfo[*iele].node2_localface[*iface];
				*d3 = eleminfo[*iele].node3_localface[*iface];
				normal_eqn(&xycoord[*d1].xstart, &xycoord[*d1].ystart,
						&xycoord[*d1].zstart, &xycoord[*d2].xstart,
						&xycoord[*d2].ystart, &xycoord[*d2].zstart,
						&xycoord[*d3].xstart, &xycoord[*d3].ystart,
						&xycoord[*d3].zstart, &centroid[*iele].xc,
						&centroid[*iele].yc, &centroid[*iele].zc,
						&nx_boun_face[i][*j], &ny_boun_face[i][*j],
						&nz_boun_face[i][*j], deno, z1, z2, z3, sum, &aa[i][*j],
						&bb[i][*j], &cc[i][*j], &dd[i][*j],
						&trans_coord_dx[i][*j], &trans_coord_dy[i][*j],
						&trans_coord_dz[i][*j], &trans_coord_mx[i][*j],
						&trans_coord_my[i][*j], &trans_coord_mz[i][*j]);
				*t = (-dd[i][*j] - aa[i][*j] * centroid[*iele].xc
						- bb[i][*j] * centroid[*iele].yc
						- cc[i][*j] * centroid[*iele].zc)
						/ (aa[i][*j] * nx_boun_face[i][*j]
								+ bb[i][*j] * ny_boun_face[i][*j]
								+ cc[i][*j] * nz_boun_face[i][*j]);
				*xp = centroid[*iele].xc + (nx_boun_face[i][*j] * *t);
				*yp = centroid[*iele].yc + (ny_boun_face[i][*j] * *t);
				*zp = centroid[*iele].zc + (nz_boun_face[i][*j] * *t);
				trans_matrix(&xycoord[*d1].xstart, &xycoord[*d1].ystart,
						&xycoord[*d1].zstart, &nx_boun_face[i][*j],
						&ny_boun_face[i][*j], &nz_boun_face[i][*j],
						&trans_coord_dx[i][*j], &trans_coord_dy[i][*j],
						&trans_coord_dz[i][*j], &trans_coord_mx[i][*j],
						&trans_coord_my[i][*j], &trans_coord_mz[i][*j], deno,
						z1, z2, z3, &centroid[*iele].xc, &centroid[*iele].yc,
						&centroid[*iele].zc, xp, yp, zp);
				centroid[*nele].xc = *xp;
				centroid[*nele].yc = *yp;
				centroid[*nele].zc = *zp;
			}
			break;
		}
	}
	return i;
}

int boun_sol_pt(int i, int b_cnt, int iface, int d2, int d3, int* boun_type,
		int* j, int* cnt_boun_type, int* iele, int** boun_elem, int* nele,
		double** ghost_elem, double* t, double** dd, double** aa,
		struct node4* sol, double** bb, double** cc, double** nx_boun_face,
		double** ny_boun_face, double** nz_boun_face, int** boun_face, int* d1,
		struct node2* eleminfo, double* xp, double* yp, double* zp,
		struct node1* xycoord, double** trans_coord_dx, double** trans_coord_dy,
		double** trans_coord_dz, double** trans_coord_mx,
		double** trans_coord_my, double** trans_coord_mz, double* deno,
		double* z1, double* z2, double* z3, int** periodic_elem,
		int* l_periodic_mapping)
{

	//assigning solution points for boundaries
	for (i = 1; i <= b_cnt; i++)
	{
		switch (boun_type[i])
		{
		//1: curved wall
		case 1:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*nele = ghost_elem[i][*j];
				*t = (-dd[i][*j] - aa[i][*j] * sol[*iele].spx
						- bb[i][*j] * sol[*iele].spy
						- cc[i][*j] * sol[*iele].spz)
						/ (aa[i][*j] * nx_boun_face[i][*j]
								+ bb[i][*j] * ny_boun_face[i][*j]
								+ cc[i][*j] * nz_boun_face[i][*j]);
				sol[*nele].spx = sol[*iele].spx + (nx_boun_face[i][*j] * *t);
				sol[*nele].spy = sol[*iele].spy + (ny_boun_face[i][*j] * *t);
				sol[*nele].spz = sol[*iele].spz + (nz_boun_face[i][*j] * *t);
			}
			break;
			//2: curved_inflow
		case 2:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				iface = boun_face[i][*j];
				*nele = ghost_elem[i][*j];
				*d1 = eleminfo[*iele].node1_localface[iface];
				d2 = eleminfo[*iele].node2_localface[iface];
				d3 = eleminfo[*iele].node3_localface[iface];
				*t = (-dd[i][*j] - aa[i][*j] * sol[*iele].spx
						- bb[i][*j] * sol[*iele].spy
						- cc[i][*j] * sol[*iele].spz)
						/ (aa[i][*j] * nx_boun_face[i][*j]
								+ bb[i][*j] * ny_boun_face[i][*j]
								+ cc[i][*j] * nz_boun_face[i][*j]);
				*xp = sol[*iele].spx + (nx_boun_face[i][*j] * *t);
				*yp = sol[*iele].spy + (ny_boun_face[i][*j] * *t);
				*zp = sol[*iele].spz + (nz_boun_face[i][*j] * *t);
				trans_matrix(&xycoord[*d1].xstart, &xycoord[*d1].ystart,
						&xycoord[*d1].zstart, &nx_boun_face[i][*j],
						&ny_boun_face[i][*j], &nz_boun_face[i][*j],
						&trans_coord_dx[i][*j], &trans_coord_dy[i][*j],
						&trans_coord_dz[i][*j], &trans_coord_mx[i][*j],
						&trans_coord_my[i][*j], &trans_coord_mz[i][*j], deno,
						z1, z2, z3, &sol[*iele].spx, &sol[*iele].spy,
						&sol[*iele].spz, xp, yp, zp);
				sol[*nele].spx = *xp;
				sol[*nele].spy = *yp;
				sol[*nele].spz = *zp;
			}
			break;
			//3: curved_outflow
		case 3:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				iface = boun_face[i][*j];
				*nele = ghost_elem[i][*j];
				*d1 = eleminfo[*iele].node1_localface[iface];
				d2 = eleminfo[*iele].node2_localface[iface];
				d3 = eleminfo[*iele].node3_localface[iface];
				*t = (-dd[i][*j] - aa[i][*j] * sol[*iele].spx
						- bb[i][*j] * sol[*iele].spy
						- cc[i][*j] * sol[*iele].spz)
						/ (aa[i][*j] * nx_boun_face[i][*j]
								+ bb[i][*j] * ny_boun_face[i][*j]
								+ cc[i][*j] * nz_boun_face[i][*j]);
				*xp = sol[*iele].spx + (nx_boun_face[i][*j] * *t);
				*yp = sol[*iele].spy + (ny_boun_face[i][*j] * *t);
				*zp = sol[*iele].spz + (nz_boun_face[i][*j] * *t);
				trans_matrix(&xycoord[*d1].xstart, &xycoord[*d1].ystart,
						&xycoord[*d1].zstart, &nx_boun_face[i][*j],
						&ny_boun_face[i][*j], &nz_boun_face[i][*j],
						&trans_coord_dx[i][*j], &trans_coord_dy[i][*j],
						&trans_coord_dz[i][*j], &trans_coord_mx[i][*j],
						&trans_coord_my[i][*j], &trans_coord_mz[i][*j], deno,
						z1, z2, z3, &sol[*iele].spx, &sol[*iele].spy,
						&sol[*iele].spz, xp, yp, zp);
				sol[*nele].spx = *xp;
				sol[*nele].spy = *yp;
				sol[*nele].spz = *zp;
			}
			break;
			//4: curved_upper
		case 4:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				iface = boun_face[i][*j];
				*nele = ghost_elem[i][*j];
				*d1 = eleminfo[*iele].node1_localface[iface];
				d2 = eleminfo[*iele].node2_localface[iface];
				d3 = eleminfo[*iele].node3_localface[iface];
				*t = (-dd[i][*j] - aa[i][*j] * sol[*iele].spx
						- bb[i][*j] * sol[*iele].spy
						- cc[i][*j] * sol[*iele].spz)
						/ (aa[i][*j] * nx_boun_face[i][*j]
								+ bb[i][*j] * ny_boun_face[i][*j]
								+ cc[i][*j] * nz_boun_face[i][*j]);
				*xp = sol[*iele].spx + (nx_boun_face[i][*j] * *t);
				*yp = sol[*iele].spy + (ny_boun_face[i][*j] * *t);
				*zp = sol[*iele].spz + (nz_boun_face[i][*j] * *t);
				trans_matrix(&xycoord[*d1].xstart, &xycoord[*d1].ystart,
						&xycoord[*d1].zstart, &nx_boun_face[i][*j],
						&ny_boun_face[i][*j], &nz_boun_face[i][*j],
						&trans_coord_dx[i][*j], &trans_coord_dy[i][*j],
						&trans_coord_dz[i][*j], &trans_coord_mx[i][*j],
						&trans_coord_my[i][*j], &trans_coord_mz[i][*j], deno,
						z1, z2, z3, &sol[*iele].spx, &sol[*iele].spy,
						&sol[*iele].spz, xp, yp, zp);
				sol[*nele].spx = *xp;
				sol[*nele].spy = *yp;
				sol[*nele].spz = *zp;
			}
			break;
			//5:FLAT_WALL perpendicular to x direction
		case 5:
			*iele = boun_elem[i][1];
			iface = boun_face[i][1];
			*d1 = eleminfo[*iele].node1_localface[iface];
			*xp = xycoord[*d1].xstart;
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*nele = ghost_elem[i][*j];
				sol[*nele].spx = *xp;
				sol[*nele].spy = sol[*iele].spy;
				sol[*nele].spz = sol[*iele].spz;
			}
			break;
			//6:FLAT_WALL perpendicular to y direction
		case 6:
			*iele = boun_elem[i][1];
			iface = boun_face[i][1];
			*d1 = eleminfo[*iele].node1_localface[iface];
			*yp = xycoord[*d1].ystart;
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*nele = ghost_elem[i][*j];
				sol[*nele].spx = sol[*iele].spx;
				sol[*nele].spy = *yp;
				sol[*nele].spz = sol[*iele].spz;
			}
			break;
			//7:FLAT_WALL perpendicular to z direction
		case 7:
			*iele = boun_elem[i][1];
			iface = boun_face[i][1];
			*d1 = eleminfo[*iele].node1_localface[iface];
			*zp = xycoord[*d1].zstart;
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*nele = ghost_elem[i][*j];
				sol[*nele].spx = sol[*iele].spx;
				sol[*nele].spy = sol[*iele].spy;
				sol[*nele].spz = *zp;
			}
			break;
			//8:FLAT_INFLOW perpendicular to x direction
		case 8:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*nele = ghost_elem[i][*j];
				iface = boun_face[i][*j];
				*d1 = eleminfo[*iele].node1_localface[iface];
				*xp = 2.0 * xycoord[*d1].xstart - sol[*iele].spx;
				sol[*nele].spx = *xp;
				sol[*nele].spy = sol[*iele].spy;
				sol[*nele].spz = sol[*iele].spz;
			}
			break;
			//9:FLAT_INFLOW perpendicular to y direction
		case 9:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*nele = ghost_elem[i][*j];
				iface = boun_face[i][*j];
				*d1 = eleminfo[*iele].node1_localface[iface];
				*yp = 2.0 * xycoord[*d1].ystart - sol[*iele].spy;
				sol[*nele].spx = sol[*iele].spx;
				sol[*nele].spy = *yp;
				sol[*nele].spz = sol[*iele].spz;
			}
			break;
			//10:FLAT_INFLOW perpendicular to z direction
		case 10:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*nele = ghost_elem[i][*j];
				iface = boun_face[i][*j];
				*d1 = eleminfo[*iele].node1_localface[iface];
				*zp = 2.0 * xycoord[*d1].zstart - sol[*iele].spz;
				sol[*nele].spx = sol[*iele].spx;
				sol[*nele].spy = sol[*iele].spy;
				sol[*nele].spz = *zp;
			}
			break;
			//11:FLAT_OUTFLOW perpendicular to x direction
		case 11:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*nele = ghost_elem[i][*j];
				iface = boun_face[i][*j];
				*d1 = eleminfo[*iele].node1_localface[iface];
				*xp = 2.0 * xycoord[*d1].xstart - sol[*iele].spx;
				sol[*nele].spx = *xp;
				sol[*nele].spy = sol[*iele].spy;
				sol[*nele].spz = sol[*iele].spz;
			}
			break;
			//12:FLAT_OUTFLOW perpendicular to y direction
		case 12:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*nele = ghost_elem[i][*j];
				iface = boun_face[i][*j];
				*d1 = eleminfo[*iele].node1_localface[iface];
				*yp = 2.0 * xycoord[*d1].ystart - sol[*iele].spy;
				sol[*nele].spx = sol[*iele].spx;
				sol[*nele].spy = *yp;
				sol[*nele].spz = sol[*iele].spz;
			}
			break;
			//13:FLAT_OUTFLOW perpendicular to z direction
		case 13:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*nele = ghost_elem[i][*j];
				iface = boun_face[i][*j];
				*d1 = eleminfo[*iele].node1_localface[iface];
				*zp = 2.0 * xycoord[*d1].zstart - sol[*iele].spz;
				sol[*nele].spx = sol[*iele].spx;
				sol[*nele].spy = sol[*iele].spy;
				sol[*nele].spz = *zp;
			}
			break;
			//14:FLAT_UPPER perpendicular to x direction
		case 14:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*nele = ghost_elem[i][*j];
				iface = boun_face[i][*j];
				*d1 = eleminfo[*iele].node1_localface[iface];
				*xp = 2.0 * xycoord[*d1].xstart - sol[*iele].spx;
				sol[*nele].spx = *xp;
				sol[*nele].spy = sol[*iele].spy;
				sol[*nele].spz = sol[*iele].spz;
			}
			break;
			//15:FLAT_UPPER perpendicular to y direction
		case 15:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*nele = ghost_elem[i][*j];
				iface = boun_face[i][*j];
				*d1 = eleminfo[*iele].node1_localface[iface];
				*yp = 2.0 * xycoord[*d1].ystart - sol[*iele].spy;
				sol[*nele].spx = sol[*iele].spx;
				sol[*nele].spy = *yp;
				sol[*nele].spz = sol[*iele].spz;
			}
			break;
			//16:FLAT_INFLOW perpendicular to z direction
		case 16:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				*nele = ghost_elem[i][*j];
				iface = boun_face[i][*j];
				*d1 = eleminfo[*iele].node1_localface[iface];
				*zp = 2.0 * xycoord[*d1].zstart - sol[*iele].spz;
				sol[*nele].spx = sol[*iele].spx;
				sol[*nele].spy = sol[*iele].spy;
				sol[*nele].spz = *zp;
			}
			break;
			//17:periodic
		case 17:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				iface = boun_face[i][*j];
				*nele = ghost_elem[i][*j];
				periodic_elem[i][*j] = l_periodic_mapping[*iele];
				*d1 = eleminfo[*iele].node1_localface[iface];
				d2 = eleminfo[*iele].node2_localface[iface];
				d3 = eleminfo[*iele].node3_localface[iface];
				*t = (-dd[i][*j] - aa[i][*j] * sol[*iele].spx
						- bb[i][*j] * sol[*iele].spy
						- cc[i][*j] * sol[*iele].spz)
						/ (aa[i][*j] * nx_boun_face[i][*j]
								+ bb[i][*j] * ny_boun_face[i][*j]
								+ cc[i][*j] * nz_boun_face[i][*j]);
				*xp = sol[*iele].spx + (nx_boun_face[i][*j] * *t);
				*yp = sol[*iele].spy + (ny_boun_face[i][*j] * *t);
				*zp = sol[*iele].spz + (nz_boun_face[i][*j] * *t);
				trans_matrix(&xycoord[*d1].xstart, &xycoord[*d1].ystart,
						&xycoord[*d1].zstart, &nx_boun_face[i][*j],
						&ny_boun_face[i][*j], &nz_boun_face[i][*j],
						&trans_coord_dx[i][*j], &trans_coord_dy[i][*j],
						&trans_coord_dz[i][*j], &trans_coord_mx[i][*j],
						&trans_coord_my[i][*j], &trans_coord_mz[i][*j], deno,
						z1, z2, z3, &sol[*iele].spx, &sol[*iele].spy,
						&sol[*iele].spz, xp, yp, zp);
				sol[*nele].spx = *xp;
				sol[*nele].spy = *yp;
				sol[*nele].spz = *zp;
			}
			break;
		}
	}
	return i;
}

void curved_wall(double *t_ux1, double *t_uy1, double *t_uz1, double *t_ux2,
		double *t_uy2, double *t_uz2, double *t_ux3, double *t_uy3,
		double *t_uz3, double *t_ux4, double *t_uy4, double *t_uz4,
		double *t_ux5, double *t_uy5, double *t_uz5, double *t_ux1_iele,
		double *t_uy1_iele, double *t_uz1_iele, double *t_ux2_iele,
		double *t_uy2_iele, double *t_uz2_iele, double *t_ux3_iele,
		double *t_uy3_iele, double *t_uz3_iele, double *t_ux4_iele,
		double *t_uy4_iele, double *t_uz4_iele, double *t_ux5_iele,
		double *t_uy5_iele, double *t_uz5_iele, double *dx, double *dy,
		double *dz, double *mx, double *my, double *mz, double *nx, double *ny,
		double *nz, double *u1_nele, double *u2_nele, double *u3_nele,
		double *u4_nele, double *u5_nele, double *u1_iele, double *u2_iele,
		double *u3_iele, double *u4_iele, double *u5_iele, double *ux1_iele,
		double *uy1_iele, double *uz1_iele, double *ux2_iele, double *uy2_iele,
		double *uz2_iele, double *ux3_iele, double *uy3_iele, double *uz3_iele,
		double *ux4_iele, double *uy4_iele, double *uz4_iele, double *ux5_iele,
		double *uy5_iele, double *uz5_iele)
{
	*u1_nele = *u1_iele;
	*u2_nele = 0.0;
	*u3_nele = 0.0;
	*u4_nele = 0.0;
	*u5_nele = *u5_iele
			- (0.5
					* (*u2_iele * *u2_iele + *u3_iele * *u3_iele
							+ *u4_iele * *u4_iele) / (*u1_iele));

	*t_ux1_iele = *dx * *ux1_iele + *dy * *uy1_iele + *dz * *uz1_iele;
	*t_uy1_iele = *mx * *ux1_iele + *my * *uy1_iele + *mz * *uz1_iele;
	*t_uz1_iele = *nx * *ux1_iele + *ny * *uy1_iele + *nz * *uz1_iele;
	*t_ux2_iele = *dx * *ux2_iele + *dy * *uy2_iele + *dz * *uz2_iele;
	*t_uy2_iele = *mx * *ux2_iele + *my * *uy2_iele + *mz * *uz2_iele;
	*t_uz2_iele = *nx * *ux2_iele + *ny * *uy2_iele + *nz * *uz2_iele;
	*t_ux3_iele = *dx * *ux3_iele + *dy * *uy3_iele + *dz * *uz3_iele;
	*t_uy3_iele = *mx * *ux3_iele + *my * *uy3_iele + *mz * *uz3_iele;
	*t_uz3_iele = *nx * *ux3_iele + *ny * *uy3_iele + *nz * *uz3_iele;
	*t_ux4_iele = *dx * *ux4_iele + *dy * *uy4_iele + *dz * *uz4_iele;
	*t_uy4_iele = *mx * *ux4_iele + *my * *uy4_iele + *mz * *uz4_iele;
	*t_uz4_iele = *nx * *ux4_iele + *ny * *uy4_iele + *nz * *uz4_iele;
	*t_ux5_iele = *dx * *ux5_iele + *dy * *uy5_iele + *dz * *uz5_iele;
	*t_uy5_iele = *mx * *ux5_iele + *my * *uy5_iele + *mz * *uz5_iele;
	*t_uz5_iele = *nx * *ux5_iele + *ny * *uy5_iele + *nz * *uz5_iele;

	*t_ux1 = *t_ux1_iele;
	*t_uy1 = *t_uy1_iele;
	*t_uz1 = 0.0;
	*t_ux2 = 0.0;
	*t_uy2 = 0.0;
	*t_uz2 = *t_uz2_iele;
	*t_ux3 = 0.0;
	*t_uy3 = 0.0;
	*t_uz3 = *t_uz3_iele;
	*t_ux4 = 0.0;
	*t_uy4 = 0.0;
	*t_uz4 = *t_uz4_iele;
	*t_ux5 = *t_ux5_iele;
	*t_uy5 = *t_uy5_iele;
	*t_uz5 = 0.0;
}

void derv_trans(double *ux1_nele, double *uy1_nele, double *uz1_nele,
		double *ux2_nele, double *uy2_nele, double *uz2_nele, double *ux3_nele,
		double *uy3_nele, double *uz3_nele, double *ux4_nele, double *uy4_nele,
		double *uz4_nele, double *ux5_nele, double *uy5_nele, double *uz5_nele,
		double *dx, double *dy, double *dz, double *mx, double *my, double *mz,
		double *nx, double *ny, double *nz, double *t_ux1, double *t_uy1,
		double *t_uz1, double *t_ux2, double *t_uy2, double *t_uz2,
		double *t_ux3, double *t_uy3, double *t_uz3, double *t_ux4,
		double *t_uy4, double *t_uz4, double *t_ux5, double *t_uy5,
		double *t_uz5)
{
	*ux1_nele = *dx * *t_ux1 + *mx * *t_uy1 + *nx * *t_uz1;
	*uy1_nele = *dy * *t_ux1 + *my * *t_uy1 + *ny * *t_uz1;
	*uz1_nele = *dz * *t_ux1 + *mz * *t_uy1 + *nz * *t_uz1;

	*ux2_nele = *dx * *t_ux2 + *mx * *t_uy2 + *nx * *t_uz2;
	*uy2_nele = *dy * *t_ux2 + *my * *t_uy2 + *ny * *t_uz2;
	*uz2_nele = *dz * *t_ux2 + *mz * *t_uy2 + *nz * *t_uz2;

	*ux3_nele = *dx * *t_ux3 + *mx * *t_uy3 + *nx * *t_uz3;
	*uy3_nele = *dy * *t_ux3 + *my * *t_uy3 + *ny * *t_uz3;
	*uz3_nele = *dz * *t_ux3 + *mz * *t_uy3 + *nz * *t_uz3;

	*ux4_nele = *dx * *t_ux4 + *mx * *t_uy4 + *nx * *t_uz4;
	*uy4_nele = *dy * *t_ux4 + *my * *t_uy4 + *ny * *t_uz4;
	*uz4_nele = *dz * *t_ux4 + *mz * *t_uy4 + *nz * *t_uz4;

	*ux5_nele = *dx * *t_ux5 + *mx * *t_uy5 + *nx * *t_uz5;
	*uy5_nele = *dy * *t_ux5 + *my * *t_uy5 + *ny * *t_uz5;
	*uz5_nele = *dz * *t_ux5 + *mz * *t_uy5 + *nz * *t_uz5;
}

void curved_inflow(double *u1_nele, double *u2_nele, double *u3_nele,
		double *u4_nele, double *u5_nele, double *ux1_nele, double *uy1_nele,
		double *uz1_nele, double *ux2_nele, double *uy2_nele, double *uz2_nele,
		double *ux3_nele, double *uy3_nele, double *uz3_nele, double *ux4_nele,
		double *uy4_nele, double *uz4_nele, double *ux5_nele, double *uy5_nele,
		double *uz5_nele, double *u1initial, double *u2initial,
		double *u3initial, double *u4initial, double *u5initial)
{
	*u1_nele = *u1initial;
	*u2_nele = *u2initial;
	*u3_nele = *u3initial;
	*u4_nele = *u4initial;
	*u5_nele = *u5initial;

	*ux1_nele = 0.0;
	*uy1_nele = 0.0;
	*uz1_nele = 0.0;

	*ux2_nele = 0.0;
	*uy2_nele = 0.0;
	*uz2_nele = 0.0;

	*ux3_nele = 0.0;
	*uy3_nele = 0.0;
	*uz3_nele = 0.0;

	*ux4_nele = 0.0;
	*uy4_nele = 0.0;
	*uz4_nele = 0.0;

	*ux5_nele = 0.0;
	*uy5_nele = 0.0;
	*uz5_nele = 0.0;
}

void curved_outflow(double *t_ux1, double *t_uy1, double *t_uz1, double *t_ux2,
		double *t_uy2, double *t_uz2, double *t_ux3, double *t_uy3,
		double *t_uz3, double *t_ux4, double *t_uy4, double *t_uz4,
		double *t_ux5, double *t_uy5, double *t_uz5, double *t_ux1_iele,
		double *t_uy1_iele, double *t_uz1_iele, double *t_ux2_iele,
		double *t_uy2_iele, double *t_uz2_iele, double *t_ux3_iele,
		double *t_uy3_iele, double *t_uz3_iele, double *t_ux4_iele,
		double *t_uy4_iele, double *t_uz4_iele, double *t_ux5_iele,
		double *t_uy5_iele, double *t_uz5_iele, double *dx, double *dy,
		double *dz, double *mx, double *my, double *mz, double *nx, double *ny,
		double *nz, double *u1_nele, double *u2_nele, double *u3_nele,
		double *u4_nele, double *u5_nele, double *u1_iele, double *u2_iele,
		double *u3_iele, double *u4_iele, double *u5_iele, double *ux1_iele,
		double *uy1_iele, double *uz1_iele, double *ux2_iele, double *uy2_iele,
		double *uz2_iele, double *ux3_iele, double *uy3_iele, double *uz3_iele,
		double *ux4_iele, double *uy4_iele, double *uz4_iele, double *ux5_iele,
		double *uy5_iele, double *uz5_iele)
{
	*u1_nele = *u1_iele;
	*u2_nele = *u2_iele;
	*u3_nele = *u3_iele;
	*u4_nele = *u4_iele;
	*u5_nele = *u5_iele;

	*t_ux1_iele = *dx * *ux1_iele + *dy * *uy1_iele + *dz * *uz1_iele;
	*t_uy1_iele = *mx * *ux1_iele + *my * *uy1_iele + *mz * *uz1_iele;
	*t_uz1_iele = *nx * *ux1_iele + *ny * *uy1_iele + *nz * *uz1_iele;
	*t_ux2_iele = *dx * *ux2_iele + *dy * *uy2_iele + *dz * *uz2_iele;
	*t_uy2_iele = *mx * *ux2_iele + *my * *uy2_iele + *mz * *uz2_iele;
	*t_uz2_iele = *nx * *ux2_iele + *ny * *uy2_iele + *nz * *uz2_iele;
	*t_ux3_iele = *dx * *ux3_iele + *dy * *uy3_iele + *dz * *uz3_iele;
	*t_uy3_iele = *mx * *ux3_iele + *my * *uy3_iele + *mz * *uz3_iele;
	*t_uz3_iele = *nx * *ux3_iele + *ny * *uy3_iele + *nz * *uz3_iele;
	*t_ux4_iele = *dx * *ux4_iele + *dy * *uy4_iele + *dz * *uz4_iele;
	*t_uy4_iele = *mx * *ux4_iele + *my * *uy4_iele + *mz * *uz4_iele;
	*t_uz4_iele = *nx * *ux4_iele + *ny * *uy4_iele + *nz * *uz4_iele;
	*t_ux5_iele = *dx * *ux5_iele + *dy * *uy5_iele + *dz * *uz5_iele;
	*t_uy5_iele = *mx * *ux5_iele + *my * *uy5_iele + *mz * *uz5_iele;
	*t_uz5_iele = *nx * *ux5_iele + *ny * *uy5_iele + *nz * *uz5_iele;

	*t_ux1 = *t_ux1_iele;
	*t_uy1 = *t_uy1_iele;
	*t_uz1 = 0.0;
	*t_ux2 = *t_ux1_iele;
	*t_uy2 = *t_uy1_iele;
	*t_uz2 = 0.0;
	*t_ux3 = *t_ux1_iele;
	*t_uy3 = *t_uy1_iele;
	*t_uz3 = 0.0;
	*t_ux4 = *t_ux1_iele;
	*t_uy4 = *t_uy1_iele;
	*t_uz4 = 0.0;
	*t_ux5 = *t_ux1_iele;
	*t_uy5 = *t_uy1_iele;
	*t_uz5 = 0.0;
}

void flat_wall_nx(double *u1_nele, double *u2_nele, double *u3_nele,
		double *u4_nele, double *u5_nele, double *ux1_nele, double *uy1_nele,
		double *uz1_nele, double *ux2_nele, double *uy2_nele, double *uz2_nele,
		double *ux3_nele, double *uy3_nele, double *uz3_nele, double *ux4_nele,
		double *uy4_nele, double *uz4_nele, double *ux5_nele, double *uy5_nele,
		double *uz5_nele, double *u1_iele, double *u2_iele, double *u3_iele,
		double *u4_iele, double *u5_iele, double *ux1_iele, double *uy1_iele,
		double *uz1_iele, double *ux2_iele, double *uy2_iele, double *uz2_iele,
		double *ux3_iele, double *uy3_iele, double *uz3_iele, double *ux4_iele,
		double *uy4_iele, double *uz4_iele, double *ux5_iele, double *uy5_iele,
		double *uz5_iele)
{
	*u1_nele = *u1_iele;
	*u2_nele = 0.0;
	*u3_nele = 0.0;
	*u4_nele = 0.0;
	*u5_nele = *u5_iele
			- (0.5
					* (*u2_iele * *u2_iele + *u3_iele * *u3_iele
							+ *u4_iele * *u4_iele) / (*u1_iele));
	*ux1_nele = 0;
	*uy1_nele = *uy1_iele;
	*uz1_nele = *uz1_iele;
	*ux2_nele = *ux2_iele;
	*uy2_nele = 0.0;
	*uz2_nele = 0.0;
	*ux3_nele = *ux3_iele;
	*uy3_nele = 0.0;
	*uz3_nele = 0.0;
	*ux4_nele = *ux4_iele;
	*uy4_nele = 0.0;
	*uz4_nele = 0.0;
	*ux5_nele = 0.0;
	*uy5_nele = *uy5_iele;
	*uz5_nele = *uz5_iele;
}

void flat_wall_ny(double *u1_nele, double *u2_nele, double *u3_nele,
		double *u4_nele, double *u5_nele, double *ux1_nele, double *uy1_nele,
		double *uz1_nele, double *ux2_nele, double *uy2_nele, double *uz2_nele,
		double *ux3_nele, double *uy3_nele, double *uz3_nele, double *ux4_nele,
		double *uy4_nele, double *uz4_nele, double *ux5_nele, double *uy5_nele,
		double *uz5_nele, double *u1_iele, double *u2_iele, double *u3_iele,
		double *u4_iele, double *u5_iele, double *ux1_iele, double *uy1_iele,
		double *uz1_iele, double *ux2_iele, double *uy2_iele, double *uz2_iele,
		double *ux3_iele, double *uy3_iele, double *uz3_iele, double *ux4_iele,
		double *uy4_iele, double *uz4_iele, double *ux5_iele, double *uy5_iele,
		double *uz5_iele)
{
	*u1_nele = *u1_iele;
	*u2_nele = 0.0;
	*u3_nele = 0.0;
	*u4_nele = 0.0;
	*u5_nele = *u5_iele
			- (0.5
					* (*u2_iele * *u2_iele + *u3_iele * *u3_iele
							+ *u4_iele * *u4_iele) / (*u1_iele));
	*ux1_nele = *ux1_iele;
	*uy1_nele = 0.0;
	*uz1_nele = *uz1_iele;
	*ux2_nele = 0.0;
	*uy2_nele = *uy2_iele;
	*uz2_nele = 0.0;
	*ux3_nele = 0.0;
	*uy3_nele = *uy3_iele;
	*uz3_nele = 0.0;
	*ux4_nele = 0.0;
	*uy4_nele = *uy4_iele;
	*uz4_nele = 0.0;
	*ux5_nele = *ux5_iele;
	*uy5_nele = 0;
	*uz5_nele = *uz5_iele;
}

void flat_wall_nz(double *u1_nele, double *u2_nele, double *u3_nele,
		double *u4_nele, double *u5_nele, double *ux1_nele, double *uy1_nele,
		double *uz1_nele, double *ux2_nele, double *uy2_nele, double *uz2_nele,
		double *ux3_nele, double *uy3_nele, double *uz3_nele, double *ux4_nele,
		double *uy4_nele, double *uz4_nele, double *ux5_nele, double *uy5_nele,
		double *uz5_nele, double *u1_iele, double *u2_iele, double *u3_iele,
		double *u4_iele, double *u5_iele, double *ux1_iele, double *uy1_iele,
		double *uz1_iele, double *ux2_iele, double *uy2_iele, double *uz2_iele,
		double *ux3_iele, double *uy3_iele, double *uz3_iele, double *ux4_iele,
		double *uy4_iele, double *uz4_iele, double *ux5_iele, double *uy5_iele,
		double *uz5_iele)
{
	*u1_nele = *u1_iele;
	*u2_nele = 0.0;
	*u3_nele = 0.0;
	*u4_nele = 0.0;
	*u5_nele = *u5_iele
			- (0.5
					* (*u2_iele * *u2_iele + *u3_iele * *u3_iele
							+ *u4_iele * *u4_iele) / (*u1_iele));
	*ux1_nele = *ux1_iele;
	*uy1_nele = *uy1_iele;
	*uz1_nele = 0.0;
	*ux2_nele = 0.0;
	*uy2_nele = 0.0;
	*uz2_nele = *uz2_iele;
	*ux3_nele = 0.0;
	*uy3_nele = 0.0;
	*uz3_nele = *uz3_iele;
	*ux4_nele = 0.0;
	*uy4_nele = 0.0;
	*uz4_nele = *uz4_iele;
	*ux5_nele = *ux5_iele;
	*uy5_nele = *uy5_iele;
	*uz5_nele = 0;
}

void flat_outflow_nx(double *u1_nele, double *u2_nele, double *u3_nele,
		double *u4_nele, double *u5_nele, double *ux1_nele, double *uy1_nele,
		double *uz1_nele, double *ux2_nele, double *uy2_nele, double *uz2_nele,
		double *ux3_nele, double *uy3_nele, double *uz3_nele, double *ux4_nele,
		double *uy4_nele, double *uz4_nele, double *ux5_nele, double *uy5_nele,
		double *uz5_nele, double *u1_iele, double *u2_iele, double *u3_iele,
		double *u4_iele, double *u5_iele, double *ux1_iele, double *uy1_iele,
		double *uz1_iele, double *ux2_iele, double *uy2_iele, double *uz2_iele,
		double *ux3_iele, double *uy3_iele, double *uz3_iele, double *ux4_iele,
		double *uy4_iele, double *uz4_iele, double *ux5_iele, double *uy5_iele,
		double *uz5_iele)
{
	*u1_nele = *u1_iele;
	*u2_nele = *u2_iele;
	*u3_nele = *u3_iele;
	*u4_nele = *u4_iele;
	*u5_nele = *u5_iele;
	*ux1_nele = 0.0;
	*uy1_nele = *uy1_iele;
	*uz1_nele = *uz1_iele;
	*ux2_nele = 0.0;
	*uy2_nele = *uy2_iele;
	*uz2_nele = *uz2_iele;
	*ux3_nele = 0.0;
	*uy3_nele = *uy3_iele;
	*uz3_nele = *uz3_iele;
	*ux4_nele = 0.0;
	*uy4_nele = *uy4_iele;
	*uz4_nele = *uz4_iele;
	*ux5_nele = 0.0;
	*uy5_nele = *uy5_iele;
	*uz5_nele = *uz5_iele;
}

void flat_outflow_ny(double *u1_nele, double *u2_nele, double *u3_nele,
		double *u4_nele, double *u5_nele, double *ux1_nele, double *uy1_nele,
		double *uz1_nele, double *ux2_nele, double *uy2_nele, double *uz2_nele,
		double *ux3_nele, double *uy3_nele, double *uz3_nele, double *ux4_nele,
		double *uy4_nele, double *uz4_nele, double *ux5_nele, double *uy5_nele,
		double *uz5_nele, double *u1_iele, double *u2_iele, double *u3_iele,
		double *u4_iele, double *u5_iele, double *ux1_iele, double *uy1_iele,
		double *uz1_iele, double *ux2_iele, double *uy2_iele, double *uz2_iele,
		double *ux3_iele, double *uy3_iele, double *uz3_iele, double *ux4_iele,
		double *uy4_iele, double *uz4_iele, double *ux5_iele, double *uy5_iele,
		double *uz5_iele)
{
	*u1_nele = *u1_iele;
	*u2_nele = *u2_iele;
	*u3_nele = *u3_iele;
	*u4_nele = *u4_iele;
	*u5_nele = *u5_iele;
	*ux1_nele = *ux1_iele;
	*uy1_nele = 0.0;
	*uz1_nele = *uz1_iele;
	*ux2_nele = *ux2_iele;
	*uy2_nele = 0.0;
	*uz2_nele = *uz2_iele;
	*ux3_nele = *ux3_iele;
	*uy3_nele = 0.0;
	*uz3_nele = *uz3_iele;
	*ux4_nele = *ux4_iele;
	*uy4_nele = 0.0;
	*uz4_nele = *uz4_iele;
	*ux5_nele = *ux5_iele;
	*uy5_nele = 0.0;
	*uz5_nele = *uz5_iele;
}

void flat_outflow_nz(double *u1_nele, double *u2_nele, double *u3_nele,
		double *u4_nele, double *u5_nele, double *ux1_nele, double *uy1_nele,
		double *uz1_nele, double *ux2_nele, double *uy2_nele, double *uz2_nele,
		double *ux3_nele, double *uy3_nele, double *uz3_nele, double *ux4_nele,
		double *uy4_nele, double *uz4_nele, double *ux5_nele, double *uy5_nele,
		double *uz5_nele, double *u1_iele, double *u2_iele, double *u3_iele,
		double *u4_iele, double *u5_iele, double *ux1_iele, double *uy1_iele,
		double *uz1_iele, double *ux2_iele, double *uy2_iele, double *uz2_iele,
		double *ux3_iele, double *uy3_iele, double *uz3_iele, double *ux4_iele,
		double *uy4_iele, double *uz4_iele, double *ux5_iele, double *uy5_iele,
		double *uz5_iele)
{
	*u1_nele = *u1_iele;
	*u2_nele = *u2_iele;
	*u3_nele = *u3_iele;
	*u4_nele = *u4_iele;
	*u5_nele = *u5_iele;
	*ux1_nele = *ux1_iele;
	*uy1_nele = *uy1_iele;
	*uz1_nele = 0.0;
	*ux2_nele = *ux2_iele;
	*uy2_nele = *uy2_iele;
	*uz2_nele = 0.0;
	*ux3_nele = *ux3_iele;
	*uy3_nele = *uy3_iele;
	*uz3_nele = 0.0;
	*ux4_nele = *ux4_iele;
	*uy4_nele = *uy4_iele;
	*uz4_nele = 0.0;
	*ux5_nele = *ux5_iele;
	*uy5_nele = *uy5_iele;
	*uz5_nele = 0.0;
}

void periodic(double *u1_nele, double *u2_nele, double *u3_nele,
		double *u4_nele, double *u5_nele, double *ux1_nele, double *uy1_nele,
		double *uz1_nele, double *ux2_nele, double *uy2_nele, double *uz2_nele,
		double *ux3_nele, double *uy3_nele, double *uz3_nele, double *ux4_nele,
		double *uy4_nele, double *uz4_nele, double *ux5_nele, double *uy5_nele,
		double *uz5_nele, double *u1_iele, double *u2_iele, double *u3_iele,
		double *u4_iele, double *u5_iele, double *ux1_iele, double *uy1_iele,
		double *uz1_iele, double *ux2_iele, double *uy2_iele, double *uz2_iele,
		double *ux3_iele, double *uy3_iele, double *uz3_iele, double *ux4_iele,
		double *uy4_iele, double *uz4_iele, double *ux5_iele, double *uy5_iele,
		double *uz5_iele)
{
	*u1_nele = *u1_iele;
	*u2_nele = *u2_iele;
	*u3_nele = *u3_iele;
	*u4_nele = *u4_iele;
	*u5_nele = *u5_iele;
	*ux1_nele = *ux1_iele;
	*uy1_nele = *uy1_iele;
	*uz1_nele = *uz1_iele;
	*ux2_nele = *ux2_iele;
	*uy2_nele = *uy2_iele;
	*uz2_nele = *uz2_iele;
	*ux3_nele = *ux3_iele;
	*uy3_nele = *uy3_iele;
	*uz3_nele = *uz3_iele;
	*ux4_nele = *ux4_iele;
	*uy4_nele = *uy4_iele;
	*uz4_nele = *uz4_iele;
	*ux5_nele = *ux5_iele;
	*uy5_nele = *uy5_iele;
	*uz5_nele = *uz5_iele;
}

//bound_cond => function for allocating boundary condition

int boun_cond(int i, int b_cnt, int nele, int* boun_type, int* j,
		int* cnt_boun_type, int* iele, int** boun_elem, double** ghost_elem,
		double* t_ux1, double* t_uy1, double* t_uz1, double* t_ux2,
		double* t_uy2, double* t_uz2, double* t_ux3, double* t_uy3,
		double* t_uz3, double* t_ux4, double* t_uy4, double* t_uz4,
		double* t_ux5, double* t_uy5, double* t_uz5, double* t_ux1_iele,
		double* t_uy1_iele, double* t_uz1_iele, double* t_ux2_iele,
		double* t_uy2_iele, double* t_uz2_iele, double* t_ux3_iele,
		double* t_uy3_iele, double* t_uz3_iele, double* t_ux4_iele,
		double* t_uy4_iele, double* t_uz4_iele, double* t_ux5_iele,
		double* t_uy5_iele, double* t_uz5_iele, double** trans_coord_dx,
		double** trans_coord_dy, double** trans_coord_dz,
		double** trans_coord_mx, double** trans_coord_my,
		double** trans_coord_mz, double** nx_boun_face, double** ny_boun_face,
		double** nz_boun_face, double* u1_p, double* u2_p, double* u3_p,
		double* u4_p, double* u5_p, double* ux1_p, double* uy1_p, double* uz1_p,
		double* ux2_p, double* uy2_p, double* uz2_p, double* ux3_p,
		double* uy3_p, double* uz3_p, double* ux4_p, double* uy4_p,
		double* uz4_p, double* ux5_p, double* uy5_p, double* uz5_p,
		double* u1initial, double* u2initial, double* u3initial,
		double* u4initial, double* u5initial, int** periodic_elem,
		struct node5* primitive, struct node6* derv)
{
	//boundary conditions
	for (i = 1; i <= b_cnt; i++)
	{
		switch (boun_type[i])
		{
		//1: curved wall
		case 1:
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				nele = ghost_elem[i][*j];
				curved_wall(t_ux1, t_uy1, t_uz1, t_ux2, t_uy2, t_uz2, t_ux3,
						t_uy3, t_uz3, t_ux4, t_uy4, t_uz4, t_ux5, t_uy5, t_uz5,
						t_ux1_iele, t_uy1_iele, t_uz1_iele, t_ux2_iele,
						t_uy2_iele, t_uz2_iele, t_ux3_iele, t_uy3_iele,
						t_uz3_iele, t_ux4_iele, t_uy4_iele, t_uz4_iele,
						t_ux5_iele, t_uy5_iele, t_uz5_iele,
						&trans_coord_dx[i][*j], &trans_coord_dy[i][*j],
						&trans_coord_dz[i][*j], &trans_coord_mx[i][*j],
						&trans_coord_my[i][*j], &trans_coord_mz[i][*j],
						&nx_boun_face[i][*j], &ny_boun_face[i][*j],
						&nz_boun_face[i][*j], &u1_p[nele], &u2_p[nele],
						&u3_p[nele], &u4_p[nele], &u5_p[nele], &u1_p[*iele],
						&u2_p[*iele], &u3_p[*iele], &u4_p[*iele], &u5_p[*iele],
						&ux1_p[*iele], &uy1_p[*iele], &uz1_p[*iele],
						&ux2_p[*iele], &uy2_p[*iele], &uz2_p[*iele],
						&ux3_p[*iele], &uy3_p[*iele], &uz3_p[*iele],
						&ux4_p[*iele], &uy4_p[*iele], &uz4_p[*iele],
						&ux5_p[*iele], &uy5_p[*iele], &uz5_p[*iele]);
				derv_trans(&ux1_p[nele], &uy1_p[nele], &uz1_p[nele],
						&ux2_p[nele], &uy2_p[nele], &uz2_p[nele], &ux3_p[nele],
						&uy3_p[nele], &uz3_p[nele], &ux4_p[nele], &uy4_p[nele],
						&uz4_p[nele], &ux5_p[nele], &uy5_p[nele], &uz5_p[nele],
						&trans_coord_dx[i][*j], &trans_coord_dy[i][*j],
						&trans_coord_dz[i][*j], &trans_coord_mx[i][*j],
						&trans_coord_my[i][*j], &trans_coord_mz[i][*j],
						&nx_boun_face[i][*j], &ny_boun_face[i][*j],
						&nz_boun_face[i][*j], t_ux1, t_uy1, t_uz1, t_ux2, t_uy2,
						t_uz2, t_ux3, t_uy3, t_uz3, t_ux4, t_uy4, t_uz4, t_ux5,
						t_uy5, t_uz5);
			}
			break;
		case 2:
			//curved inflow
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				nele = ghost_elem[i][*j];
				curved_inflow(&u1_p[nele], &u2_p[nele], &u3_p[nele],
						&u4_p[nele], &u5_p[nele], &ux1_p[nele], &uy1_p[nele],
						&uz1_p[nele], &ux2_p[nele], &uy2_p[nele], &uz2_p[nele],
						&ux3_p[nele], &uy3_p[nele], &uz3_p[nele], &ux4_p[nele],
						&uy4_p[nele], &uz4_p[nele], &ux5_p[nele], &uy5_p[nele],
						&uz5_p[nele], u1initial, u2initial, u3initial,
						u4initial, u5initial);
			}
			break;
		case 3:
			//curved outflow
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				nele = ghost_elem[i][*j];
				curved_outflow(t_ux1, t_uy1, t_uz1, t_ux2, t_uy2, t_uz2, t_ux3,
						t_uy3, t_uz3, t_ux4, t_uy4, t_uz4, t_ux5, t_uy5, t_uz5,
						t_ux1_iele, t_uy1_iele, t_uz1_iele, t_ux2_iele,
						t_uy2_iele, t_uz2_iele, t_ux3_iele, t_uy3_iele,
						t_uz3_iele, t_ux4_iele, t_uy4_iele, t_uz4_iele,
						t_ux5_iele, t_uy5_iele, t_uz5_iele,
						&trans_coord_dx[i][*j], &trans_coord_dy[i][*j],
						&trans_coord_dz[i][*j], &trans_coord_mx[i][*j],
						&trans_coord_my[i][*j], &trans_coord_mz[i][*j],
						&nx_boun_face[i][*j], &ny_boun_face[i][*j],
						&nz_boun_face[i][*j], &u1_p[nele], &u2_p[nele],
						&u3_p[nele], &u4_p[nele], &u5_p[nele], &u1_p[*iele],
						&u2_p[*iele], &u3_p[*iele], &u4_p[*iele], &u5_p[*iele],
						&ux1_p[*iele], &uy1_p[*iele], &uz1_p[*iele],
						&ux2_p[*iele], &uy2_p[*iele], &uz2_p[*iele],
						&ux3_p[*iele], &uy3_p[*iele], &uz3_p[*iele],
						&ux4_p[*iele], &uy4_p[*iele], &uz4_p[*iele],
						&ux5_p[*iele], &uy5_p[*iele], &uz5_p[*iele]);
				derv_trans(&ux1_p[nele], &uy1_p[nele], &uz1_p[nele],
						&ux2_p[nele], &uy2_p[nele], &uz2_p[nele], &ux3_p[nele],
						&uy3_p[nele], &uz3_p[nele], &ux4_p[nele], &uy4_p[nele],
						&uz4_p[nele], &ux5_p[nele], &uy5_p[nele], &uz5_p[nele],
						&trans_coord_dx[i][*j], &trans_coord_dy[i][*j],
						&trans_coord_dz[i][*j], &trans_coord_mx[i][*j],
						&trans_coord_my[i][*j], &trans_coord_mz[i][*j],
						&nx_boun_face[i][*j], &ny_boun_face[i][*j],
						&nz_boun_face[i][*j], t_ux1, t_uy1, t_uz1, t_ux2, t_uy2,
						t_uz2, t_ux3, t_uy3, t_uz3, t_ux4, t_uy4, t_uz4, t_ux5,
						t_uy5, t_uz5);
			}
			break;
		case 4:
			//curved upper
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				nele = ghost_elem[i][*j];
				if ((nx_boun_face[i][*j] * u2_p[*iele]
						+ ny_boun_face[i][*j] * u3_p[*iele]
						+ nz_boun_face[i][*j] * u4_p[*iele]) > 0.0)
				{
					curved_outflow(t_ux1, t_uy1, t_uz1, t_ux2, t_uy2, t_uz2,
							t_ux3, t_uy3, t_uz3, t_ux4, t_uy4, t_uz4, t_ux5,
							t_uy5, t_uz5, t_ux1_iele, t_uy1_iele, t_uz1_iele,
							t_ux2_iele, t_uy2_iele, t_uz2_iele, t_ux3_iele,
							t_uy3_iele, t_uz3_iele, t_ux4_iele, t_uy4_iele,
							t_uz4_iele, t_ux5_iele, t_uy5_iele, t_uz5_iele,
							&trans_coord_dx[i][*j], &trans_coord_dy[i][*j],
							&trans_coord_dz[i][*j], &trans_coord_mx[i][*j],
							&trans_coord_my[i][*j], &trans_coord_mz[i][*j],
							&nx_boun_face[i][*j], &ny_boun_face[i][*j],
							&nz_boun_face[i][*j], &u1_p[nele], &u2_p[nele],
							&u3_p[nele], &u4_p[nele], &u5_p[nele], &u1_p[*iele],
							&u2_p[*iele], &u3_p[*iele], &u4_p[*iele],
							&u5_p[*iele], &ux1_p[*iele], &uy1_p[*iele],
							&uz1_p[*iele], &ux2_p[*iele], &uy2_p[*iele],
							&uz2_p[*iele], &ux3_p[*iele], &uy3_p[*iele],
							&uz3_p[*iele], &ux4_p[*iele], &uy4_p[*iele],
							&uz4_p[*iele], &ux5_p[*iele], &uy5_p[*iele],
							&uz5_p[*iele]);
					derv_trans(&ux1_p[nele], &uy1_p[nele], &uz1_p[nele],
							&ux2_p[nele], &uy2_p[nele], &uz2_p[nele],
							&ux3_p[nele], &uy3_p[nele], &uz3_p[nele],
							&ux4_p[nele], &uy4_p[nele], &uz4_p[nele],
							&ux5_p[nele], &uy5_p[nele], &uz5_p[nele],
							&trans_coord_dx[i][*j], &trans_coord_dy[i][*j],
							&trans_coord_dz[i][*j], &trans_coord_mx[i][*j],
							&trans_coord_my[i][*j], &trans_coord_mz[i][*j],
							&nx_boun_face[i][*j], &ny_boun_face[i][*j],
							&nz_boun_face[i][*j], t_ux1, t_uy1, t_uz1, t_ux2,
							t_uy2, t_uz2, t_ux3, t_uy3, t_uz3, t_ux4, t_uy4,
							t_uz4, t_ux5, t_uy5, t_uz5);
				}
				else
				{
					curved_inflow(&u1_p[nele], &u2_p[nele], &u3_p[nele],
							&u4_p[nele], &u5_p[nele], &ux1_p[nele],
							&uy1_p[nele], &uz1_p[nele], &ux2_p[nele],
							&uy2_p[nele], &uz2_p[nele], &ux3_p[nele],
							&uy3_p[nele], &uz3_p[nele], &ux4_p[nele],
							&uy4_p[nele], &uz4_p[nele], &ux5_p[nele],
							&uy5_p[nele], &uz5_p[nele], u1initial, u2initial,
							u3initial, u4initial, u5initial);
				}
			}
			break;
		case 5:
			//flat_wall_nx
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				nele = ghost_elem[i][*j];
				flat_wall_nx(&u1_p[nele], &u2_p[nele], &u3_p[nele], &u4_p[nele],
						&u5_p[nele], &ux1_p[nele], &uy1_p[nele], &uz1_p[nele],
						&ux2_p[nele], &uy2_p[nele], &uz2_p[nele], &ux3_p[nele],
						&uy3_p[nele], &uz3_p[nele], &ux4_p[nele], &uy4_p[nele],
						&uz4_p[nele], &ux5_p[nele], &uy5_p[nele], &uz5_p[nele],
						&u1_p[*iele], &u2_p[*iele], &u3_p[*iele], &u4_p[*iele],
						&u5_p[*iele], &ux1_p[*iele], &uy1_p[*iele],
						&uz1_p[*iele], &ux2_p[*iele], &uy2_p[*iele],
						&uz2_p[*iele], &ux3_p[*iele], &uy3_p[*iele],
						&uz3_p[*iele], &ux4_p[*iele], &uy4_p[*iele],
						&uz4_p[*iele], &ux5_p[*iele], &uy5_p[*iele],
						&uz5_p[*iele]);
			}
			break;
		case 6:
			//flat_wall_ny
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				nele = ghost_elem[i][*j];
				flat_wall_ny(&u1_p[nele], &u2_p[nele], &u3_p[nele], &u4_p[nele],
						&u5_p[nele], &ux1_p[nele], &uy1_p[nele], &uz1_p[nele],
						&ux2_p[nele], &uy2_p[nele], &uz2_p[nele], &ux3_p[nele],
						&uy3_p[nele], &uz3_p[nele], &ux4_p[nele], &uy4_p[nele],
						&uz4_p[nele], &ux5_p[nele], &uy5_p[nele], &uz5_p[nele],
						&u1_p[*iele], &u2_p[*iele], &u3_p[*iele], &u4_p[*iele],
						&u5_p[*iele], &ux1_p[*iele], &uy1_p[*iele],
						&uz1_p[*iele], &ux2_p[*iele], &uy2_p[*iele],
						&uz2_p[*iele], &ux3_p[*iele], &uy3_p[*iele],
						&uz3_p[*iele], &ux4_p[*iele], &uy4_p[*iele],
						&uz4_p[*iele], &ux5_p[*iele], &uy5_p[*iele],
						&uz5_p[*iele]);
			}
			break;
		case 7:
			//flat_wall_nz
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				nele = ghost_elem[i][*j];
				flat_wall_nz(&u1_p[nele], &u2_p[nele], &u3_p[nele], &u4_p[nele],
						&u5_p[nele], &ux1_p[nele], &uy1_p[nele], &uz1_p[nele],
						&ux2_p[nele], &uy2_p[nele], &uz2_p[nele], &ux3_p[nele],
						&uy3_p[nele], &uz3_p[nele], &ux4_p[nele], &uy4_p[nele],
						&uz4_p[nele], &ux5_p[nele], &uy5_p[nele], &uz5_p[nele],
						&u1_p[*iele], &u2_p[*iele], &u3_p[*iele], &u4_p[*iele],
						&u5_p[*iele], &ux1_p[*iele], &uy1_p[*iele],
						&uz1_p[*iele], &ux2_p[*iele], &uy2_p[*iele],
						&uz2_p[*iele], &ux3_p[*iele], &uy3_p[*iele],
						&uz3_p[*iele], &ux4_p[*iele], &uy4_p[*iele],
						&uz4_p[*iele], &ux5_p[*iele], &uy5_p[*iele],
						&uz5_p[*iele]);
			}
			break;
		case 8:
			//flat inflow is same as curved inflow:flat inflow perpendicular to x axis
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				nele = ghost_elem[i][*j];
				curved_inflow(&u1_p[nele], &u2_p[nele], &u3_p[nele],
						&u4_p[nele], &u5_p[nele], &ux1_p[nele], &uy1_p[nele],
						&uz1_p[nele], &ux2_p[nele], &uy2_p[nele], &uz2_p[nele],
						&ux3_p[nele], &uy3_p[nele], &uz3_p[nele], &ux4_p[nele],
						&uy4_p[nele], &uz4_p[nele], &ux5_p[nele], &uy5_p[nele],
						&uz5_p[nele], u1initial, u2initial, u3initial,
						u4initial, u5initial);
			}
			break;
		case 9:
			//flat inflow is same as curved inflow:flat inflow perpendicular to y axis
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				nele = ghost_elem[i][*j];
				curved_inflow(&u1_p[nele], &u2_p[nele], &u3_p[nele],
						&u4_p[nele], &u5_p[nele], &ux1_p[nele], &uy1_p[nele],
						&uz1_p[nele], &ux2_p[nele], &uy2_p[nele], &uz2_p[nele],
						&ux3_p[nele], &uy3_p[nele], &uz3_p[nele], &ux4_p[nele],
						&uy4_p[nele], &uz4_p[nele], &ux5_p[nele], &uy5_p[nele],
						&uz5_p[nele], u1initial, u2initial, u3initial,
						u4initial, u5initial);
			}
			break;
		case 10:
			//flat inflow is same as curved inflow:flat inflow perpendicular to z axis
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				nele = ghost_elem[i][*j];
				curved_inflow(&u1_p[nele], &u2_p[nele], &u3_p[nele],
						&u4_p[nele], &u5_p[nele], &ux1_p[nele], &uy1_p[nele],
						&uz1_p[nele], &ux2_p[nele], &uy2_p[nele], &uz2_p[nele],
						&ux3_p[nele], &uy3_p[nele], &uz3_p[nele], &ux4_p[nele],
						&uy4_p[nele], &uz4_p[nele], &ux5_p[nele], &uy5_p[nele],
						&uz5_p[nele], u1initial, u2initial, u3initial,
						u4initial, u5initial);
			}
			break;
		case 11:
			//flow_outflow_nx
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				nele = ghost_elem[i][*j];
				flat_outflow_nx(&u1_p[nele], &u2_p[nele], &u3_p[nele],
						&u4_p[nele], &u5_p[nele], &ux1_p[nele], &uy1_p[nele],
						&uz1_p[nele], &ux2_p[nele], &uy2_p[nele], &uz2_p[nele],
						&ux3_p[nele], &uy3_p[nele], &uz3_p[nele], &ux4_p[nele],
						&uy4_p[nele], &uz4_p[nele], &ux5_p[nele], &uy5_p[nele],
						&uz5_p[nele], &u1_p[*iele], &u2_p[*iele], &u3_p[*iele],
						&u4_p[*iele], &u5_p[*iele], &ux1_p[*iele],
						&uy1_p[*iele], &uz1_p[*iele], &ux2_p[*iele],
						&uy2_p[*iele], &uz2_p[*iele], &ux3_p[*iele],
						&uy3_p[*iele], &uz3_p[*iele], &ux4_p[*iele],
						&uy4_p[*iele], &uz4_p[*iele], &ux5_p[*iele],
						&uy5_p[*iele], &uz5_p[*iele]);
			}
			break;
		case 12:
			//flow_outflow_ny
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				nele = ghost_elem[i][*j];
				flat_outflow_ny(&u1_p[nele], &u2_p[nele], &u3_p[nele],
						&u4_p[nele], &u5_p[nele], &ux1_p[nele], &uy1_p[nele],
						&uz1_p[nele], &ux2_p[nele], &uy2_p[nele], &uz2_p[nele],
						&ux3_p[nele], &uy3_p[nele], &uz3_p[nele], &ux4_p[nele],
						&uy4_p[nele], &uz4_p[nele], &ux5_p[nele], &uy5_p[nele],
						&uz5_p[nele], &u1_p[*iele], &u2_p[*iele], &u3_p[*iele],
						&u4_p[*iele], &u5_p[*iele], &ux1_p[*iele],
						&uy1_p[*iele], &uz1_p[*iele], &ux2_p[*iele],
						&uy2_p[*iele], &uz2_p[*iele], &ux3_p[*iele],
						&uy3_p[*iele], &uz3_p[*iele], &ux4_p[*iele],
						&uy4_p[*iele], &uz4_p[*iele], &ux5_p[*iele],
						&uy5_p[*iele], &uz5_p[*iele]);
			}
			break;
		case 13:
			//flow_outflow_nz
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				nele = ghost_elem[i][*j];
				flat_outflow_nz(&u1_p[nele], &u2_p[nele], &u3_p[nele],
						&u4_p[nele], &u5_p[nele], &ux1_p[nele], &uy1_p[nele],
						&uz1_p[nele], &ux2_p[nele], &uy2_p[nele], &uz2_p[nele],
						&ux3_p[nele], &uy3_p[nele], &uz3_p[nele], &ux4_p[nele],
						&uy4_p[nele], &uz4_p[nele], &ux5_p[nele], &uy5_p[nele],
						&uz5_p[nele], &u1_p[*iele], &u2_p[*iele], &u3_p[*iele],
						&u4_p[*iele], &u5_p[*iele], &ux1_p[*iele],
						&uy1_p[*iele], &uz1_p[*iele], &ux2_p[*iele],
						&uy2_p[*iele], &uz2_p[*iele], &ux3_p[*iele],
						&uy3_p[*iele], &uz3_p[*iele], &ux4_p[*iele],
						&uy4_p[*iele], &uz4_p[*iele], &ux5_p[*iele],
						&uy5_p[*iele], &uz5_p[*iele]);
			}
			break;
		case 14:
			//flow_upper_nx
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				nele = ghost_elem[i][*j];
				if ((nx_boun_face[i][*j] * u2_p[*iele]
						+ ny_boun_face[i][*j] * u3_p[*iele]
						+ nz_boun_face[i][*j] * u4_p[*iele]) > 0.0)
				{
					flat_outflow_nx(&u1_p[nele], &u2_p[nele], &u3_p[nele],
							&u4_p[nele], &u5_p[nele], &ux1_p[nele],
							&uy1_p[nele], &uz1_p[nele], &ux2_p[nele],
							&uy2_p[nele], &uz2_p[nele], &ux3_p[nele],
							&uy3_p[nele], &uz3_p[nele], &ux4_p[nele],
							&uy4_p[nele], &uz4_p[nele], &ux5_p[nele],
							&uy5_p[nele], &uz5_p[nele], &u1_p[*iele],
							&u2_p[*iele], &u3_p[*iele], &u4_p[*iele],
							&u5_p[*iele], &ux1_p[*iele], &uy1_p[*iele],
							&uz1_p[*iele], &ux2_p[*iele], &uy2_p[*iele],
							&uz2_p[*iele], &ux3_p[*iele], &uy3_p[*iele],
							&uz3_p[*iele], &ux4_p[*iele], &uy4_p[*iele],
							&uz4_p[*iele], &ux5_p[*iele], &uy5_p[*iele],
							&uz5_p[*iele]);
				}
				else
				{
					curved_inflow(&u1_p[nele], &u2_p[nele], &u3_p[nele],
							&u4_p[nele], &u5_p[nele], &ux1_p[nele],
							&uy1_p[nele], &uz1_p[nele], &ux2_p[nele],
							&uy2_p[nele], &uz2_p[nele], &ux3_p[nele],
							&uy3_p[nele], &uz3_p[nele], &ux4_p[nele],
							&uy4_p[nele], &uz4_p[nele], &ux5_p[nele],
							&uy5_p[nele], &uz5_p[nele], u1initial, u2initial,
							u3initial, u4initial, u5initial);
				}
			}
			break;
		case 15:
			//flow_upper_ny
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				nele = ghost_elem[i][*j];
				if ((nx_boun_face[i][*j] * u2_p[*iele]
						+ ny_boun_face[i][*j] * u3_p[*iele]
						+ nz_boun_face[i][*j] * u4_p[*iele]) > 0.0)
				{
					flat_outflow_ny(&u1_p[nele], &u2_p[nele], &u3_p[nele],
							&u4_p[nele], &u5_p[nele], &ux1_p[nele],
							&uy1_p[nele], &uz1_p[nele], &ux2_p[nele],
							&uy2_p[nele], &uz2_p[nele], &ux3_p[nele],
							&uy3_p[nele], &uz3_p[nele], &ux4_p[nele],
							&uy4_p[nele], &uz4_p[nele], &ux5_p[nele],
							&uy5_p[nele], &uz5_p[nele], &u1_p[*iele],
							&u2_p[*iele], &u3_p[*iele], &u4_p[*iele],
							&u5_p[*iele], &ux1_p[*iele], &uy1_p[*iele],
							&uz1_p[*iele], &ux2_p[*iele], &uy2_p[*iele],
							&uz2_p[*iele], &ux3_p[*iele], &uy3_p[*iele],
							&uz3_p[*iele], &ux4_p[*iele], &uy4_p[*iele],
							&uz4_p[*iele], &ux5_p[*iele], &uy5_p[*iele],
							&uz5_p[*iele]);
				}
				else
				{
					curved_inflow(&u1_p[nele], &u2_p[nele], &u3_p[nele],
							&u4_p[nele], &u5_p[nele], &ux1_p[nele],
							&uy1_p[nele], &uz1_p[nele], &ux2_p[nele],
							&uy2_p[nele], &uz2_p[nele], &ux3_p[nele],
							&uy3_p[nele], &uz3_p[nele], &ux4_p[nele],
							&uy4_p[nele], &uz4_p[nele], &ux5_p[nele],
							&uy5_p[nele], &uz5_p[nele], u1initial, u2initial,
							u3initial, u4initial, u5initial);
				}
			}
			break;
		case 16:
			//flow_upper_nz
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = boun_elem[i][*j];
				nele = ghost_elem[i][*j];
				if ((nx_boun_face[i][*j] * u2_p[*iele]
						+ ny_boun_face[i][*j] * u3_p[*iele]
						+ nz_boun_face[i][*j] * u4_p[*iele]) > 0.0)
				{
					flat_outflow_nz(&u1_p[nele], &u2_p[nele], &u3_p[nele],
							&u4_p[nele], &u5_p[nele], &ux1_p[nele],
							&uy1_p[nele], &uz1_p[nele], &ux2_p[nele],
							&uy2_p[nele], &uz2_p[nele], &ux3_p[nele],
							&uy3_p[nele], &uz3_p[nele], &ux4_p[nele],
							&uy4_p[nele], &uz4_p[nele], &ux5_p[nele],
							&uy5_p[nele], &uz5_p[nele], &u1_p[*iele],
							&u2_p[*iele], &u3_p[*iele], &u4_p[*iele],
							&u5_p[*iele], &ux1_p[*iele], &uy1_p[*iele],
							&uz1_p[*iele], &ux2_p[*iele], &uy2_p[*iele],
							&uz2_p[*iele], &ux3_p[*iele], &uy3_p[*iele],
							&uz3_p[*iele], &ux4_p[*iele], &uy4_p[*iele],
							&uz4_p[*iele], &ux5_p[*iele], &uy5_p[*iele],
							&uz5_p[*iele]);
				}
				else
				{
					curved_inflow(&u1_p[nele], &u2_p[nele], &u3_p[nele],
							&u4_p[nele], &u5_p[nele], &ux1_p[nele],
							&uy1_p[nele], &uz1_p[nele], &ux2_p[nele],
							&uy2_p[nele], &uz2_p[nele], &ux3_p[nele],
							&uy3_p[nele], &uz3_p[nele], &ux4_p[nele],
							&uy4_p[nele], &uz4_p[nele], &ux5_p[nele],
							&uy5_p[nele], &uz5_p[nele], u1initial, u2initial,
							u3initial, u4initial, u5initial);
				}
			}
			break;
		case 17:
			//periodic
			for (*j = 1; *j <= cnt_boun_type[i]; *j++)
			{
				*iele = periodic_elem[i][*j];
				nele = ghost_elem[i][*j];
				periodic(&u1_p[nele], &u2_p[nele], &u3_p[nele], &u4_p[nele],
						&u5_p[nele], &ux1_p[nele], &uy1_p[nele], &uz1_p[nele],
						&ux2_p[nele], &uy2_p[nele], &uz2_p[nele], &ux3_p[nele],
						&uy3_p[nele], &uz3_p[nele], &ux4_p[nele], &uy4_p[nele],
						&uz4_p[nele], &ux5_p[nele], &uy5_p[nele], &uz5_p[nele],
						&primitive[*iele].u1, &primitive[*iele].u2,
						&primitive[*iele].u3, &primitive[*iele].u4,
						&primitive[*iele].u5, &derv[*iele].ux1,
						&derv[*iele].uy1, &derv[*iele].uz1, &derv[*iele].ux2,
						&derv[*iele].uy2, &derv[*iele].uz2, &derv[*iele].ux3,
						&derv[*iele].uy3, &derv[*iele].uz3, &derv[*iele].ux4,
						&derv[*iele].uy4, &derv[*iele].uz4, &derv[*iele].ux5,
						&derv[*iele].uy5, &derv[*iele].uz5);
			}
			break;
		}
		//switch case
	}	//for loopreturn i;
}

int previous_time_step_alloc(int iele, int subdomain_cnt, double* u1_p,
		struct node5* primitive, double* ux1_p, struct node6* derv,
		double* uy1_p, double* uz1_p, double* ut1_p, struct node7* secondary,
		double* f1_p, double* fx1_p, double* fy1_p, double* fz1_p,
		double* ft1_p, double* g1_p, double* gx1_p, double* gy1_p,
		double* gz1_p, double* gt1_p, double* e1_p, double* ex1_p,
		double* ey1_p, double* ez1_p, double* et1_p, double* u2_p,
		double* ux2_p, double* uy2_p, double* uz2_p, double* ut2_p,
		double* f2_p, double* fx2_p, double* fy2_p, double* fz2_p,
		double* ft2_p, double* g2_p, double* gx2_p, double* gy2_p,
		double* gz2_p, double* gt2_p, double* e2_p, double* ex2_p,
		double* ey2_p, double* ez2_p, double* et2_p, double* u3_p,
		double* ux3_p, double* uy3_p, double* uz3_p, double* ut3_p,
		double* f3_p, double* fx3_p, double* fy3_p, double* fz3_p,
		double* ft3_p, double* g3_p, double* gx3_p, double* gy3_p,
		double* gz3_p, double* gt3_p, double* e3_p, double* ex3_p,
		double* ey3_p, double* ez3_p, double* et3_p, double* u4_p,
		double* ux4_p, double* uy4_p, double* uz4_p, double* ut4_p,
		double* f4_p, double* fx4_p, double* fy4_p, double* fz4_p,
		double* ft4_p, double* g4_p, double* gx4_p, double* gy4_p,
		double* gz4_p, double* gt4_p, double* e4_p, double* ex4_p,
		double* ey4_p, double* ez4_p, double* et4_p, double* u5_p,
		double* ux5_p, double* uy5_p, double* uz5_p, double* ut5_p,
		double* f5_p, double* fx5_p, double* fy5_p, double* fz5_p,
		double* ft5_p, double* g5_p, double* gx5_p, double* gy5_p,
		double* gz5_p, double* gt5_p, double* e5_p, double* ex5_p,
		double* ey5_p, double* ez5_p, double* et5_p)
{
	//allocating values for u[]_p:previous values
	for (iele = 1; iele <= subdomain_cnt; iele++)
	{
		u1_p[iele] = primitive[iele].u1; //p stands for previous time step
		ux1_p[iele] = derv[iele].ux1;
		uy1_p[iele] = derv[iele].uy1;
		uz1_p[iele] = derv[iele].uz1;
		ut1_p[iele] = secondary[iele].ut1;
		f1_p[iele] = secondary[iele].f1;
		fx1_p[iele] = secondary[iele].fx1;
		fy1_p[iele] = secondary[iele].fy1;
		fz1_p[iele] = secondary[iele].fz1;
		ft1_p[iele] = secondary[iele].ft1;
		g1_p[iele] = secondary[iele].g1;
		gx1_p[iele] = secondary[iele].gx1;
		gy1_p[iele] = secondary[iele].gy1;
		gz1_p[iele] = secondary[iele].gz1;
		gt1_p[iele] = secondary[iele].gt1;
		e1_p[iele] = secondary[iele].e1;
		ex1_p[iele] = secondary[iele].ex1;
		ey1_p[iele] = secondary[iele].ey1;
		ez1_p[iele] = secondary[iele].ez1;
		et1_p[iele] = secondary[iele].et1;
		u2_p[iele] = primitive[iele].u2; //p stands for previous time step
		ux2_p[iele] = derv[iele].ux2;
		uy2_p[iele] = derv[iele].uy2;
		uz2_p[iele] = derv[iele].uz2;
		ut2_p[iele] = secondary[iele].ut2;
		f2_p[iele] = secondary[iele].f2;
		fx2_p[iele] = secondary[iele].fx2;
		fy2_p[iele] = secondary[iele].fy2;
		fz2_p[iele] = secondary[iele].fz2;
		ft2_p[iele] = secondary[iele].ft2;
		g2_p[iele] = secondary[iele].g2;
		gx2_p[iele] = secondary[iele].gx2;
		gy2_p[iele] = secondary[iele].gy2;
		gz2_p[iele] = secondary[iele].gz2;
		gt2_p[iele] = secondary[iele].gt2;
		e2_p[iele] = secondary[iele].e2;
		ex2_p[iele] = secondary[iele].ex2;
		ey2_p[iele] = secondary[iele].ey2;
		ez2_p[iele] = secondary[iele].ez2;
		et2_p[iele] = secondary[iele].et2;
		u3_p[iele] = primitive[iele].u3; //p stands for previous time step
		ux3_p[iele] = derv[iele].ux3;
		uy3_p[iele] = derv[iele].uy3;
		uz3_p[iele] = derv[iele].uz3;
		ut3_p[iele] = secondary[iele].ut3;
		f3_p[iele] = secondary[iele].f3;
		fx3_p[iele] = secondary[iele].fx3;
		fy3_p[iele] = secondary[iele].fy3;
		fz3_p[iele] = secondary[iele].fz3;
		ft3_p[iele] = secondary[iele].ft3;
		g3_p[iele] = secondary[iele].g3;
		gx3_p[iele] = secondary[iele].gx3;
		gy3_p[iele] = secondary[iele].gy3;
		gz3_p[iele] = secondary[iele].gz3;
		gt3_p[iele] = secondary[iele].gt3;
		e3_p[iele] = secondary[iele].e3;
		ex3_p[iele] = secondary[iele].ex3;
		ey3_p[iele] = secondary[iele].ey3;
		ez3_p[iele] = secondary[iele].ez3;
		et3_p[iele] = secondary[iele].et3;
		u4_p[iele] = primitive[iele].u4; //p stands for previous time step
		ux4_p[iele] = derv[iele].ux4;
		uy4_p[iele] = derv[iele].uy4;
		uz4_p[iele] = derv[iele].uz4;
		ut4_p[iele] = secondary[iele].ut4;
		f4_p[iele] = secondary[iele].f4;
		fx4_p[iele] = secondary[iele].fx4;
		fy4_p[iele] = secondary[iele].fy4;
		fz4_p[iele] = secondary[iele].fz4;
		ft4_p[iele] = secondary[iele].ft4;
		g4_p[iele] = secondary[iele].g4;
		gx4_p[iele] = secondary[iele].gx4;
		gy4_p[iele] = secondary[iele].gy4;
		gz4_p[iele] = secondary[iele].gz4;
		gt4_p[iele] = secondary[iele].gt4;
		e4_p[iele] = secondary[iele].e4;
		ex4_p[iele] = secondary[iele].ex4;
		ey4_p[iele] = secondary[iele].ey4;
		ez4_p[iele] = secondary[iele].ez4;
		et4_p[iele] = secondary[iele].et4;
		u5_p[iele] = primitive[iele].u5; //p stands for previous time step
		ux5_p[iele] = derv[iele].ux5;
		uy5_p[iele] = derv[iele].uy5;
		uz5_p[iele] = derv[iele].uz5;
		ut5_p[iele] = secondary[iele].ut5;
		f5_p[iele] = secondary[iele].f5;
		fx5_p[iele] = secondary[iele].fx5;
		fy5_p[iele] = secondary[iele].fy5;
		fz5_p[iele] = secondary[iele].fz5;
		ft5_p[iele] = secondary[iele].ft5;
		g5_p[iele] = secondary[iele].g5;
		gx5_p[iele] = secondary[iele].gx5;
		gy5_p[iele] = secondary[iele].gy5;
		gz5_p[iele] = secondary[iele].gz5;
		gt5_p[iele] = secondary[iele].gt5;
		e5_p[iele] = secondary[iele].e5;
		ex5_p[iele] = secondary[iele].ex5;
		ey5_p[iele] = secondary[iele].ey5;
		ez5_p[iele] = secondary[iele].ez5;
		et5_p[iele] = secondary[iele].et5;
	}
	return iele;
}

#endif /* SLAVE_BC_IC_H_ */

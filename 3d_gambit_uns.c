//3d CESE program for complex geometry

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/unistd.h>

#define master 0

//declaration "struct" are placed outside main() to avoid error while struct is called in function

struct node1
{
	double xstart, ystart, zstart;
};
struct node1 *xycoord;

struct node2
{
	int gele_nodes[5], node1_localface[5], node2_localface[5],
			node3_localface[5];
};
struct node2 *eleminfo;

struct node5
{
	double u1, u2, u3, u4, u5;
};
struct node5 *primitive;

struct node6
{
	double ux1, uy1, uz1, ux2, uy2, uz2, ux3, uy3, uz3, ux4, uy4, uz4, ux5, uy5,
			uz5;
};
struct node6 *derv;

struct node3
{
	double xc, yc, zc;
};
struct node3 *centroid;

struct node4
{
	double spx, spy, spz;
};
struct node4 *sol;

struct node7
{
	double ut1, ut2, ut3, ut4, ut5, f1, fx1, fy1, fz1, ft1, g1, gx1, gy1, gz1,
			gt1, e1, ex1, ey1, ez1, et1, f2, fx2, fy2, fz2, ft2, g2, gx2, gy2,
			gz2, gt2, e2, ex2, ey2, ez2, et2, f3, fx3, fy3, fz3, ft3, g3, gx3,
			gy3, gz3, gt3, e3, ex3, ey3, ez3, et3, f4, fx4, fy4, fz4, ft4, g4,
			gx4, gy4, gz4, gt4, e4, ex4, ey4, ez4, et4, f5, fx5, fy5, fz5, ft5,
			g5, gx5, gy5, gz5, gt5, e5, ex5, ey5, ez5, et5;
};
struct node7 *secondary;

double a1, a2, a3, a4, a5, a6, a7, a8, a9;
double ga = 1.4;

double f21, f22, f23, f24, f25, g21, g22, g23, e21, e22, e24, g31, g32, g33,
		g34, g35, e31, e33, e34, e41, e42, e43, e44, e45, f51, f52, f53, f54,
		f55, g51, g52, g53, g54, g55, e51, e52, e53, e54, e55;

#include "master_result_writing.h"
#include "master_preproc.h"
#include "slave_solver.h"
#include "slave_bc_ic.h"
#include "slave_preproc.h"

int main(int argc, char *argv[])
{
	FILE *fp, *fp1, *fp2, *fp3, *fp4, *fp5;

	char *buffer, *buffer1, *buffer2, *buffer3, *buffer4;

	int iele, iface, jele, d_node[4], point, dstore, dstore1, check;
	int temp, d_cnt[3], d_ele, proc;
	int position, position1, position2, position3;
	int maxsize, maxsize1, maxsize2, maxsize3, maxsize4, maxsize5, maxsize6,
			maxsize7, maxsize8;

	int slaves[1];
	slaves[0] = 0;
	MPI_Group group_a, group_b;
	MPI_Comm comm_slaves;

	int node, r1, l;

//no_of_nodes=Total number of nodes (without boundary)
//no_of_ele = Total number of elements including boundary: here 2 for 2 sides
//front & back;left &right
//int_ele=only interior elements
	int n, no_of_iter, start_iter, dw_iter,ns;
	int no_of_nodes, int_ele, no_of_ele;

	int i, j, k, ii, jj, kk, m, jface, ineigh, iside, nele, fg, pg;

	double *re, *pr;
	re = malloc(sizeof(double));
	pr = malloc(sizeof(double));

	int tot_ele;
	double *dummy, *maximum;
	dummy = malloc(sizeof(double));
	maximum = malloc(sizeof(double));

	double *u1initial, *u2initial, *u3initial, *u4initial, *u5initial;
	u1initial = malloc(sizeof(double));
	u2initial = malloc(sizeof(double));
	u3initial = malloc(sizeof(double));
	u4initial = malloc(sizeof(double));
	u5initial = malloc(sizeof(double));

	double *ma;
	ma = malloc(sizeof(double));

	double *delt, *d_delt;
	delt = malloc(sizeof(double));
	d_delt = malloc(sizeof(double));

	fp = fopen("userinput_data.txt", "r");
	 fscanf(fp, "%d", &start_iter);
	 fscanf(fp, "%d", &no_of_iter);
	 fscanf(fp, "%d", &ns);
	 fscanf(fp, "%d", u1initial);
	 fscanf(fp, "%d", u2initial);
	 fscanf(fp, "%d", u3initial);
	 fscanf(fp, "%d", u4initial);
	 fscanf(fp, "%d", u5initial);
	 fscanf(fp, "%d", re);
	 fscanf(fp, "%d", pr);
	 fscanf(fp, "%d", ma);
	 fscanf(fp, "%d", delt);
	 fscanf(fp, "%d", &pg);
   fclose(fp);

	double *scont, *tdash, *tcont, *mudash;
	scont = malloc(sizeof(double));
	tdash = malloc(sizeof(double));
	tcont = malloc(sizeof(double));
	mudash = malloc(sizeof(double));
//tref is taken as 293; scont is taken as 110.4
	*scont = 110.4 / 293.0;
	*tcont = 1.0 / (*ma * *ma * ga);

	char filename[100], filename1[100];

	int d1, d2, d3, d4, d5;
//ele indicates interior element count
	int ele, l1, l2, l3, l4, m1, m2, m3, m4, n1, n2, n3, n4, s1, s2, s3;

	int id, prs;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size( MPI_COMM_WORLD, &prs);
	MPI_Status status, stat;
	MPI_Request request;
//communicator creation for slaves[1-n]
	MPI_Comm_group(MPI_COMM_WORLD, &group_a);
	MPI_Group_excl(group_a, 1, slaves, &group_b);
	MPI_Comm_create(MPI_COMM_WORLD, group_b, &comm_slaves);

//printf("t1\n");
//MPI_Datatype
	MPI_Datatype trans_xycoord;
	MPI_Datatype type1[3] =
	{ MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };
	int blocklen1[3] =
	{ 1, 1, 1 };
	MPI_Aint disp1[3] =
	{ 0, sizeof(double), 2 * sizeof(double) };
	MPI_Type_struct(3, blocklen1, disp1, type1, &trans_xycoord);
	MPI_Type_commit(&trans_xycoord);

	MPI_Datatype trans_eleminfo;
	MPI_Datatype type2[4] =
	{ MPI_INT, MPI_INT, MPI_INT, MPI_INT };
	int blocklen2[4] =
	{ 5, 5, 5, 5 };
	MPI_Aint disp2[4] =
	{ 0, 5 * sizeof(int), 10 * sizeof(int), 15 * sizeof(int) };
	MPI_Type_struct(4, blocklen2, disp2, type2, &trans_eleminfo);
	MPI_Type_commit(&trans_eleminfo);

	MPI_Datatype trans_primitive;
	MPI_Datatype type5[5] =
	{ MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };
	int blocklen5[5] =
	{ 1, 1, 1, 1, 1 };
	MPI_Aint disp5[5] =
	{ 0, sizeof(double), 2 * sizeof(double), 3 * sizeof(double), 4
			* sizeof(double) };
	MPI_Type_struct(5, blocklen5, disp5, type5, &trans_primitive);
	MPI_Type_commit(&trans_primitive);

	MPI_Datatype trans_derv;
	MPI_Datatype type6[15] =
	{ MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
			MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
			MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };
	int blocklen6[15] =
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	MPI_Aint disp6[15] =
	{ 0, sizeof(double), 2 * sizeof(double), 3 * sizeof(double), 4
			* sizeof(double),\
 5 * sizeof(double), 6 * sizeof(double), 7
			* sizeof(double), 8 * sizeof(double), 9 * sizeof(double), 10
			* sizeof(double),\
 11 * sizeof(double), 12 * sizeof(double), 13
			* sizeof(double), 14 * sizeof(double) };
	MPI_Type_struct(15, blocklen6, disp6, type6, &trans_derv);
	MPI_Type_commit(&trans_derv);

	char oneword[1000];
	int bc_set, result;

	MPI_Barrier(MPI_COMM_WORLD);

	if (id == 0)
	{

		int periodic_cnt, period, b1, b2;
		periodic_cnt = 0;
		period = 0;
		b1 = 1;
		b2 = 2;

		if (start_iter==1)
		{
			fg = 1;
		}
		else
		{
			fg = 0;
		}

//Reading gambit file

		fp = pre_meshfile_read(oneword, temp, result, prs, filename, fp, &i,
				&no_of_nodes, &int_ele, &bc_set, &ele, &node, xycoord, eleminfo,
				fp1, &d1, &periodic_cnt, &j, &iele, &jface, &k);

//memory allocation

		int **elem_neighbour;
		elem_neighbour = (int **) malloc(5 * sizeof(int *));
		for (i = 1; i <= 5; i++)
		{
			elem_neighbour[i] = (int *) malloc(int_ele * sizeof(int));
			if (elem_neighbour[i] == NULL)
				printf("memory failed\n");
		}

		primitive = (struct node5*) malloc(int_ele * sizeof(struct node5));
		if (primitive == NULL)
			printf("memory failed4\n");

		derv = (struct node6*) malloc(int_ele * sizeof(struct node6));
		if (derv == NULL)
			printf("memory failed4\n");

		double *max;
		int *loc_ele_surr_pts;
		max = (double *) malloc(5 * sizeof(double));
		loc_ele_surr_pts = (int *) malloc((no_of_nodes + 2) * sizeof(int));

		if (loc_ele_surr_pts == NULL)
			printf("loc_ele_surr_pts memory failed4\n");

		int *zonetype, t1, t2, t3, t4;
		zonetype = (int *) malloc(101 * sizeof(int));
		char c1;
		int temp1, temp2;

		int *periodic_pair1, *periodic_pair2, *periodic_face,
				*periodic_elem_pair1, *periodic_elem_pair2, *periodic_mapping;
		periodic_pair1 = (int *) malloc((periodic_cnt + 1) * sizeof(int));
		periodic_pair2 = (int *) malloc((periodic_cnt + 1) * sizeof(int));
		periodic_elem_pair1 = (int *) malloc((periodic_cnt + 1) * sizeof(int));
		periodic_elem_pair2 = (int *) malloc((periodic_cnt + 1) * sizeof(int));
		periodic_mapping = (int *) malloc((int_ele + 1) * sizeof(int));

		if (periodic_pair1 == NULL)
			printf("periodic_pair1 memory failed\n");

		if (periodic_pair2 == NULL)
			printf("periodic_pair2 memory failed\n");

		if (periodic_elem_pair1 == NULL)
			printf("periodic_elem_pair1 memory failed\n");

		if (periodic_elem_pair2 == NULL)
			printf("periodic_elem_pair2 memory failed\n");

		if (periodic_mapping == NULL)
			printf("periodic_mapping memory failed\n");

		for (i = 1; i <= int_ele; i++)
		{
			periodic_mapping[i] = 0;
		}


		if (pg == 1)
		{

			fp = pre_meshdata_periodic(periodic_cnt, oneword, c1, t1, b1, b2,
					t2, period, temp1, temp2, t3, t4, check, prs, filename, ele,
					fp, zonetype, periodic_pair1, periodic_pair2,
					&periodic_face, &i, &j, periodic_elem_pair1,
					periodic_elem_pair2, periodic_mapping, &iele);
		}

//comment the next six lines if there r no periodic boundary conditions

		free(periodic_pair1);
		free(periodic_pair2);
		free(periodic_face);
		free(periodic_elem_pair1);
		free(periodic_elem_pair2);
		free(periodic_mapping);

//loc_ele_surr_pts[ele_no+2]

		iele = neigh_iden(iele, ele, jj, node, dstore, dstore1, point, jele,
				d_node, ii, prs, filename, &s1, &s2, &s3, eleminfo, &i,
				loc_ele_surr_pts, &iface, elem_neighbour, fp);

		iele = initial_cond_read_write(fg, iele, ele, fp, dummy, u1initial,
				u2initial, u3initial, u4initial, u5initial, primitive, derv);
	} //id==0

	MPI_Bcast(&fg, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&pg, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ele, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&node, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&bc_set, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	if (id == 0)
	{
		int *count;
		count = (int *) malloc((prs + 1) * sizeof(int));

		double *max_id;
		max_id = (double *) malloc(prs * sizeof(double));

		double *mxm;
		mxm = (double *) malloc(sizeof(double));

		for (i = 1; i < prs; i++)
		{
			MPI_Recv(&count[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
		}

		int **glelem;

		glelem = (int **) malloc((prs) * sizeof(int *));
		if (glelem == NULL)
			printf("memory failed5\n");
		for (j = 1; j < (prs); j++)
		{
			glelem[j] = (int *) malloc((count[j] + 1) * sizeof(int));
			if (glelem == NULL)
				printf("memory failed5\n");
		}

		for (i = 1; i < prs; i++)
		{
			position = 0;
			position1 = 0;
			position2 = 0;
			position3 = 0;

			MPI_Recv(&d_cnt[1], 2, MPI_INT, i, i, MPI_COMM_WORLD, &status);

			MPI_Pack_size(d_cnt[1], MPI_INT, MPI_COMM_WORLD, &maxsize);
			buffer = malloc(maxsize);

			MPI_Pack_size(d_cnt[2], MPI_INT, MPI_COMM_WORLD, &maxsize1);
			buffer1 = malloc(maxsize1);

			MPI_Recv(buffer, maxsize, MPI_PACKED, i, i, MPI_COMM_WORLD,
					&status);
			MPI_Recv(buffer1, maxsize1, MPI_PACKED, i, i, MPI_COMM_WORLD,
					&status);

			MPI_Pack_size(d_cnt[1], trans_eleminfo, MPI_COMM_WORLD, &maxsize4);
			MPI_Pack_size(d_cnt[1], trans_primitive, MPI_COMM_WORLD, &maxsize5);
			MPI_Pack_size(d_cnt[1], trans_derv, MPI_COMM_WORLD, &maxsize6);
			maxsize2 = maxsize4 + maxsize5 + maxsize6;
			buffer2 = malloc(maxsize2);

			if (buffer2 == NULL)
				printf("buffer2 memory failed\n");

			MPI_Pack_size(d_cnt[2], trans_xycoord, MPI_COMM_WORLD, &maxsize3);
			buffer3 = malloc(maxsize3);

			if (buffer3 == NULL)
				printf("buffer3 memory failed\n");

			for (j = 1; j <= d_cnt[1]; j++)
			{
				MPI_Unpack(buffer, maxsize, &position, &glelem[i][j], 1,
						MPI_INT, MPI_COMM_WORLD);
				d_ele = glelem[i][j];
				MPI_Pack(&eleminfo[d_ele], 1, trans_eleminfo, buffer2, maxsize2,
						&position2, MPI_COMM_WORLD);
				MPI_Pack(&primitive[d_ele], 1, trans_primitive, buffer2,
						maxsize2, &position2, MPI_COMM_WORLD);
				MPI_Pack(&derv[d_ele], 1, trans_derv, buffer2, maxsize2,
						&position2, MPI_COMM_WORLD);
			}

			for (j = 1; j <= d_cnt[2]; j++)
			{
				MPI_Unpack(buffer1, maxsize1, &position1, &d_ele, 1, MPI_INT,
						MPI_COMM_WORLD);
				MPI_Pack(&xycoord[d_ele], 1, trans_xycoord, buffer3, maxsize3,
						&position3, MPI_COMM_WORLD);
			}

			MPI_Send(buffer2, position2, MPI_PACKED, i, i, MPI_COMM_WORLD);
			MPI_Send(buffer3, position3, MPI_PACKED, i, i, MPI_COMM_WORLD);

			free(buffer);
			free(buffer1);
			free(buffer2);
			free(buffer3);

			buffer = NULL;
			buffer1 = NULL;
			buffer2 = NULL;
			buffer3 = NULL;

		}

		n = result_file_write(n, start_iter, no_of_iter, prs, ns, trans_primitive,
				position1, d_ele, filename, ele, node, trans_derv, fp4, &i,
				max_id, &status, mxm, maximum, delt, count, &maxsize, buffer,
				&j, glelem, primitive, fp2, fp5, xycoord, eleminfo, derv, fp3);
	} //id==0

	if (id > 0)
	{

		int *gb, *sd_elem, *elem2proc;
		gb = (int *) malloc((ele + 1) * sizeof(int));
		elem2proc = (int *) malloc((ele + 1) * sizeof(int));

		if (gb == NULL)
			printf("gb_slave memory failed\n");

		if (elem2proc == NULL)
			printf("elem2proc_slave memory failed\n");

//initializing gb to zero

		for (i = 1; i <= ele; i++)
		{
			gb[i] = 0;
		}

		double *z1, *z2, *z3, *sum, *deno, *d_delt, *CFL_no, *aa1, *bb1, *cc1,
				*dd1, *t, *normal_x, *normal_y, *normal_z, *maxprime, *del, *nu,
				*d_nu, *cfl_iele;
		z1 = malloc(sizeof(double));
		z2 = malloc(sizeof(double));
		z3 = malloc(sizeof(double));
		sum = malloc(sizeof(double));
		deno = malloc(sizeof(double));
		d_delt = malloc(sizeof(double));
		CFL_no = malloc(sizeof(double));
		cfl_iele = malloc(sizeof(double));
		aa1 = malloc(sizeof(double));
		bb1 = malloc(sizeof(double));
		cc1 = malloc(sizeof(double));
		dd1 = malloc(sizeof(double));
		t = malloc(sizeof(double));
		normal_x = malloc(sizeof(double));
		normal_y = malloc(sizeof(double));
		normal_z = malloc(sizeof(double));
		maxprime = malloc(sizeof(double));
		del = malloc(sizeof(double));
		nu = malloc(sizeof(double));
		d_nu = malloc(sizeof(double));
		*CFL_no = 1.0;

		double *max;
		max = (double *) malloc(6 * sizeof(double));

		int count[3];

		int d_store, check;
		int *id_proc, *id_periodic_proc, *pp2id;
		int *sd_cnt, *send_cnt, subdomain_cnt, *pp_cnt, *p_recvcnt;
		int boun = 0;
		int sd = 0;
		int lele_cnt;
		kk = 0;
		int **lelem_send, **lelem_recv, *k2id;
		int d_store1, d_store2, d_store3;
		int n11, n22;

		double *fvm, *gvm, *evm, *flux, **umdash;
		fvm = (double *) malloc(4 * sizeof(double));
		gvm = (double *) malloc(4 * sizeof(double));
		evm = (double *) malloc(4 * sizeof(double));
		flux = (double *) malloc(5 * sizeof(double));
		umdash = (double **) malloc(4 * sizeof(double));

		for (i = 0; i < 4; i++)
		{
			umdash[i] = (double *) malloc(5 * sizeof(double));
		}

		k2id = (int *) malloc((prs + 1) * sizeof(int));
		id_proc = (int *) malloc((prs + 1) * sizeof(int));
		id_periodic_proc = (int *) malloc((prs + 1) * sizeof(int));
		pp2id = (int *) malloc((prs + 1) * sizeof(int));
		sd_cnt = (int *) malloc((prs + 1) * sizeof(int));
		send_cnt = (int *) malloc((prs + 1) * sizeof(int));
		pp_cnt = (int *) malloc((prs + 1) * sizeof(int));
		p_recvcnt = (int *) malloc((prs + 1) * sizeof(int));

		double *xc_side, *yc_side, *zc_side, *f1_side, *f2_side, *f3_side,
				*f4_side, *um, *uxm, *uym, *uzm, *utm, *fm, *fxm, *fym, *fzm,
				*ftm, *gm, *gxm, *gym, *gzm, *gtm, *em, *exm, *eym, *ezm, *etm;
		xc_side = malloc(sizeof(double));
		yc_side = malloc(sizeof(double));
		zc_side = malloc(sizeof(double));
		f1_side = malloc(sizeof(double));
		f2_side = malloc(sizeof(double));
		f3_side = malloc(sizeof(double));
		f4_side = malloc(sizeof(double));
		um = malloc(sizeof(double));
		uxm = malloc(sizeof(double));
		uym = malloc(sizeof(double));
		uzm = malloc(sizeof(double));
		utm = malloc(sizeof(double));
		fm = malloc(sizeof(double));
		fxm = malloc(sizeof(double));
		fym = malloc(sizeof(double));
		fzm = malloc(sizeof(double));
		ftm = malloc(sizeof(double));
		gm = malloc(sizeof(double));
		gxm = malloc(sizeof(double));
		gym = malloc(sizeof(double));
		gzm = malloc(sizeof(double));
		gtm = malloc(sizeof(double));
		em = malloc(sizeof(double));
		exm = malloc(sizeof(double));
		eym = malloc(sizeof(double));
		ezm = malloc(sizeof(double));
		etm = malloc(sizeof(double));
		dummy = malloc(sizeof(double));
		maximum = malloc(sizeof(double));

		double **u1_side, **u2_side, **u3_side, **u4_side, **u5_side, **ux_side,
				**vx_side, **wx_side, **uy_side, **vy_side, **wy_side,
				**uz_side, **vz_side, **wz_side, **towxx, **towyy, **towzz,
				**towxy, **towyz, **towzx, **qx, **qy, **qz;

		u1_side = (double **) malloc(4 * sizeof(double *));
		u2_side = (double **) malloc(4 * sizeof(double *));
		u3_side = (double **) malloc(4 * sizeof(double *));
		u4_side = (double **) malloc(4 * sizeof(double *));
		u5_side = (double **) malloc(4 * sizeof(double *));
		ux_side = (double **) malloc(4 * sizeof(double *));
		vx_side = (double **) malloc(4 * sizeof(double *));
		wx_side = (double **) malloc(4 * sizeof(double *));
		uy_side = (double **) malloc(4 * sizeof(double *));
		vy_side = (double **) malloc(4 * sizeof(double *));
		wy_side = (double **) malloc(4 * sizeof(double *));
		uz_side = (double **) malloc(4 * sizeof(double *));
		vz_side = (double **) malloc(4 * sizeof(double *));
		wz_side = (double **) malloc(4 * sizeof(double *));
		towxx = (double **) malloc(4 * sizeof(double *));
		towyy = (double **) malloc(4 * sizeof(double *));
		towzz = (double **) malloc(4 * sizeof(double *));
		towxy = (double **) malloc(4 * sizeof(double *));
		towyz = (double **) malloc(4 * sizeof(double *));
		towzx = (double **) malloc(4 * sizeof(double *));
		qx = (double **) malloc(4 * sizeof(double *));
		qy = (double **) malloc(4 * sizeof(double *));
		qz = (double **) malloc(4 * sizeof(double *));

		for (i = 0; i < 4; i++)
		{
			u1_side[i] = (double *) malloc(5 * sizeof(double));
			u2_side[i] = (double *) malloc(5 * sizeof(double));
			u3_side[i] = (double *) malloc(5 * sizeof(double));
			u4_side[i] = (double *) malloc(5 * sizeof(double));
			u5_side[i] = (double *) malloc(5 * sizeof(double));
			ux_side[i] = (double *) malloc(5 * sizeof(double));
			vx_side[i] = (double *) malloc(5 * sizeof(double));
			wx_side[i] = (double *) malloc(5 * sizeof(double));
			uy_side[i] = (double *) malloc(5 * sizeof(double));
			vy_side[i] = (double *) malloc(5 * sizeof(double));
			wy_side[i] = (double *) malloc(5 * sizeof(double));
			uz_side[i] = (double *) malloc(5 * sizeof(double));
			vz_side[i] = (double *) malloc(5 * sizeof(double));
			wz_side[i] = (double *) malloc(5 * sizeof(double));
			towxx[i] = (double *) malloc(5 * sizeof(double));
			towyy[i] = (double *) malloc(5 * sizeof(double));
			towzz[i] = (double *) malloc(5 * sizeof(double));
			towxy[i] = (double *) malloc(5 * sizeof(double));
			towyz[i] = (double *) malloc(5 * sizeof(double));
			towzx[i] = (double *) malloc(5 * sizeof(double));
			qx[i] = (double *) malloc(5 * sizeof(double));
			qy[i] = (double *) malloc(5 * sizeof(double));
			qz[i] = (double *) malloc(5 * sizeof(double));
		}

//element allocation to processor id
		i = 1;
		sprintf(filename, "proc2elem_%d.txt", id);
		fp = fopen(filename, "r");
		while (!feof(fp))
		{
			fscanf(fp, "%d", &proc);
			elem2proc[i] = proc;
			i = i + 1;
		}
		fclose(fp);

		lele_cnt = 0;
		;

//for determining the size of lele_cnt
		for (i = 1; i <= ele; i++)
		{
			if (elem2proc[i] == id - 1)
			{
				lele_cnt = lele_cnt + 1;
			}
		}

		dummy = 0;

//local numbering of element neighbours
		int **lelem_neighbour;
		double **xr, **yr, **zr, **num_dd, *c;

		lelem_neighbour = (int **) malloc(5 * sizeof(int *));
		xr = (double **) malloc(5 * sizeof(double *));
		yr = (double **) malloc(5 * sizeof(double *));
		zr = (double **) malloc(5 * sizeof(double *));
		num_dd = (double **) malloc(5 * sizeof(double *));
		c = (double *) malloc((lele_cnt + 1) * sizeof(double));

		if (lelem_neighbour == NULL)
			printf("memory failed\n");

		for (i = 1; i < 5; i++)
		{
			lelem_neighbour[i] = (int *) malloc((lele_cnt + 1) * sizeof(int));
			xr[i] = (double *) malloc((lele_cnt + 1) * sizeof(double));
			if (xr[i] == NULL)
				printf("xr_slave memory failed\n");
			yr[i] = (double *) malloc((lele_cnt + 1) * sizeof(double));
			if (yr[i] == NULL)
				printf("yr_slave memory failed\n");
			zr[i] = (double *) malloc((lele_cnt + 1) * sizeof(double));
			if (zr[i] == NULL)
				printf("zr_slave memory failed\n");
			num_dd[i] = (double *) malloc((lele_cnt + 1) * sizeof(double));
			if (num_dd[i] == NULL)
				printf("num_dd_slave memory failed\n");
		}

		int *lelem;
		lelem = (int *) malloc((ele + 1) * sizeof(int));
		if (lelem == NULL)
			printf("lelem_slave memory failed\n");

		lele_cnt = slave_global2local_numbering(lele_cnt, filename, id, ele, n4,
				d_store, check, fp1, &i, elem2proc, lelem, gb, &n1, &n2, &n3,
				lelem_neighbour, dummy, &sd, &kk, id_proc, k2id, sd_cnt, fp,
				&j);

//finding the subdomain neigh elements:lelem_send-local element numbers of subdomain boundary elements

		int* group = slave_subdomain_search(kk, sd, filename, id, check, d_ele,
				lele_cnt, d_store, &lelem_send, &i, sd_cnt, send_cnt, sd_elem,
				fp, &n1, &n2, &n3, id_proc, &j, gb, &temp, &subdomain_cnt,
				lelem);

//Find subdomain boundary data

		int belem_cnt[bc_set + 1], bcnt_iden[bc_set + 1], b_cnt;
		b_cnt = 0;

		i = sd_boun_data_find1(i, bc_set, belem_cnt, bcnt_iden, filename, id,
				oneword, check, fp, &d1, &j, &d2, &temp, elem2proc, &b_cnt);

		int **boun_elem, *cnt_boun_type, *boun_type, **boun_face;
		boun_elem = (int **) malloc((b_cnt + 1) * sizeof(int *));
		boun_face = (int **) malloc((b_cnt + 1) * sizeof(int *));
		cnt_boun_type = (int *) malloc((b_cnt + 1) * sizeof(int));
		boun_type = (int *) malloc((b_cnt + 1) * sizeof(int));

		b_cnt = 0;

		int l_periodic_cnt = sd_boun_data_find2(filename, id, bc_set, bcnt_iden,
				oneword, temp, belem_cnt, subdomain_cnt, pg, fp, &i, &d1, &j,
				&k, &b_cnt, cnt_boun_type, boun_elem, boun_face, boun_type,
				elem2proc, gb, &tot_ele, &n1);

		int p_tot_ele;

		p_tot_ele = tot_ele + l_periodic_cnt;

//memory allocation for primitive and dervative

		eleminfo = (struct node2*) malloc(
				(lele_cnt + 1) * sizeof(struct node2));

		if (eleminfo == NULL)
			printf("eleminfo_slave memory failed\n");

		primitive = (struct node5*) malloc(p_tot_ele * sizeof(struct node5));
		if (primitive == NULL)
			printf("primitive_slave memory failed\n");

		derv = (struct node6*) malloc(p_tot_ele * sizeof(struct node6));
		if (derv == NULL)
			printf("derv_slave memory failed\n");

		secondary = (struct node7*) malloc(p_tot_ele * sizeof(struct node7));
		if (secondary == NULL)
			printf("secondary_slave memory failed\n");

		MPI_Datatype trans_secondary;
		MPI_Datatype type7[80] =
		{ MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
				MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
				MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
				MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
				MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
				MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
				MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
				MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
				MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
				MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
				MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
				MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
				MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
				MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
				MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
				MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };
		int blocklen7[80] =
		{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
				1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
				1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
				1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
		MPI_Aint disp7[80] =
		{ 0, sizeof(double), 2 * sizeof(double), 3 * sizeof(double), 4
				* sizeof(double), 5 * sizeof(double), 6 * sizeof(double), 7
				* sizeof(double), 8 * sizeof(double), 9 * sizeof(double), 10
				* sizeof(double), 11 * sizeof(double), 12 * sizeof(double), 13
				* sizeof(double), 14 * sizeof(double), 15 * sizeof(double), 16
				* sizeof(double), 17 * sizeof(double), 18 * sizeof(double), 19
				* sizeof(double), 20 * sizeof(double), 21 * sizeof(double), 22
				* sizeof(double), 23 * sizeof(double), 24 * sizeof(double), 25
				* sizeof(double), 26 * sizeof(double), 27 * sizeof(double), 28
				* sizeof(double), 29 * sizeof(double), 30 * sizeof(double), 31
				* sizeof(double), 32 * sizeof(double), 33 * sizeof(double), 34
				* sizeof(double), 35 * sizeof(double), 36 * sizeof(double), 37
				* sizeof(double), 38 * sizeof(double), 39 * sizeof(double), 40
				* sizeof(double), 41 * sizeof(double), 42 * sizeof(double), 43
				* sizeof(double), 44 * sizeof(double), 45 * sizeof(double), 46
				* sizeof(double), 47 * sizeof(double), 48 * sizeof(double), 49
				* sizeof(double), 50 * sizeof(double), 51 * sizeof(double), 52
				* sizeof(double), 53 * sizeof(double), 54 * sizeof(double), 55
				* sizeof(double), 56 * sizeof(double), 57 * sizeof(double), 58
				* sizeof(double), 59 * sizeof(double), 60 * sizeof(double), 61
				* sizeof(double), 62 * sizeof(double), 63 * sizeof(double), 64
				* sizeof(double), 65 * sizeof(double), 66 * sizeof(double), 67
				* sizeof(double), 68 * sizeof(double), 69 * sizeof(double), 70
				* sizeof(double), 71 * sizeof(double), 72 * sizeof(double), 73
				* sizeof(double), 74 * sizeof(double), 75 * sizeof(double), 76
				* sizeof(double), 77 * sizeof(double), 78 * sizeof(double), 79
				* sizeof(double) };
		MPI_Type_struct(80, blocklen7, disp7, type7, &trans_secondary);
		MPI_Type_commit(&trans_secondary);

		double *u1_p, *ux1_p, *uy1_p, *uz1_p, *ut1_p, *f1_p, *fx1_p, *fy1_p,
				*fz1_p, *ft1_p, *g1_p, *gx1_p, *gy1_p, *gz1_p, *gt1_p, *e1_p,
				*ex1_p, *ey1_p, *ez1_p, *et1_p, *u2_p, *ux2_p, *uy2_p, *uz2_p,
				*ut2_p, *f2_p, *fx2_p, *fy2_p, *fz2_p, *ft2_p, *g2_p, *gx2_p,
				*gy2_p, *gz2_p, *gt2_p, *e2_p, *ex2_p, *ey2_p, *ez2_p, *et2_p,
				*u3_p, *ux3_p, *uy3_p, *uz3_p, *ut3_p, *f3_p, *fx3_p, *fy3_p,
				*fz3_p, *ft3_p, *g3_p, *gx3_p, *gy3_p, *gz3_p, *gt3_p, *e3_p,
				*ex3_p, *ey3_p, *ez3_p, *et3_p, *u4_p, *ux4_p, *uy4_p, *uz4_p,
				*ut4_p, *f4_p, *fx4_p, *fy4_p, *fz4_p, *ft4_p, *g4_p, *gx4_p,
				*gy4_p, *gz4_p, *gt4_p, *e4_p, *ex4_p, *ey4_p, *ez4_p, *et4_p,
				*u5_p, *ux5_p, *uy5_p, *uz5_p, *ut5_p, *f5_p, *fx5_p, *fy5_p,
				*fz5_p, *ft5_p, *g5_p, *gx5_p, *gy5_p, *gz5_p, *gt5_p, *e5_p,
				*ex5_p, *ey5_p, *ez5_p, *et5_p;

		u1_p = (double *) malloc(tot_ele * sizeof(double));
		if (u1_p == NULL)
			printf("u1_p_slave memory failed\n");
		ux1_p = (double *) malloc(tot_ele * sizeof(double));
		if (ux1_p == NULL)
			printf("ux1_p_slave memory failed\n");
		uy1_p = (double *) malloc(tot_ele * sizeof(double));
		if (uy1_p == NULL)
			printf("uy1_p_slave memory failed\n");
		uz1_p = (double *) malloc(tot_ele * sizeof(double));
		if (uz1_p == NULL)
			printf("uz1_p_slave memory failed\n");
		ut1_p = (double *) malloc(tot_ele * sizeof(double));
		if (ut1_p == NULL)
			printf("ut1_p_slave memory failed\n");
		f1_p = (double *) malloc(tot_ele * sizeof(double));
		if (f1_p == NULL)
			printf("f1_p_slave memory failed\n");
		fx1_p = (double *) malloc(tot_ele * sizeof(double));
		if (fx1_p == NULL)
			printf("fx1_p_slave memory failed\n");
		fy1_p = (double *) malloc(tot_ele * sizeof(double));
		if (fy1_p == NULL)
			printf("fy1_p_slave memory failed\n");
		fz1_p = (double *) malloc(tot_ele * sizeof(double));
		if (fz1_p == NULL)
			printf("fz1_p_slave memory failed\n");
		ft1_p = (double *) malloc(tot_ele * sizeof(double));
		if (ft1_p == NULL)
			printf("ft1_p_slave memory failed\n");
		g1_p = (double *) malloc(tot_ele * sizeof(double));
		if (g1_p == NULL)
			printf("g1_p_slave memory failed\n");
		gx1_p = (double *) malloc(tot_ele * sizeof(double));
		if (gx1_p == NULL)
			printf("gx1_p_slave memory failed\n");
		gy1_p = (double *) malloc(tot_ele * sizeof(double));
		if (gy1_p == NULL)
			printf("gy1_p_slave memory failed\n");
		gz1_p = (double *) malloc(tot_ele * sizeof(double));
		if (gz1_p == NULL)
			printf("gz1_p_slave memory failed\n");
		gt1_p = (double *) malloc(tot_ele * sizeof(double));
		if (gt1_p == NULL)
			printf("gt1_p_slave memory failed\n");
		e1_p = (double *) malloc(tot_ele * sizeof(double));
		if (e1_p == NULL)
			printf("e1_p_slave memory failed\n");
		ex1_p = (double *) malloc(tot_ele * sizeof(double));
		if (ex1_p == NULL)
			printf("ex1_p_slave memory failed\n");
		ey1_p = (double *) malloc(tot_ele * sizeof(double));
		if (ey1_p == NULL)
			printf("ey1_p_slave memory failed\n");
		ez1_p = (double *) malloc(tot_ele * sizeof(double));
		if (ez1_p == NULL)
			printf("ez1_p_slave memory failed\n");
		et1_p = (double *) malloc(tot_ele * sizeof(double));
		if (et1_p == NULL)
			printf("et1_p_slave memory failed\n");
		u2_p = (double *) malloc(tot_ele * sizeof(double));
		if (u2_p == NULL)
			printf("u2_p_slave memory failed\n");
		ux2_p = (double *) malloc(tot_ele * sizeof(double));
		if (ux2_p == NULL)
			printf("ux2_p_slave memory failed\n");
		uy2_p = (double *) malloc(tot_ele * sizeof(double));
		if (uy2_p == NULL)
			printf("uy2_p_slave memory failed\n");
		uz2_p = (double *) malloc(tot_ele * sizeof(double));
		if (uz2_p == NULL)
			printf("uz2_p_slave memory failed\n");
		ut2_p = (double *) malloc(tot_ele * sizeof(double));
		if (ut2_p == NULL)
			printf("ut2_p_slave memory failed\n");
		f2_p = (double *) malloc(tot_ele * sizeof(double));
		if (f2_p == NULL)
			printf("f2_p_slave memory failed\n");
		fx2_p = (double *) malloc(tot_ele * sizeof(double));
		if (fx2_p == NULL)
			printf("fx2_p_slave memory failed\n");
		fy2_p = (double *) malloc(tot_ele * sizeof(double));
		if (fy2_p == NULL)
			printf("fy2_p_slave memory failed\n");
		fz2_p = (double *) malloc(tot_ele * sizeof(double));
		if (fz2_p == NULL)
			printf("fz2_p_slave memory failed\n");
		ft2_p = (double *) malloc(tot_ele * sizeof(double));
		if (ft2_p == NULL)
			printf("ft2_p_slave memory failed\n");
		g2_p = (double *) malloc(tot_ele * sizeof(double));
		if (g2_p == NULL)
			printf("g2_p_slave memory failed\n");
		gx2_p = (double *) malloc(tot_ele * sizeof(double));
		if (gx2_p == NULL)
			printf("gx2_p_slave memory failed\n");
		gy2_p = (double *) malloc(tot_ele * sizeof(double));
		if (gy2_p == NULL)
			printf("gy2_p_slave memory failed\n");
		gz2_p = (double *) malloc(tot_ele * sizeof(double));
		if (gz2_p == NULL)
			printf("gz2_p_slave memory failed\n");
		gt2_p = (double *) malloc(tot_ele * sizeof(double));
		if (gt2_p == NULL)
			printf("gt2_p_slave memory failed\n");
		e2_p = (double *) malloc(tot_ele * sizeof(double));
		if (e2_p == NULL)
			printf("e2_p_slave memory failed\n");
		ex2_p = (double *) malloc(tot_ele * sizeof(double));
		if (ex2_p == NULL)
			printf("ex2_p_slave memory failed\n");
		ey2_p = (double *) malloc(tot_ele * sizeof(double));
		if (ey2_p == NULL)
			printf("ey2_p_slave memory failed\n");
		ez2_p = (double *) malloc(tot_ele * sizeof(double));
		if (ez2_p == NULL)
			printf("ez2_p_slave memory failed\n");
		et2_p = (double *) malloc(tot_ele * sizeof(double));
		if (et2_p == NULL)
			printf("et2_p_slave memory failed\n");

		u3_p = (double *) malloc(tot_ele * sizeof(double));
		if (u3_p == NULL)
			printf("u3_p_slave memory failed\n");
		ux3_p = (double *) malloc(tot_ele * sizeof(double));
		if (ux3_p == NULL)
			printf("ux3_p_slave memory failed\n");
		uy3_p = (double *) malloc(tot_ele * sizeof(double));
		if (uy3_p == NULL)
			printf("uy3_p_slave memory failed\n");
		uz3_p = (double *) malloc(tot_ele * sizeof(double));
		if (uz3_p == NULL)
			printf("uz3_p_slave memory failed\n");
		f3_p = (double *) malloc(tot_ele * sizeof(double));
		if (f3_p == NULL)
			printf("f3_p_slave memory failed\n");
		fx3_p = (double *) malloc(tot_ele * sizeof(double));
		if (fx3_p == NULL)
			printf("fx3_p_slave memory failed\n");
		fy3_p = (double *) malloc(tot_ele * sizeof(double));
		if (fy3_p == NULL)
			printf("fy3_p_slave memory failed\n");
		fz3_p = (double *) malloc(tot_ele * sizeof(double));
		if (fz3_p == NULL)
			printf("fz3_p_slave memory failed\n");
		ft3_p = (double *) malloc(tot_ele * sizeof(double));
		if (ft3_p == NULL)
			printf("ft3_p_slave memory failed\n");
		g3_p = (double *) malloc(tot_ele * sizeof(double));
		if (g3_p == NULL)
			printf("g3_p_slave memory failed\n");
		gx3_p = (double *) malloc(tot_ele * sizeof(double));
		if (gx3_p == NULL)
			printf("gx3_p_slave memory failed\n");
		gy3_p = (double *) malloc(tot_ele * sizeof(double));
		if (gy3_p == NULL)
			printf("gy3_p_slave memory failed\n");
		gz3_p = (double *) malloc(tot_ele * sizeof(double));
		if (gz3_p == NULL)
			printf("gz3_p_slave memory failed\n");
		gt3_p = (double *) malloc(tot_ele * sizeof(double));
		if (gt3_p == NULL)
			printf("gt3_p_slave memory failed\n");
		e3_p = (double *) malloc(tot_ele * sizeof(double));
		if (e3_p == NULL)
			printf("e3_p_slave memory failed\n");
		ex3_p = (double *) malloc(tot_ele * sizeof(double));
		if (ex3_p == NULL)
			printf("ex3_p_slave memory failed\n");
		ey3_p = (double *) malloc(tot_ele * sizeof(double));
		if (ey3_p == NULL)
			printf("ey3_p_slave memory failed\n");
		ez3_p = (double *) malloc(tot_ele * sizeof(double));
		if (ez3_p == NULL)
			printf("ez3_p_slave memory failed\n");
		et3_p = (double *) malloc(tot_ele * sizeof(double));
		if (et3_p == NULL)
			printf("et3_p_slave memory failed\n");

		u4_p = (double *) malloc(tot_ele * sizeof(double));
		if (u4_p == NULL)
			printf("u4_p_slave memory failed\n");
		ux4_p = (double *) malloc(tot_ele * sizeof(double));
		if (ux4_p == NULL)
			printf("ux4_p_slave memory failed\n");
		uy4_p = (double *) malloc(tot_ele * sizeof(double));
		if (uy4_p == NULL)
			printf("uy4_p_slave memory failed\n");
		uz4_p = (double *) malloc(tot_ele * sizeof(double));
		if (uz4_p == NULL)
			printf("uz4_p_slave memory failed\n");
		f4_p = (double *) malloc(tot_ele * sizeof(double));
		if (f4_p == NULL)
			printf("f4_p_slave memory failed\n");
		fx4_p = (double *) malloc(tot_ele * sizeof(double));
		if (fx4_p == NULL)
			printf("fx4_p_slave memory failed\n");
		fy4_p = (double *) malloc(tot_ele * sizeof(double));
		if (fy4_p == NULL)
			printf("fy4_p_slave memory failed\n");
		fz4_p = (double *) malloc(tot_ele * sizeof(double));
		if (fz4_p == NULL)
			printf("fz4_p_slave memory failed\n");
		ft4_p = (double *) malloc(tot_ele * sizeof(double));
		if (ft4_p == NULL)
			printf("ft4_p_slave memory failed\n");
		g4_p = (double *) malloc(tot_ele * sizeof(double));
		if (g4_p == NULL)
			printf("g4_p_slave memory failed\n");
		gx4_p = (double *) malloc(tot_ele * sizeof(double));
		if (gx4_p == NULL)
			printf("gx4_p_slave memory failed\n");
		gy4_p = (double *) malloc(tot_ele * sizeof(double));
		if (gy4_p == NULL)
			printf("gy4_p_slave memory failed\n");
		gz4_p = (double *) malloc(tot_ele * sizeof(double));
		if (gz4_p == NULL)
			printf("gz4_p_slave memory failed\n");
		gt4_p = (double *) malloc(tot_ele * sizeof(double));
		if (gt4_p == NULL)
			printf("gt4_p_slave memory failed\n");
		e4_p = (double *) malloc(tot_ele * sizeof(double));
		if (e4_p == NULL)
			printf("e4_p_slave memory failed\n");
		ex4_p = (double *) malloc(tot_ele * sizeof(double));
		if (ex4_p == NULL)
			printf("ex4_p_slave memory failed\n");
		ey4_p = (double *) malloc(tot_ele * sizeof(double));
		if (ey4_p == NULL)
			printf("ey4_p_slave memory failed\n");
		ez4_p = (double *) malloc(tot_ele * sizeof(double));
		if (ez4_p == NULL)
			printf("ez4_p_slave memory failed\n");
		et4_p = (double *) malloc(tot_ele * sizeof(double));
		if (et4_p == NULL)
			printf("et4_p_slave memory failed\n");

		u5_p = (double *) malloc(tot_ele * sizeof(double));
		if (u5_p == NULL)
			printf("u5_p_slave memory failed\n");
		ux5_p = (double *) malloc(tot_ele * sizeof(double));
		if (ux5_p == NULL)
			printf("ux5_p_slave memory failed\n");
		uy5_p = (double *) malloc(tot_ele * sizeof(double));
		if (uy5_p == NULL)
			printf("uy5_p_slave memory failed\n");
		uz5_p = (double *) malloc(tot_ele * sizeof(double));
		if (uz5_p == NULL)
			printf("uz5_p_slave memory failed\n");
		f5_p = (double *) malloc(tot_ele * sizeof(double));
		if (f5_p == NULL)
			printf("f5_p_slave memory failed\n");
		fx5_p = (double *) malloc(tot_ele * sizeof(double));
		if (fx5_p == NULL)
			printf("fx5_p_slave memory failed\n");
		fy5_p = (double *) malloc(tot_ele * sizeof(double));
		if (fy5_p == NULL)
			printf("fy5_p_slave memory failed\n");
		fz5_p = (double *) malloc(tot_ele * sizeof(double));
		if (fz5_p == NULL)
			printf("fz5_p_slave memory failed\n");
		ft5_p = (double *) malloc(tot_ele * sizeof(double));
		if (ft5_p == NULL)
			printf("ft5_p_slave memory failed\n");
		g5_p = (double *) malloc(tot_ele * sizeof(double));
		if (g5_p == NULL)
			printf("g5_p_slave memory failed\n");
		gx5_p = (double *) malloc(tot_ele * sizeof(double));
		if (gx5_p == NULL)
			printf("gx5_p_slave memory failed\n");
		gy5_p = (double *) malloc(tot_ele * sizeof(double));
		if (gy5_p == NULL)
			printf("gy5_p_slave memory failed\n");
		gz5_p = (double *) malloc(tot_ele * sizeof(double));
		if (gz5_p == NULL)
			printf("gz5_p_slave memory failed\n");
		gt5_p = (double *) malloc(tot_ele * sizeof(double));
		if (gt5_p == NULL)
			printf("gt5_p_slave memory failed\n");
		e5_p = (double *) malloc(tot_ele * sizeof(double));
		if (e5_p == NULL)
			printf("e5_p_slave memory failed\n");
		ex5_p = (double *) malloc(tot_ele * sizeof(double));
		if (ex5_p == NULL)
			printf("ex5_p_slave memory failed\n");
		ey5_p = (double *) malloc(tot_ele * sizeof(double));
		if (ey5_p == NULL)
			printf("ey5_p_slave memory failed\n");
		ez5_p = (double *) malloc(tot_ele * sizeof(double));
		if (ez5_p == NULL)
			printf("ez5_p_slave memory failed\n");
		et5_p = (double *) malloc(tot_ele * sizeof(double));
		if (et5_p == NULL)
			printf("et5_p_slave memory failed\n");

		double **xc_lateral_side1, **yc_lateral_side1, **zc_lateral_side1,
				**nx_lateral_side1, **ny_lateral_side1, **nz_lateral_side1,
				**area_lateral_side1, **xc_lateral_side2, **yc_lateral_side2,
				**zc_lateral_side2, **nx_lateral_side2, **ny_lateral_side2,
				**nz_lateral_side2, **area_lateral_side2, **xc_lateral_side3,
				**yc_lateral_side3, **zc_lateral_side3, **nx_lateral_side3,
				**ny_lateral_side3, **nz_lateral_side3, **area_lateral_side3,
				**xc_neighface, **yc_neighface, **zc_neighface, **vol_neighface,
				**xc_face, **yc_face, **zc_face, **vol_face,
				**vol_hexahedra_face, **xc_hexahedra_face, **yc_hexahedra_face,
				**zc_hexahedra_face, *xp_node, *yp_node, *zp_node, *xpstar,
				*ypstar, *zpstar, **xpstar_face, **ypstar_face, **zpstar_face,
				**theta;

		xc_lateral_side1 = (double **) malloc(5 * sizeof(double *));
		yc_lateral_side1 = (double **) malloc(5 * sizeof(double *));
		zc_lateral_side1 = (double **) malloc(5 * sizeof(double *));
		nx_lateral_side1 = (double **) malloc(5 * sizeof(double *));
		ny_lateral_side1 = (double **) malloc(5 * sizeof(double *));
		nz_lateral_side1 = (double **) malloc(5 * sizeof(double *));
		area_lateral_side1 = (double **) malloc(5 * sizeof(double *));
		xc_lateral_side2 = (double **) malloc(5 * sizeof(double *));
		yc_lateral_side2 = (double **) malloc(5 * sizeof(double *));
		zc_lateral_side2 = (double **) malloc(5 * sizeof(double *));
		nx_lateral_side2 = (double **) malloc(5 * sizeof(double *));
		ny_lateral_side2 = (double **) malloc(5 * sizeof(double *));
		nz_lateral_side2 = (double **) malloc(5 * sizeof(double *));
		area_lateral_side2 = (double **) malloc(5 * sizeof(double *));
		xc_lateral_side3 = (double **) malloc(5 * sizeof(double *));
		yc_lateral_side3 = (double **) malloc(5 * sizeof(double *));
		zc_lateral_side3 = (double **) malloc(5 * sizeof(double *));
		nx_lateral_side3 = (double **) malloc(5 * sizeof(double *));
		ny_lateral_side3 = (double **) malloc(5 * sizeof(double *));
		nz_lateral_side3 = (double **) malloc(5 * sizeof(double *));
		area_lateral_side3 = (double **) malloc(5 * sizeof(double *));
		xc_neighface = (double **) malloc(5 * sizeof(double *));
		yc_neighface = (double **) malloc(5 * sizeof(double *));
		zc_neighface = (double **) malloc(5 * sizeof(double *));
		vol_neighface = (double **) malloc(5 * sizeof(double *));
		xc_face = (double **) malloc(5 * sizeof(double *));
		yc_face = (double **) malloc(5 * sizeof(double *));
		zc_face = (double **) malloc(5 * sizeof(double *));
		vol_face = (double **) malloc(5 * sizeof(double *));
		vol_hexahedra_face = (double **) malloc(5 * sizeof(double *));
		xc_hexahedra_face = (double **) malloc(5 * sizeof(double *));
		yc_hexahedra_face = (double **) malloc(5 * sizeof(double *));
		zc_hexahedra_face = (double **) malloc(5 * sizeof(double *));
		xp_node = (double *) malloc(5 * sizeof(double));
		yp_node = (double *) malloc(5 * sizeof(double));
		zp_node = (double *) malloc(5 * sizeof(double));
		xpstar = (double *) malloc(5 * sizeof(double));
		ypstar = (double *) malloc(5 * sizeof(double));
		zpstar = (double *) malloc(5 * sizeof(double));
		xpstar_face = (double **) malloc(5 * sizeof(double *));
		ypstar_face = (double **) malloc(5 * sizeof(double *));
		zpstar_face = (double **) malloc(5 * sizeof(double *));
		theta = (double **) malloc(5 * sizeof(double *));

		for (i = 1; i < 5; i++)
		{
			xc_lateral_side1[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (xc_lateral_side1[i] == NULL)
				printf("xc_lateral_side1_slave memory failed\n");
			yc_lateral_side1[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (yc_lateral_side1[i] == NULL)
				printf("yc_lateral_side1_slave memory failed\n");
			zc_lateral_side1[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (zc_lateral_side1[i] == NULL)
				printf("zc_lateral_side1_slave memory failed\n");
			nx_lateral_side1[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (nx_lateral_side1[i] == NULL)
				printf("nx_lateral_side1_slave memory failed\n");
			ny_lateral_side1[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (ny_lateral_side1[i] == NULL)
				printf("ny_lateral_side1_slave memory failed\n");
			nz_lateral_side1[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (nz_lateral_side1[i] == NULL)
				printf("nz_lateral_side1_slave memory failed\n");
			area_lateral_side1[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (area_lateral_side1[i] == NULL)
				printf("area_lateral_side1_slave memory failed\n");
			xc_lateral_side2[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (xc_lateral_side2[i] == NULL)
				printf("xc_lateral_side2_slave memory failed\n");
			yc_lateral_side2[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (yc_lateral_side2[i] == NULL)
				printf("yc_lateral_side2_slave memory failed\n");
			zc_lateral_side2[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (zc_lateral_side2[i] == NULL)
				printf("zc_lateral_side2_slave memory failed\n");
			nx_lateral_side2[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (nx_lateral_side2[i] == NULL)
				printf("nx_lateral_side2_slave memory failed\n");
			ny_lateral_side2[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (ny_lateral_side2[i] == NULL)
				printf("ny_lateral_side2_slave memory failed\n");
			nz_lateral_side2[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (nz_lateral_side2[i] == NULL)
				printf("nz_lateral_side2_slave memory failed\n");
			area_lateral_side2[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (area_lateral_side2[i] == NULL)
				printf("area_lateral_side2_slave memory failed\n");
			xc_lateral_side3[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (xc_lateral_side3[i] == NULL)
				printf("xc_lateral_side3_slave memory failed\n");
			yc_lateral_side3[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (yc_lateral_side3[i] == NULL)
				printf("yc_lateral_side3_slave memory failed\n");
			zc_lateral_side3[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (zc_lateral_side3[i] == NULL)
				printf("zc_lateral_side3_slave memory failed\n");
			nx_lateral_side3[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (nx_lateral_side3[i] == NULL)
				printf("nx_lateral_side3_slave memory failed\n");
			ny_lateral_side3[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (ny_lateral_side3[i] == NULL)
				printf("ny_lateral_side3_slave memory failed\n");
			nz_lateral_side3[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (nz_lateral_side3[i] == NULL)
				printf("nz_lateral_side3_slave memory failed\n");
			area_lateral_side3[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (area_lateral_side3[i] == NULL)
				printf("area_lateral_side1_slave memory failed\n");
			xc_neighface[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (xc_neighface[i] == NULL)
				printf("xc_neighface_slave memory failed\n");
			yc_neighface[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (yc_neighface[i] == NULL)
				printf("yc_neighface_slave memory failed\n");
			zc_neighface[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (zc_neighface[i] == NULL)
				printf("zc_neighface_slave memory failed\n");
			vol_neighface[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (vol_neighface[i] == NULL)
				printf("vol_neighface_slave memory failed\n");
			xc_face[i] = (double *) malloc((lele_cnt + 1) * sizeof(double));
			if (xc_face[i] == NULL)
				printf("xc_face_slave memory failed\n");
			yc_face[i] = (double *) malloc((lele_cnt + 1) * sizeof(double));
			if (yc_face[i] == NULL)
				printf("yc_face_slave memory failed\n");
			zc_face[i] = (double *) malloc((lele_cnt + 1) * sizeof(double));
			if (zc_face[i] == NULL)
				printf("zc_face_slave memory failed\n");
			vol_face[i] = (double *) malloc((lele_cnt + 1) * sizeof(double));
			if (vol_face[i] == NULL)
				printf("vol_face_slave memory failed\n");
			vol_hexahedra_face[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (vol_hexahedra_face[i] == NULL)
				printf("vol_hexa_face_slave memory failed\n");
			xc_hexahedra_face[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (xc_hexahedra_face[i] == NULL)
				printf("xc_hexa_face_slave memory failed\n");
			yc_hexahedra_face[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (yc_hexahedra_face[i] == NULL)
				printf("yc_hexa_face_slave memory failed\n");
			zc_hexahedra_face[i] = (double *) malloc(
					(lele_cnt + 1) * sizeof(double));
			if (zc_hexahedra_face[i] == NULL)
				printf("zc_hexa_face_slave memory failed\n");
			xpstar_face[i] = (double *) malloc(5 * sizeof(double));
			if (xpstar_face[i] == NULL)
				printf("xpstar_face_slave memory failed\n");
			ypstar_face[i] = (double *) malloc(5 * sizeof(double));
			if (ypstar_face[i] == NULL)
				printf("ypstar_face_slave memory failed\n");
			zpstar_face[i] = (double *) malloc(5 * sizeof(double));
			if (zpstar_face[i] == NULL)
				printf("zpstar_face_slave memory failed\n");
			theta[i] = (double *) malloc((lele_cnt + 1) * sizeof(double));
			if (theta[i] == NULL)
				printf("theta_slave memory failed\n");
		}

		double *vol_polyhedra, *vol;
		vol_polyhedra = (double *) malloc((lele_cnt + 1) * sizeof(double));
		if (vol_polyhedra == NULL)
			printf("vol_polyhedra_slave memory failed\n");
		vol = (double *) malloc((lele_cnt + 1) * sizeof(double));
		if (vol == NULL)
			printf("vol_slave memory failed\n");

		int jnode;

//reprocessing for nodes: local numbering of nodes

		int lnode_cnt;

		int *lnode, *gb_node;
		lnode = (int *) malloc((lele_cnt + 1) * sizeof(int));
		if (lnode == NULL)
			printf("lnode_slave memory failed\n");
		gb_node = (int *) malloc((node + 1) * sizeof(int));
		if (gb_node == NULL)
			printf("gb_node_slave memory failed\n");

		lnode_cnt = 0;

		int r[5];

		i = slave_node_local_numbering(filename, id, i, lele_cnt, r, check, fp,
				&j, lelem, &lnode_cnt, lnode, gb_node, &k, &n1);

		xycoord = (struct node1*) malloc(
				(lnode_cnt + 1) * sizeof(struct node1));
		if (xycoord == NULL)
			printf("xycoord_slave memory failed4\n");
		MPI_Send(&lele_cnt, 1, MPI_INT, master, id, MPI_COMM_WORLD);

		count[1] = lele_cnt;
		count[2] = lnode_cnt;

		MPI_Send(&count[1], 2, MPI_INT, master, id, MPI_COMM_WORLD);

		MPI_Pack_size(count[1], MPI_INT, MPI_COMM_WORLD, &maxsize);
		buffer = malloc(maxsize);

		MPI_Pack_size(count[2], MPI_INT, MPI_COMM_WORLD, &maxsize1);
		buffer1 = malloc(maxsize1);

		position = 0;

		for (i = 1; i <= lele_cnt; i++)
		{
			MPI_Pack(&lelem[i], 1, MPI_INT, buffer, maxsize, &position,
					MPI_COMM_WORLD);
		}

		MPI_Send(buffer, position, MPI_PACKED, master, id, MPI_COMM_WORLD);

		free(buffer);
		buffer = NULL;

		position1 = 0;

		for (i = 1; i <= lnode_cnt; i++)
		{
			MPI_Pack(&lnode[i], 1, MPI_INT, buffer1, maxsize1, &position1,
					MPI_COMM_WORLD);
		}

		MPI_Send(buffer1, position1, MPI_PACKED, master, id, MPI_COMM_WORLD);

		free(buffer1);
		buffer1 = NULL;

		MPI_Pack_size(count[1], trans_eleminfo, MPI_COMM_WORLD, &maxsize1);
		MPI_Pack_size(count[1], trans_primitive, MPI_COMM_WORLD, &maxsize2);
		MPI_Pack_size(count[1], trans_derv, MPI_COMM_WORLD, &maxsize3);
		maxsize = maxsize1 + maxsize2 + maxsize3;
		buffer = malloc(maxsize);

		MPI_Recv(buffer, maxsize, MPI_PACKED, master, id, MPI_COMM_WORLD,
				&status);

		position = 0;

		for (j = 1; j <= lele_cnt; j++)
		{
			MPI_Unpack(buffer, maxsize, &position, &eleminfo[j], 1,
					trans_eleminfo, MPI_COMM_WORLD);
			MPI_Unpack(buffer, maxsize, &position, &primitive[j], 1,
					trans_primitive, MPI_COMM_WORLD);
			MPI_Unpack(buffer, maxsize, &position, &derv[j], 1, trans_derv,
					MPI_COMM_WORLD);
		}

		free(buffer);
		buffer = NULL;

		for (j = 1; j <= lele_cnt; j++)
		{
			d_ele = eleminfo[j].gele_nodes[1];
			d_store = gb_node[d_ele];
			eleminfo[j].gele_nodes[1] = d_store;

			d_ele = eleminfo[j].gele_nodes[2];
			d_store = gb_node[d_ele];
			eleminfo[j].gele_nodes[2] = d_store;

			d_ele = eleminfo[j].gele_nodes[3];
			d_store = gb_node[d_ele];
			eleminfo[j].gele_nodes[3] = d_store;

			d_ele = eleminfo[j].gele_nodes[4];
			d_store = gb_node[d_ele];
			eleminfo[j].gele_nodes[4] = d_store;

			d_ele = eleminfo[j].node1_localface[1];
			d_store = gb_node[d_ele];
			eleminfo[j].node1_localface[1] = d_store;

			d_ele = eleminfo[j].node1_localface[2];
			d_store = gb_node[d_ele];
			eleminfo[j].node1_localface[2] = d_store;

			d_ele = eleminfo[j].node1_localface[3];
			d_store = gb_node[d_ele];
			eleminfo[j].node1_localface[3] = d_store;

			d_ele = eleminfo[j].node1_localface[4];
			d_store = gb_node[d_ele];
			eleminfo[j].node1_localface[4] = d_store;

			d_ele = eleminfo[j].node2_localface[1];
			d_store = gb_node[d_ele];
			eleminfo[j].node2_localface[1] = d_store;

			d_ele = eleminfo[j].node2_localface[2];
			d_store = gb_node[d_ele];
			eleminfo[j].node2_localface[2] = d_store;

			d_ele = eleminfo[j].node2_localface[3];
			d_store = gb_node[d_ele];
			eleminfo[j].node2_localface[3] = d_store;

			d_ele = eleminfo[j].node2_localface[4];
			d_store = gb_node[d_ele];
			eleminfo[j].node2_localface[4] = d_store;

			d_ele = eleminfo[j].node3_localface[1];
			d_store = gb_node[d_ele];
			eleminfo[j].node3_localface[1] = d_store;

			d_ele = eleminfo[j].node3_localface[2];
			d_store = gb_node[d_ele];
			eleminfo[j].node3_localface[2] = d_store;

			d_ele = eleminfo[j].node3_localface[3];
			d_store = gb_node[d_ele];
			eleminfo[j].node3_localface[3] = d_store;

			d_ele = eleminfo[j].node3_localface[4];
			d_store = gb_node[d_ele];
			eleminfo[j].node3_localface[4] = d_store;

		}

		MPI_Pack_size(count[2], trans_xycoord, MPI_COMM_WORLD, &maxsize);
		buffer = malloc(maxsize);

		MPI_Recv(buffer, maxsize, MPI_PACKED, master, id, MPI_COMM_WORLD,
				&status);

		position = 0;

		for (j = 1; j <= lnode_cnt; j++)
		{
			MPI_Unpack(buffer, maxsize, &position, &xycoord[j], 1,
					trans_xycoord, MPI_COMM_WORLD);
		}

		free(buffer);
		buffer = NULL;

		MPI_Barrier(comm_slaves);

//converting global numbers into local numbering in lelem_neighbour

		i = slave_global_2_local_no(i, lele_cnt, d_store, &j, lelem_neighbour,
				gb);

//calculate centroids

		centroid = (struct node3*) malloc(tot_ele * sizeof(struct node3));

		if (centroid == NULL)
			printf("centroid_slave memory failed4\n");

		MPI_Datatype trans_centroid;
		MPI_Datatype type3[3] =
		{ MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };
		int blocklen3[3] =
		{ 1, 1, 1 };
		MPI_Aint disp3[3] =
		{ 0, sizeof(double), 2 * sizeof(double) };
		MPI_Type_struct(3, blocklen3, disp3, type3, &trans_centroid);
		MPI_Type_commit(&trans_centroid);

		sol = (struct node4*) malloc(tot_ele * sizeof(struct node4));
		if (sol == NULL)
			printf("sol_slave memory failed4\n");

		MPI_Datatype trans_sol;
		MPI_Datatype type4[3] =
		{ MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE };
		int blocklen4[3] =
		{ 1, 1, 1 };
		MPI_Aint disp4[3] =
		{ 0, sizeof(double), 2 * sizeof(double) };
		MPI_Type_struct(3, blocklen4, disp4, type4, &trans_sol);
		MPI_Type_commit(&trans_sol);

//centroid of tetrahedrons

		iele = centroid_cal(iele, lele_cnt, &d1, eleminfo, &d2, &d3, &d4,
				centroid, xycoord);

//periodic_elem[i][j]:maps with the periodic boundary  cell
		int **periodic_elem;

		double **nx_boun_face, **ny_boun_face, **nz_boun_face, **trans_coord_dx,
				**trans_coord_dy, **trans_coord_dz, **trans_coord_mx,
				**trans_coord_my, **trans_coord_mz, **aa, **bb, **cc, **dd;

		nx_boun_face = (double **) malloc((b_cnt + 1) * sizeof(double *));
		ny_boun_face = (double **) malloc((b_cnt + 1) * sizeof(double *));
		nz_boun_face = (double **) malloc((b_cnt + 1) * sizeof(double *));
		trans_coord_dx = (double **) malloc((b_cnt + 1) * sizeof(double *));
		trans_coord_dy = (double **) malloc((b_cnt + 1) * sizeof(double *));
		trans_coord_dz = (double **) malloc((b_cnt + 1) * sizeof(double *));
		trans_coord_mx = (double **) malloc((b_cnt + 1) * sizeof(double *));
		trans_coord_my = (double **) malloc((b_cnt + 1) * sizeof(double *));
		trans_coord_mz = (double **) malloc((b_cnt + 1) * sizeof(double *));
		aa = (double **) malloc((b_cnt + 1) * sizeof(double *));
		bb = (double **) malloc((b_cnt + 1) * sizeof(double *));
		cc = (double **) malloc((b_cnt + 1) * sizeof(double *));
		dd = (double **) malloc((b_cnt + 1) * sizeof(double *));

		periodic_elem = (int **) malloc((bc_set + 1) * sizeof(int *));

//for curved boundary sets

		for (i = 1; i <= b_cnt; i++)
		{
			nx_boun_face[i] = (double *) malloc(
					(cnt_boun_type[i] + 1) * sizeof(double));
			if (nx_boun_face[i] == NULL)
				printf("nx_boun_face[i]_slave memory failed4\n");
			ny_boun_face[i] = (double *) malloc(
					(cnt_boun_type[i] + 1) * sizeof(double));
			if (ny_boun_face[i] == NULL)
				printf("ny_boun_face[i]_slave memory failed4\n");
			nz_boun_face[i] = (double *) malloc(
					(cnt_boun_type[i] + 1) * sizeof(double));
			if (nz_boun_face[i] == NULL)
				printf("nz_boun_face[i]_slave memory failed4\n");
			trans_coord_dx[i] = (double *) malloc(
					(cnt_boun_type[i] + 1) * sizeof(double));
			if (trans_coord_dx[i] == NULL)
				printf("trans_coord_dx[i]_slave memory failed4\n");
			trans_coord_dy[i] = (double *) malloc(
					(cnt_boun_type[i] + 1) * sizeof(double));
			if (trans_coord_dy[i] == NULL)
				printf("trans_coord_dy[i]_slave memory failed4\n");
			trans_coord_dz[i] = (double *) malloc(
					(cnt_boun_type[i] + 1) * sizeof(double));
			if (trans_coord_dz[i] == NULL)
				printf("trans_coord_dz[i]_slave memory failed4\n");
			trans_coord_mx[i] = (double *) malloc(
					(cnt_boun_type[i] + 1) * sizeof(double));
			if (trans_coord_mx[i] == NULL)
				printf("trans_coord_mx[i]_slave memory failed4\n");
			trans_coord_my[i] = (double *) malloc(
					(cnt_boun_type[i] + 1) * sizeof(double));
			if (trans_coord_my[i] == NULL)
				printf("trans_coord_my[i]_slave memory failed4\n");
			trans_coord_mz[i] = (double *) malloc(
					(cnt_boun_type[i] + 1) * sizeof(double));
			if (trans_coord_mz[i] == NULL)
				printf("trans_coord_mz[i]_slave memory failed4\n");
			aa[i] = (double *) malloc((cnt_boun_type[i] + 1) * sizeof(double));
			if (aa[i] == NULL)
				printf("aa[i]_slave memory failed4\n");
			bb[i] = (double *) malloc((cnt_boun_type[i] + 1) * sizeof(double));
			if (bb[i] == NULL)
				printf("bb[i]_slave memory failed4\n");
			cc[i] = (double *) malloc((cnt_boun_type[i] + 1) * sizeof(double));
			if (cc[i] == NULL)
				printf("cc[i]_slave memory failed4\n");
			dd[i] = (double *) malloc((cnt_boun_type[i] + 1) * sizeof(double));
			if (dd[i] == NULL)
				printf("dd[i]_slave memory failed4\n");

		}

		double **ghost_elem, *xp, *yp, *zp;
		ghost_elem = (double **) malloc((b_cnt + 1) * sizeof(double *));
		xp = malloc(sizeof(double));
		yp = malloc(sizeof(double));
		zp = malloc(sizeof(double));

		nele = subdomain_cnt;

		int boun_cnt_start = nele + 1;

//boundary element calculation

		i = boun_ele_calc(i, b_cnt, ghost_elem, cnt_boun_type, boun_type, &j,
				&iele, boun_elem, &iface, boun_face, &nele, lelem_neighbour,
				&d1, eleminfo, &d2, &d3, xycoord, centroid, nx_boun_face,
				ny_boun_face, nz_boun_face, deno, z1, z2, z3, sum, aa, bb, cc,
				dd, trans_coord_dx, trans_coord_dy, trans_coord_dz,
				trans_coord_mx, trans_coord_my, trans_coord_mz, t, xp, yp, zp,
				periodic_elem);

		int boun_cnt_end;
		boun_cnt_end = nele;

		MPI_Barrier(comm_slaves);

		int periodic_proc, pp;
		periodic_proc = 0;
		pp = 0;
		int *g_period, *l_period, *l_periodic_mapping;
		g_period = (int *) malloc((l_periodic_cnt + 1) * sizeof(int));
		if (g_period == NULL)
			printf("g_period_slave memory failed4\n");
		l_period = (int *) malloc((l_periodic_cnt + 1) * sizeof(int));
		if (l_period == NULL)
			printf("l_period_slave memory failed4\n");
		l_periodic_mapping = (int *) malloc((lele_cnt + 1) * sizeof(int));
		if (l_periodic_mapping == NULL)
			printf("l_periodic_mappinng_slave memory failed4\n");

		if (pg == 1)
		{

			i = slave_periodic_boun_iden(filename, id, i, l_periodic_cnt,
					subdomain_cnt, periodic_proc, d_store, check, fp, &n1,
					elem2proc, g_period, l_period, gb, fp1, &d1,
					l_periodic_mapping, &j, lelem_neighbour, &nele, lelem, &pp,
					id_periodic_proc, pp2id, pp_cnt);
		}

//p_elem_recv: stores the elements which r to be received while exchanging periodic elements
//The algorithm has been written in a different way as compared to subdomain exchange since some
//nodes have to send only and some nodes have to receive only

		int **p_elem_recv;

		p_elem_recv = (int **) malloc((pp + 1) * sizeof(int *));

		for (i = 1; i <= pp; i++)
		{
			p_elem_recv[i] = (int *) malloc((pp_cnt[i] + 1) * sizeof(int));
		}

		for (i = 1; i <= pp; i++)
		{
			p_recvcnt[i] = 0;
		}

		int periodic_id;
		periodic_id = 0;
		if (pg == 1)
		{

			sprintf(filename, "periodic_boun_%d.txt", id);

			if (access(filename, F_OK) != -1)
			{
				periodic_id = 1;
				fp = fopen(filename, "r");
				while (!feof(fp))
				{
					fscanf(fp, "%d %d %d\n", &n1, &n2, &n3);
					for (i = 1; i <= pp; i++)
					{
						if (id_periodic_proc[i] == n3)
						{
							if (p_recvcnt[i] == 0)
							{
								p_recvcnt[i] = 1;
								p_elem_recv[i][p_recvcnt[i]] = n2;
							}
							else
							{
								check = 0;
								for (j = 1; j <= p_recvcnt[i]; j++)
								{
//Most likely the elements wont be repeated as in identification of subdomain neighbouring elements.
									if (n2 == p_elem_recv[i][j])
									{
										check = 1;
										break;
									}
								}
								if (check == 0)
								{
									p_recvcnt[i] = p_recvcnt[i] + 1;
									p_elem_recv[i][p_recvcnt[i]] = n2;
								}
							}
						}
					}
				}

				fclose(fp);
			}
		}

//converting global numbering of p_elem_recv to local numbering
		if (pg == 1)
		{
			if (periodic_id == 1)
			{
				for (i = 1; i <= pp; i++)
				{
					for (j = 1; j <= p_recvcnt[i]; j++)
					{
						d_ele = p_elem_recv[i][j];
						p_elem_recv[i][j] = gb[d_ele];
					}
				}
			}
		}

		int *periodic_send, *pp_send, *pp_recv, *pp_sendrecv, m, pp1,
				*int_p_send;
		periodic_send = (int *) malloc((prs + 1) * sizeof(int));
		pp_send = (int *) malloc((prs + 1) * sizeof(int));
		pp_recv = (int *) malloc((prs + 1) * sizeof(int));
		pp_sendrecv = (int *) malloc((prs + 1) * sizeof(int));
		int_p_send = (int *) malloc((prs + 1) * sizeof(int));

		for (i = 1; i <= prs; i++)
		{
			pp_send[i] = 0;
			pp_recv[i] = 0;
			pp_sendrecv[i] = 0;
		}

		if (pg == 1)
		{
			if (periodic_id == 1)
			{
				pp1 = pp;

				for (i = 1; i <= pp; i++)
				{
					pp_recv[i] = 1;
				}

			}
		}

		m = 0;
		if (pg == 1)
		{
			for (i = 1; i < prs; i++)
			{
				if (i != id)
				{
					MPI_Recv(&d2, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
					if (d2 == 1)
					{
						MPI_Recv(&d1, 1, MPI_INT, i, i, MPI_COMM_WORLD,
								&status);
						MPI_Recv(&int_p_send[1], d1, MPI_INT, i, i,
								MPI_COMM_WORLD, &status);

						for (k = 1; k <= d1; k++)
						{
							if (int_p_send[k] == (id - 1))
							{
								m = m + 1;
								periodic_send[m] = i - 1;
								break;
							}
						}
					}
				}
				else
				{
					for (j = 1; j < id; j++)
					{
						MPI_Send(&periodic_id, 1, MPI_INT, j, i,
								MPI_COMM_WORLD);
					}
					for (j = id + 1; j < prs; j++)
					{
						MPI_Send(&periodic_id, 1, MPI_INT, j, i,
								MPI_COMM_WORLD);
					}

					if (periodic_id == 1)
					{
						for (j = 1; j < id; j++)
						{
							MPI_Send(&pp, 1, MPI_INT, j, i, MPI_COMM_WORLD);
							MPI_Send(&id_periodic_proc[1], pp, MPI_INT, j, i,
									MPI_COMM_WORLD);
						}
						for (j = id + 1; j < prs; j++)
						{
							MPI_Send(&pp, 1, MPI_INT, j, i, MPI_COMM_WORLD);
							MPI_Send(&id_periodic_proc[1], pp, MPI_INT, j, i,
									MPI_COMM_WORLD);
						}
					}
				}
			}
		}

		if (pg == 1)
		{
			if (periodic_id == 1)
			{
				for (k = 1; k <= m; k++)
				{
					check = 0;
					for (i = 1; i <= pp1; i++)
					{
						if (id_periodic_proc[i] == periodic_send[k])
						{
							pp_recv[i] = 0;
							pp_sendrecv[i] = 1;
							check = 1;
							break;
						}
					}
					if (check == 0)
					{
						pp = pp + 1;
						id_periodic_proc[pp] = periodic_send[k];
						pp2id[periodic_send[k]] = pp;
						pp_send[pp] = 1;
					}
				}
			}
		}

		int *p_group;
		p_group = (int *) malloc(sizeof(int) * (pp + 1));

		if (pg == 1)
		{
			i = periodic_proc_rank_arrangement(periodic_id, i, pp, temp,
					p_group, id_periodic_proc, &j);
		}
//preprocessing for exchange of data between nodes
		int recvcnt[kk + 1];

		for (i = 1; i <= kk; i++)
		{
			d_store = group[i];
			d_ele = k2id[d_store];
			MPI_Sendrecv(&send_cnt[d_ele], 1, MPI_INT, d_store + 1, d_store + 1,
					&recvcnt[d_ele], 1, MPI_INT, d_store + 1, id,
					MPI_COMM_WORLD, &status);

		}

		MPI_Barrier(comm_slaves);
//memory allocation for lele_recv: variable used for identifying local_element number that is being received from another node.
//lelem_send:sending elements ; lelem_recv:receiving elements

		lelem_recv = (int **) malloc(sizeof(int *) * (kk + 1));

		for (i = 1; i <= kk; i++)
		{
			lelem_recv[i] = (int *) malloc(sizeof(int) * (recvcnt[i] + 1));
		}

		for (i = 1; i <= kk; i++)
		{
			d_store = group[i];
			d_ele = k2id[d_store];

			MPI_Pack_size(send_cnt[d_ele], MPI_INT, MPI_COMM_WORLD, &maxsize1);
			buffer1 = malloc(maxsize1);

			position1 = 0;
			for (j = 1; j <= send_cnt[d_ele]; j++)
			{
				d_store1 = lelem_send[d_ele][j];
				MPI_Pack(&lelem[d_store1], 1, MPI_INT, buffer1, maxsize1,
						&position1, MPI_COMM_WORLD);
			}

			MPI_Pack_size(recvcnt[d_ele], MPI_INT, MPI_COMM_WORLD, &maxsize2);
			buffer2 = malloc(maxsize2);

			MPI_Sendrecv(buffer1, position1, MPI_PACKED, d_store + 1,
					d_store + 1, buffer2, maxsize2, MPI_PACKED, d_store + 1, id,
					MPI_COMM_WORLD, &status);

			position1 = 0;

			for (j = 1; j <= recvcnt[d_ele]; j++)
			{
				MPI_Unpack(buffer2, maxsize2, &position1, &k, 1, MPI_INT,
						MPI_COMM_WORLD);
				lelem_recv[d_ele][j] = gb[k];
			}

			free(buffer1);
			buffer1 = NULL;
			free(buffer2);
			buffer2 = NULL;

		}

//exchange of centroids

		for (i = 1; i <= kk; i++)
		{
			d_store = group[i];
			d_ele = k2id[d_store];
			MPI_Pack_size(send_cnt[d_ele], trans_centroid, MPI_COMM_WORLD,
					&maxsize1);
			buffer1 = malloc(maxsize1);
			position1 = 0;

			for (j = 1; j <= send_cnt[d_ele]; j++)
			{
				d_store1 = lelem_send[d_ele][j];
				MPI_Pack(&centroid[d_store1], 1, trans_centroid, buffer1,
						maxsize1, &position1, MPI_COMM_WORLD);
			}

			MPI_Pack_size(recvcnt[d_ele], trans_centroid, MPI_COMM_WORLD,
					&maxsize2);
			buffer2 = malloc(maxsize2);
			MPI_Sendrecv(buffer1, position1, MPI_PACKED, d_store + 1,
					d_store + 1, buffer2, maxsize2, MPI_PACKED, d_store + 1, id,
					MPI_COMM_WORLD, &status);

			position1 = 0;

			for (j = 1; j <= recvcnt[d_ele]; j++)
			{
				d_store2 = lelem_recv[d_ele][j];
				MPI_Unpack(buffer2, maxsize2, &position1, &centroid[d_store2],
						1, trans_centroid, MPI_COMM_WORLD);
			}

			free(buffer1);
			buffer1 = NULL;
			free(buffer2);
			buffer2 = NULL;

		} //i bracket

		int *p_send_cnt;
		p_send_cnt = (int *) malloc(sizeof(int) * (pp + 1));

		MPI_Barrier(comm_slaves);
//memory allocation for lele_recv: variable used for identifying local_element number that is being received from another node.
//lelem_send:sending elements ; lelem_recv:receiving elements
		int **p_elem_send;
		p_elem_send = (int **) malloc(sizeof(int *) * (pp + 1));

		if (pg == 1)
		{
			if (periodic_id == 1)
			{

				for (i = 1; i <= pp; i++)
				{
					d_store = p_group[i];
					d_ele = pp2id[d_store];
					if (pp_send[d_ele] == 1)
					{
						MPI_Recv(&p_send_cnt[d_ele], 1, MPI_INT, d_store + 1,
								d_store + 1, MPI_COMM_WORLD, &status);
					}
					else if (pp_recv[d_ele] == 1)
					{
						MPI_Send(&p_recvcnt[d_ele], 1, MPI_INT, d_store + 1, id,
								MPI_COMM_WORLD);
					}
					else if (pp_sendrecv[d_ele] == 1)
					{
						MPI_Sendrecv(&p_recvcnt[d_ele], 1, MPI_INT, d_store + 1,
								d_store + 1, &p_send_cnt[d_ele], 1, MPI_INT,
								d_store + 1, id, MPI_COMM_WORLD, &status);
					}
				}
			}
		}

		for (i = 1; i <= pp; i++)
		{
			if ((pp_send[i] == 1) || (pp_sendrecv[i] == 1))
			{
				p_elem_send[i] = (int *) malloc(
						sizeof(int) * (p_send_cnt[i] + 1));
			}
		}

		if (pg == 1)
		{
			if (periodic_id == 1)
			{
				for (i = 1; i <= pp; i++)
				{
					d_store = p_group[i];
					d_ele = pp2id[d_store];
					if (pp_recv[d_ele] == 1)
					{
						MPI_Pack_size(p_recvcnt[d_ele], MPI_INT, MPI_COMM_WORLD,
								&maxsize1);
						buffer1 = malloc(maxsize1);

						position1 = 0;
						for (j = 1; j <= p_recvcnt[d_ele]; j++)
						{
							d_store1 = p_elem_recv[d_ele][j];
							MPI_Pack(&lelem[d_store1], 1, MPI_INT, buffer1,
									maxsize1, &position1, MPI_COMM_WORLD);
						}

						MPI_Send(buffer1, position1, MPI_PACKED, d_store + 1,
								d_store + 1, MPI_COMM_WORLD);

						free(buffer1);
						buffer1 = NULL;
					}
//receiving end
					else if (pp_send[d_ele] == 1)
					{

						MPI_Pack_size(p_send_cnt[d_ele], MPI_INT,
								MPI_COMM_WORLD, &maxsize2);
						buffer2 = malloc(maxsize2);

						MPI_Recv(buffer2, maxsize2, MPI_PACKED, d_store + 1, id,
								MPI_COMM_WORLD, &status);

						position1 = 0;

						for (j = 1; j <= p_send_cnt[d_ele]; j++)
						{
							MPI_Unpack(buffer2, maxsize2, &position1, &k, 1,
									MPI_INT, MPI_COMM_WORLD);
							p_elem_send[d_ele][j] = gb[k];
						}

						free(buffer2);
						buffer2 = NULL;
					}

					else if (pp_sendrecv[d_ele] == 1)
					{

						MPI_Pack_size(p_recvcnt[d_ele], MPI_INT, MPI_COMM_WORLD,
								&maxsize1);
						buffer1 = malloc(maxsize1);

						position1 = 0;
						for (j = 1; j <= p_recvcnt[d_ele]; j++)
						{
							d_store1 = p_elem_recv[d_ele][j];
							MPI_Pack(&lelem[d_store1], 1, MPI_INT, buffer1,
									maxsize1, &position1, MPI_COMM_WORLD);
						}

						MPI_Pack_size(p_send_cnt[d_ele], MPI_INT,
								MPI_COMM_WORLD, &maxsize2);
						buffer2 = malloc(maxsize2);

						MPI_Sendrecv(buffer1, position1, MPI_PACKED,
								d_store + 1, d_store + 1, buffer2, maxsize2,
								MPI_PACKED, d_store + 1, id, MPI_COMM_WORLD,
								&status);

						position1 = 0;

						for (j = 1; j <= p_send_cnt[d_ele]; j++)
						{
							MPI_Unpack(buffer2, maxsize2, &position1, &k, 1,
									MPI_INT, MPI_COMM_WORLD);
							p_elem_send[d_ele][j] = gb[k];
						}

						free(buffer1);
						buffer1 = NULL;
						free(buffer2);
						buffer2 = NULL;
					}
				}

			}
		}

//calculation of solution point

		iele = sol_pt_cal(iele, lele_cnt, iface, d1, d2, d3, d4, eleminfo,
				lelem_neighbour, xc_neighface, xycoord, centroid, yc_neighface,
				zc_neighface, vol_neighface, xc_face, yc_face, zc_face,
				vol_face, vol_hexahedra_face, xc_hexahedra_face,
				yc_hexahedra_face, zc_hexahedra_face, xc_lateral_side1,
				yc_lateral_side1, zc_lateral_side1, xc_lateral_side2,
				yc_lateral_side2, zc_lateral_side2, xc_lateral_side3,
				yc_lateral_side3, zc_lateral_side3, nx_lateral_side1,
				ny_lateral_side1, nz_lateral_side1, area_lateral_side1, deno,
				z1, z2, z3, sum, nx_lateral_side2, ny_lateral_side2,
				nz_lateral_side2, area_lateral_side2, nx_lateral_side3,
				ny_lateral_side3, nz_lateral_side3, area_lateral_side3,
				vol_polyhedra, sol);

		MPI_Barrier(comm_slaves);
//exchange of solution_point

		for (i = 1; i <= kk; i++)
		{

			d_store = group[i];
			d_ele = k2id[d_store];
			MPI_Pack_size(send_cnt[d_ele], trans_sol, MPI_COMM_WORLD,
					&maxsize1);
			buffer1 = malloc(maxsize1);
			position1 = 0;

			for (j = 1; j <= send_cnt[d_ele]; j++)
			{
				d_store1 = lelem_send[d_ele][j];
				MPI_Pack(&sol[d_store1], 1, trans_sol, buffer1, maxsize1,
						&position1, MPI_COMM_WORLD);
			}

			MPI_Pack_size(recvcnt[d_ele], trans_sol, MPI_COMM_WORLD, &maxsize2);
			buffer2 = malloc(maxsize2);
			MPI_Sendrecv(buffer1, position1, MPI_PACKED, d_store + 1,
					d_store + 1, buffer2, maxsize2, MPI_PACKED, d_store + 1, id,
					MPI_COMM_WORLD, &status);

			position1 = 0;

			for (j = 1; j <= recvcnt[d_ele]; j++)
			{
				d_store2 = lelem_recv[d_ele][j];
				MPI_Unpack(buffer2, maxsize2, &position1, &sol[d_store2], 1,
						trans_sol, MPI_COMM_WORLD);
			}

			free(buffer1);
			buffer1 = NULL;
			free(buffer2);
			buffer2 = NULL;

		} //i bracket


//assigning solution points for boundaries

		i = boun_sol_pt(i, b_cnt, iface, d2, d3, boun_type, &j, cnt_boun_type,
				&iele, boun_elem, &nele, ghost_elem, t, dd, aa, sol, bb, cc,
				nx_boun_face, ny_boun_face, nz_boun_face, boun_face, &d1,
				eleminfo, xp, yp, zp, xycoord, trans_coord_dx, trans_coord_dy,
				trans_coord_dz, trans_coord_mx, trans_coord_my, trans_coord_mz,
				deno, z1, z2, z3, periodic_elem, l_periodic_mapping);

//calculation of secondary variables

		for (iele = 1; iele <= lele_cnt; iele++)
		{
			secondary_variables(&primitive[iele].u1, &derv[iele].ux1,
					&derv[iele].uy1, &derv[iele].uz1, &primitive[iele].u2,
					&derv[iele].ux2, &derv[iele].uy2, &derv[iele].uz2,
					&primitive[iele].u3, &derv[iele].ux3, &derv[iele].uy3,
					&derv[iele].uz3, &primitive[iele].u4, &derv[iele].ux4,
					&derv[iele].uy4, &derv[iele].uz4, &primitive[iele].u5,
					&derv[iele].ux5, &derv[iele].uy5, &derv[iele].uz5,
					&secondary[iele].ut1, &secondary[iele].ut2,
					&secondary[iele].ut3, &secondary[iele].ut4,
					&secondary[iele].ut5, &secondary[iele].f1,
					&secondary[iele].fx1, &secondary[iele].fy1,
					&secondary[iele].fz1, &secondary[iele].g1,
					&secondary[iele].gx1, &secondary[iele].gy1,
					&secondary[iele].gz1, &secondary[iele].e1,
					&secondary[iele].ex1, &secondary[iele].ey1,
					&secondary[iele].ez1, &secondary[iele].ut1,
					&secondary[iele].ft1, &secondary[iele].gt1,
					&secondary[iele].et1, &secondary[iele].f2,
					&secondary[iele].fx2, &secondary[iele].fy2,
					&secondary[iele].fz2, &secondary[iele].g2,
					&secondary[iele].gx2, &secondary[iele].gy2,
					&secondary[iele].gz2, &secondary[iele].e2,
					&secondary[iele].ex2, &secondary[iele].ey2,
					&secondary[iele].ez2, &secondary[iele].ut2,
					&secondary[iele].ft2, &secondary[iele].gt2,
					&secondary[iele].et2, &secondary[iele].f3,
					&secondary[iele].fx3, &secondary[iele].fy3,
					&secondary[iele].fz3, &secondary[iele].g3,
					&secondary[iele].gx3, &secondary[iele].gy3,
					&secondary[iele].gz3, &secondary[iele].e3,
					&secondary[iele].ex3, &secondary[iele].ey3,
					&secondary[iele].ez3, &secondary[iele].ut3,
					&secondary[iele].ft3, &secondary[iele].gt3,
					&secondary[iele].et3, &secondary[iele].f4,
					&secondary[iele].fx4, &secondary[iele].fy4,
					&secondary[iele].fz4, &secondary[iele].g4,
					&secondary[iele].gx4, &secondary[iele].gy4,
					&secondary[iele].gz4, &secondary[iele].e4,
					&secondary[iele].ex4, &secondary[iele].ey4,
					&secondary[iele].ez4, &secondary[iele].ut4,
					&secondary[iele].ft4, &secondary[iele].gt4,
					&secondary[iele].et4, &secondary[iele].f5,
					&secondary[iele].fx5, &secondary[iele].fy5,
					&secondary[iele].fz5, &secondary[iele].g5,
					&secondary[iele].gx5, &secondary[iele].gy5,
					&secondary[iele].gz5, &secondary[iele].e5,
					&secondary[iele].ex5, &secondary[iele].ey5,
					&secondary[iele].ez5, &secondary[iele].ut5,
					&secondary[iele].ft5, &secondary[iele].gt5,
					&secondary[iele].et5);
		}

//exchange of data for primitive,derv and secondary variable betweeen nodes

		for (i = 1; i <= kk; i++)
		{
			d_store = group[i];
			d_ele = k2id[d_store];
			MPI_Pack_size(send_cnt[d_ele], trans_primitive, MPI_COMM_WORLD,
					&maxsize3);
			MPI_Pack_size(send_cnt[d_ele], trans_derv, MPI_COMM_WORLD,
					&maxsize4);
			MPI_Pack_size(send_cnt[d_ele], trans_secondary, MPI_COMM_WORLD,
					&maxsize5);
			maxsize1 = maxsize3 + maxsize4 + maxsize5;
			buffer1 = malloc(maxsize1);

			position1 = 0;

			for (j = 1; j <= send_cnt[d_ele]; j++)
			{
				d_store1 = lelem_send[d_ele][j];
				MPI_Pack(&primitive[d_store1], 1, trans_primitive, buffer1,
						maxsize1, &position1, MPI_COMM_WORLD);
				MPI_Pack(&derv[d_store1], 1, trans_derv, buffer1, maxsize1,
						&position1, MPI_COMM_WORLD);
				MPI_Pack(&secondary[d_store1], 1, trans_secondary, buffer1,
						maxsize1, &position1, MPI_COMM_WORLD);
			}

			MPI_Pack_size(recvcnt[d_ele], trans_primitive, MPI_COMM_WORLD,
					&maxsize6);
			MPI_Pack_size(recvcnt[d_ele], trans_derv, MPI_COMM_WORLD,
					&maxsize7);
			MPI_Pack_size(recvcnt[d_ele], trans_secondary, MPI_COMM_WORLD,
					&maxsize8);
			maxsize2 = maxsize6 + maxsize7 + maxsize8;
			buffer2 = malloc(maxsize2);
			MPI_Sendrecv(buffer1, position1, MPI_PACKED, d_store + 1,
					d_store + 1, buffer2, maxsize2, MPI_PACKED, d_store + 1, id,
					MPI_COMM_WORLD, &status);

			position1 = 0;

			for (j = 1; j <= recvcnt[d_ele]; j++)
			{
				d_store2 = lelem_recv[d_ele][j];
				MPI_Unpack(buffer2, maxsize2, &position1, &primitive[d_store2],
						1, trans_primitive, MPI_COMM_WORLD);
				MPI_Unpack(buffer2, maxsize2, &position1, &derv[d_store2], 1,
						trans_derv, MPI_COMM_WORLD);
				MPI_Unpack(buffer2, maxsize2, &position1, &secondary[d_store2],
						1, trans_secondary, MPI_COMM_WORLD);
			}

			free(buffer1);
			buffer1 = NULL;
			free(buffer2);
			buffer2 = NULL;

		} //i bracket

		if (pg == 1)
		{
//exchange of data for periodic elements
			if (periodic_id == 1)
			{
				for (i = 1; i <= pp; i++)
				{
					d_store = p_group[i];
					d_ele = pp2id[d_store];

					if (pp_send[d_ele] == 1)
					{

						MPI_Pack_size(p_send_cnt[d_ele], trans_primitive,
								MPI_COMM_WORLD, &maxsize3);
						MPI_Pack_size(p_send_cnt[d_ele], trans_derv,
								MPI_COMM_WORLD, &maxsize4);
						maxsize1 = maxsize3 + maxsize4;
						buffer1 = malloc(maxsize1);

						position1 = 0;

						for (j = 1; j <= p_send_cnt[d_ele]; j++)
						{
							d_store1 = p_elem_send[d_ele][j];
							MPI_Pack(&primitive[d_store1], 1, trans_primitive,
									buffer1, maxsize1, &position1,
									MPI_COMM_WORLD);
							MPI_Pack(&derv[d_store1], 1, trans_derv, buffer1,
									maxsize1, &position1, MPI_COMM_WORLD);
						}
						MPI_Send(buffer1, position1, MPI_PACKED, d_store + 1,
								d_store + 1, MPI_COMM_WORLD);

						free(buffer1);
						buffer1 = NULL;
					}

					else if (pp_recv[d_ele] == 1)
					{

						MPI_Pack_size(p_recvcnt[d_ele], trans_primitive,
								MPI_COMM_WORLD, &maxsize6);
						MPI_Pack_size(p_recvcnt[d_ele], trans_derv,
								MPI_COMM_WORLD, &maxsize7);
//MPI_Pack_size(p_recvcnt[d_ele],trans_secondary,MPI_COMM_WORLD,&maxsize8);
						maxsize2 = maxsize6 + maxsize7;
						buffer2 = malloc(maxsize2);

						MPI_Recv(buffer2, maxsize2, MPI_PACKED, d_store + 1, id,
								MPI_COMM_WORLD, &status);

						position1 = 0;

						for (j = 1; j <= p_recvcnt[d_ele]; j++)
						{
							d_store2 = p_elem_recv[d_ele][j];
							MPI_Unpack(buffer2, maxsize2, &position1,
									&primitive[d_store2], 1, trans_primitive,
									MPI_COMM_WORLD);
							MPI_Unpack(buffer2, maxsize2, &position1,
									&derv[d_store2], 1, trans_derv,
									MPI_COMM_WORLD);

						}

						free(buffer2);
						buffer2 = NULL;
					}

					else if (pp_sendrecv[d_ele] == 1)
					{
						MPI_Pack_size(p_send_cnt[d_ele], trans_primitive,
								MPI_COMM_WORLD, &maxsize3);
						MPI_Pack_size(p_send_cnt[d_ele], trans_derv,
								MPI_COMM_WORLD, &maxsize4);

						maxsize1 = maxsize3 + maxsize4;
						buffer1 = malloc(maxsize1);

						position1 = 0;

						for (j = 1; j <= p_send_cnt[d_ele]; j++)
						{
							d_store1 = p_elem_send[d_ele][j];
							MPI_Pack(&primitive[d_store1], 1, trans_primitive,
									buffer1, maxsize1, &position1,
									MPI_COMM_WORLD);
							MPI_Pack(&derv[d_store1], 1, trans_derv, buffer1,
									maxsize1, &position1, MPI_COMM_WORLD);

						}

						MPI_Pack_size(p_recvcnt[d_ele], trans_primitive,
								MPI_COMM_WORLD, &maxsize6);
						MPI_Pack_size(p_recvcnt[d_ele], trans_derv,
								MPI_COMM_WORLD, &maxsize7);

						maxsize2 = maxsize6 + maxsize7;
						buffer2 = malloc(maxsize2);
						MPI_Sendrecv(buffer1, position1, MPI_PACKED,
								d_store + 1, d_store + 1, buffer2, maxsize2,
								MPI_PACKED, d_store + 1, id, MPI_COMM_WORLD,
								&status);

						position1 = 0;

						for (j = 1; j <= p_recvcnt[d_ele]; j++)
						{
							d_store2 = p_elem_recv[d_ele][j];
							MPI_Unpack(buffer2, maxsize2, &position1,
									&primitive[d_store2], 1, trans_primitive,
									MPI_COMM_WORLD);
							MPI_Unpack(buffer2, maxsize2, &position1,
									&derv[d_store2], 1, trans_derv,
									MPI_COMM_WORLD);

						}

						free(buffer1);
						buffer1 = NULL;
						free(buffer2);
						buffer2 = NULL;
					}

				} //i bracket
			}
		}

//allocating values for u[]_p:previous values

		iele = previous_time_step_alloc(iele, subdomain_cnt, u1_p, primitive,
				ux1_p, derv, uy1_p, uz1_p, ut1_p, secondary, f1_p, fx1_p, fy1_p,
				fz1_p, ft1_p, g1_p, gx1_p, gy1_p, gz1_p, gt1_p, e1_p, ex1_p,
				ey1_p, ez1_p, et1_p, u2_p, ux2_p, uy2_p, uz2_p, ut2_p, f2_p,
				fx2_p, fy2_p, fz2_p, ft2_p, g2_p, gx2_p, gy2_p, gz2_p, gt2_p,
				e2_p, ex2_p, ey2_p, ez2_p, et2_p, u3_p, ux3_p, uy3_p, uz3_p,
				ut3_p, f3_p, fx3_p, fy3_p, fz3_p, ft3_p, g3_p, gx3_p, gy3_p,
				gz3_p, gt3_p, e3_p, ex3_p, ey3_p, ez3_p, et3_p, u4_p, ux4_p,
				uy4_p, uz4_p, ut4_p, f4_p, fx4_p, fy4_p, fz4_p, ft4_p, g4_p,
				gx4_p, gy4_p, gz4_p, gt4_p, e4_p, ex4_p, ey4_p, ez4_p, et4_p,
				u5_p, ux5_p, uy5_p, uz5_p, ut5_p, f5_p, fx5_p, fy5_p, fz5_p,
				ft5_p, g5_p, gx5_p, gy5_p, gz5_p, gt5_p, e5_p, ex5_p, ey5_p,
				ez5_p, et5_p);

		double nu1;

//initial boundary condition

		double *t_ux1, *t_uy1, *t_uz1, *t_ux2, *t_uy2, *t_uz2, *t_ux3, *t_uy3,
				*t_uz3, *t_ux4, *t_uy4, *t_uz4, *t_ux5, *t_uy5, *t_uz5,
				*t_ux1_iele, *t_uy1_iele, *t_uz1_iele, *t_ux2_iele, *t_uy2_iele,
				*t_uz2_iele, *t_ux3_iele, *t_uy3_iele, *t_uz3_iele, *t_ux4_iele,
				*t_uy4_iele, *t_uz4_iele, *t_ux5_iele, *t_uy5_iele, *t_uz5_iele;

		t_ux1 = malloc(sizeof(double));
		t_uy1 = malloc(sizeof(double));
		t_uz1 = malloc(sizeof(double));
		t_ux2 = malloc(sizeof(double));
		t_uy2 = malloc(sizeof(double));
		t_uz2 = malloc(sizeof(double));
		t_ux3 = malloc(sizeof(double));
		t_uy3 = malloc(sizeof(double));
		t_uz3 = malloc(sizeof(double));
		t_ux4 = malloc(sizeof(double));
		t_uy4 = malloc(sizeof(double));
		t_uz4 = malloc(sizeof(double));
		t_ux5 = malloc(sizeof(double));
		t_uy5 = malloc(sizeof(double));
		t_uz5 = malloc(sizeof(double));
		t_ux1_iele = malloc(sizeof(double));
		t_uy1_iele = malloc(sizeof(double));
		t_uz1_iele = malloc(sizeof(double));
		t_ux2_iele = malloc(sizeof(double));
		t_uy2_iele = malloc(sizeof(double));
		t_uz2_iele = malloc(sizeof(double));
		t_ux3_iele = malloc(sizeof(double));
		t_uy3_iele = malloc(sizeof(double));
		t_uz3_iele = malloc(sizeof(double));
		t_ux4_iele = malloc(sizeof(double));
		t_uy4_iele = malloc(sizeof(double));
		t_uz4_iele = malloc(sizeof(double));
		t_ux5_iele = malloc(sizeof(double));
		t_uy5_iele = malloc(sizeof(double));
		t_uz5_iele = malloc(sizeof(double));

		i = boun_cond(i, b_cnt, nele, boun_type, &j, cnt_boun_type, &iele,
				boun_elem, ghost_elem, t_ux1, t_uy1, t_uz1, t_ux2, t_uy2, t_uz2,
				t_ux3, t_uy3, t_uz3, t_ux4, t_uy4, t_uz4, t_ux5, t_uy5, t_uz5,
				t_ux1_iele, t_uy1_iele, t_uz1_iele, t_ux2_iele, t_uy2_iele,
				t_uz2_iele, t_ux3_iele, t_uy3_iele, t_uz3_iele, t_ux4_iele,
				t_uy4_iele, t_uz4_iele, t_ux5_iele, t_uy5_iele, t_uz5_iele,
				trans_coord_dx, trans_coord_dy, trans_coord_dz, trans_coord_mx,
				trans_coord_my, trans_coord_mz, nx_boun_face, ny_boun_face,
				nz_boun_face, u1_p, u2_p, u3_p, u4_p, u5_p, ux1_p, uy1_p, uz1_p,
				ux2_p, uy2_p, uz2_p, ux3_p, uy3_p, uz3_p, ux4_p, uy4_p, uz4_p,
				ux5_p, uy5_p, uz5_p, u1initial, u2initial, u3initial, u4initial,
				u5initial, periodic_elem, primitive, derv);

		sprintf(filename, "test1_%d.txt", id);
		fp = fopen(filename, "w");
		fprintf(fp, "there\n");
		fclose(fp);

		for (iele = boun_cnt_start; iele < boun_cnt_end; iele++)
		{

			sprintf(filename, "test2_%d.txt", id);
			fp = fopen(filename, "a");
			fprintf(fp, "%d %d %d\n", iele, boun_cnt_end, tot_ele);
			fclose(fp);
			secondary_variables(&u1_p[iele], &ux1_p[iele], &uy1_p[iele],
					&uz1_p[iele], &u2_p[iele], &ux2_p[iele], &uy2_p[iele],
					&uz2_p[iele], &u3_p[iele], &ux3_p[iele], &uy3_p[iele],
					&uz3_p[iele], &u4_p[iele], &ux4_p[iele], &uy4_p[iele],
					&uz4_p[iele], &u5_p[iele], &ux5_p[iele], &uy5_p[iele],
					&uz5_p[iele], &ut1_p[iele], &ut2_p[iele], &ut3_p[iele],
					&ut4_p[iele], &ut5_p[iele], &f1_p[iele], &fx1_p[iele],
					&fy1_p[iele], &fz1_p[iele], &g1_p[iele], &gx1_p[iele],
					&gy1_p[iele], &gz1_p[iele], &e1_p[iele], &ex1_p[iele],
					&ey1_p[iele], &ez1_p[iele], &ut1_p[iele], &ft1_p[iele],
					&gt1_p[iele], &et1_p[iele], &f2_p[iele], &fx2_p[iele],
					&fy2_p[iele], &fz2_p[iele], &g2_p[iele], &gx2_p[iele],
					&gy2_p[iele], &gz2_p[iele], &e2_p[iele], &ex2_p[iele],
					&ey2_p[iele], &ez2_p[iele], &ut2_p[iele], &ft2_p[iele],
					&gt2_p[iele], &et2_p[iele], &f3_p[iele], &fx3_p[iele],
					&fy3_p[iele], &fz3_p[iele], &g3_p[iele], &gx3_p[iele],
					&gy3_p[iele], &gz3_p[iele], &e3_p[iele], &ex3_p[iele],
					&ey3_p[iele], &ez3_p[iele], &ut3_p[iele], &ft3_p[iele],
					&gt3_p[iele], &et3_p[iele], &f4_p[iele], &fx4_p[iele],
					&fy4_p[iele], &fz4_p[iele], &g4_p[iele], &gx4_p[iele],
					&gy4_p[iele], &gz4_p[iele], &e4_p[iele], &ex4_p[iele],
					&ey4_p[iele], &ez4_p[iele], &ut4_p[iele], &ft4_p[iele],
					&gt4_p[iele], &et4_p[iele], &f5_p[iele], &fx5_p[iele],
					&fy5_p[iele], &fz5_p[iele], &g5_p[iele], &gx5_p[iele],
					&gy5_p[iele], &gz5_p[iele], &e5_p[iele], &ex5_p[iele],
					&ey5_p[iele], &ez5_p[iele], &ut5_p[iele], &ft5_p[iele],
					&gt5_p[iele], &et5_p[iele]);

		}

		iele = face_normal_calc(iele, lele_cnt, d2, d3, &jface, &s1, &s2, &s3,
				&d1, lelem_neighbour, aa1, sol, bb1, cc1, dd1, deno, num_dd,
				normal_x, normal_y, normal_z, t, xr, yr, zr);

//exit(0);

//flux and derivative calculation
		for (n = start_iter; n <= no_of_iter; n++)
		{

			iele = main_solver(iele, lele_cnt, ga, d1, nu1, r1, c, u5_p, u2_p,
					u3_p, u4_p, u1_p, maxprime, &jface, theta, xr, sol, yr, zr,
					maximum, &ineigh, lelem_neighbour, &iside, xc_side,
					xc_lateral_side1, yc_side, yc_lateral_side1, zc_side,
					zc_lateral_side1, xc_lateral_side2, yc_lateral_side2,
					zc_lateral_side2, xc_lateral_side3, yc_lateral_side3,
					zc_lateral_side3, u1_side, ux1_p, uy1_p, uz1_p, u2_side,
					ux2_p, uy2_p, uz2_p, u3_side, ux3_p, uy3_p, uz3_p, u4_side,
					ux4_p, uy4_p, uz4_p, u5_side, ux5_p, uy5_p, uz5_p, ux_side,
					uy_side, uz_side, vx_side, vy_side, vz_side, wx_side,
					wy_side, wz_side, tdash, tcont, mudash, scont, towxx, re,
					towyy, towzz, towxy, towyz, towzx, qx, pr, qy, qz, d_nu,
					delt, num_dd, cfl_iele, nu, &jnode, eleminfo, xp_node,
					xycoord, yp_node, zp_node, &s1, &s2, &s3, xpstar, ypstar,
					zpstar, &m, um, uxm, uym, uzm, utm, ut1_p, fm, f1_p, fxm,
					fx1_p, fym, fy1_p, fzm, fz1_p, ftm, ft1_p, gm, g1_p, gxm,
					gx1_p, gym, gy1_p, gzm, gz1_p, gtm, gt1_p, em, e1_p, exm,
					ex1_p, eym, ey1_p, ezm, ez1_p, etm, et1_p, fvm, gvm, evm,
					ut2_p, f2_p, fx2_p, fy2_p, fz2_p, ft2_p, g2_p, gx2_p, gy2_p,
					gz2_p, gt2_p, e2_p, ex2_p, ey2_p, ez2_p, et2_p, ut3_p, f3_p,
					fx3_p, fy3_p, fz3_p, ft3_p, g3_p, gx3_p, gy3_p, gz3_p,
					gt3_p, e3_p, ex3_p, ey3_p, ez3_p, et3_p, ut4_p, f4_p, fx4_p,
					fy4_p, fz4_p, ft4_p, g4_p, gx4_p, gy4_p, gz4_p, gt4_p, e4_p,
					ex4_p, ey4_p, ez4_p, et4_p, ut5_p, f5_p, fx5_p, fy5_p,
					fz5_p, ft5_p, g5_p, gx5_p, gy5_p, gz5_p, gt5_p, e5_p, ex5_p,
					ey5_p, ez5_p, et5_p, f1_side, area_lateral_side1,
					nx_lateral_side1, ny_lateral_side1, nz_lateral_side1,
					f2_side, area_lateral_side2, nx_lateral_side2,
					ny_lateral_side2, nz_lateral_side2, f3_side,
					area_lateral_side3, nx_lateral_side3, ny_lateral_side3,
					nz_lateral_side3, f4_side, vol_hexahedra_face,
					xc_hexahedra_face, yc_hexahedra_face, zc_hexahedra_face,
					flux, xpstar_face, ypstar_face, zpstar_face, umdash,
					primitive, vol_polyhedra, derv, secondary, max, &l);

			MPI_Send(maximum, 1, MPI_DOUBLE, master, id, MPI_COMM_WORLD);

//exchange of data for primitive,derv and secondary variable betweeen nodes

			for (i = 1; i <= kk; i++)
			{
				d_store = group[i];
				d_ele = k2id[d_store];
				MPI_Pack_size(send_cnt[d_ele], trans_primitive, MPI_COMM_WORLD,
						&maxsize3);
				MPI_Pack_size(send_cnt[d_ele], trans_derv, MPI_COMM_WORLD,
						&maxsize4);
				MPI_Pack_size(send_cnt[d_ele], trans_secondary, MPI_COMM_WORLD,
						&maxsize5);
				maxsize1 = maxsize3 + maxsize4 + maxsize5;
				buffer1 = malloc(maxsize1);

				position1 = 0;

				for (j = 1; j <= send_cnt[d_ele]; j++)
				{
					d_store1 = lelem_send[d_ele][j];
					MPI_Pack(&primitive[d_store1], 1, trans_primitive, buffer1,
							maxsize1, &position1, MPI_COMM_WORLD);
					MPI_Pack(&derv[d_store1], 1, trans_derv, buffer1, maxsize1,
							&position1, MPI_COMM_WORLD);
					MPI_Pack(&secondary[d_store1], 1, trans_secondary, buffer1,
							maxsize1, &position1, MPI_COMM_WORLD);
				}

				MPI_Pack_size(recvcnt[d_ele], trans_primitive, MPI_COMM_WORLD,
						&maxsize6);
				MPI_Pack_size(recvcnt[d_ele], trans_derv, MPI_COMM_WORLD,
						&maxsize7);
				MPI_Pack_size(recvcnt[d_ele], trans_secondary, MPI_COMM_WORLD,
						&maxsize8);
				maxsize2 = maxsize6 + maxsize7 + maxsize8;
				buffer2 = malloc(maxsize2);
				MPI_Sendrecv(buffer1, position1, MPI_PACKED, d_store + 1,
						d_store + 1, buffer2, maxsize2, MPI_PACKED, d_store + 1,
						id, MPI_COMM_WORLD, &status);

				position1 = 0;

				for (j = 1; j <= recvcnt[d_ele]; j++)
				{
					d_store2 = lelem_recv[d_ele][j];
					MPI_Unpack(buffer2, maxsize2, &position1,
							&primitive[d_store2], 1, trans_primitive,
							MPI_COMM_WORLD);
					MPI_Unpack(buffer2, maxsize2, &position1, &derv[d_store2],
							1, trans_derv, MPI_COMM_WORLD);
					MPI_Unpack(buffer2, maxsize2, &position1,
							&secondary[d_store2], 1, trans_secondary,
							MPI_COMM_WORLD);
				}

				free(buffer1);
				buffer1 = NULL;
				free(buffer2);
				buffer2 = NULL;

			} //i bracket

			if (pg == 1)
			{
//exchange of data for periodic elements
				if (periodic_id == 1)
				{
					for (i = 1; i <= pp; i++)
					{
						d_store = p_group[i];
						d_ele = pp2id[d_store];

						if (pp_send[d_ele] == 1)
						{

							MPI_Pack_size(p_send_cnt[d_ele], trans_primitive,
									MPI_COMM_WORLD, &maxsize3);
							MPI_Pack_size(p_send_cnt[d_ele], trans_derv,
									MPI_COMM_WORLD, &maxsize4);

							maxsize1 = maxsize3 + maxsize4;
							buffer1 = malloc(maxsize1);

							position1 = 0;

							for (j = 1; j <= p_send_cnt[d_ele]; j++)
							{
								d_store1 = p_elem_send[d_ele][j];
								MPI_Pack(&primitive[d_store1], 1,
										trans_primitive, buffer1, maxsize1,
										&position1, MPI_COMM_WORLD);
								MPI_Pack(&derv[d_store1], 1, trans_derv,
										buffer1, maxsize1, &position1,
										MPI_COMM_WORLD);

							}
							MPI_Send(buffer1, position1, MPI_PACKED,
									d_store + 1, d_store + 1, MPI_COMM_WORLD);

							free(buffer1);
							buffer1 = NULL;
						}

						else if (pp_recv[d_ele] == 1)
						{

							MPI_Pack_size(p_recvcnt[d_ele], trans_primitive,
									MPI_COMM_WORLD, &maxsize6);
							MPI_Pack_size(p_recvcnt[d_ele], trans_derv,
									MPI_COMM_WORLD, &maxsize7);

							maxsize2 = maxsize6 + maxsize7;
							buffer2 = malloc(maxsize2);

							MPI_Recv(buffer2, maxsize2, MPI_PACKED, d_store + 1,
									id, MPI_COMM_WORLD, &status);

							position1 = 0;

							for (j = 1; j <= p_recvcnt[d_ele]; j++)
							{
								d_store2 = p_elem_recv[d_ele][j];
								MPI_Unpack(buffer2, maxsize2, &position1,
										&primitive[d_store2], 1,
										trans_primitive, MPI_COMM_WORLD);
								MPI_Unpack(buffer2, maxsize2, &position1,
										&derv[d_store2], 1, trans_derv,
										MPI_COMM_WORLD);

							}

							free(buffer2);
							buffer2 = NULL;
						}

						else if (pp_sendrecv[d_ele] == 1)
						{
							MPI_Pack_size(p_send_cnt[d_ele], trans_primitive,
									MPI_COMM_WORLD, &maxsize3);
							MPI_Pack_size(p_send_cnt[d_ele], trans_derv,
									MPI_COMM_WORLD, &maxsize4);

							maxsize1 = maxsize3 + maxsize4;
							buffer1 = malloc(maxsize1);

							position1 = 0;

							for (j = 1; j <= p_send_cnt[d_ele]; j++)
							{
								d_store1 = p_elem_send[d_ele][j];
								MPI_Pack(&primitive[d_store1], 1,
										trans_primitive, buffer1, maxsize1,
										&position1, MPI_COMM_WORLD);
								MPI_Pack(&derv[d_store1], 1, trans_derv,
										buffer1, maxsize1, &position1,
										MPI_COMM_WORLD);

							}

							MPI_Pack_size(p_recvcnt[d_ele], trans_primitive,
									MPI_COMM_WORLD, &maxsize6);
							MPI_Pack_size(p_recvcnt[d_ele], trans_derv,
									MPI_COMM_WORLD, &maxsize7);

							maxsize2 = maxsize6 + maxsize7;
							buffer2 = malloc(maxsize2);
							MPI_Sendrecv(buffer1, position1, MPI_PACKED,
									d_store + 1, d_store + 1, buffer2, maxsize2,
									MPI_PACKED, d_store + 1, id, MPI_COMM_WORLD,
									&status);

							position1 = 0;

							for (j = 1; j <= p_recvcnt[d_ele]; j++)
							{
								d_store2 = p_elem_recv[d_ele][j];
								MPI_Unpack(buffer2, maxsize2, &position1,
										&primitive[d_store2], 1,
										trans_primitive, MPI_COMM_WORLD);
								MPI_Unpack(buffer2, maxsize2, &position1,
										&derv[d_store2], 1, trans_derv,
										MPI_COMM_WORLD);

							}

							free(buffer1);
							buffer1 = NULL;
							free(buffer2);
							buffer2 = NULL;
						}

					} //i bracket
				}
			}

//allocating values for u[]_p:previous values

			iele = previous_time_step_alloc(iele, subdomain_cnt, u1_p,
					primitive, ux1_p, derv, uy1_p, uz1_p, ut1_p, secondary,
					f1_p, fx1_p, fy1_p, fz1_p, ft1_p, g1_p, gx1_p, gy1_p, gz1_p,
					gt1_p, e1_p, ex1_p, ey1_p, ez1_p, et1_p, u2_p, ux2_p, uy2_p,
					uz2_p, ut2_p, f2_p, fx2_p, fy2_p, fz2_p, ft2_p, g2_p, gx2_p,
					gy2_p, gz2_p, gt2_p, e2_p, ex2_p, ey2_p, ez2_p, et2_p, u3_p,
					ux3_p, uy3_p, uz3_p, ut3_p, f3_p, fx3_p, fy3_p, fz3_p,
					ft3_p, g3_p, gx3_p, gy3_p, gz3_p, gt3_p, e3_p, ex3_p, ey3_p,
					ez3_p, et3_p, u4_p, ux4_p, uy4_p, uz4_p, ut4_p, f4_p, fx4_p,
					fy4_p, fz4_p, ft4_p, g4_p, gx4_p, gy4_p, gz4_p, gt4_p, e4_p,
					ex4_p, ey4_p, ez4_p, et4_p, u5_p, ux5_p, uy5_p, uz5_p,
					ut5_p, f5_p, fx5_p, fy5_p, fz5_p, ft5_p, g5_p, gx5_p, gy5_p,
					gz5_p, gt5_p, e5_p, ex5_p, ey5_p, ez5_p, et5_p);

//boundary conditions

			i = boun_cond(i, b_cnt, nele, boun_type, &j, cnt_boun_type, &iele,
					boun_elem, ghost_elem, t_ux1, t_uy1, t_uz1, t_ux2, t_uy2,
					t_uz2, t_ux3, t_uy3, t_uz3, t_ux4, t_uy4, t_uz4, t_ux5,
					t_uy5, t_uz5, t_ux1_iele, t_uy1_iele, t_uz1_iele,
					t_ux2_iele, t_uy2_iele, t_uz2_iele, t_ux3_iele, t_uy3_iele,
					t_uz3_iele, t_ux4_iele, t_uy4_iele, t_uz4_iele, t_ux5_iele,
					t_uy5_iele, t_uz5_iele, trans_coord_dx, trans_coord_dy,
					trans_coord_dz, trans_coord_mx, trans_coord_my,
					trans_coord_mz, nx_boun_face, ny_boun_face, nz_boun_face,
					u1_p, u2_p, u3_p, u4_p, u5_p, ux1_p, uy1_p, uz1_p, ux2_p,
					uy2_p, uz2_p, ux3_p, uy3_p, uz3_p, ux4_p, uy4_p, uz4_p,
					ux5_p, uy5_p, uz5_p, u1initial, u2initial, u3initial,
					u4initial, u5initial, periodic_elem, primitive, derv);

			for (iele = boun_cnt_start; iele < boun_cnt_end; iele++)
			{

				secondary_variables(&u1_p[iele], &ux1_p[iele], &uy1_p[iele],
						&uz1_p[iele], &u2_p[iele], &ux2_p[iele], &uy2_p[iele],
						&uz2_p[iele], &u3_p[iele], &ux3_p[iele], &uy3_p[iele],
						&uz3_p[iele], &u4_p[iele], &ux4_p[iele], &uy4_p[iele],
						&uz4_p[iele], &u5_p[iele], &ux5_p[iele], &uy5_p[iele],
						&uz5_p[iele], &ut1_p[iele], &ut2_p[iele], &ut3_p[iele],
						&ut4_p[iele], &ut5_p[iele], &f1_p[iele], &fx1_p[iele],
						&fy1_p[iele], &fz1_p[iele], &g1_p[iele], &gx1_p[iele],
						&gy1_p[iele], &gz1_p[iele], &e1_p[iele], &ex1_p[iele],
						&ey1_p[iele], &ez1_p[iele], &ut1_p[iele], &ft1_p[iele],
						&gt1_p[iele], &et1_p[iele], &f2_p[iele], &fx2_p[iele],
						&fy2_p[iele], &fz2_p[iele], &g2_p[iele], &gx2_p[iele],
						&gy2_p[iele], &gz2_p[iele], &e2_p[iele], &ex2_p[iele],
						&ey2_p[iele], &ez2_p[iele], &ut2_p[iele], &ft2_p[iele],
						&gt2_p[iele], &et2_p[iele], &f3_p[iele], &fx3_p[iele],
						&fy3_p[iele], &fz3_p[iele], &g3_p[iele], &gx3_p[iele],
						&gy3_p[iele], &gz3_p[iele], &e3_p[iele], &ex3_p[iele],
						&ey3_p[iele], &ez3_p[iele], &ut3_p[iele], &ft3_p[iele],
						&gt3_p[iele], &et3_p[iele], &f4_p[iele], &fx4_p[iele],
						&fy4_p[iele], &fz4_p[iele], &g4_p[iele], &gx4_p[iele],
						&gy4_p[iele], &gz4_p[iele], &e4_p[iele], &ex4_p[iele],
						&ey4_p[iele], &ez4_p[iele], &ut4_p[iele], &ft4_p[iele],
						&gt4_p[iele], &et4_p[iele], &f5_p[iele], &fx5_p[iele],
						&fy5_p[iele], &fz5_p[iele], &g5_p[iele], &gx5_p[iele],
						&gy5_p[iele], &gz5_p[iele], &e5_p[iele], &ex5_p[iele],
						&ey5_p[iele], &ez5_p[iele], &ut5_p[iele], &ft5_p[iele],
						&gt5_p[iele], &et5_p[iele]);

			}

//data to be sent to master for writing

			if (n % ns == 0)
			{
				MPI_Pack_size(lele_cnt, trans_primitive, MPI_COMM_WORLD,
						&maxsize1);
				buffer1 = malloc(maxsize1);
				position1 = 0;
				for (j = 1; j <= lele_cnt; j++)
				{
					MPI_Pack(&primitive[j], 1, trans_primitive, buffer1,
							maxsize1, &position1, MPI_COMM_WORLD);
				}
				MPI_Send(buffer1, position1, MPI_PACKED, master, id,
						MPI_COMM_WORLD);
				free(buffer1);
				buffer1 = NULL;
			}

			if (n % 10000 == 0)
			{
				MPI_Pack_size(lele_cnt, trans_derv, MPI_COMM_WORLD, &maxsize1);
				buffer1 = malloc(maxsize1);
				position1 = 0;
				for (j = 1; j <= lele_cnt; j++)
				{
					MPI_Pack(&derv[j], 1, trans_derv, buffer1, maxsize1,
							&position1, MPI_COMM_WORLD);
				}
				MPI_Send(buffer1, position1, MPI_PACKED, master, id,
						MPI_COMM_WORLD);
				free(buffer1);
				buffer1 = NULL;
			}

			MPI_Barrier(comm_slaves);

		} //n bracket

	} //id >0 bracket

	MPI_Finalize();

	return 0;

}


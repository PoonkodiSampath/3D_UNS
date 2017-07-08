/*
 * slave_solver.h
 *
 *  Created on: 12-Feb-2017
 *      Author: poonkodi
 */

#ifndef SLAVE_SOLVER_H_
#define SLAVE_SOLVER_H_

void sigma_cal(double *xstart_d1, double *ystart_d1, double *zstart_d1,
		double *xstart_d2, double *ystart_d2, double *zstart_d2,
		double *xstart_d3, double *ystart_d3, double *zstart_d3,
		double *xstart_d4, double *ystart_d4, double *zstart_d4, double *sigma)
{
	a1 = *xstart_d1 - *xstart_d4;
	a2 = *ystart_d1 - *ystart_d4;
	a3 = *zstart_d1 - *zstart_d4;
	a4 = *xstart_d2 - *xstart_d4;
	a5 = *ystart_d2 - *ystart_d4;
	a6 = *zstart_d2 - *zstart_d4;
	a7 = *xstart_d3 - *xstart_d4;
	a8 = *ystart_d3 - *ystart_d4;
	a9 = *zstart_d3 - *zstart_d4;

	*sigma = ((a1 * (a5 * a9 - a8 * a6)) - (a2 * (a4 * a9 - a7 * a6))
			+ (a3 * (a4 * a8 - a7 * a5)));
}

void derv_cal(double *um_dash11, double *um_dash21, double *um_dash31,
		double *um_dash12, double *um_dash22, double *um_dash32,
		double *um_dash13, double *um_dash23, double *um_dash33,
		double *um_dash14, double *um_dash24, double *um_dash34,
		double *x_star11, double *x_star21, double *x_star31, double *x_star12,
		double *x_star22, double *x_star32, double *x_star13, double *x_star23,
		double *x_star33, double *x_star14, double *x_star24, double *x_star34,
		double *y_star11, double *y_star21, double *y_star31, double *y_star12,
		double *y_star22, double *y_star32, double *y_star13, double *y_star23,
		double *y_star33, double *y_star14, double *y_star24, double *y_star34,
		double *z_star11, double *z_star21, double *z_star31, double *z_star12,
		double *z_star22, double *z_star32, double *z_star13, double *z_star23,
		double *z_star33, double *z_star14, double *z_star24, double *z_star34,
		double *u_m, double *sp_x, double *sp_y, double *sp_z, double *ux_m,
		double *uy_m, double *uz_m, double *cfl)
{
	int iderv, s1, s2, s3, p;
	double sigma_x, sigma_y, sigma_z, e;
	double umx[5], umy[5], umz[5];
	double o[5], eta[5], w1, w2, w3, w4, sm_dash[5], f_cfl;
	double osq[5];
	double omega;
	double umdash[4][5], xstar[4][5], ystar[4][5], zstar[4][5], sigmadeno;
	umdash[1][1] = *um_dash11;
	umdash[2][1] = *um_dash21;
	umdash[3][1] = *um_dash31;
	umdash[1][2] = *um_dash12;
	umdash[2][2] = *um_dash22;
	umdash[3][2] = *um_dash32;
	umdash[1][3] = *um_dash13;
	umdash[2][3] = *um_dash23;
	umdash[3][3] = *um_dash33;
	umdash[1][4] = *um_dash14;
	umdash[2][4] = *um_dash24;
	umdash[3][4] = *um_dash34;

	xstar[1][1] = *x_star11;
	xstar[2][1] = *x_star21;
	xstar[3][1] = *x_star31;
	xstar[1][2] = *x_star12;
	xstar[2][2] = *x_star22;
	xstar[3][2] = *x_star32;
	xstar[1][3] = *x_star13;
	xstar[2][3] = *x_star23;
	xstar[3][3] = *x_star33;
	xstar[1][4] = *x_star14;
	xstar[2][4] = *x_star24;
	xstar[3][4] = *x_star34;

	ystar[1][1] = *y_star11;
	ystar[2][1] = *y_star21;
	ystar[3][1] = *y_star31;
	ystar[1][2] = *y_star12;
	ystar[2][2] = *y_star22;
	ystar[3][2] = *y_star32;
	ystar[1][3] = *y_star13;
	ystar[2][3] = *y_star23;
	ystar[3][3] = *y_star33;
	ystar[1][4] = *y_star14;
	ystar[2][4] = *y_star24;
	ystar[3][4] = *y_star34;

	zstar[1][1] = *z_star11;
	zstar[2][1] = *z_star21;
	zstar[3][1] = *z_star31;
	zstar[1][2] = *z_star12;
	zstar[2][2] = *z_star22;
	zstar[3][2] = *z_star32;
	zstar[1][3] = *z_star13;
	zstar[2][3] = *z_star23;
	zstar[3][3] = *z_star33;
	zstar[1][4] = *z_star14;
	zstar[2][4] = *z_star24;
	zstar[3][4] = *z_star34;

	for (iderv = 1; iderv <= 4; iderv++)
	{

		sigma_cal(&umdash[1][iderv], &ystar[1][iderv], &zstar[1][iderv],
				&umdash[2][iderv], &ystar[2][iderv], &zstar[2][iderv],
				&umdash[3][iderv], &ystar[3][iderv], &zstar[3][iderv], u_m,
				sp_y, sp_z, &sigma_x);

		sigma_cal(&xstar[1][iderv], &umdash[1][iderv], &zstar[1][iderv],
				&xstar[2][iderv], &umdash[2][iderv], &zstar[2][iderv],
			&xstar[3][iderv], &umdash[3][iderv], &zstar[3][iderv], sp_x,
				u_m, sp_z, &sigma_y);

		sigma_cal(&xstar[1][iderv], &ystar[1][iderv], &umdash[1][iderv],
				&xstar[2][iderv], &ystar[2][iderv], &umdash[2][iderv],
				&xstar[3][iderv], &ystar[3][iderv], &umdash[3][iderv], sp_x,
				sp_y, u_m, &sigma_z);

		sigma_cal(&xstar[1][iderv], &ystar[1][iderv], &zstar[1][iderv],
				&xstar[2][iderv], &ystar[2][iderv], &zstar[2][iderv],
				&xstar[3][iderv], &ystar[3][iderv], &zstar[3][iderv], sp_x,
				sp_y, sp_z, &sigmadeno);

		umx[iderv] = sigma_x / (sigmadeno);
		umy[iderv] = sigma_y / (sigmadeno);
		umz[iderv] = sigma_z / (sigmadeno);
	}
	for (p = 1; p <= 4; p++)
	{
		osq[p] = sqrt(
				(umx[p] * umx[p]) + (umy[p] * umy[p]) + (umz[p] * umz[p]));
	}
	o[1] = osq[2] * osq[3] * osq[4];
	o[2] = osq[1] * osq[3] * osq[4];
	o[3] = osq[1] * osq[2] * osq[4];
	o[4] = osq[1] * osq[2] * osq[3];

	omega = o[1] + o[2] + o[3] + o[4] + 10e-20;

	o[1] = o[1] / omega;
	o[2] = o[2] / omega;
	o[3] = o[3] / omega;
	o[4] = o[4] / omega;

	for (s1 = 1; s1 < 4; s1++)
	{
		for (s2 = 4; s2 > s1; s2--)
		{
			if (o[s1] > o[s2])
			{
				omega = o[s1];
				o[s1] = o[s2];
				o[s2] = omega;
				omega = umx[s1];
				umx[s1] = umx[s2];
				umx[s2] = omega;
				omega = umy[s1];
				umy[s1] = umy[s2];
				umy[s2] = omega;
				omega = umz[s1];
				umz[s1] = umz[s2];
				umz[s2] = omega;
			}
		}
	}
	if (o[1] > 0.0)
	{
		for (p = 1; p <= 3; p++)
		{
			if (fabs(o[p]) < 0.00005)
			{
				eta[p] = -1.0;
			}
			else
			{
				eta[p] = (o[p + 1] / o[p]) - 1.0;
			}
		}

		if (*cfl > 0.3)
		{
			f_cfl = 1 / (*cfl);
		}
		else
		{
			*cfl = 0.3;
			f_cfl = 1 / (*cfl);
		}

		sm_dash[1] = o[1];

		for (p = 1; p <= 3; p++)
		{
			sm_dash[p + 1] = (1 + f_cfl * eta[p]) * sm_dash[p];
		}

		e = sm_dash[1] + sm_dash[2] + sm_dash[3] + sm_dash[4];

		if (fabs(e) < 0.000005)
		{
			for (p = 1; p <= 4; p++)
			{
				*ux_m = 0.0;
				*uy_m = 0.0;
				*uz_m = 0.0;
			}
		}
		else
		{

			w1 = sm_dash[1] / e;
			w2 = sm_dash[2] / e;
			w3 = sm_dash[3] / e;
			w4 = sm_dash[4] / e;

			*ux_m = (w1 * umx[1] + w2 * umx[2] + w3 * umx[3] + w4 * umx[4]);
			*uy_m = (w1 * umy[1] + w2 * umy[2] + w3 * umy[3] + w4 * umy[4]);
			*uz_m = (w1 * umz[1] + w2 * umz[2] + w3 * umz[3] + w4 * umz[4]);

		}
	}
	else
	{
		*ux_m = 0.0;
		*uy_m = 0.0;
		*uz_m = 0.0;
	}

}

void secondary_variables(double *p_u1, double *p_u1x, double *p_u1y,
		double *p_u1z, double *p_u2, double *p_u2x, double *p_u2y,
		double *p_u2z, double *p_u3, double *p_u3x, double *p_u3y,
		double *p_u3z, double *p_u4, double *p_u4x, double *p_u4y,
		double *p_u4z, double *p_u5, double *p_u5x, double *p_u5y,
		double *p_u5z, double *p_u1t, double *p_u2t, double *p_u3t,
		double *p_u4t, double *p_u5t, double *s_f1, double *s_f1x,
		double *s_f1y, double *s_f1z, double *s_g1, double *s_g1x,
		double *s_g1y, double *s_g1z, double *s_e1, double *s_e1x,
		double *s_e1y, double *s_e1z, double *s_u1t, double *s_f1t,
		double *s_g1t, double *s_e1t, double *s_f2, double *s_f2x,
		double *s_f2y, double *s_f2z, double *s_g2, double *s_g2x,
		double *s_g2y, double *s_g2z, double *s_e2, double *s_e2x,
		double *s_e2y, double *s_e2z, double *s_u2t, double *s_f2t,
		double *s_g2t, double *s_e2t, double *s_f3, double *s_f3x,
		double *s_f3y, double *s_f3z, double *s_g3, double *s_g3x,
		double *s_g3y, double *s_g3z, double *s_e3, double *s_e3x,
		double *s_e3y, double *s_e3z, double *s_u3t, double *s_f3t,
		double *s_g3t, double *s_e3t, double *s_f4, double *s_f4x,
		double *s_f4y, double *s_f4z, double *s_g4, double *s_g4x,
		double *s_g4y, double *s_g4z, double *s_e4, double *s_e4x,
		double *s_e4y, double *s_e4z, double *s_u4t, double *s_f4t,
		double *s_g4t, double *s_e4t, double *s_f5, double *s_f5x,
		double *s_f5y, double *s_f5z, double *s_g5, double *s_g5x,
		double *s_g5y, double *s_g5z, double *s_e5, double *s_e5x,
		double *s_e5y, double *s_e5z, double *s_u5t, double *s_f5t,
		double *s_g5t, double *s_e5t)
{
	*s_f1 = *p_u2;
	*s_f1x = *p_u2x;
	*s_f1y = *p_u2y;
	*s_f1z = *p_u2z;

	*s_g1 = *p_u3;
	*s_g1x = *p_u3x;
	*s_g1y = *p_u3y;
	*s_g1z = *p_u3z;

	*s_e1 = *p_u4;
	*s_e1x = *p_u4x;
	*s_e1y = *p_u4y;
	*s_e1z = *p_u4z;

	*s_u1t = -*s_f1x - *s_g1y - *s_e1z;

	*s_f2 = ((*p_u5) * (ga - 1.0))
			+ ((3.0 - ga) * (*p_u2 * *p_u2) / (2.0 * *p_u1))
			- ((ga - 1.0) * (*p_u3 * *p_u3 + *p_u4 * *p_u4) / (2 * *p_u1));

	f21 = (((ga - 1) * (*p_u3 * *p_u3 + *p_u4 * *p_u4) / 2.0)
			- ((3.0 - ga) * (*p_u2 * *p_u2) / 2.0)) / (*p_u1 * *p_u1);

	f22 = (3.0 - ga) * (*p_u2) / (*p_u1);

	f23 = -(ga - 1.0) * *p_u3 / (*p_u1);

	f24 = -(ga - 1.0) * *p_u4 / (*p_u1);

	f25 = ga - 1.0;

	*s_f2x = f21 * *p_u1x + f22 * *p_u2x + f23 * *p_u3x + f24 * *p_u4x
			+ f25 * *p_u5x;
	*s_f2y = f21 * *p_u1y + f22 * *p_u2y + f23 * *p_u3y + f24 * *p_u4y
			+ f25 * *p_u5y;
	*s_f2z = f21 * *p_u1z + f22 * *p_u2z + f23 * *p_u3z + f24 * *p_u4z
			+ f25 * *p_u5z;

	*s_g2 = (*p_u2 * *p_u3) / (*p_u1);

	g21 = -(*p_u2 * *p_u3) / (*p_u1 * *p_u1);

	g22 = *p_u3 / (*p_u1);

	g23 = *p_u2 / (*p_u1);

	*s_g2x = g21 * *p_u1x + g22 * *p_u2x + g23 * *p_u3x;
	*s_g2y = g21 * *p_u1y + g22 * *p_u2y + g23 * *p_u3y;
	*s_g2z = g21 * *p_u1z + g22 * *p_u2z + g23 * *p_u3z;

	*s_e2 = *p_u2 * *p_u4 / (*p_u1);

	e21 = -(*p_u2 * *p_u4) / (*p_u1 * *p_u1);

	e22 = *p_u4 / (*p_u1);

	e24 = *p_u2 / (*p_u1);

	*s_e2x = e21 * *p_u1x + e22 * *p_u2x + e24 * *p_u4x;
	*s_e2y = e21 * *p_u1y + e22 * *p_u2y + e24 * *p_u4y;
	*s_e2z = e21 * *p_u1z + e22 * *p_u2z + e24 * *p_u4z;

	*s_u2t = -*s_f2x - *s_g2y - *s_e2z;

	*s_f3 = *s_g2;

	*s_f3x = *s_g2x;
	*s_f3y = *s_g2y;
	*s_f3z = *s_g2z;

	*s_g3 = *p_u5 * (ga - 1.0) + ((3.0 - ga) * (*p_u3 * *p_u3) / (2.0 * *p_u1))
			- ((ga - 1.0) * (*p_u2 * *p_u2 + *p_u4 * *p_u4) / (2 * *p_u1));

	g31 = (((ga - 1.0) * (*p_u2 * *p_u2 + *p_u4 * *p_u4) / 2.0)
			- ((3.0 - ga) * (*p_u3 * *p_u3) / 2.0)) / (*p_u1 * *p_u1);

	g32 = -(ga - 1) * *p_u2 / (*p_u1);

	g33 = (3.0 - ga) * *p_u3 / (*p_u1);

	g34 = f24;

	g35 = f25;

	*s_g3x = g31 * *p_u1x + g32 * *p_u2x + g33 * *p_u3x + g34 * *p_u4x
			+ g35 * *p_u5x;
	*s_g3y = g31 * *p_u1y + g32 * *p_u2y + g33 * *p_u3y + g34 * *p_u4y
			+ g35 * *p_u5y;
	*s_g3z = g31 * *p_u1z + g32 * *p_u2z + g33 * *p_u3z + g34 * *p_u4z
			+ g35 * *p_u5z;

	*s_e3 = *p_u3 * *p_u4 / (*p_u1);

	e31 = -(*p_u3 * *p_u4) / (*p_u1 * *p_u1);

	e33 = e22;

	e34 = g22;

	*s_e3x = e31 * *p_u1x + e33 * *p_u3x + e34 * *p_u4x;
	*s_e3y = e31 * *p_u1y + e33 * *p_u3y + e34 * *p_u4y;
	*s_e3z = e31 * *p_u1z + e33 * *p_u3z + e34 * *p_u4z;

	*s_u3t = -*s_f3x - *s_g3y - *s_e3z;

	*s_f4 = *s_e2;

	*s_f4x = *s_e2x;
	*s_f4y = *s_e2y;
	*s_f4z = *s_e2z;

	*s_g4 = *s_e3;

	*s_g4x = *s_e3x;
	*s_g4y = *s_e3y;
	*s_g4z = *s_e3z;

	*s_e4 = (*p_u5 * (ga - 1.0))
			+ ((3.0 - ga) * (*p_u4 * *p_u4) / (2.0 * *p_u1))
			- ((ga - 1.0) * (*p_u2 * *p_u2 + *p_u3 * *p_u3) / (2.0 * *p_u1));

	e41 = (((ga - 1.0) * (*p_u2 * *p_u2 + *p_u3 * *p_u3) / 2.0)
			- ((3.0 - ga) * (*p_u4 * *p_u4) / 2.0)) / (*p_u1 * *p_u1);

	e42 = g32;

	e43 = f23;

	e44 = (3.0 - ga) * *p_u4 / (*p_u1);

	e45 = f25;

	*s_e4x = e41 * *p_u1x + e42 * *p_u2x + e43 * *p_u3x + e44 * *p_u4x
			+ e45 * *p_u5x;
	*s_e4y = e41 * *p_u1y + e42 * *p_u2y + e43 * *p_u3y + e44 * *p_u4y
			+ e45 * *p_u5y;
	*s_e4z = e41 * *p_u1z + e42 * *p_u2z + e43 * *p_u3z + e44 * *p_u4z
			+ e45 * *p_u5z;

	*s_u4t = -*s_f4x - *s_g4y - *s_e4z;

	*s_f5 = (*p_u2 * *p_u5 * ga / (*p_u1))
			- ((ga - 1.0)
					* (*p_u2 * *p_u2 * *p_u2 + *p_u3 * *p_u3 * *p_u2
							+ *p_u4 * *p_u4 * *p_u2) / (2.0 * *p_u1 * *p_u1));

	f51 = -(*p_u2 * *p_u5 * ga / (*p_u1 * *p_u1))
			+ ((ga - 1.0)
					* (*p_u2 * *p_u2 * *p_u2 + *p_u3 * *p_u3 * *p_u2
							+ *p_u4 * *p_u4 * *p_u2) / (*p_u1 * *p_u1 * *p_u1));

	f52 = (*p_u5 * ga / (*p_u1))
			- ((ga - 1.0)
					* (3.0 * *p_u2 * *p_u2 + *p_u3 * *p_u3 + *p_u4 * *p_u4)
					/ (2.0 * *p_u1 * *p_u1));

	f53 = -(ga - 1.0) * (*p_u3 * *p_u2) / (*p_u1 * *p_u1);

	f54 = -(ga - 1.0) * (*p_u4 * *p_u2) / (*p_u1 * *p_u1);

	f55 = *p_u2 * ga / (*p_u1);

	*s_f5x = f51 * *p_u1x + f52 * *p_u2x + f53 * *p_u3x + f54 * *p_u4x
			+ f55 * *p_u5x;
	*s_f5y = f51 * *p_u1y + f52 * *p_u2y + f53 * *p_u3y + f54 * *p_u4y
			+ f55 * *p_u5y;
	*s_f5z = f51 * *p_u1z + f52 * *p_u2z + f53 * *p_u3z + f54 * *p_u4z
			+ f55 * *p_u5z;

	*s_g5 = (*p_u3 * *p_u5 * ga / (*p_u1))
			- ((ga - 1.0)
					* (*p_u2 * *p_u2 * *p_u3 + *p_u3 * *p_u3 * *p_u3
							+ *p_u4 * *p_u4 * *p_u3) / (2.0 * *p_u1 * *p_u1));

	g51 = -(*p_u3 * *p_u5 * ga / (*p_u1 * *p_u1))
			+ ((ga - 1.0)
					* (*p_u2 * *p_u2 * *p_u3 + *p_u3 * *p_u3 * *p_u3
							+ *p_u4 * *p_u4 * *p_u3) / (*p_u1 * *p_u1 * *p_u1));

	g52 = f53;

	g53 = (*p_u5 * ga / (*p_u1))
			- ((ga - 1.0)
					* (*p_u2 * *p_u2 + 3.0 * *p_u3 * *p_u3 + *p_u4 * *p_u4)
					/ (2.0 * *p_u1 * *p_u1));

	g54 = -(ga - 1.0) * (*p_u4 * *p_u3) / (*p_u1 * *p_u1);

	g55 = *p_u3 * ga / (*p_u1);

	*s_g5x = g51 * *p_u1x + g52 * *p_u2x + g53 * *p_u3x + g54 * *p_u4x
			+ g55 * *p_u5x;
	*s_g5y = g51 * *p_u1y + g52 * *p_u2y + g53 * *p_u3y + g54 * *p_u4y
			+ g55 * *p_u5y;
	*s_g5z = g51 * *p_u1z + g52 * *p_u2z + g53 * *p_u3z + g54 * *p_u4z
			+ g55 * *p_u5z;

	*s_e5 = (*p_u4 * *p_u5 * ga / (*p_u1))
			- ((ga - 1.0)
					* (*p_u2 * *p_u2 * *p_u4 + *p_u3 * *p_u3 * *p_u4
							+ *p_u4 * *p_u4 * *p_u4) / (2.0 * *p_u1 * *p_u1));

	e51 = -(*p_u4 * *p_u5 * ga / (*p_u1 * *p_u1))
			+ ((ga - 1.0)
					* (*p_u2 * *p_u2 * *p_u4 + *p_u3 * *p_u3 * *p_u4
							+ *p_u4 * *p_u4 * *p_u4) / (*p_u1 * *p_u1 * *p_u1));

	e52 = -(ga - 1.0) * (*p_u4 * *p_u2) / (*p_u1 * *p_u1);

	e53 = g54;

	e54 = (*p_u5 * ga / (*p_u1))
			- ((ga - 1.0)
					* (*p_u2 * *p_u2 + *p_u3 * *p_u3 + 3.0 * *p_u4 * *p_u4)
					/ (2.0 * *p_u1 * *p_u1));

	e55 = *p_u4 * ga / (*p_u1);

	*s_e5x = e51 * *p_u1x + e52 * *p_u2x + e53 * *p_u3x + e54 * *p_u4x
			+ e55 * *p_u5x;
	*s_e5y = e51 * *p_u1y + e52 * *p_u2y + e53 * *p_u3y + e54 * *p_u4y
			+ e55 * *p_u5y;
	*s_e5z = e51 * *p_u1z + e52 * *p_u2z + e53 * *p_u3z + e54 * *p_u4z
			+ e55 * *p_u5z;

	*s_u5t = -*s_f5x - *s_g5y - *s_e5z;

	*s_f1t = *s_u2t;
	*s_g1t = *s_u3t;
	*s_e1t = *s_u4t;

	*s_f2t = f21 * *p_u1t + f22 * *p_u2t + f23 * *p_u3t + f24 * *p_u4t
			+ f25 * *p_u5t;
	*s_g2t = g21 * *p_u1t + g22 * *p_u2t + g23 * *p_u3t;
	*s_e2t = e21 * *p_u1t + e22 * *p_u2t + e24 * *p_u4t;

	*s_f3t = *s_g2t;
	*s_g3t = g31 * *p_u1t + g32 * *p_u2t + g33 * *p_u3t + g34 * *p_u4t
			+ g35 * *p_u5t;
	*s_e3t = e31 * *p_u1t + e33 * *p_u3t + e34 * *p_u4t;

	*s_f4t = *s_e2t;
	*s_g4t = *s_e3t;
	*s_e4t = e41 * *p_u1t + e42 * *p_u2t + e43 * *p_u3t + e44 * *p_u4t
			+ e45 * *p_u5t;

	*s_f5t = f51 * *p_u1t + f52 * *p_u2t + f53 * *p_u3t + f54 * *p_u4t
			+ f55 * *p_u5t;
	*s_g5t = g51 * *p_u1t + g52 * *p_u2t + g53 * *p_u3t + g54 * *p_u4t
			+ g55 * *p_u5t;
	*s_e5t = e51 * *p_u1t + e52 * *p_u2t + e53 * *p_u3t + e54 * *p_u4t
			+ e55 * *p_u5t;
}

int main_solver(int iele, int lele_cnt, double ga, int d1, double nu1, int r1,
		double* c, double* u5_p, double* u2_p, double* u3_p, double* u4_p,
		double* u1_p, double* maxprime, int* jface, double** theta, double** xr,
		struct node4* sol, double** yr, double** zr, double* maximum,
		int* ineigh, int** lelem_neighbour, int* iside, double* xc_side,
		double** xc_lateral_side1, double* yc_side, double** yc_lateral_side1,
		double* zc_side, double** zc_lateral_side1, double** xc_lateral_side2,
		double** yc_lateral_side2, double** zc_lateral_side2,
		double** xc_lateral_side3, double** yc_lateral_side3,
		double** zc_lateral_side3, double** u1_side, double* ux1_p,
		double* uy1_p, double* uz1_p, double** u2_side, double* ux2_p,
		double* uy2_p, double* uz2_p, double** u3_side, double* ux3_p,
		double* uy3_p, double* uz3_p, double** u4_side, double* ux4_p,
		double* uy4_p, double* uz4_p, double** u5_side, double* ux5_p,
		double* uy5_p, double* uz5_p, double** ux_side, double** uy_side,
		double** uz_side, double** vx_side, double** vy_side, double** vz_side,
		double** wx_side, double** wy_side, double** wz_side, double* tdash,
		double* tcont, double* mudash, double* scont, double** towxx,
		double* re, double** towyy, double** towzz, double** towxy,
		double** towyz, double** towzx, double** qx, double* pr, double** qy,
		double** qz, double* d_nu, double* delt, double** num_dd,
		double* cfl_iele, double* nu, int* jnode, struct node2* eleminfo,
		double* xp_node, struct node1* xycoord, double* yp_node,
		double* zp_node, int* s1, int* s2, int* s3, double* xpstar,
		double* ypstar, double* zpstar, int* m, double* um, double* uxm,
		double* uym, double* uzm, double* utm, double* ut1_p, double* fm,
		double* f1_p, double* fxm, double* fx1_p, double* fym, double* fy1_p,
		double* fzm, double* fz1_p, double* ftm, double* ft1_p, double* gm,
		double* g1_p, double* gxm, double* gx1_p, double* gym, double* gy1_p,
		double* gzm, double* gz1_p, double* gtm, double* gt1_p, double* em,
		double* e1_p, double* exm, double* ex1_p, double* eym, double* ey1_p,
		double* ezm, double* ez1_p, double* etm, double* et1_p, double* fvm,
		double* gvm, double* evm, double* ut2_p, double* f2_p, double* fx2_p,
		double* fy2_p, double* fz2_p, double* ft2_p, double* g2_p,
		double* gx2_p, double* gy2_p, double* gz2_p, double* gt2_p,
		double* e2_p, double* ex2_p, double* ey2_p, double* ez2_p,
		double* et2_p, double* ut3_p, double* f3_p, double* fx3_p,
		double* fy3_p, double* fz3_p, double* ft3_p, double* g3_p,
		double* gx3_p, double* gy3_p, double* gz3_p, double* gt3_p,
		double* e3_p, double* ex3_p, double* ey3_p, double* ez3_p,
		double* et3_p, double* ut4_p, double* f4_p, double* fx4_p,
		double* fy4_p, double* fz4_p, double* ft4_p, double* g4_p,
		double* gx4_p, double* gy4_p, double* gz4_p, double* gt4_p,
		double* e4_p, double* ex4_p, double* ey4_p, double* ez4_p,
		double* et4_p, double* ut5_p, double* f5_p, double* fx5_p,
		double* fy5_p, double* fz5_p, double* ft5_p, double* g5_p,
		double* gx5_p, double* gy5_p, double* gz5_p, double* gt5_p,
		double* e5_p, double* ex5_p, double* ey5_p, double* ez5_p,
		double* et5_p, double* f1_side, double** area_lateral_side1,
		double** nx_lateral_side1, double** ny_lateral_side1,
		double** nz_lateral_side1, double* f2_side, double** area_lateral_side2,
		double** nx_lateral_side2, double** ny_lateral_side2,
		double** nz_lateral_side2, double* f3_side, double** area_lateral_side3,
		double** nx_lateral_side3, double** ny_lateral_side3,
		double** nz_lateral_side3, double* f4_side, double** vol_hexahedra_face,
		double** xc_hexahedra_face, double** yc_hexahedra_face,
		double** zc_hexahedra_face, double* flux, double** xpstar_face,
		double** ypstar_face, double** zpstar_face, double** umdash,
		struct node5* primitive, double* vol_polyhedra, struct node6* derv,
		struct node7* secondary, double* max, int* l)
{

	for (iele = 1; iele <= lele_cnt; iele++)
	{
		c[iele] = pow(
				((ga) * (ga - 1)
						* (u5_p[iele]
								- ((u2_p[iele] * u2_p[iele]
										+ u3_p[iele] * u3_p[iele]
										+ u4_p[iele] * u4_p[iele])
										/ (2.0 * u1_p[iele]))) / (u1_p[iele])),
				0.5);
	}
	*maxprime = 10000;
	for (iele = 1; iele <= lele_cnt; iele++)
	{
		for (*jface = 1; *jface <= 4; *jface++)
		{
			if ((fabs(u2_p[iele]) < 0.000001) && (fabs(u3_p[iele]) < 0.000001)
					&& (fabs(u4_p[iele]) < 0.000001))
			{
				theta[*jface][iele] = 0.0;
			}
			else
			{
				theta[*jface][iele] =
						(u2_p[iele] * (xr[*jface][iele] - sol[iele].spx)
								+ u3_p[iele]
										* (yr[*jface][iele] - sol[iele].spy)
								+ u4_p[iele]
										* (zr[*jface][iele] - sol[iele].spz))
								/ (pow(
										(u2_p[iele] * u2_p[iele]
												+ u3_p[iele] * u3_p[iele]
												+ u4_p[iele] * u4_p[iele]), 0.5)
										* pow(
												((xr[*jface][iele]
														- sol[iele].spx)
														* (xr[*jface][iele]
																- sol[iele].spx)
														+ (yr[*jface][iele]
																- sol[iele].spy)
																* (yr[*jface][iele]
																		- sol[iele].spy)
														+ (zr[*jface][iele]
																- sol[iele].spz)
																* (zr[*jface][iele]
																		- sol[iele].spz)),
												0.5));
			}

		}

	}
	*maximum = 0.0;
	for (iele = 1; iele <= lele_cnt; iele++)
	{
		for (*ineigh = 1; *ineigh <= 4; *ineigh++)
		{
			d1 = lelem_neighbour[*ineigh][iele];
			for (*iside = 1; *iside <= 3; *iside++)
			{
				if (*iside == 1)
				{
					*xc_side = xc_lateral_side1[*ineigh][iele];
					*yc_side = yc_lateral_side1[*ineigh][iele];
					*zc_side = zc_lateral_side1[*ineigh][iele];
				}
				else if (*iside == 2)
				{
					*xc_side = xc_lateral_side2[*ineigh][iele];
					*yc_side = yc_lateral_side2[*ineigh][iele];
					*zc_side = zc_lateral_side2[*ineigh][iele];
				}
				else if (*iside == 3)
				{
					*xc_side = xc_lateral_side3[*ineigh][iele];
					*yc_side = yc_lateral_side3[*ineigh][iele];
					*zc_side = zc_lateral_side3[*ineigh][iele];
				}

				//u1_p=previous time step value
				u1_side[*iside][*ineigh] = u1_p[d1]
						+ ((*xc_side - sol[d1].spx) * ux1_p[d1])
						+ ((*yc_side - sol[d1].spy) * uy1_p[d1])
						+ ((*zc_side - sol[d1].spz) * uz1_p[d1]);
				u2_side[*iside][*ineigh] = u2_p[d1]
						+ ((*xc_side - sol[d1].spx) * ux2_p[d1])
						+ ((*yc_side - sol[d1].spy) * uy2_p[d1])
						+ ((*zc_side - sol[d1].spz) * uz2_p[d1]);
				u3_side[*iside][*ineigh] = u3_p[d1]
						+ ((*xc_side - sol[d1].spx) * ux3_p[d1])
						+ ((*yc_side - sol[d1].spy) * uy3_p[d1])
						+ ((*zc_side - sol[d1].spz) * uz3_p[d1]);
				u4_side[*iside][*ineigh] = u4_p[d1]
						+ ((*xc_side - sol[d1].spx) * ux4_p[d1])
						+ ((*yc_side - sol[d1].spy) * uy4_p[d1])
						+ ((*zc_side - sol[d1].spz) * uz4_p[d1]);
				u5_side[*iside][*ineigh] = u5_p[d1]
						+ ((*xc_side - sol[d1].spx) * ux5_p[d1])
						+ ((*yc_side - sol[d1].spy) * uy5_p[d1])
						+ ((*zc_side - sol[d1].spz) * uz5_p[d1]);
				ux_side[*iside][*ineigh] =
						(ux2_p[d1] / u1_side[*iside][*ineigh])
								- ((u2_side[*iside][*ineigh] * ux1_p[d1])
										/ (u1_side[*iside][*ineigh]
												* u1_side[*iside][*ineigh]));
				uy_side[*iside][*ineigh] =
						(uy2_p[d1] / u1_side[*iside][*ineigh])
								- ((u2_side[*iside][*ineigh] * uy1_p[d1])
										/ (u1_side[*iside][*ineigh]
												* u1_side[*iside][*ineigh]));
				uz_side[*iside][*ineigh] =
						(uz2_p[d1] / u1_side[*iside][*ineigh])
								- ((u2_side[*iside][*ineigh] * uz1_p[d1])
										/ (u1_side[*iside][*ineigh]
												* u1_side[*iside][*ineigh]));
				vx_side[*iside][*ineigh] =
						(ux3_p[d1] / u1_side[*iside][*ineigh])
								- ((u3_side[*iside][*ineigh] * ux1_p[d1])
										/ (u1_side[*iside][*ineigh]
												* u1_side[*iside][*ineigh]));
				vy_side[*iside][*ineigh] =
						(uy3_p[d1] / u1_side[*iside][*ineigh])
								- ((u3_side[*iside][*ineigh] * uy1_p[d1])
										/ (u1_side[*iside][*ineigh]
												* u1_side[*iside][*ineigh]));
				vz_side[*iside][*ineigh] =
						(uz3_p[d1] / u1_side[*iside][*ineigh])
								- ((u3_side[*iside][*ineigh] * uz1_p[d1])
										/ (u1_side[*iside][*ineigh]
												* u1_side[*iside][*ineigh]));
				wx_side[*iside][*ineigh] =
						(ux4_p[d1] / u1_side[*iside][*ineigh])
								- ((u4_side[*iside][*ineigh] * ux1_p[d1])
										/ (u1_side[*iside][*ineigh]
												* u1_side[*iside][*ineigh]));
				wy_side[*iside][*ineigh] =
						(uy4_p[d1] / u1_side[*iside][*ineigh])
								- ((u4_side[*iside][*ineigh] * uy1_p[d1])
										/ (u1_side[*iside][*ineigh]
												* u1_side[*iside][*ineigh]));
				wz_side[*iside][*ineigh] =
						(uz4_p[d1] / u1_side[*iside][*ineigh])
								- ((u4_side[*iside][*ineigh] * uz1_p[d1])
										/ (u1_side[*iside][*ineigh]
												* u1_side[*iside][*ineigh]));
				//effect of viscosity has been taken into account
				*tdash = ((ga - 1.0) / u1_side[*iside][*ineigh])
						* ((u5_side[*iside][*ineigh]
								- (((u2_side[*iside][*ineigh]
										* u2_side[*iside][*ineigh])
										+ (u3_side[*iside][*ineigh]
												* u3_side[*iside][*ineigh])
										+ (u4_side[*iside][*ineigh]
												* u4_side[*iside][*ineigh]))
										/ (2.0 * u1_side[*iside][*ineigh])))
								/ (*tcont));
				*mudash = pow(*tdash, 1.5) * ((1 + *scont) / (*tdash + *scont));
				towxx[*iside][*ineigh] = (*mudash / (3.0 * *re))
						* (4.0 * ux_side[*iside][*ineigh]
								- 2.0 * vy_side[*iside][*ineigh]
								- 2.0 * wz_side[*iside][*ineigh]);
				towyy[*iside][*ineigh] = (*mudash / (3.0 * *re))
						* (4.0 * vy_side[*iside][*ineigh]
								- 2.0 * ux_side[*iside][*ineigh]
								- 2.0 * wz_side[*iside][*ineigh]);
				towzz[*iside][*ineigh] = (*mudash / (3.0 * *re))
						* (4.0 * wz_side[*iside][*ineigh]
								- 2.0 * ux_side[*iside][*ineigh]
								- 2.0 * vy_side[*iside][*ineigh]);
				towxy[*iside][*ineigh] = (*mudash / (*re))
						* (uy_side[*iside][*ineigh] + vx_side[*iside][*ineigh]);
				towyz[*iside][*ineigh] = (*mudash / (*re))
						* (wy_side[*iside][*ineigh] + vz_side[*iside][*ineigh]);
				towzx[*iside][*ineigh] = (*mudash / (*re))
						* (uz_side[*iside][*ineigh] + wx_side[*iside][*ineigh]);
				qx[*iside][*ineigh] = -(ga * *mudash
						/ (*re * *pr * u1_side[*iside][*ineigh]))
						* (ux5_p[d1]
								- (u2_side[*iside][*ineigh]
										* ux_side[*iside][*ineigh])
								- (u3_side[*iside][*ineigh]
										* vx_side[*iside][*ineigh])
								- (u4_side[*iside][*ineigh]
										* wx_side[*iside][*ineigh])
								- (u5_side[*iside][*ineigh] * ux1_p[d1]
										/ u1_side[*iside][*ineigh]));
				qy[*iside][*ineigh] = -(ga * *mudash
						/ (*re * *pr * u1_side[*iside][*ineigh]))
						* (uy5_p[d1]
								- (u2_side[*iside][*ineigh]
										* uy_side[*iside][*ineigh])
								- (u3_side[*iside][*ineigh]
										* vy_side[*iside][*ineigh])
								- (u4_side[*iside][*ineigh]
										* wy_side[*iside][*ineigh])
								- (u5_side[*iside][*ineigh] * uy1_p[d1]
										/ u1_side[*iside][*ineigh]));
				qz[*iside][*ineigh] = -(ga * *mudash
						/ (*re * *pr * u1_side[*iside][*ineigh]))
						* (uz5_p[d1]
								- (u2_side[*iside][*ineigh]
										* uz_side[*iside][*ineigh])
								- (u3_side[*iside][*ineigh]
										* vz_side[*iside][*ineigh])
								- (u4_side[*iside][*ineigh]
										* wz_side[*iside][*ineigh])
								- (u5_side[*iside][*ineigh] * uz1_p[d1]
										/ u1_side[*iside][*ineigh]));
			}
		}

		//determination of nu
		*maxprime = 0.0;
		for (*jface = 1; *jface <= 4; *jface++)
		{
			*d_nu = *delt
					* ((cos(theta[*jface][iele])
							* pow(
									(u2_p[iele] * u2_p[iele]
											+ u3_p[iele] * u3_p[iele]
											+ u4_p[iele] * u4_p[iele]), 0.5)
							* u1_p[iele]) + c[iele]) / (num_dd[*jface][iele]);
			if (*d_nu > *maxprime)
			{
				*maxprime = *d_nu;
			}
		}
		*cfl_iele = *maxprime;
		*cfl_iele = *maxprime / 3.0;
		nu1 = *cfl_iele;
		if (nu1 < 0.3)
		{
			*nu = 1.25;
		}
		else if ((nu1 > 0.3) && (nu1 < 0.8))
		{
			*nu = 1.3;
		}
		else if ((nu1 > 0.8) && (nu1 < 2.0))
		{
			*nu = 1.4;
		}
		else if (nu1 > 2.0)
		{
			printf("try to redefine the mesh\n");
			exit(0);
		}

		for (*jnode = 1; *jnode <= 4; *jnode++)
		{
			d1 = eleminfo[iele].gele_nodes[*jnode];
			xp_node[*jnode] = ((*nu) * (xycoord[d1].xstart))
					+ (sol[iele].spx * (1 - *nu));
			yp_node[*jnode] = ((*nu) * (xycoord[d1].ystart))
					+ (sol[iele].spy * (1 - *nu));
			zp_node[*jnode] = ((*nu) * (xycoord[d1].zstart))
					+ (sol[iele].spz * (1 - *nu));
		}
		//xpstar calculation:Parallel Translation of Tetrahedron
		for (*jnode = 1; *jnode <= 4; *jnode++)
		{
			if (*jnode == 1)
			{
				*s1 = 2;
				*s2 = 3;
				*s3 = 4;
			}
			else if (*jnode == 2)
			{
				*s1 = 1;
				*s2 = 3;
				*s3 = 4;
			}
			else if (*jnode == 3)
			{
				*s1 = 1;
				*s2 = 2;
				*s3 = 4;
			}
			else if (*jnode == 4)
			{
				*s1 = 1;
				*s2 = 2;
				*s3 = 3;
			}

			xpstar[*jnode] = (4 * sol[iele].spx + 3 * xp_node[*jnode]
					- xp_node[*s1] - xp_node[*s2] - xp_node[*s3]) / 4.0;
			ypstar[*jnode] = (4 * sol[iele].spy + 3 * yp_node[*jnode]
					- yp_node[*s1] - yp_node[*s2] - yp_node[*s3]) / 4.0;
			zpstar[*jnode] = (4 * sol[iele].spz + 3 * zp_node[*jnode]
					- zp_node[*s1] - zp_node[*s2] - zp_node[*s3]) / 4.0;
		}
		for (*m = 1; *m <= 5; *m++)
		{
			/*m is used for representing um->u1,u2,u3,u4,u5*/
			for (*jface = 1; *jface <= 4; *jface++)
			{
				//identifying neighbour
				d1 = lelem_neighbour[*jface][iele];
				if (*m == 1)
				{
					*um = u1_p[d1]; //p stands for previous time step
					*uxm = ux1_p[d1];
					*uym = uy1_p[d1];
					*uzm = uz1_p[d1];
					*utm = ut1_p[d1];
					*fm = f1_p[d1];
					*fxm = fx1_p[d1];
					*fym = fy1_p[d1];
					*fzm = fz1_p[d1];
					*ftm = ft1_p[d1];
					*gm = g1_p[d1];
					*gxm = gx1_p[d1];
					*gym = gy1_p[d1];
					*gzm = gz1_p[d1];
					*gtm = gt1_p[d1];
					*em = e1_p[d1];
					*exm = ex1_p[d1];
					*eym = ey1_p[d1];
					*ezm = ez1_p[d1];
					*etm = et1_p[d1];
					for (*iside = 1; *iside <= 3; *iside++)
					{
						fvm[*iside] = 0;
						gvm[*iside] = 0;
						evm[*iside] = 0;
					}
				}
				else if (*m == 2)
				{
					*um = u2_p[d1];
					*uxm = ux2_p[d1];
					*uym = uy2_p[d1];
					*uzm = uz2_p[d1];
					*utm = ut2_p[d1];
					*fm = f2_p[d1];
					*fxm = fx2_p[d1];
					*fym = fy2_p[d1];
					*fzm = fz2_p[d1];
					*ftm = ft2_p[d1];
					*gm = g2_p[d1];
					*gxm = gx2_p[d1];
					*gym = gy2_p[d1];
					*gzm = gz2_p[d1];
					*gtm = gt2_p[d1];
					*em = e2_p[d1];
					*exm = ex2_p[d1];
					*eym = ey2_p[d1];
					*ezm = ez2_p[d1];
					*etm = et2_p[d1];
					for (*iside = 1; *iside <= 3; *iside++)
					{
						fvm[*iside] = towxx[*iside][*jface];
						gvm[*iside] = towxy[*iside][*jface];
						evm[*iside] = towzx[*iside][*jface];
					}
				}
				else if (*m == 3)
				{
					*um = u3_p[d1];
					*uxm = ux3_p[d1];
					*uym = uy3_p[d1];
					*uzm = uz3_p[d1];
					*utm = ut3_p[d1];
					*fm = f3_p[d1];
					*fxm = fx3_p[d1];
					*fym = fy3_p[d1];
					*fzm = fz3_p[d1];
					*ftm = ft3_p[d1];
					*gm = g3_p[d1];
					*gxm = gx3_p[d1];
					*gym = gy3_p[d1];
					*gzm = gz3_p[d1];
					*gtm = gt3_p[d1];
					*em = e3_p[d1];
					*exm = ex3_p[d1];
					*eym = ey3_p[d1];
					*ezm = ez3_p[d1];
					*etm = et3_p[d1];
					for (*iside = 1; *iside <= 3; *iside++)
					{
						fvm[*iside] = towxy[*iside][*jface];
						gvm[*iside] = towyy[*iside][*jface];
						evm[*iside] = towyz[*iside][*jface];
					}
				}
				else if (*m == 4)
				{
					*um = u4_p[d1];
					*uxm = ux4_p[d1];
					*uym = uy4_p[d1];
					*uzm = uz4_p[d1];
					*utm = ut4_p[d1];
					*fm = f4_p[d1];
					*fxm = fx4_p[d1];
					*fym = fy4_p[d1];
					*fzm = fz4_p[d1];
					*ftm = ft4_p[d1];
					*gm = g4_p[d1];
					*gxm = gx4_p[d1];
					*gym = gy4_p[d1];
					*gzm = gz4_p[d1];
					*gtm = gt4_p[d1];
					*em = e4_p[d1];
					*exm = ex4_p[d1];
					*eym = ey4_p[d1];
					*ezm = ez4_p[d1];
					*etm = et4_p[d1];
					for (*iside = 1; *iside <= 3; *iside++)
					{
						fvm[*iside] = towzx[*iside][*jface];
						gvm[*iside] = towyz[*iside][*jface];
						evm[*iside] = towzz[*iside][*jface];
					}
				}
				else if (*m == 5)
				{
					*um = u5_p[d1];
					*uxm = ux5_p[d1];
					*uym = uy5_p[d1];
					*uzm = uz5_p[d1];
					*utm = ut5_p[d1];
					*fm = f5_p[d1];
					*fxm = fx5_p[d1];
					*fym = fy5_p[d1];
					*fzm = fz5_p[d1];
					*ftm = ft5_p[d1];
					*gm = g5_p[d1];
					*gxm = gx5_p[d1];
					*gym = gy5_p[d1];
					*gzm = gz5_p[d1];
					*gtm = gt5_p[d1];
					*em = e5_p[d1];
					*exm = ex5_p[d1];
					*eym = ey5_p[d1];
					*ezm = ez5_p[d1];
					*etm = et5_p[d1];
					for (*iside = 1; *iside <= 3; *iside++)
					{
						fvm[*iside] = ((u2_side[*iside][*jface]
								/ u1_side[*iside][*jface])
								* towxx[*iside][*jface])
								+ ((u3_side[*iside][*jface]
										/ u1_side[*iside][*jface])
										* towxy[*iside][*jface])
								+ ((u4_side[*iside][*jface]
										/ u1_side[*iside][*jface])
										* towzx[*iside][*jface])
								- qx[*iside][*jface];
						gvm[*iside] = ((u2_side[*iside][*jface]
								/ u1_side[*iside][*jface])
								* towxy[*iside][*jface])
								+ ((u3_side[*iside][*jface]
										/ u1_side[*iside][*jface])
										* towyy[*iside][*jface])
								+ ((u4_side[*iside][*jface]
										/ u1_side[*iside][*jface])
										* towyz[*iside][*jface])
								- qy[*iside][*jface];
						evm[*iside] = ((u2_side[*iside][*jface]
								/ u1_side[*iside][*jface])
								* towzx[*iside][*jface])
								+ ((u3_side[*iside][*jface]
										/ u1_side[*iside][*jface])
										* towyz[*iside][*jface])
								+ ((u4_side[*iside][*jface]
										/ u1_side[*iside][*jface])
										* towzz[*iside][*jface])
								- qz[*iside][*jface];
					}
				}

				//f1,f2,f3 r flux calculated for side1,side2 and side3 of a particular neighbour facing "jface"
				*f1_side =
						((*delt / 2.0)
								* (area_lateral_side1[*jface][iele]
										* nx_lateral_side1[*jface][iele])
								* (*fm
										+ ((xc_lateral_side1[*jface][iele]
												- sol[d1].spx) * (*fxm))
										+ ((yc_lateral_side1[*jface][iele]
												- sol[d1].spy) * (*fym))
										+ ((zc_lateral_side1[*jface][iele]
												- sol[d1].spz) * (*fzm))
										+ (*delt * *ftm / 4.0) - fvm[1]))
								+ ((*delt / 2.0)
										* (area_lateral_side1[*jface][iele]
												* ny_lateral_side1[*jface][iele])
										* (*gm
												+ ((xc_lateral_side1[*jface][iele]
														- sol[d1].spx) * (*gxm))
												+ ((yc_lateral_side1[*jface][iele]
														- sol[d1].spy) * (*gym))
												+ ((zc_lateral_side1[*jface][iele]
														- sol[d1].spz) * (*gzm))
												+ (*delt * *gtm / 4.0) - gvm[1]))
								+ ((*delt / 2.0)
										* (area_lateral_side1[*jface][iele]
												* nz_lateral_side1[*jface][iele])
										* (*em
												+ ((xc_lateral_side1[*jface][iele]
														- sol[d1].spx) * (*exm))
												+ ((yc_lateral_side1[*jface][iele]
														- sol[d1].spy) * (*eym))
												+ ((zc_lateral_side1[*jface][iele]
														- sol[d1].spz) * (*ezm))
												+ (*delt * *etm / 4.0) - evm[1]));
				*f2_side =
						((*delt / 2.0)
								* (area_lateral_side2[*jface][iele]
										* nx_lateral_side2[*jface][iele])
								* (*fm
										+ ((xc_lateral_side2[*jface][iele]
												- sol[d1].spx) * (*fxm))
										+ ((yc_lateral_side2[*jface][iele]
												- sol[d1].spy) * (*fym))
										+ ((zc_lateral_side2[*jface][iele]
												- sol[d1].spz) * (*fzm))
										+ (*delt * *ftm / 4.0) - fvm[2]))
								+ ((*delt / 2.0)
										* (area_lateral_side2[*jface][iele]
												* ny_lateral_side2[*jface][iele])
										* (*gm
												+ ((xc_lateral_side2[*jface][iele]
														- sol[d1].spx) * (*gxm))
												+ ((yc_lateral_side2[*jface][iele]
														- sol[d1].spy) * (*gym))
												+ ((zc_lateral_side2[*jface][iele]
														- sol[d1].spz) * (*gzm))
												+ (*delt * *gtm / 4.0) - gvm[2]))
								+ ((*delt / 2.0)
										* (area_lateral_side2[*jface][iele]
												* nz_lateral_side2[*jface][iele])
										* (*em
												+ ((xc_lateral_side2[*jface][iele]
														- sol[d1].spx) * (*exm))
												+ ((yc_lateral_side2[*jface][iele]
														- sol[d1].spy) * (*eym))
												+ ((zc_lateral_side2[*jface][iele]
														- sol[d1].spz) * (*ezm))
												+ (*delt * *etm / 4.0) - evm[2]));
				*f3_side =
						((*delt / 2.0)
								* (area_lateral_side3[*jface][iele]
										* nx_lateral_side3[*jface][iele])
								* (*fm
										+ ((xc_lateral_side3[*jface][iele]
												- sol[d1].spx) * (*fxm))
										+ ((yc_lateral_side3[*jface][iele]
												- sol[d1].spy) * (*fym))
										+ ((zc_lateral_side3[*jface][iele]
												- sol[d1].spz) * (*fzm))
										+ (*delt * *ftm / 4.0) - fvm[3]))
								+ ((*delt / 2.0)
										* (area_lateral_side3[*jface][iele]
												* ny_lateral_side3[*jface][iele])
										* (*gm
												+ ((xc_lateral_side3[*jface][iele]
														- sol[d1].spx) * (*gxm))
												+ ((yc_lateral_side3[*jface][iele]
														- sol[d1].spy) * (*gym))
												+ ((zc_lateral_side3[*jface][iele]
														- sol[d1].spz) * (*gzm))
												+ (*delt * *gtm / 4.0) - gvm[3]))
								+ ((*delt / 2.0)
										* (area_lateral_side3[*jface][iele]
												* nz_lateral_side3[*jface][iele])
										* (*em
												+ ((xc_lateral_side3[*jface][iele]
														- sol[d1].spx) * (*exm))
												+ ((yc_lateral_side3[*jface][iele]
														- sol[d1].spy) * (*eym))
												+ ((zc_lateral_side3[*jface][iele]
														- sol[d1].spz) * (*ezm))
												+ (*delt * *etm / 4.0) - evm[3]));
				*f4_side = -vol_hexahedra_face[*jface][iele]
						* (*um
								+ ((xc_hexahedra_face[*jface][iele]
										- sol[d1].spx) * *uxm)
								+ ((yc_hexahedra_face[*jface][iele]
										- sol[d1].spy) * *uym)
								+ ((zc_hexahedra_face[*jface][iele]
										- sol[d1].spz) * *uzm));
				flux[*jface] = *f1_side + *f2_side + *f3_side + *f4_side;
				//The order in which nodes r arranged for each face: for gambit mesh input it differs
				if (*jface == 1)
				{
					*s1 = 2;
					*s2 = 1;
					*s3 = 3;
				}
				else if (*jface == 2)
				{
					*s1 = 1;
					*s2 = 2;
					*s3 = 4;
				}
				else if (*jface == 3)
				{
					*s1 = 2;
					*s2 = 3;
					*s3 = 4;
				}
				else if (*jface == 4)
				{
					*s1 = 3;
					*s2 = 1;
					*s3 = 4;
				}

				xpstar_face[1][*jface] = xpstar[*s1];
				xpstar_face[2][*jface] = xpstar[*s2];
				xpstar_face[3][*jface] = xpstar[*s3];
				ypstar_face[1][*jface] = ypstar[*s1];
				ypstar_face[2][*jface] = ypstar[*s2];
				ypstar_face[3][*jface] = ypstar[*s3];
				zpstar_face[1][*jface] = zpstar[*s1];
				zpstar_face[2][*jface] = zpstar[*s2];
				zpstar_face[3][*jface] = zpstar[*s3];
				umdash[1][*jface] = *um + ((*delt / 2.0) * *utm)
						+ ((xpstar_face[1][*jface] - sol[d1].spx) * *uxm)
						+ ((ypstar_face[1][*jface] - sol[d1].spy) * *uym)
						+ ((zpstar_face[1][*jface] - sol[d1].spz) * *uzm);
				umdash[2][*jface] = *um + ((*delt / 2.0) * *utm)
						+ ((xpstar_face[2][*jface] - sol[d1].spx) * *uxm)
						+ ((ypstar_face[2][*jface] - sol[d1].spy) * *uym)
						+ ((zpstar_face[2][*jface] - sol[d1].spz) * *uzm);
				umdash[3][*jface] = *um + ((*delt / 2.0) * *utm)
						+ ((xpstar_face[3][*jface] - sol[d1].spx) * *uxm)
						+ ((ypstar_face[3][*jface] - sol[d1].spy) * *uym)
						+ ((zpstar_face[3][*jface] - sol[d1].spz) * *uzm);
			}
			if (*m == 1)
			{
				primitive[iele].u1 = -(flux[1] + flux[2] + flux[3] + flux[4])
						/ vol_polyhedra[iele];
				derv_cal(&umdash[1][1], &umdash[2][1], &umdash[3][1],
						&umdash[1][2], &umdash[2][2], &umdash[3][2],
						&umdash[1][3], &umdash[2][3], &umdash[3][3],
						&umdash[1][4], &umdash[2][4], &umdash[3][4],
						&xpstar_face[1][1], &xpstar_face[2][1],
						&xpstar_face[3][1], &xpstar_face[1][2],
						&xpstar_face[2][2], &xpstar_face[3][2],
						&xpstar_face[1][3], &xpstar_face[2][3],
						&xpstar_face[3][3], &xpstar_face[1][4],
						&xpstar_face[2][4], &xpstar_face[3][4],
						&ypstar_face[1][1], &ypstar_face[2][1],
						&ypstar_face[3][1], &ypstar_face[1][2],
						&ypstar_face[2][2], &ypstar_face[3][2],
						&ypstar_face[1][3], &ypstar_face[2][3],
						&ypstar_face[3][3], &ypstar_face[1][4],
						&ypstar_face[2][4], &ypstar_face[3][4],
						&zpstar_face[1][1], &zpstar_face[2][1],
						&zpstar_face[3][1], &zpstar_face[1][2],
						&zpstar_face[2][2], &zpstar_face[3][2],
						&zpstar_face[1][3], &zpstar_face[2][3],
						&zpstar_face[3][3], &zpstar_face[1][4],
						&zpstar_face[2][4], &zpstar_face[3][4],
						&primitive[iele].u1, &sol[iele].spx, &sol[iele].spy,
						&sol[iele].spz, &derv[iele].ux1, &derv[iele].uy1,
						&derv[iele].uz1, cfl_iele);
			}
			else if (*m == 2)
			{
				primitive[iele].u2 = -(flux[1] + flux[2] + flux[3] + flux[4])
						/ vol_polyhedra[iele];
				derv_cal(&umdash[1][1], &umdash[2][1], &umdash[3][1],
						&umdash[1][2], &umdash[2][2], &umdash[3][2],
						&umdash[1][3], &umdash[2][3], &umdash[3][3],
						&umdash[1][4], &umdash[2][4], &umdash[3][4],
						&xpstar_face[1][1], &xpstar_face[2][1],
						&xpstar_face[3][1], &xpstar_face[1][2],
						&xpstar_face[2][2], &xpstar_face[3][2],
						&xpstar_face[1][3], &xpstar_face[2][3],
						&xpstar_face[3][3], &xpstar_face[1][4],
						&xpstar_face[2][4], &xpstar_face[3][4],
						&ypstar_face[1][1], &ypstar_face[2][1],
						&ypstar_face[3][1], &ypstar_face[1][2],
						&ypstar_face[2][2], &ypstar_face[3][2],
						&ypstar_face[1][3], &ypstar_face[2][3],
						&ypstar_face[3][3], &ypstar_face[1][4],
						&ypstar_face[2][4], &ypstar_face[3][4],
						&zpstar_face[1][1], &zpstar_face[2][1],
						&zpstar_face[3][1], &zpstar_face[1][2],
						&zpstar_face[2][2], &zpstar_face[3][2],
						&zpstar_face[1][3], &zpstar_face[2][3],
						&zpstar_face[3][3], &zpstar_face[1][4],
						&zpstar_face[2][4], &zpstar_face[3][4],
						&primitive[iele].u2, &sol[iele].spx, &sol[iele].spy,
						&sol[iele].spz, &derv[iele].ux2, &derv[iele].uy2,
						&derv[iele].uz2, cfl_iele);
			}
			else if (*m == 3)
			{
				primitive[iele].u3 = -(flux[1] + flux[2] + flux[3] + flux[4])
						/ vol_polyhedra[iele];
				derv_cal(&umdash[1][1], &umdash[2][1], &umdash[3][1],
						&umdash[1][2], &umdash[2][2], &umdash[3][2],
						&umdash[1][3], &umdash[2][3], &umdash[3][3],
						&umdash[1][4], &umdash[2][4], &umdash[3][4],
						&xpstar_face[1][1], &xpstar_face[2][1],
						&xpstar_face[3][1], &xpstar_face[1][2],
						&xpstar_face[2][2], &xpstar_face[3][2],
						&xpstar_face[1][3], &xpstar_face[2][3],
						&xpstar_face[3][3], &xpstar_face[1][4],
						&xpstar_face[2][4], &xpstar_face[3][4],
						&ypstar_face[1][1], &ypstar_face[2][1],
						&ypstar_face[3][1], &ypstar_face[1][2],
						&ypstar_face[2][2], &ypstar_face[3][2],
						&ypstar_face[1][3], &ypstar_face[2][3],
						&ypstar_face[3][3], &ypstar_face[1][4],
						&ypstar_face[2][4], &ypstar_face[3][4],
						&zpstar_face[1][1], &zpstar_face[2][1],
						&zpstar_face[3][1], &zpstar_face[1][2],
						&zpstar_face[2][2], &zpstar_face[3][2],
						&zpstar_face[1][3], &zpstar_face[2][3],
						&zpstar_face[3][3], &zpstar_face[1][4],
						&zpstar_face[2][4], &zpstar_face[3][4],
						&primitive[iele].u3, &sol[iele].spx, &sol[iele].spy,
						&sol[iele].spz, &derv[iele].ux3, &derv[iele].uy3,
						&derv[iele].uz3, cfl_iele);
			}
			else if (*m == 4)
			{
				primitive[iele].u4 = -(flux[1] + flux[2] + flux[3] + flux[4])
						/ vol_polyhedra[iele];
				derv_cal(&umdash[1][1], &umdash[2][1], &umdash[3][1],
						&umdash[1][2], &umdash[2][2], &umdash[3][2],
						&umdash[1][3], &umdash[2][3], &umdash[3][3],
						&umdash[1][4], &umdash[2][4], &umdash[3][4],
						&xpstar_face[1][1], &xpstar_face[2][1],
						&xpstar_face[3][1], &xpstar_face[1][2],
						&xpstar_face[2][2], &xpstar_face[3][2],
						&xpstar_face[1][3], &xpstar_face[2][3],
						&xpstar_face[3][3], &xpstar_face[1][4],
						&xpstar_face[2][4], &xpstar_face[3][4],
						&ypstar_face[1][1], &ypstar_face[2][1],
						&ypstar_face[3][1], &ypstar_face[1][2],
						&ypstar_face[2][2], &ypstar_face[3][2],
						&ypstar_face[1][3], &ypstar_face[2][3],
						&ypstar_face[3][3], &ypstar_face[1][4],
						&ypstar_face[2][4], &ypstar_face[3][4],
						&zpstar_face[1][1], &zpstar_face[2][1],
						&zpstar_face[3][1], &zpstar_face[1][2],
						&zpstar_face[2][2], &zpstar_face[3][2],
						&zpstar_face[1][3], &zpstar_face[2][3],
						&zpstar_face[3][3], &zpstar_face[1][4],
						&zpstar_face[2][4], &zpstar_face[3][4],
						&primitive[iele].u4, &sol[iele].spx, &sol[iele].spy,
						&sol[iele].spz, &derv[iele].ux4, &derv[iele].uy4,
						&derv[iele].uz4, cfl_iele);
			}
			else if (*m == 5)
			{
				primitive[iele].u5 = -(flux[1] + flux[2] + flux[3] + flux[4])
						/ vol_polyhedra[iele];
				derv_cal(&umdash[1][1], &umdash[2][1], &umdash[3][1],
						&umdash[1][2], &umdash[2][2], &umdash[3][2],
						&umdash[1][3], &umdash[2][3], &umdash[3][3],
						&umdash[1][4], &umdash[2][4], &umdash[3][4],
						&xpstar_face[1][1], &xpstar_face[2][1],
						&xpstar_face[3][1], &xpstar_face[1][2],
						&xpstar_face[2][2], &xpstar_face[3][2],
						&xpstar_face[1][3], &xpstar_face[2][3],
						&xpstar_face[3][3], &xpstar_face[1][4],
						&xpstar_face[2][4], &xpstar_face[3][4],
						&ypstar_face[1][1], &ypstar_face[2][1],
						&ypstar_face[3][1], &ypstar_face[1][2],
						&ypstar_face[2][2], &ypstar_face[3][2],
						&ypstar_face[1][3], &ypstar_face[2][3],
						&ypstar_face[3][3], &ypstar_face[1][4],
						&ypstar_face[2][4], &ypstar_face[3][4],
						&zpstar_face[1][1], &zpstar_face[2][1],
						&zpstar_face[3][1], &zpstar_face[1][2],
						&zpstar_face[2][2], &zpstar_face[3][2],
						&zpstar_face[1][3], &zpstar_face[2][3],
						&zpstar_face[3][3], &zpstar_face[1][4],
						&zpstar_face[2][4], &zpstar_face[3][4],
						&primitive[iele].u5, &sol[iele].spx, &sol[iele].spy,
						&sol[iele].spz, &derv[iele].ux5, &derv[iele].uy5,
						&derv[iele].uz5, cfl_iele);
			}
		}
		secondary_variables(&primitive[iele].u1, &derv[iele].ux1,
				&derv[iele].uy1, &derv[iele].uz1, &primitive[iele].u2,
				&derv[iele].ux2, &derv[iele].uy2, &derv[iele].uz2,
				&primitive[iele].u3, &derv[iele].ux3, &derv[iele].uy3,
				&derv[iele].uz3, &primitive[iele].u4, &derv[iele].ux4,
				&derv[iele].uy4, &derv[iele].uz4, &primitive[iele].u5,
				&derv[iele].ux5, &derv[iele].uy5, &derv[iele].uz5,
				&secondary[iele].ut1, &secondary[iele].ut2,
				&secondary[iele].ut3, &secondary[iele].ut4,
				&secondary[iele].ut5, &secondary[iele].f1, &secondary[iele].fx1,
				&secondary[iele].fy1, &secondary[iele].fz1, &secondary[iele].g1,
				&secondary[iele].gx1, &secondary[iele].gy1,
				&secondary[iele].gz1, &secondary[iele].e1, &secondary[iele].ex1,
				&secondary[iele].ey1, &secondary[iele].ez1,
				&secondary[iele].ut1, &secondary[iele].ft1,
				&secondary[iele].gt1, &secondary[iele].et1, &secondary[iele].f2,
				&secondary[iele].fx2, &secondary[iele].fy2,
				&secondary[iele].fz2, &secondary[iele].g2, &secondary[iele].gx2,
				&secondary[iele].gy2, &secondary[iele].gz2, &secondary[iele].e2,
				&secondary[iele].ex2, &secondary[iele].ey2,
				&secondary[iele].ez2, &secondary[iele].ut2,
				&secondary[iele].ft2, &secondary[iele].gt2,
				&secondary[iele].et2, &secondary[iele].f3, &secondary[iele].fx3,
				&secondary[iele].fy3, &secondary[iele].fz3, &secondary[iele].g3,
				&secondary[iele].gx3, &secondary[iele].gy3,
				&secondary[iele].gz3, &secondary[iele].e3, &secondary[iele].ex3,
				&secondary[iele].ey3, &secondary[iele].ez3,
				&secondary[iele].ut3, &secondary[iele].ft3,
				&secondary[iele].gt3, &secondary[iele].et3, &secondary[iele].f4,
				&secondary[iele].fx4, &secondary[iele].fy4,
				&secondary[iele].fz4, &secondary[iele].g4, &secondary[iele].gx4,
				&secondary[iele].gy4, &secondary[iele].gz4, &secondary[iele].e4,
				&secondary[iele].ex4, &secondary[iele].ey4,
				&secondary[iele].ez4, &secondary[iele].ut4,
				&secondary[iele].ft4, &secondary[iele].gt4,
				&secondary[iele].et4, &secondary[iele].f5, &secondary[iele].fx5,
				&secondary[iele].fy5, &secondary[iele].fz5, &secondary[iele].g5,
				&secondary[iele].gx5, &secondary[iele].gy5,
				&secondary[iele].gz5, &secondary[iele].e5, &secondary[iele].ex5,
				&secondary[iele].ey5, &secondary[iele].ez5,
				&secondary[iele].ut5, &secondary[iele].ft5,
				&secondary[iele].gt5, &secondary[iele].et5);
		/*to find the max error*/
		max[1] = primitive[iele].u1 - u1_p[iele];
		max[2] = primitive[iele].u2 - u2_p[iele];
		max[3] = primitive[iele].u3 - u3_p[iele];
		max[4] = primitive[iele].u4 - u4_p[iele];
		max[5] = primitive[iele].u5 - u5_p[iele];
		for (*l = 1; *l <= 5; *l++)
		{
			if ((fabs((*maximum))) < (fabs(max[*l])))
			{
				r1 = iele;
				(*maximum) = fabs(max[*l]);
			}
		}
	}
	return iele;
}

#endif /* SLAVE_SOLVER_H_ */

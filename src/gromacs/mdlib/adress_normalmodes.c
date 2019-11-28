#include "math.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/legacyheaders/types/simple.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/vec.h"
//#include "adress_normalmodes.h"
#include "gromacs/legacyheaders/update.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/legacyheaders/types/state.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/legacyheaders/nrnb.h"

//#define PARTDECOMP(cr)     ((cr)->pd != NULL)

#include "gromacs/math/vec.h"
#include "adress_normalmodes.h"
#include "gromacs/legacyheaders/update.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/legacyheaders/types/state.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/legacyheaders/nrnb.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/random/random.h"

//#define PARTDECOMP(cr)     ((cr)->pd != NULL)

typedef struct {
	double gdt;
	double eph;
	double emh;
	double em;
	double b;
	double c;
	double d;
} gmx_sd_const_t;

typedef struct {
	real V;
	real X;
	real Yv;
	real Yx;
} gmx_sd_sigma_t;


typedef struct { 
    /* BD stuff */
    real           *bd_rf;
    /* SD stuff */
    gmx_sd_const_t *sdc;
    gmx_sd_sigma_t *sdsig;
    rvec           *sd_V;
    int             sd_V_nalloc;
    /* andersen temperature control stuff */
    gmx_bool       *randomize_group;
    real           *boltzfac;
} gmx_stochd_t;

typedef struct gmx_update
{
	gmx_stochd_t *sd;
	/* xprime for constraint algorithms */
	rvec *xp;
	int  xp_nalloc;

	/* Variables for the deform algorithm */
	gmx_int64_t deformref_step;
	matrix     deformref_box;
} t_gmx_update;

static gmx_stochd_t *init_stochd(t_inputrec *ir)
{
    gmx_stochd_t   *sd;
    gmx_sd_const_t *sdc;
    int             ngtc, n;
    real            y;

    snew(sd, 1);

    ngtc = ir->opts.ngtc;

    if (EI_SD(ir->eI))
    {
        snew(sd->sdc, ngtc);
        snew(sd->sdsig, ngtc);

        sdc = sd->sdc;
        for (n = 0; n < ngtc; n++)
        {
            if (ir->opts.tau_t[n] > 0)
            {
                sdc[n].gdt = ir->delta_t/ir->opts.tau_t[n];
                sdc[n].eph = exp(sdc[n].gdt/2);
                sdc[n].emh = exp(-sdc[n].gdt/2);
                sdc[n].em  = exp(-sdc[n].gdt);
            }
            else
            {
                /* No friction and noise on this group */
                sdc[n].gdt = 0;
                sdc[n].eph = 1;
                sdc[n].emh = 1;
                sdc[n].em  = 1;
            }
        }
    }
    return sd;
}

/*gmx_update_t init_update(t_inputrec *ir)
{
    t_gmx_update *upd;

    snew(upd, 1);

    if (ir->eI == eiBD || EI_SD(ir->eI) || ir->etc == etcVRESCALE || ETC_ANDERSEN(ir->etc))
    {
        upd->sd    = init_stochd(ir);
    }

    upd->xp        = NULL;
    upd->xp_nalloc = 0;

    return upd;
}
*/

static rvec *get_xprime(const t_state *state,gmx_update_t upd)
{
	if (state->nalloc > upd->xp_nalloc)
	{
	    upd->xp_nalloc = state->nalloc;
	    srenew(upd->xp,upd->xp_nalloc);
	}
	return upd->xp;
}


void calc_force_on_cg(int cg0, int cg1, t_block * cgs, rvec x[], t_forcerec * fr, t_mdatoms * mdatoms, rvec f[])
{
    int            i, icg,k,k0,k1,d, nrcg;
    atom_id *      cgindex=cgs->index;;
    rvec fsum;

    for (i = 0; i < mdatoms->nr; i++)
    {
        mdatoms->ring_force[i][0] = 0.0;
	mdatoms->ring_force[i][1] = 0.0;
	mdatoms->ring_force[i][2] = 0.0;
    }

    for(icg=cg0; (icg<cg1); icg++)
    {
        k0      = cgindex[icg];
        k1      = cgindex[icg+1];
        nrcg    = k1-k0;
        fsum[0] = 0.0;
        fsum[1] = 0.0;
        fsum[2] = 0.0;

        /* calculate sum all forces on a cg*/
        for(k=k0; (k<k1); k++)
        {
            fsum[0]+=f[k][0];
            fsum[1]+=f[k][1];
	    fsum[2]+=f[k][2];
	}

	    /* store sum in mdatos struct*/
	    for(k=k0; (k<k1); k++)
	    {
	        mdatoms->ring_force[k][0]=fsum[0];
	        mdatoms->ring_force[k][1]=fsum[1];
	        mdatoms->ring_force[k][2]=fsum[2];
	    }
    }
}


static void nm_transform(rvec *x, rvec *u, int P, t_forcerec * fr)
{
    int i, j, d, k;
    for (d = 0; d < DIM; d++)
    {
        for (k = 0; k < P; k++)
       	{
            u[k][d] = 0;

            for (j = 0; j < P; j++) 
	    {
                u[k][d] += fr->adress_NM_M[j][k] * x[j][d];
            }
        }
    }
}

static void inverse_nm_transform(rvec *x, rvec *u, int P, t_forcerec * fr)
{
    int i, j, d, k;
    for (d = 0; d < DIM; d++)
    {
        for (j = 0; j < P; j++)
       	{
            u[j][d] = 0;

            for (k = 0; k < P; k++)
	    {
                u[j][d] += fr->adress_NM_M[j][k] * x[k][d];
	    }
        }
    }
}

static void do_update_sd1_nm(gmx_stochd_t *sd, 
		          int start, int nrend, double dt,
			  rvec accel[], ivec nFreeze[],
			  real invmass[], unsigned short ptype[],
			  unsigned short cFREEZE[], unsigned short cACC[],
			  unsigned short cTC[],
			  rvec x[], rvec xprime[], rvec v[], rvec f[],
			  rvec sd_X[], int ngtc, real ref_t[],
			  int cg0, int cg1, t_block * cgs, t_commrec *cr, 
			  t_forcerec *fr, gmx_ekindata_t *ekind, t_nrnb *nrnb,
			  t_mdatoms *md, t_inputrec *inputrec,
			  gmx_int64_t step, int seed, int* gatindex)
{
    gmx_sd_const_t *sdc;
    gmx_sd_sigma_t *sig;
    real kT;
    int gf = 0, ga = 0, gt = 0;
    real ism, sd_V, fdd, bdd;
    int d;
    atom_id * cgindex=cgs->index;;
    int icg, n0, n1, n, nrcg, k, m, g, a, b, fp, bp;
    rvec * qv; /*velocity in nm space*/
    rvec * qx; /*coordintates in nm space*/
    rvec * fh; /*forces without harmonic ring interactions*/
    rvec * qvnew; 
    t_grpopts *opts = &inputrec->opts;
    t_grp_tcstat *tcstat = ekind->tcstat;
    t_grp_acc *grpstat = ekind->grpstat;
    real dekindl;
    real hm, kk;
    dekindl = 0;
    real omg = (BOLTZ*ref_t[0]*2*M_PI)/PLANCK;

//    printf("cg0 %d, cg1 %d\n", cg0, cg1);
    
    //printf( "do_update_sd1_nm omg %f \n", omg);
    
    snew(qv, fr->n_pi_grps);
    snew(qx, fr->n_pi_grps);
    snew(fh, fr->n_pi_grps);
    snew(qvnew, fr->n_pi_grps);
    sdc = sd->sdc;
    sig = sd->sdsig;

/*	printf("hibefore\n");
	printf("sdc.em = %f\n", sdc[0].em);
	printf("hiafter\n");
*/

    if (nrend > sd->sd_V_nalloc)
    {
        sd->sd_V_nalloc = over_alloc_dd(nrend);
    	srenew(sd->sd_V,sd->sd_V_nalloc);
    }

    for(n = 0; n < ngtc; n++)
    {
        //printf("n = %d\n", n);
    	kT = BOLTZ*ref_t[n];
    	/* The mass is accounted for later, since this differs per atom */
    	sig[n].V  = sqrt(kT*(1 - sdc[n].em * sdc[n].em));
    }

    /* for kinetic energy*/
    for (g = 0; (g < opts->ngtc); g++)
    {
        copy_mat(tcstat[g].ekinh,tcstat[g].ekinh_old);
    	clear_mat(tcstat[g].ekinh);
    }


    /******************NEW IMPLEMENTATION****************/
    //printf("cg0 %d cg1 %d n1 %d n2 %d \n", cg0, cg1, cgindex[cg0], cgindex[cg1]);
    for (icg = cg0; (icg < cg1); icg++) 
    {
        n0 = cgindex[icg];
    	n1 = cgindex[icg + 1];
//	printf("n0 %d, n1 %d\n", n0, n1);
    	/* loop over all atoms in cg*/
        nrcg = n1 - n0;
        if (nrcg == 1)
       	{
            for (d = 0; d < DIM; d++)
	    {
                v[n0][d] = 0.0;
                xprime[n0][d] = x[n0][d];
            }
            continue;
        }
    
        for (n = n0; (n < n1); n++) 
        {
    
    	    if (cFREEZE)
	    {
                gf = cFREEZE[n];
    	    }
    	    if (cACC)
	    {
    		ga = cACC[n];
    	    }
    	    if (cTC)
	    {
    		gt = cTC[n];
    	    }
    
    	    real rnd[3];
    	    int ng = gatindex ? gatindex[n] : n;
    	    ism = sqrt(invmass[n]);
            kk = md->massT[n] * omg * omg;                                              
    //      printf("kk %f \n", kk);
    
    	    fp = n + 1;
            if (fp == n1)
	    {
                fp = n0;
	    }
    
    	    bp = n - 1;
            if (bp < n0)
	    {
                bp = n1 - 1;
	    }
    
//    	    printf("n0 %d n1 %d n %d fp %d bp %d \n", n0, n1, n, fp, bp);
    		
    	    gmx_rng_cycle_3gaussian_table(step, ng, seed, RND_SEED_UPDATE, rnd);
    
    	    /* velocity update step in real space*/
    	    for (d = 0; d < DIM; d++)
    	    {
    	        if ((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d])
    		{
             	    real vn;
    
                    fdd = -kk * (x[n][d] - x[fp][d]);
                    bdd = -kk * (x[bp][d] - x[n][d]);
                    fh[n-n0][d] = f[n][d] - fdd + bdd;
                    // printf("n %d dim %d force %f \n", n, d, fh[n-n0][d]);
                    sd_V = ism*sig[gt].V*rnd[d];
    		    vn   = v[n][d] + (invmass[n]*fh[n-n0][d] + accel[ga][d])*dt;
    		    /* sdc.em is defined as exp(-delta_t/tau) in update.c */ // ???
    		    v[n][d] = vn*sdc[gt].em + sd_V;   
    		}
                else
                {
                    v[n][d] = 0.0; 
                } 
    	    }
        }   
        
   	//   printf("before %d %f %f %f \n", n0, v[n0][0], v[n0][1], v[n0][2]);
    	/* transform velocities and positions to nm space*/
   
	nm_transform(&v[n0], qv, fr->n_pi_grps, fr);
	nm_transform(&x[n0], qx, fr->n_pi_grps, fr);
    	
	/* evolution of a free ring polymer*/
    	for (n = n0; n < n1; n++)
       	{
            k = n - n0;
            hm = 0.5 * md->massT[n];
          
   //       real omg = (fr->n_pi_grps*BOLTZ*ref_t[0]*2*M_PI*md->massT[n])/PLANCK;
           
	    if (cFREEZE)
	    {
                gf = cFREEZE[n];
            }
            if (cACC)
	    {
                ga = cACC[n];
            }
            if (cTC)
	    {
                gt = cTC[n];
            }
    
            for (d = 0; d < DIM; d++) 
	    {
    	        if ((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze[gf][d])
	       	{
                    if (k == 0)
                    {
                        qvnew[k][d] = qv[k][d];
                        qx[k][d] = qv[k][d]*dt + qx[k][d]; 
                    }
                    else
                    {  
                        qvnew[k][d] = cos(fr->adress_NM_mu[k]*omg*dt) * qv[k][d] - omg * fr->adress_NM_mu[k] * 
				      sin(omg * fr->adress_NM_mu[k]*dt) *qx[k][d];
                        
			qx[k][d] = (1/(fr->adress_NM_mu[k]*omg)) * sin(fr->adress_NM_mu[k] * omg * dt) *
                                   qv[k][d] + cos(fr->adress_NM_mu[k] * omg * dt) * qx[k][d];
                    }
                    
		    qv[k][d] = qvnew[k][d];                                      
    		} 
    	    }
                   
                /*    for (d = 0; (d < DIM); d++) {
                            for (m = 0; (m < DIM); m++) {
                                    tcstat[gt].ekinh[m][d] += hm * qv[k][m] * qv[k][d];// + hm * fr->adress_NM_mu[k]* omg * fr->adress_NM_mu[k] * omg * qx[k][m] * qx[k][d]);
                            }
                    }*/
                  
              
        }
    
    /* transform velocities and coordinates back to real space*/
    inverse_nm_transform(qv, &v[n0], fr->n_pi_grps, fr);
    inverse_nm_transform(qx, &xprime[n0], fr->n_pi_grps, fr);
               
 // printf("after %d %f %f %f \n", n0, v[n0][0], v[n0][1], v[n0][2]);
    
    } 
    
    
    ekind->dekindl = dekindl;
    inc_nrnb(nrnb, eNR_EKIN, nrend);
    
    sfree(qv);
    sfree(qvnew);
    sfree(qx);
    sfree(fh);
}

void update_coords_nm(
		  gmx_int64_t  step,
		  t_inputrec   *inputrec,      /* input record and box stuff  */
		  t_mdatoms    *md,
		  t_state      *state,
		  rvec         *f,        /* forces on home particles */
		  gmx_update_t upd,
		  gmx_localtop_t *top,
		  t_commrec    *cr,
		  t_forcerec   *fr,
		  gmx_ekindata_t *ekind,
		  t_nrnb *nrnb)
{
    double           dt, alpha;
    rvec             *force;
    real             dt_1;
    int              start, homenr, nrend, i, j, d, n, m, g;
    int              cg0, cg1, th, nth;
    rvec             *xprime;
    
    start  = md->start;
    homenr = md->homenr;
    nrend = start+homenr;
    
    dt   = inputrec->delta_t;
    dt_1 = 1.0/dt;
    
    force = f;
    
//    printf("update_coords_nm\n");
//    printf("cr nodeid %d\n", cr->nodeid);
    
    cg0 = 0;
    if (DOMAINDECOMP(cr))
    {
        cg1 = cr->dd->ncg_home;
    }
    else
    {
        cg1 = top->cgs.nr;
    }
     
    
    // printf("cg0 %d cg1 %d \n", cg0, cg1);
     xprime = get_xprime(state,upd);
    
     nth = gmx_omp_nthreads_get(emntUpdate);
	 
#pragma omp parallel for num_threads(nth) schedule(static) private(alpha)
         for (th = 0; th < nth; th++)  
         {
         int start_th, end_th;
         start_th = start + ((nrend-start)* th   )/nth;
         end_th   = start + ((nrend-start)*(th+1))/nth;
	 do_update_sd1_nm(upd->sd,
			  start_th, end_th, dt,
			  inputrec->opts.acc, inputrec->opts.nFreeze,
			  md->invmass, md->ptype,
			  md->cFREEZE, md->cACC, md->cTC,
			  state->x, xprime, state->v, force,
			  state->sd_X,
			  inputrec->opts.ngtc, inputrec->opts.ref_t,
			  cg0, cg1, &(top->cgs), cr, fr, ekind, nrnb, 
			  md, inputrec, step, inputrec->ld_seed,
			  DOMAINDECOMP(cr) ? cr->dd->gatindex : NULL);
         }
}

void InitNMMatrix(int P, t_forcerec *fr)
{
    int i, j, n, k;

    real invP= 1.0/P;
    real twopi=M_PI*2.0;
//   real omg = P*BOLTZ*ref_[t]*twopi/PLANK; 
    snew(fr->adress_NM_M, P);

    for(i=0; i< P;i++)
    {
        snew(fr->adress_NM_M[i], P);
    }

	for (j = 0; j < P; j++)
       	{
            for (k = 0; k < P; k++)
	    {
                        
                if (k == 0)
                {
                    fr->adress_NM_M[j][k] = sqrt(invP);
                }
                else if (k >= 1 && k <= (P/2 - 1))
                {
                    fr->adress_NM_M[j][k] = sqrt(2*invP)*cos(twopi*k*(j+1)*invP);
                }   
                else if (k == P/2)
                { 
                    fr->adress_NM_M[j][k] = sqrt(invP)*pow(-1, j+1); 
                }
                else 
                {
                    fr->adress_NM_M[j][k] = sqrt(2*invP)*sin(twopi*k*(j+1)*invP);
                }
		//	fr->adress_NM_M[i][j]=(cos(twopi* i * j * invP) - sin(twopi*  i * j * invP))*sqrt(invP);
                    printf("j = %d k = %d fr->adress_NM_M[j][k] = %f \n", j, k, fr->adress_NM_M[j][k]); 
	    }
	}
    
    snew(fr->adress_NM_mu, P);

    for (i = 0; i < P; i++) 
    {
	fr->adress_NM_mu[i]=2.0*(sin(M_PI*i/P));
 //       printf("fr->adress_NM_mu[i] = %f time = %f \n", fr->adress_NM_mu[i], 2*M_PI/(fr->adress_NM_mu[i]));
    }

}

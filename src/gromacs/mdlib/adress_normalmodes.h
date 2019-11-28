#ifndef _adress_normalmodes_h_
#define _adress_normalmodes_h_

#include "math.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/legacyheaders/types/simple.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/vec.h"
//#include "adress_normalmodes.h"
#include "gromacs/legacyheaders/update.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/legacyheaders/types/state.h"
//#include "partdec.h"
#include "gromacs/legacyheaders/types/commrec.h"

#ifdef __cplusplus
extern "C" {
#endif

void calc_force_on_cg (int                  cg0,
		         int                  cg1,
			 t_block *            cgs,
			 rvec                 x[],
			 t_forcerec *         fr,
			 t_mdatoms *          mdatoms,
			 rvec       f[]);

void update_coords_nm(
		  gmx_int64_t  step,
		  t_inputrec   *inputrec,      /* input record and box stuff  */
		  t_mdatoms    *md,
		  t_state      *state,
		  rvec         *f,        /* forces on home particles */
		  gmx_update_t upd,
		  gmx_localtop_t *top,
		  t_commrec    *cr,
		  t_forcerec  *fr,
		  gmx_ekindata_t *ekind,
		  t_nrnb *nrnb);
void InitNMMatrix(int P, t_forcerec *fr);

#ifdef __cplusplus
}
#endif

#endif

/****************************************************************
Copyright (C) 1997, 2001 Lucent Technologies
All Rights Reserved

Permission to use, copy, modify, and distribute this software and
its documentation for any purpose and without fee is hereby
granted, provided that the above copyright notice appear in all
copies and that both that the copyright notice and this
permission notice and warranty disclaimer appear in supporting
documentation, and that the name of Lucent or any of its entities
not be used in advertising or publicity pertaining to
distribution of the software without specific, written prior
permission.

LUCENT DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
IN NO EVENT SHALL LUCENT OR ANY OF ITS ENTITIES BE LIABLE FOR ANY
SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.
****************************************************************/

#include "asl.h"

/* Fortran-callable routines: */
/* objval_, objgrd_, conval_, jacval_, hvcomp_, cnival_, congrd_ */

 static void
#ifdef KR_headers
bad_N(asl, N, who) ASL *asl; fint *N; char *who;
#else
bad_N(ASL *asl, fint *N, char *who)
#endif
{
	what_prog();
	fprintf(Stderr, "%s: got N = %ld; expected %d\n", who, (long)*N, n_var);
	exit(1);
	}

 real
#ifdef KR_headers
objval_(N, X, NOBJ, nerror) fint *N; real *X; fint *NOBJ; fint *nerror;
#else
objval_(fint *N, real *X, fint *NOBJ, fint *nerror)
#endif
{
	ASL *a;
	static char who[] = "objval_";
	if (!(a = cur_ASL))
		badasl_ASL(a,0,who);
	if (a->i.n_var_ != *N)
		bad_N(a, N, who);
	return (*a->p.Objval)(a, (int)*NOBJ, X, nerror);
	}

 void
#ifdef KR_headers
objgrd_(N, X, NOBJ, G, nerror) fint *N; real *X; fint *NOBJ; real *G; fint *nerror;
#else
objgrd_(fint *N, real *X, fint *NOBJ, real *G, fint *nerror)
#endif
{
	ASL *a;
	static char who[] = "objgrd_";
	if (!(a = cur_ASL))
		badasl_ASL(a, 0, who);
	if (a->i.n_var_ != *N)
		bad_N(a, N, who);
	(*a->p.Objgrd)(a, (int)*NOBJ, X, G, nerror);
	}

 void
#ifdef KR_headers
conval_(M, N, X, F, nerror) fint *M; fint *N; real *X; real *F; fint *nerror;
#else
conval_(fint *M, fint *N, real *X, real *F, fint *nerror)
#endif
{
	ASL *a;
	static char who[] = "conval_";
	if (!(a = cur_ASL))
		badasl_ASL(a, 0, who);
	if (*M != a->i.n_con_ || *N != a->i.n_var_) {
		what_prog();
		fprintf(Stderr, "%s: got M = %ld, N = %ld; expected %d, %d\n",
			who, *M, *N, a->i.n_con_, a->i.n_var_);
		exit(1);
		}
	(*a->p.Conval)(a, X, F, nerror);
	}

 void
#ifdef KR_headers
jacval_(M, N, NZ, X, JAC, nerror) fint *M; fint *N; fint *NZ; real *X; real *JAC; fint *nerror;
#else
jacval_(fint *M, fint *N, fint *NZ, real *X, real *JAC, fint *nerror)
#endif
{
	ASL *a = cur_ASL;
	mnnzchk_ASL(a, M, N, NZ, "jacval_");
	(*a->p.Jacval)(a, X, JAC, nerror);
	}

 void
#ifdef KR_headers
hvcomp_(hv, p, nobj, ow, y) real *hv, *p, *ow, *y; fint *nobj;
#else
hvcomp_(real *hv, real *p, fint *nobj, real *ow, real *y)
#endif
{
	ASL *a;
	if (!(a = cur_ASL))
		badasl_ASL(a,0,"objval");
	(*a->p.Hvcomp)(a, hv, p, (int)*nobj, ow, y);
	}

 void
#ifdef KR_headers
hvinit_(nobj, ow, y) fint *nobj; real *ow, *y;
#else
hvinit_(fint *nobj, real *ow, real *y)
#endif
{
	ASL *a;
	if (!(a = cur_ASL))
		badasl_ASL(a,0,"hvinit");
	(*a->p.Hvinit)(a, a->p.ihd_limit_, (int)*nobj, ow, y);
	}

 static ASL *
#ifdef KR_headers
NI_check(I, N, who) fint *I, *N; char *who;
#else
NI_check(fint *I, fint *N, char *who)
#endif
{
	ASL *asl;
	fint i, m;

	if (!(asl = cur_ASL))
		badasl_ASL(asl,0,who);
	if (*N != n_var)
		bad_N(asl, N, who);
	i = *I;
	m = n_con;
	if (i < 1 || i > m) {
		what_prog();
		fprintf(Stderr, "%s: got I = %ld; expected 1 <= I <= %ld\n",
			who, (long)i, (long)m);
		exit(1);
		}
	return asl;
	}

 void
#ifdef KR_headers
congrd_(N, I, X, G, nerror) fint *N, *I, *nerror; real *X, *G;
#else
congrd_(fint *N, fint *I, real *X, real *G, fint *nerror)
#endif
{
	ASL *a = NI_check(I, N, "congrd_");
	(*a->p.Congrd)(a, (int)*I, X, G, nerror);
	}

 real
#ifdef KR_headers
cnival_(N, I, X, nerror) fint *N; fint *I; real *X; fint *nerror;
#else
cnival_(fint *N, fint *I, real *X, fint *nerror)
#endif
{
	ASL *a = NI_check(I, N, "cnival_");
	return (*a->p.Conival)(a, (int)*I, X, nerror);
	}

 void
delprb_(VOID)
{
	ASL_free(&cur_ASL);
	}
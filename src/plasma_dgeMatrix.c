#include "plasma_dgeMatrix.h"
#include "P.h"
#include "plasma_init.h"
#include "hiplar_dbg.h"

static
double get_norm(SEXP obj, const char *typstr)
{
#ifdef HIPLAR_WITH_PLASMA
    if(any_NA_in_x(obj)) {
    	return NA_REAL;
    } else {
    	char typnm[] = {'\0', '\0'};
    	int *dims = INTEGER(GET_SLOT(obj, Matrix_DimSym));
    	double *work = (double *) NULL;

    	typnm[0] = La_norm_type(typstr);

        if ((*typnm == 'F') && (CHECK_VVERSION_BEQ(2,4,5))) {
#ifdef HIBLAR_DBG
		  R_ShowMessage("not implemented");
#endif
		return F77_CALL(dlange)(typnm, dims, dims+1,
					      REAL(GET_SLOT(obj, Matrix_xSym)),
						dims, work);
	  }

    	if (*typnm == 'F') {
    	  	work = (double *) R_alloc(2*R_PLASMA_NUM_THREADS, sizeof(double));
    	} else {
    		work = (double *) R_alloc(R_PLASMA_NUM_THREADS, sizeof(double));
        }
        return P_dlange(typnm, dims[0], dims[1], REAL(GET_SLOT(obj, Matrix_xSym)), dims[0], work);

    }
#endif
	return 0.0;
}

SEXP plasma_dgeMatrix_norm(SEXP obj, SEXP type)
{
#ifdef HIPLAR_WITH_PLASMA

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dgeMatrix_norm");
#endif

    return ScalarReal(get_norm(obj, CHAR(asChar(type))));
#endif
    return R_NilValue;
}

SEXP plasma_dgeMatrix_rcond(SEXP obj, SEXP type)
{
#ifdef HIPLAR_WITH_PLASMA

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dgeMatrix_rcond");
#endif

    SEXP LU = PROTECT(plasma_dgeMatrix_LU_(obj, FALSE));/* <- not warning about singularity */
    char typnm[] = {'\0', '\0'};
    int *dims = INTEGER(GET_SLOT(LU, Matrix_DimSym)), info;
    double anorm, rcond;

    if (dims[0] != dims[1] || dims[0] < 1) {
	UNPROTECT(1);
	error(_("rcond requires a square, non-empty matrix"));
    }
    typnm[0] = La_rcond_type(CHAR(asChar(type)));
    anorm = get_norm(obj, typnm);
    F77_CALL(dgecon)(typnm,
		     dims, REAL(GET_SLOT(LU, Matrix_xSym)),
		     dims, &anorm, &rcond,
		     (double *) R_alloc(4*dims[0], sizeof(double)),
		     (int *) R_alloc(dims[0], sizeof(int)), &info);
    UNPROTECT(1);
    return ScalarReal(rcond);
#endif
    return R_NilValue;
}

SEXP plasma_dgeMatrix_crossprod(SEXP x, SEXP trans)
{
#ifdef HIPLAR_WITH_PLASMA
    int tr = asLogical(trans);/* trans=TRUE: tcrossprod(x) */
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dpoMatrix"))),
	nms = VECTOR_ELT(GET_SLOT(x, Matrix_DimNamesSym), tr ? 0 : 1),
	vDnms = ALLOC_SLOT(val, Matrix_DimNamesSym, VECSXP, 2);
    int *Dims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	*vDims = INTEGER(ALLOC_SLOT(val, Matrix_DimSym, INTSXP, 2));
    int k = tr ? Dims[1] : Dims[0], n = tr ? Dims[0] : Dims[1];
    double *vx = REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, n * n)),
	one = 1.0, zero = 0.0;
    int info;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dgeMatrix_crossprod");
#endif

    AZERO(vx, n * n);
    SET_SLOT(val, Matrix_uploSym, mkString("U"));
    ALLOC_SLOT(val, Matrix_factorSym, VECSXP, 0);
    vDims[0] = vDims[1] = n;
    SET_VECTOR_ELT(vDnms, 0, duplicate(nms));
    SET_VECTOR_ELT(vDnms, 1, duplicate(nms));
    if(n) {

    info = P_dsyrk("U", tr ? "N" : "T", n, k,
		    one, REAL(GET_SLOT(x, Matrix_xSym)), Dims[0],
		    zero, vx, n);
   	if (info) {
        error(_("PLASMA routine %s returned error code %d"), "PLASMA_dsyrk", info);
    }
    }

    SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));
    UNPROTECT(1);
    return val;
#endif
    return R_NilValue;
}

SEXP plasma_dgeMatrix_dgeMatrix_crossprod(SEXP x, SEXP y, SEXP trans)
{
#ifdef HIPLAR_WITH_PLASMA
    int tr = asLogical(trans);/* trans=TRUE: tcrossprod(x,y) */
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
    int *xDims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	*yDims = INTEGER(GET_SLOT(y, Matrix_DimSym)),
	*vDims;
    int m  = xDims[!tr],  n = yDims[!tr];/* -> result dim */
    int xd = xDims[ tr], yd = yDims[ tr];/* the conformable dims */
    double one = 1.0, zero = 0.0;
    int info;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dgeMatrix_dgeMatrix_crossprod");
#endif

    SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));
    SET_SLOT(val, Matrix_DimSym, allocVector(INTSXP, 2));
    vDims = INTEGER(GET_SLOT(val, Matrix_DimSym));
    if (xd > 0 && yd > 0 && n > 0 && m > 0) {
	if (xd != yd)
	    error(_("Dimensions of x and y are not compatible for %s"),
		  tr ? "tcrossprod" : "crossprod");
	vDims[0] = m; vDims[1] = n;
	SET_SLOT(val, Matrix_xSym, allocVector(REALSXP, m * n));

    info = P_dgemm(tr ? "N" : "T", tr ? "T" : "N", m, n, xd, one,
			REAL(GET_SLOT(x, Matrix_xSym)), xDims[0],
			REAL(GET_SLOT(y, Matrix_xSym)), yDims[0],
			zero, REAL(GET_SLOT(val, Matrix_xSym)), m);
   	if (info) {
        error(_("PLASMA routine %s returned error code %d"), "PLASMA_dgemm", info);
    }

    }
    UNPROTECT(1);
    return val;
#endif
    return R_NilValue;
}

SEXP plasma_dgeMatrix_matrix_crossprod(SEXP x, SEXP y, SEXP trans)
{
#ifdef HIPLAR_WITH_PLASMA
    int tr = asLogical(trans);/* trans=TRUE: tcrossprod(x,y) */
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
    int *xDims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	*yDims = INTEGER(getAttrib(y, R_DimSymbol)),
	*vDims, nprot = 1;
    int m  = xDims[!tr],  n = yDims[!tr];/* -> result dim */
    int xd = xDims[ tr], yd = yDims[ tr];/* the conformable dims */
    double one = 1.0, zero = 0.0;
    int info;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dgeMatrix_matrix_crossprod");
#endif

    if (isInteger(y)) {
	y = PROTECT(coerceVector(y, REALSXP));
	nprot++;
    }
    if (!(isMatrix(y) && isReal(y)))
	error(_("Argument y must be a numeric matrix"));
    SET_SLOT(val, Matrix_factorSym, allocVector(VECSXP, 0));
    SET_SLOT(val, Matrix_DimSym, allocVector(INTSXP, 2));
    vDims = INTEGER(GET_SLOT(val, Matrix_DimSym));
    if (xd > 0 && yd > 0 && n > 0 && m > 0) {
	if (xd != yd)
	    error(_("Dimensions of x and y are not compatible for %s"),
		  tr ? "tcrossprod" : "crossprod");
	vDims[0] = m; vDims[1] = n;
	SET_SLOT(val, Matrix_xSym, allocVector(REALSXP, m * n));

    info = P_dgemm(tr ? "N" : "T", tr ? "T" : "N", m, n, xd, one,
			REAL(GET_SLOT(x, Matrix_xSym)), xDims[0],
			REAL(y), yDims[0],
			zero, REAL(GET_SLOT(val, Matrix_xSym)), m);
   	if (info) {
        error(_("PLASMA routine %s returned error code %d"), "PLASMA_dgemm", info);
    }

    }
    UNPROTECT(nprot);
    return val;
#endif
    return R_NilValue;
}


SEXP plasma_dgeMatrix_LU_(SEXP x, Rboolean warn_sing)
{
#ifdef HIPLAR_WITH_PLASMA
    SEXP val = get_factors(x, "LU");
    int *dims, npiv, info;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dgeMatrix_LU_");
#endif

    if (val != R_NilValue) { /* nothing to do if it's there in 'factors' slot */
#ifdef HIPLAR_DBG
R_ShowMessage("DBG: plasma_dgeMatrix_LU_ using cached value");
#endif
    	return val;
    }

    dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    if (dims[0] < 1 || dims[1] < 1)
	error(_("Cannot factor a matrix with zero extents"));
    npiv = (dims[0] < dims[1]) ? dims[0] : dims[1];
    val = PROTECT(NEW_OBJECT(MAKE_CLASS("denseLU")));
    slot_dup(val, x, Matrix_xSym);
    slot_dup(val, x, Matrix_DimSym);

    info = P_dgetrf(dims[0], dims[1], REAL(GET_SLOT(val, Matrix_xSym)), dims[0],
        INTEGER(ALLOC_SLOT(val, Matrix_permSym, INTSXP, npiv)));

    if (info < 0)
	error(_("PLASMA routine %s returned error code %d"), "PLASMA_dgetrf", info);
    else if (info > 0 && warn_sing)
	warning(_("Exact singularity detected during LU decomposition: %s, i=%d."),
		"U[i,i]=0", info);

    UNPROTECT(1);
    return set_factors(x, val, "LU");
#endif
    return R_NilValue;
}

SEXP plasma_dgeMatrix_LU(SEXP x, SEXP warn_singularity)
{
#ifdef HIPLAR_WITH_PLASMA

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dgeMatrix_LU");
#endif

    return plasma_dgeMatrix_LU_(x, asLogical(warn_singularity));
#endif
    return R_NilValue;
}

SEXP plasma_dgeMatrix_determinant(SEXP x, SEXP logarithm)
{
#ifdef HIPLAR_WITH_PLASMA

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dgeMatrix_determinant");
#endif

    int lg = asLogical(logarithm);
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym)),
	n = dims[0], sign = 1;
    double modulus = lg ? 0. : 1; /* initialize; = result for n == 0 */

    if (n != dims[1])
	error(_("Determinant requires a square matrix"));
    if (n > 0) {
	SEXP lu = plasma_dgeMatrix_LU_(x, /* do not warn about singular LU: */ FALSE);
	int i, *jpvt = INTEGER(GET_SLOT(lu, Matrix_permSym));
	double *luvals = REAL(GET_SLOT(lu, Matrix_xSym));

	for (i = 0; i < n; i++) if (jpvt[i] != (i + 1)) sign = -sign;
	if (lg) {
	    for (i = 0; i < n; i++) {
		double dii = luvals[i*(n + 1)]; /* ith diagonal element */
		modulus += log(dii < 0 ? -dii : dii);
		if (dii < 0) sign = -sign;
	    }
	} else {
	    for (i = 0; i < n; i++)
		modulus *= luvals[i*(n + 1)];
	    if (modulus < 0) {
		modulus = -modulus;
		sign = -sign;
	    }
	}
    }
    return as_det_obj(modulus, lg, sign);
#endif
    return R_NilValue;
}

SEXP plasma_dgeMatrix_solve(SEXP a)
{
#ifdef HIPLAR_WITH_PLASMA

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dgeMatrix_solve");
#endif

    /*  compute the 1-norm of the matrix, which is needed
	later for the computation of the reciprocal condition number. */
    double aNorm = get_norm(a, "1");

    /* the LU decomposition : */
    SEXP val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix"))),
	lu = plasma_dgeMatrix_LU_(a, TRUE);
    int *dims = INTEGER(GET_SLOT(lu, Matrix_DimSym)),
	*pivot = INTEGER(GET_SLOT(lu, Matrix_permSym));

    /* prepare variables for the dgetri calls */
    double *x;//, tmp;
    int	info;//, lwork = -1;


    if (dims[0] != dims[1]) error(_("Solve requires a square matrix"));
    slot_dup(val, lu, Matrix_xSym);
    x = REAL(GET_SLOT(val, Matrix_xSym));
    slot_dup(val, lu, Matrix_DimSym);

    if(dims[0]) /* the dimension is not zero */
    {
        /* is the matrix is *computationally* singular ? */
        double rcond;
        F77_CALL(dgecon)("1", dims, x, dims, &aNorm, &rcond,
                         (double *) R_alloc(4*dims[0], sizeof(double)),
                         (int *) R_alloc(dims[0], sizeof(int)), &info);
        if (info)
            error(_("error [%d] from Lapack 'dgecon()'"), info);
        if(rcond < DOUBLE_EPS)
            error(_("Lapack dgecon(): system computationally singular, reciprocal condition number = %g"),
		  rcond);

        /* only now try the inversion and check if the matrix is *exactly* singular: */
        info = P_dgetri(dims[0], x, dims[0], pivot);
        if (info < 0) {
            error(_("PLASMA routine %s returned error code %d"), "PLASMA_dgetri", info);
        }
    	if (info > 0) {
    	    error(_("PLASMA routine PLASMA_dgetri: system is exactly singular"));
        }
    }
    UNPROTECT(1);
    return val;
#endif
    return R_NilValue;
}

SEXP plasma_dgeMatrix_matrix_solve(SEXP a, SEXP b)
{
#ifdef HIPLAR_WITH_PLASMA

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dgeMatrix_matrix_solve");
#endif

    SEXP val = PROTECT(dup_mMatrix_as_dgeMatrix(b)),
	lu = PROTECT(plasma_dgeMatrix_LU_(a, TRUE));
    int *adims = INTEGER(GET_SLOT(lu, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(val, Matrix_DimSym));
    int info, n = bdims[0], nrhs = bdims[1];

    if (*adims != *bdims || bdims[1] < 1 || *adims < 1 || *adims != adims[1])
	error(_("Dimensions of system to be solved are inconsistent"));

    info = P_dgetrs("N", n, nrhs, REAL(GET_SLOT(lu, Matrix_xSym)), n,
        INTEGER(GET_SLOT(lu, Matrix_permSym)), REAL(GET_SLOT(val, Matrix_xSym)), n);
    if (info < 0) {
        error(_("PLASMA routine %s returned error code %d"), "PLASMA_dgetrs", info);
    }
    if (info > 0) {
	    error(_("PLASMA routine PLASMA_dgetrs: system is exactly singular"));
    }

    UNPROTECT(2);
    return val;
#endif
    return R_NilValue;
}

SEXP plasma_dgeMatrix_matrix_mm(SEXP a, SEXP bP, SEXP right)
{
#ifdef HIPLAR_WITH_PLASMA
    SEXP b = PROTECT(mMatrix_as_dgeMatrix(bP)),
	val = PROTECT(NEW_OBJECT(MAKE_CLASS("dgeMatrix")));
    int *adims = INTEGER(GET_SLOT(a, Matrix_DimSym)),
	*bdims = INTEGER(GET_SLOT(b, Matrix_DimSym)),
	*cdims = INTEGER(ALLOC_SLOT(val, Matrix_DimSym, INTSXP, 2));
    double one = 1., zero = 0.;
    int info;

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering plasma_dgeMatrix_matrix_mm");
#endif

    if (asLogical(right)) {
	int m = bdims[0], n = adims[1], k = bdims[1];
	if (adims[0] != k)
	    error(_("Matrices are not conformable for multiplication"));
	cdims[0] = m; cdims[1] = n;
	if (m < 1 || n < 1 || k < 1) {
/* 	    error(_("Matrices with zero extents cannot be multiplied")); */
	    ALLOC_SLOT(val, Matrix_xSym, REALSXP, m * n);
	} else {

        info = P_dgemm("N", "N",  m, n, k, one,
			     REAL(GET_SLOT(b, Matrix_xSym)), m,
			     REAL(GET_SLOT(a, Matrix_xSym)), k, zero,
			     REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, m * n)),
			     m);
   	    if (info) {
            error(_("PLASMA routine %s returned error code %d"), "PLASMA_dgemm", info);
        }
    }

    } else {
	int m = adims[0], n = bdims[1], k = adims[1];

	if (bdims[0] != k)
	    error(_("Matrices are not conformable for multiplication"));
	cdims[0] = m; cdims[1] = n;
	if (m < 1 || n < 1 || k < 1) {
/* 	    error(_("Matrices with zero extents cannot be multiplied")); */
	    ALLOC_SLOT(val, Matrix_xSym, REALSXP, m * n);
	} else {

        info = P_dgemm("N", "N", m, n, k, one,
			     REAL(GET_SLOT(a, Matrix_xSym)), m,
			     REAL(GET_SLOT(b, Matrix_xSym)), k, zero,
			     REAL(ALLOC_SLOT(val, Matrix_xSym, REALSXP, m * n)),
			     m);
   	    if (info) {
            error(_("PLASMA routine %s returned error code %d"), "PLASMA_dgemm", info);
        }
    }

    }
    ALLOC_SLOT(val, Matrix_DimNamesSym, VECSXP, 2);
    UNPROTECT(2);
    return val;
#endif
    return R_NilValue;
}

#if(0)
SEXP dgeMatrix_svd(SEXP x, SEXP nnu, SEXP nnv)
{
    int /* nu = asInteger(nnu),
	   nv = asInteger(nnv), */
	*dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    double *xx = REAL(GET_SLOT(x, Matrix_xSym));
    SEXP val = PROTECT(allocVector(VECSXP, 3));

    if (dims[0] && dims[1]) {
	int m = dims[0], n = dims[1], mm = (m < n)?m:n,
	    lwork = -1, info;
	double tmp, *work;
	int *iwork = Alloca(8 * mm, int);
	R_CheckStack();

	SET_VECTOR_ELT(val, 0, allocVector(REALSXP, mm));
	SET_VECTOR_ELT(val, 1, allocMatrix(REALSXP, m, mm));
	SET_VECTOR_ELT(val, 2, allocMatrix(REALSXP, mm, n));
	F77_CALL(dgesdd)("S", &m, &n, xx, &m,
			 REAL(VECTOR_ELT(val, 0)),
			 REAL(VECTOR_ELT(val, 1)), &m,
			 REAL(VECTOR_ELT(val, 2)), &mm,
			 &tmp, &lwork, iwork, &info);
	lwork = (int) tmp;
	work = Alloca(lwork, double);
	R_CheckStack();
	F77_CALL(dgesdd)("S", &m, &n, xx, &m,
			 REAL(VECTOR_ELT(val, 0)),
			 REAL(VECTOR_ELT(val, 1)), &m,
			 REAL(VECTOR_ELT(val, 2)), &mm,
			 work, &lwork, iwork, &info);

    }
    UNPROTECT(1);
    return val;
}


const static double padec [] = /* for matrix exponential calculation. */
{
  5.0000000000000000e-1,
  1.1666666666666667e-1,
  1.6666666666666667e-2,
  1.6025641025641026e-3,
  1.0683760683760684e-4,
  4.8562548562548563e-6,
  1.3875013875013875e-7,
  1.9270852604185938e-9,
};

/**
 * Matrix exponential - based on the _corrected_ code for Octave's expm function.
 *
 * @param x real square matrix to exponentiate
 *
 * @return matrix exponential of x
 */
SEXP dgeMatrix_exp(SEXP x)
{
    const double one = 1.0, zero = 0.0;
    const int i1 = 1;
    int *Dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    const int n = Dims[1], nsqr = n * n, np1 = n + 1;

    SEXP val = PROTECT(duplicate(x));
    int i, ilo, ilos, ihi, ihis, j, sqpow;
    int *pivot = Calloc(n, int);
    double *dpp = Calloc(nsqr, double), /* denominator power Pade' */
	*npp = Calloc(nsqr, double), /* numerator power Pade' */
	*perm = Calloc(n, double),
	*scale = Calloc(n, double),
	*v = REAL(GET_SLOT(val, Matrix_xSym)),
	*work = Calloc(nsqr, double), inf_norm, m1_j/*= (-1)^j */, trshift;
    R_CheckStack();

    if (n < 1 || Dims[0] != n)
	error(_("Matrix exponential requires square, non-null matrix"));
    if(n == 1) {
	v[0] = exp(v[0]);
	UNPROTECT(1);
	return val;
    }

    /* Preconditioning 1.  Shift diagonal by average diagonal if positive. */
    trshift = 0;		/* determine average diagonal element */
    for (i = 0; i < n; i++) trshift += v[i * np1];
    trshift /= n;
    if (trshift > 0.) {		/* shift diagonal by -trshift */
	for (i = 0; i < n; i++) v[i * np1] -= trshift;
    }

    /* Preconditioning 2. Balancing with dgebal. */
    F77_CALL(dgebal)("P", &n, v, &n, &ilo, &ihi, perm, &j);
    if (j) error(_("dgeMatrix_exp: LAPACK routine dgebal returned %d"), j);
    F77_CALL(dgebal)("S", &n, v, &n, &ilos, &ihis, scale, &j);
    if (j) error(_("dgeMatrix_exp: LAPACK routine dgebal returned %d"), j);

    /* Preconditioning 3. Scaling according to infinity norm */
    inf_norm = F77_CALL(dlange)("I", &n, &n, v, &n, work);
    sqpow = (inf_norm > 0) ? (int) (1 + log(inf_norm)/log(2.)) : 0;
    if (sqpow < 0) sqpow = 0;
    if (sqpow > 0) {
	double scale_factor = 1.0;
	for (i = 0; i < sqpow; i++) scale_factor *= 2.;
	for (i = 0; i < nsqr; i++) v[i] /= scale_factor;
    }

    /* Pade' approximation. Powers v^8, v^7, ..., v^1 */
    AZERO(npp, nsqr);
    AZERO(dpp, nsqr);
    m1_j = -1;
    for (j = 7; j >=0; j--) {
	double mult = padec[j];
	/* npp = m * npp + padec[j] *m */
	F77_CALL(dgemm)("N", "N", &n, &n, &n, &one, v, &n, npp, &n,
			&zero, work, &n);
	for (i = 0; i < nsqr; i++) npp[i] = work[i] + mult * v[i];
	/* dpp = m * dpp + (m1_j * padec[j]) * m */
	mult *= m1_j;
	F77_CALL(dgemm)("N", "N", &n, &n, &n, &one, v, &n, dpp, &n,
			&zero, work, &n);
	for (i = 0; i < nsqr; i++) dpp[i] = work[i] + mult * v[i];
	m1_j *= -1;
    }
    /* Zero power */
    for (i = 0; i < nsqr; i++) dpp[i] *= -1.;
    for (j = 0; j < n; j++) {
	npp[j * np1] += 1.;
	dpp[j * np1] += 1.;
    }

    /* Pade' approximation is solve(dpp, npp) */
    F77_CALL(dgetrf)(&n, &n, dpp, &n, pivot, &j);
    if (j) error(_("dgeMatrix_exp: dgetrf returned error code %d"), j);
    F77_CALL(dgetrs)("N", &n, &n, dpp, &n, pivot, npp, &n, &j);
    if (j) error(_("dgeMatrix_exp: dgetrs returned error code %d"), j);
    Memcpy(v, npp, nsqr);

    /* Now undo all of the preconditioning */
    /* Preconditioning 3: square the result for every power of 2 */
    while (sqpow--) {
	F77_CALL(dgemm)("N", "N", &n, &n, &n, &one, v, &n, v, &n,
			&zero, work, &n);
	Memcpy(v, work, nsqr);
    }
    /* Preconditioning 2: apply inverse scaling */
    for (j = 0; j < n; j++)
	for (i = 0; i < n; i++)
	    v[i + j * n] *= scale[i]/scale[j];


    /* 2 b) Inverse permutation  (if not the identity permutation) */
    if (ilo != 1 || ihi != n) {
	/* Martin Maechler's code */

#define SWAP_ROW(I,J) F77_CALL(dswap)(&n, &v[(I)], &n, &v[(J)], &n)

#define SWAP_COL(I,J) F77_CALL(dswap)(&n, &v[(I)*n], &i1, &v[(J)*n], &i1)

#define RE_PERMUTE(I)				\
	int p_I = (int) (perm[I]) - 1;		\
	SWAP_COL(I, p_I);			\
	SWAP_ROW(I, p_I)

	/* reversion of "leading permutations" : in reverse order */
	for (i = (ilo - 1) - 1; i >= 0; i--) {
	    RE_PERMUTE(i);
	}

	/* reversion of "trailing permutations" : applied in forward order */
	for (i = (ihi + 1) - 1; i < n; i++) {
	    RE_PERMUTE(i);
	}
    }

    /* Preconditioning 1: Trace normalization */
    if (trshift > 0.) {
	double mult = exp(trshift);
	for (i = 0; i < nsqr; i++) v[i] *= mult;
    }

    /* Clean up */
    Free(work); Free(scale); Free(perm); Free(npp); Free(dpp); Free(pivot);
    UNPROTECT(1);
    return val;
}

SEXP dgeMatrix_Schur(SEXP x, SEXP vectors)
{
    int *dims = INTEGER(GET_SLOT(x, Matrix_DimSym));
    int vecs = asLogical(vectors), info, izero = 0, lwork = -1, n = dims[0];
    double *work, tmp;
    const char *nms[] = {"WR", "WI", "T", "Z", ""};
    SEXP val = PROTECT(Rf_mkNamed(VECSXP, nms));

    if (n != dims[1] || n < 1)
	error(_("dgeMatrix_Schur: argument x must be a non-null square matrix"));
    SET_VECTOR_ELT(val, 0, allocVector(REALSXP, n));
    SET_VECTOR_ELT(val, 1, allocVector(REALSXP, n));
    SET_VECTOR_ELT(val, 2, allocMatrix(REALSXP, n, n));
    Memcpy(REAL(VECTOR_ELT(val, 2)), REAL(GET_SLOT(x, Matrix_xSym)), n * n);
    SET_VECTOR_ELT(val, 3, allocMatrix(REALSXP, vecs ? n : 0, vecs ? n : 0));
    F77_CALL(dgees)(vecs ? "V" : "N", "N", NULL, dims, (double *) NULL, dims, &izero,
		    (double *) NULL, (double *) NULL, (double *) NULL, dims,
		    &tmp, &lwork, (int *) NULL, &info);
    if (info) error(_("dgeMatrix_Schur: first call to dgees failed"));
    lwork = (int) tmp;
    work = Alloca(lwork, double);
    R_CheckStack();
    F77_CALL(dgees)(vecs ? "V" : "N", "N", NULL, dims, REAL(VECTOR_ELT(val, 2)), dims,
		    &izero, REAL(VECTOR_ELT(val, 0)), REAL(VECTOR_ELT(val, 1)),
		    REAL(VECTOR_ELT(val, 3)), dims, work, &lwork,
		    (int *) NULL, &info);
    if (info) error(_("dgeMatrix_Schur: dgees returned code %d"), info);
    UNPROTECT(1);
    return val;
}
#endif

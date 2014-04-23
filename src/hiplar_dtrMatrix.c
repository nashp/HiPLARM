#include "hiplar_dtrMatrix.h"
#include "plasma_dtrMatrix.h"
#include "magma_dtrMatrix.h"
#include "hiplar_at.h"
#include "hiplar_dbg.h"


SEXP hiplar_dtrMatrix_solve(SEXP a) {


#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dtrMatrix_solve");
#endif

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    
    int *dims = INTEGER(GET_SLOT(a, Matrix_DimSym));
	int size = dims[0];
	if  ((hiplar_library == HIPLAR_USE_PLASMA) ||
         ((hiplar_library == HIPLAR_USE_AUTO) && (size < xover_dtrMatrix_solve))) {
#endif
#ifdef HIPLAR_WITH_PLASMA
        return plasma_dtrMatrix_solve(a);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    } else {
#endif
#ifdef HIPLAR_WITH_MAGMA
        return magma_dtrMatrix_solve(a);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    }
#endif

    return R_NilValue;
}


SEXP hiplar_dtrMatrix_chol2inv(SEXP a) {

#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dtrMatrix_chol2inv");
#endif

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    int *dims = INTEGER(GET_SLOT(a, Matrix_DimSym));
    int size = dims[0];
	if  ((hiplar_library == HIPLAR_USE_PLASMA) ||
         ((hiplar_library == HIPLAR_USE_AUTO) && (size < xover_dtrMatrix_chol2inv))) {
#endif
#ifdef HIPLAR_WITH_PLASMA
        return plasma_dtrMatrix_chol2inv(a);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    } else {
#endif
#ifdef HIPLAR_WITH_MAGMA
        return magma_dtrMatrix_chol2inv(a);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    }
#endif

    return R_NilValue;
}


SEXP hiplar_dtrMatrix_matrix_solve(SEXP a, SEXP b) {


#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dtrMatrix_matrix_solve");
#endif


#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    int *dims = INTEGER(GET_SLOT(a, Matrix_DimSym)); 
	int size = dims[0];
	if  ((hiplar_library == HIPLAR_USE_PLASMA) ||
         ((hiplar_library == HIPLAR_USE_AUTO) && (size < xover_dtrMatrix_matrix_solve))) {
#endif
#ifdef HIPLAR_WITH_PLASMA
        return plasma_dtrMatrix_matrix_solve(a, b);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    } else {
#endif
#ifdef HIPLAR_WITH_MAGMA
        return magma_dtrMatrix_matrix_solve(a, b);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    }
#endif

    return R_NilValue;
}


SEXP hiplar_dtrMatrix_matrix_mm(SEXP a, SEXP b, SEXP right, SEXP trans) {


#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dtrMatrix_matrix_mm");
#endif

#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    
    int *dims = INTEGER(GET_SLOT(a, Matrix_DimSym));
	int size = dims[0];
	if  ((hiplar_library == HIPLAR_USE_PLASMA) ||
         ((hiplar_library == HIPLAR_USE_AUTO) && (size < xover_dtrMatrix_matrix_mm))) {
#endif
#ifdef HIPLAR_WITH_PLASMA
        return plasma_dtrMatrix_matrix_mm(a, b, right, trans);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    } else {
#endif
#ifdef HIPLAR_WITH_MAGMA
        return magma_dtrMatrix_matrix_mm(a, b, right, trans);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    }
#endif

    return R_NilValue;
}


SEXP hiplar_dtrMatrix_dtrMatrix_mm(SEXP a, SEXP b, SEXP right, SEXP trans) {


#ifdef HIPLAR_DBG
R_ShowMessage("DBG: Entering hiplar_dtrMatrix_dtrMatrix_mm");
#endif


#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    int *dims = INTEGER(GET_SLOT(a, Matrix_DimSym)); 
	int size = dims[0];
	if  ((hiplar_library == HIPLAR_USE_PLASMA) ||
         ((hiplar_library == HIPLAR_USE_AUTO) && (size < xover_dtrMatrix_dtrMatrix_mm))) {
#endif
#ifdef HIPLAR_WITH_PLASMA
        return plasma_dtrMatrix_dtrMatrix_mm(a, b, right, trans);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    } else {
#endif
#ifdef HIPLAR_WITH_MAGMA
        return magma_dtrMatrix_dtrMatrix_mm(a, b, right, trans);
#endif
#if defined(HIPLAR_WITH_PLASMA) && defined(HIPLAR_WITH_MAGMA)
    }
#endif

    return R_NilValue;
}

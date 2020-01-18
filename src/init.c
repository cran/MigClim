#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void genClust(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void mcMigrate(void *, void *);
extern void validate(void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"genClust",  (DL_FUNC) &genClust,  9},
    {"mcMigrate", (DL_FUNC) &mcMigrate, 2},
    {"validate",  (DL_FUNC) &validate,  5},
    {NULL, NULL, 0}
};

void R_init_MigClim(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

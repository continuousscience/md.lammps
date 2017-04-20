#include <stdio.h>
#include <string.h>

#define lmp_str \
    LAMMPS *lmp = (LAMMPS *)sil_topointer(S, 1); \
    if(lmp == NULL) { \
        return sil_err(S, "lammps.%s: bad args", __func__); \
    } \
    char *name = sil_toastring(S, 2); \
    if(name == NULL) { \
        return sil_err(S, "lammps.%s: bad args", __func__); \
    }

#define ck_err(cleanup) \
    if(lammps_has_error(lmp->lmp)) { \
        char err[128]; \
        lammps_get_last_error_message(lmp->lmp, err, 128); \
        sil_err(S, "lammps.%s: %s", __func__, err); \
        cleanup; \
        return 0; \
    }

#define run_cmd(...) { \
    char *_cmd; \
    int _clen = asprintf(&_cmd, __VA_ARGS__); \
    int _n=0; for(; _cmd[_n] == ' ' || _cmd[_n] == '\t'; _n++); \
    if(_clen-_n >= 4 && memcmp(_cmd+_n, "run ", 4)) { \
      lammps_command(lmp->lmp, _cmd); \
      _cmd[_clen] = '\n'; \
      write(lmp->dat->fd, _cmd, _clen+1); \
      ck_err(free(_cmd)); \
      _cmd[_clen] = 0; \
      if(_clen-_n >= 11 && !memcmp(_cmd+_n, "create_box ", 11)) \
          lmp->initialized = 1; \
    } \
    free(_cmd); \
    sil_settop(S, 0); }


/* copied from library.h */
void  lammps_close(void *);
int lammps_extract_setting(void *, char *);
void *lammps_extract_global(void *, char *);
void lammps_extract_box(void *, double *, double *,
                                double *, double *, double *, int *, int *);
void *lammps_extract_atom(void *, char *);
void *lammps_extract_compute(void *, char *, int, int);
void *lammps_extract_fix(void *, char *, int, int, int, int);
void *lammps_extract_variable(void *, char *, char *);
void lammps_reset_box(void *, double *, double *, double, double, double);
int lammps_set_variable(void *, char *, char *);
double lammps_get_thermo(void *, char *);

int lammps_get_natoms(void *);
void lammps_gather_atoms(void *, char *, int, int, void *);
void lammps_scatter_atoms(void *, char *, int, int, void *);
// assumes not LAMMPS_BIGBIG
void lammps_create_atoms(void *, int, int *, int *,
                                 double *, double *, int *, int);
int lammps_has_error(void *ptr);
int lammps_get_last_error_message(void *ptr, char * buffer, int buffer_size);

// String -> ST(LAMMPS, Nil)
int command(sil_State *S) {
    size_t len;
    LAMMPS *lmp = (LAMMPS *)sil_getST(S, &len);
    if(lmp == NULL) {
        return sil_err(S, "Invalid LAMMPS ST");
    }
    const char *cmd = sil_tobinary(S, 1, &len);
    if(cmd == NULL) {
        return sil_err(S, "lammps.command: bad arg");
    }

    run_cmd("%.*s", (int)len, cmd);

    sil_pushnil(S);
    return 0;
}

// String -> ST(LAMMPS, Int)
int newDatum(sil_State *S) {
    size_t len;
    LAMMPS *lmp = (LAMMPS *)sil_getST(S, &len);
    if(lmp == NULL) {
        return sil_err(S, "Invalid LAMMPS ST");
    }
    const char *buf = sil_tobinary(S, 1, &len);
    if(buf == NULL) {
        return sil_err(S, "newDatum: String required.", __func__);
    }

    int n = new_datum(&lmp->dat, buf, len);
    declare_datum(lmp, n);
 
    sil_settop(S, 0);
    sil_pushinteger(S, n);
    return 0;
}

// Int -> ST(LAMMPS, String)
int getDatum(sil_State *S) {
    size_t len;
    LAMMPS *lmp = (LAMMPS *)sil_getST(S, &len);
    if(lmp == NULL) {
        return sil_err(S, "Invalid LAMMPS ST");
    }
    int n = sil_tointeger(S, 1);
    if(n < 0) {
        return sil_err(S, "Invalid datum n.\n");
    }
    char *buf = read_datum(lmp->dat, n, &len);

    sil_settop(S, 0);
    sil_pushbinary(S, buf, len, 1);
    return 0;
}

// Int -> String -> ST(LAMMPS, Bool)
/* Can't use this because we replay commands to restart.
 * (and over-written input files prevent this).
int putDatum(sil_State *S) {
    size_t len;
    LAMMPS *lmp = (LAMMPS *)sil_getST(S, &len);
    if(lmp == NULL) {
        return sil_err(S, "Invalid LAMMPS ST");
    }
    int n = sil_tointeger(S, 1);
    const char *buf = sil_tobinary(S, 2, &len);
    if(n < 0 || buf == NULL) {
        return sil_err(S, "lammps.putDatum: invalid args");
    }
    int ret = put_datum(lmp->dat, n, buf, len);

    sil_settop(S, 0);
    sil_pushboolean(S, !ret);
    return 0;
}*/


// String -> Int -> String -> ST(LAMMPS, Nil)
int molecule(sil_State *S) {
    size_t len, l1, l2;
    LAMMPS *lmp = (LAMMPS *)sil_getST(S, &len);
    if(lmp == NULL) {
        return sil_err(S, "Invalid LAMMPS ST");
    }
    const char *s1 = sil_tobinary(S, 1, &l1);
    int n = sil_tointeger(S, 2);
    const char *s2 = sil_tobinary(S, 3, &l2);
    if(s1 == NULL || s2 == NULL) {
        return sil_err(S, "lammps.molecule: String -> Int -> String -> ST(LAMMPS,Nil)");
    }
    /*LmpDatum *dat = *get_datum(&lmp->dat, n);
    if(dat == NULL) {
        return sil_err(S, "lammps.molecule: Invalid datum.");
    }*/

    run_cmd("molecule %.*s ${datum%03d} %.*s", (int)l1, s1, n, (int)l2, s2);

    sil_pushnil(S);
    return 0;
}

// Int -> ST(LAMMPS, Int)
int run(sil_State *S) {
    size_t len;
    LAMMPS *lmp = (LAMMPS *)sil_getST(S, &len);
    if(lmp == NULL) {
        return sil_err(S, "Invalid LAMMPS ST");
    }
    int n = sil_tointeger(S, 1);
    if(n < 0) {
        return sil_err(S, "lammps.run: Negative steps are invalid.");
    }

    char *cmd;
    int clen = asprintf(&cmd, "run %d", n);
    lammps_command(lmp->lmp, cmd);
    free(cmd);
    ck_err();
    sil_settop(S, 0);

    lmp->steps += n;

    sil_pushinteger(S, lmp->steps);
    return 0;
}

// Int -> ST(LAMMPS, Nil)
int read_data(sil_State *S) {
    size_t len;
    LAMMPS *lmp = (LAMMPS *)sil_getST(S, &len);
    if(lmp == NULL) {
        return sil_err(S, "Invalid LAMMPS ST");
    }
    int n = sil_tointeger(S, 1);
    /*LmpDatum *dat = *get_datum(&lmp->dat, n);
    if(dat == NULL) {
        return sil_err(S, "lammps.read_data: Invalid datum.");
    }*/

    run_cmd("read_data ${datum%03d}", n);
    lmp->initialized = 1;

    sil_pushnil(S);
    return 0;
}

// LAMMPS -> String -> a:* -> a
int global(sil_State *S) {
    lmp_str;

    int r = sil_tointeger(S, 3);
    if(!(r == 6 || r == 7)) {
        return sil_err(S, "lammps.global: invalid return type");
    }

    void *ret = lammps_extract_global(lmp->lmp, name);
    free(name);

    sil_settop(S, 0);
    if(r == 7) {
        sil_pushdouble(S, *(double *)ret);
    } else {
        sil_pushinteger(S, *(int *)ret);
    }
    return 0;
}

// String prop, Kind t, Vector v -> ST(LAMMPS, Nil)
int scatter_atoms(sil_State *S) {
    size_t len;
    LAMMPS *lmp = (LAMMPS *)sil_getST(S, &len);
    if(lmp == NULL) {
        return sil_err(S, "Invalid LAMMPS ST");
    }
    int t      = sil_tointeger(S, 2);
    const vector *v  = (const vector *)sil_topointer(S, 3);
    int atoms = lammps_get_natoms(lmp->lmp);
    int count  = v->n / atoms;
    if(!(t == 6 || t == 7)) {
        return sil_err(S, "scatter_atoms: invalid type");
    }
    if(count*atoms != v->n) {
        return sil_err(S, "scatter_atoms: vector length is not a multiple"
                          " of %d atoms", atoms);
    }
    char *name = sil_toastring(S, 1);
    if(name == NULL) {
        return sil_err(S, "scatter_atoms: invalid property name");
    }

    if(type == 6) {
        int *y = (int *)malloc(sizeof(int)*v->n);
        for(int i=0; i<v->n; i++) {
            y[i] = v->x[i];
        }
        lammps_scatter_atoms(lmp->lmp, name, 0, count, y);
        free(y);
    } else {
        lammps_scatter_atoms(lmp->lmp, name, 1, count, v->x);
    }

    free(name);
    sil_settop(S, 0);
    sil_pushnil();
    return 0;
}



/* TODO: add scatter_atoms
 *       implement a mechanism for inferring '+' between vectors
 *       with different implementation types...
 *       implement set_global?
 */
// LAMMPS lmp, String name, Kind t, Int count -> Vector
// count is number of values per atom (usu. 1 or 3)
int gather_atoms(sil_State *S) {
    lmp_str;
    int t = sil_tointeger(S, 3);
    int count = sil_tointeger(S, 4);
    int n = lammps_get_natoms(lmp->lmp)*count;

    if(count < 1 || count > 100) {
        free(name);
        return sil_err(S, "lammps.gather_atoms: invalid count");
    }
    if(!(t == 6 || t == 7)) { // 6 ~> 0, 7 ~> 1
        free(name);
        return sil_err(S, "lammps.gather_atoms: invalid return type");
    }

    vector *v = (vector *)malloc(sizeof(vector) + 8*n);

    lammps_gather_atoms(lmp->lmp, name, t-6, count, v->x);
    free(name);
    if(t == 6) { // convert to doubles in-place
        // assuming sizeof(double) >= sizeof(int)
        int *y = (int *)v->x;
        for(int i=n-1; i >= 0; i--) {
            v->x[i] = y[i];
        }
    }

    sil_settop(S, 0);
    sil_pushvector(S, v);
    return 0;
}

// TODO: extract_fix, extract_variable

// LAMMPS -> id : String -> style : Int -> t : Int -> If t == 0: Float else: Vector
// extract compute values
int extract_compute(sil_State *S) {
    lmp_str;

    int style = sil_tointeger(S, 3);
    int t = sil_tointeger(S, 4);
    if(!(t == 0 || t == 1)) {
        free(name);
        return sil_err(S, "lammps.extract_compute: Unsupported compute type.");
    }
    double *ret = (double *)lammps_extract_compute(lmp->lmp, name, style, t);
    free(name);
    sil_settop(S, 0);

    if(t == 1) {
        int n = 1;
        if(style == 1) n = lammps_get_natoms(lmp->lmp);
        vector *v = sil_copyvector(S, n, ret);
        if(v == NULL)
            return sil_err(S, "lammps.extract_compute: Allocation error.");
        sil_pushvector(S, v);
    } else {
        sil_pushdouble(S, *ret);
    }
    return 0;
}

// LAMMPS -> String -> t : * -> count : Int -> Vector
int gather_atoms(sil_State *S) {
    lmp_str;

    int t = sil_tointeger(S, 3);
    if(!(t == 6 || t == 7)) {
        free(name);
        return sil_err(S, "lammps.gather_atoms: Unsupported output type.");
    }

    int count = sil_tointeger(S, 4);
    int n = count * lammps_get_natoms(lmp->lmp);

    vector *v = (vector *)malloc(sizeof(vector) + 8*n);
    if(v == NULL) {
        return sil_err(S, "lammps.gather_atoms: Memory error");
    }
    v->n = n;
    lammps_gather_atoms(lmp->lmp, name, t == 7, count, v->x);

    free(name);
    sil_settop(S, 0);
    sil_pushvector(S, v);
    return 0;
}

// LAMMPS -> Int
int get_natoms(sil_State *S) {
    LAMMPS *lmp = (LAMMPS *)sil_topointer(S, 1);
    if(lmp == NULL) {
        return sil_err(S, "lammps.get_natoms: bad arg");
    }
    sil_settop(S, 0);

    sil_pushinteger(S, lammps_get_natoms(lmp->lmp));
    return 0;
}

// Int ndim -> LAMMPS
int open(sil_State *S) {
    int n = sil_tointeger(S, 1);
    if(n < 1 || n > 3)
        return sil_err(S, "lammps.open: invalid dimension");

    LAMMPS *lmp = open_lammps();
    if(lmp == NULL)
        return sil_err(S, "lammps.open: initialization error");
    startup_lammps(lmp, 0, 0);

    run_cmd("dimension %d", n);
 
    if(sil_pushlammps(S, lmp)) {
        lammps_close(lmp->lmp);
        free(lmp);
        return sil_err(S, "lammps.open: sil-lammps setup error");
    }
    return 0;
}


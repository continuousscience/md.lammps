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

#define ck_err \
    if(lammps_has_error(lmp->lmp)) { \
        char err[128]; \
        lammps_get_last_error_message(lmp->lmp, err, 128); \
        sil_err(S, "lammps.%s: %s", __func__, err); \
        return 0; \
    }

#define run_cmd(...) { \
    char *_cmd; \
    int _clen = asprintf(&_cmd, __VA_ARGS__); \
    if(strstr(_cmd, "create_box")) lmp->initialized = 1; \
    lammps_command(lmp->lmp, _cmd); \
    if(!lmp->initialized) { \
        _cmd[_clen] = '\n'; \
        write(lmp->dat->fd, _cmd, _clen+1); \
    } \
    free(_cmd); \
    ck_err; \
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
}


// Int -> ST(LAMMPS, Nil)
int molecule(sil_State *S) {
    size_t len;
    LAMMPS *lmp = (LAMMPS *)sil_getST(S, &len);
    if(lmp == NULL) {
        return sil_err(S, "Invalid LAMMPS ST");
    }
    int n = sil_tointeger(S, 1);
    LmpDatum *dat = *get_datum(&lmp->dat, n);
    if(dat == NULL) {
        return sil_err(S, "lammps.molecule: Invalid datum.");
    }

    run_cmd("molecule %s", dat->name);

    sil_pushnil(S);
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
    LmpDatum *dat = *get_datum(&lmp->dat, n);
    if(dat == NULL) {
        return sil_err(S, "lammps.read_data: Invalid datum.");
    }

    lmp->initialized = 1;
    run_cmd("read_data %s", dat->name);

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

// LAMMPS -> String -> a:* -> a
// extract atom-based quantities.
int extract_atom(sil_State *S) {
    lmp_str;

    int r = sil_tointeger(S, 3);
    if(!(r == 6 || r == 7)) {
        return sil_err(S, "lammps.extract_atom: invalid return type");
    }

    void *ret = lammps_extract_atom(lmp->lmp, name);
    free(name);

    sil_settop(S, 0);
    if(r == 7) {
        sil_pushdouble(S, *(double *)ret);
    } else {
        sil_pushinteger(S, *(int *)ret);
    }
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
    startup_lammps(lmp, 0);

    run_cmd("dimension %d", n);
 
    if(sil_pushlammps(S, lmp)) {
        lammps_close(lmp->lmp);
        free(lmp);
        return sil_err(S, "lammps.open: sil-lammps setup error");
    }
    return 0;
}


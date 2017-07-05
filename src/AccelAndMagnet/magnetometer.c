#include "magnetometer.h"
#include "fixarray.h"
#include <string.h>
#include <stdlib.h>

//if you don't undrestand this code, see ST's DT0059 (rotated ellipsoid)

fix16_t vHIcomp[3]; // Hard-Iron compensation vector (3x1)
fix16_t mSIcomp[3][3]; // Soft-Iron compensation matrix (rotation + scale) (3x3)
uvector_t vBp; //Uncalibrated measurments vector (Bpx, Bpy, Bpz)
uvector_t vBf; //Calibrated measurements vector (3x1)

mf16 *D, *DtD, *DtOnes, *Ones, *v, *A3, *A4, *vghi, *offset, *T, *B4, *B3,
        *luB3, *R, *eigenvalues = NULL; //LMS values matrix (50x50) max
//uvector_t BpCollect[MAX_SAMPLES];
uint8_t sampleCounter, *p;

//private functions
void fill_row(mf16 *mat, const uvector_t *B);

void create_auxiliary_matrix_A3_vghi(mf16 *v, mf16 *A3, mf16 *vghi);

void create_auxiliary_matrix_A4_B4_T(mf16 *A4, mf16 *B4, mf16 *T, mf16 *v, mf16 *offset);

void create_auxiliary_matrix_B3_from_B4(mf16 *B3, mf16 *B4);

void calculate_eigenvalues(mf16 *eigenvalues, mf16 *luB3, uint8_t *p, mf16 *B3);

void calculate_eigenvectors(mf16 *R, mf16 *eigenvalues, mf16 *luB3, uint8_t *p);

void calculate_axis_gains(mf16 *gains, mf16 *eigenvalues);

void copy_matrix_to_array(fix16_t *vHIcomp, fix16_t (*mSIcomp)[], mf16 *offset, mf16 *R);

enum {
    CALIB_INIT = 0,
    CALIB_COLLECT_DATA,
    CALIB_COMPUTE,
    CALIB_FINISH,
    CALIB_ERROR,
    CALIB_STOP,
} calibrationState = CALIB_INIT;

void calibrationStateMachine(void) {
    switch (calibrationState) {
        case CALIB_INIT: {
            sampleCounter = 0;
#ifndef TESTDATA
            for (uint8_t i = 0; i < MAX_SAMPLES; i++) {
                BpCollect[i].x = F16(0.0);
                BpCollect[i].y = F16(0.0);
                BpCollect[i].z = F16(0.0);
            }
#endif
            calibrationState = CALIB_COLLECT_DATA;

            break;
        }
        case CALIB_COLLECT_DATA:
#ifndef TESTDATA
            if (sampleCounter < MAX_SAMPLES) {
                //TODO call real functions
                BpCollect[sampleCounter].x = F16(0.0);//getBx();
                BpCollect[sampleCounter].y = F16(0.0);//getBy();
                BpCollect[sampleCounter].z = F16(0.0);//getBz();
                sampleCounter++;
            } else {
                calibrationState = CALIB_COMPUTE;
            }

#else
            sampleCounter = MAX_SAMPLES - 1;
            calibrationState = CALIB_COMPUTE;
#endif
            break;
        case CALIB_COMPUTE:

            D = (mf16 *) malloc(sizeof(mf16));
            DtD = (mf16 *) malloc(sizeof(mf16));
            DtOnes = (mf16 *) malloc(sizeof(mf16));
            Ones = (mf16 *) malloc(sizeof(mf16));

            D->rows = MAX_SAMPLES;
            D->columns = 9;
            mf16_fill(D, F16(0.0));
            D->rows = 1; // trick to fill array column by column

            DtD->rows = 9;
            DtD->columns = 9;
            mf16_fill(DtD, F16(0.0));

            DtOnes->rows = 1;
            DtOnes->columns = 9;
            mf16_fill(DtOnes, F16(0.0));

            Ones->rows = MAX_SAMPLES;
            Ones->columns = 1;
            mf16_fill(Ones, F16(1.0));

            //fill D matrix
            for (uint8_t i = 0; i < sampleCounter; i++) {
                fill_row(D, BpCollect);
                D->rows++;
            }

            mf16_mul_at(DtD, D, D);
            mf16_mul_at(DtOnes, D, Ones);

            //free unused memory
            free(D);
            free(Ones);

            // initalize v matrix
            v = (mf16 *) malloc(sizeof(mf16));
            v->rows = 1;
            v->columns = 9;
            mf16_fill(v, F16(0.0));

            //calculate v matrix
            mf16_invert_gj(DtD, DtD);
            mf16_mul(v, DtD, DtOnes);

            // free unused matrix
            free(DtD);
            free(DtOnes);

            //--------------------------------------------------------------------------

            A3 = (mf16 *) malloc(sizeof(mf16));
            vghi = (mf16 *) malloc(sizeof(mf16));
            offset = (mf16 *) malloc(sizeof(mf16));

            create_auxiliary_matrix_A3_vghi(v, A3, vghi);

            mf16_invert_gj(A3, A3);
            mf16_mul(offset, A3, vghi);

            free(A3);
            free(vghi);

            //--------------------------------------------------------------------------

            A4 = (mf16 *) malloc(sizeof(mf16));
            T = (mf16 *) malloc(sizeof(mf16));
            B4 = (mf16 *) malloc(sizeof(mf16));

            create_auxiliary_matrix_A4_B4_T(A4, B4, T, v, offset);

            mf16_mul(B4, T, B4);
            mf16_mul_bt(B4, B4, T);

            free(A4);
            free(T);
            free(v);

            //--------------------------------------------------------------------------

            B3 = (mf16 *) malloc(sizeof(mf16));

            create_auxiliary_matrix_B3_from_B4(B3, B4);

            free(B4);

            //--------------------------------------------------------------------------

            luB3 = (mf16 *) malloc(sizeof(mf16));
            eigenvalues = (mf16 *) malloc(sizeof(mf16));
            p = (uint8_t *) malloc(sizeof(uint8_t) * B3->rows);

            memcpy(luB3, B3, sizeof(mf16));
            eigenvalues->rows = 1;
            eigenvalues->columns = 3;
            mf16_fill(eigenvalues, F16(0.0));

            mf16_lup_decompose(luB3, p);
            calculate_eigenvalues(eigenvalues, luB3, p, B3);

            free(B3);

            R = (mf16 *) malloc(sizeof(mf16));

            calculate_eigenvectors(R, eigenvalues, luB3, p);

            free(luB3);
            free(p);

            //--------------------------------------------------------------------------
            mf16_invert_gj(R, R);

            calculate_axis_gains(eigenvalues, eigenvalues);
            //at this point eigenvalues become gains

            mf16_mul(R, R, eigenvalues);
            //at this point we have combined rotation + gain matrix

            calibrationState = CALIB_FINISH;

            break;
        case CALIB_FINISH:

            copy_matrix_to_array(vHIcomp, mSIcomp, offset, R);
            free(offset);
            free(R);

            calibrationState = CALIB_STOP;

            break;
        case CALIB_ERROR:
            break;
        case CALIB_STOP:

            calibrationState = CALIB_STOP;

            break;
        default:
            break;
    }
}

//------------------------------------------------------------------------------

void fill_row(mf16 *mat, const uvector_t *B) {
    if (mat->rows != 9) {
        return;
    }
    mat->data[mat->rows - 1][0] = fix16_sq(B->x);
    mat->data[mat->rows - 1][1] = fix16_sq(B->y);
    mat->data[mat->rows - 1][2] = fix16_sq(B->z);
    mat->data[mat->rows - 1][3] = fix16_mul(F16(2.0), fix16_mul(B->x, B->y));
    mat->data[mat->rows - 1][4] = fix16_mul(F16(2.0), fix16_mul(B->x, B->z));
    mat->data[mat->rows - 1][5] = fix16_mul(F16(2.0), fix16_mul(B->y, B->z));
    mat->data[mat->rows - 1][6] = fix16_mul(F16(2.0), B->x);
    mat->data[mat->rows - 1][7] = fix16_mul(F16(2.0), B->y);
    mat->data[mat->rows - 1][8] = fix16_mul(F16(2.0), B->z);
}

void create_auxiliary_matrix_A3_vghi(mf16 *v, mf16 *A3, mf16 *vghi) {
    uint_fast8_t i;

    A3->rows = 3;
    A3->columns = 3;
    vghi->rows = 3;
    vghi->columns = 1;

    mf16_fill(A3, F16(0.0));
    mf16_fill(vghi, F16(0.0));

    for (i = 0; i < 3; i++) {
        vghi->data[i][0] = v->data[i + 6][0];
    }

    //fill main diagonal A3
    for (i = 0; i < 4; i++) {
        A3->data[i][i] = v->data[i][0];
    }

    //fill additional diagonals A3
    A3->data[0][1] = A3->data[1][0] = v->data[3][0];
    A3->data[0][2] = A3->data[2][0] = v->data[4][0];
    A3->data[1][2] = A3->data[2][1] = v->data[5][0];

}

void create_auxiliary_matrix_A4_B4_T(mf16 *A4, mf16 *B4, mf16 *T, mf16 *v, mf16 *offset) {
    uint_fast8_t i;

    A4->rows = 4;
    A4->columns = 4;
    T->rows = 4;
    T->columns = 4;

    B4->rows = 4;
    B4->columns = 4;

    mf16_fill(A4, F16(0.0));
    mf16_fill(B4, F16(0.0));

    //fill T main diagonal
    mf16_fill(T, F16(0.0));
    mf16_fill_diagonal(T, F16(1.0));

    //fill main diagonal A4
    for (i = 0; i < 3; i++) {
        A4->data[i][i] = v->data[i][0];
    }
    A4->data[3][3] = F16(-1.0);

    //fill additional diagonals A4
    A4->data[0][1] = A4->data[1][0] = v->data[3][0];
    A4->data[0][2] = A4->data[2][0] = v->data[4][0];
    A4->data[1][2] = A4->data[2][1] = v->data[5][0];
    A4->data[0][3] = A4->data[3][0] = v->data[6][0];
    A4->data[1][3] = A4->data[3][1] = v->data[7][0];
    A4->data[2][3] = A4->data[3][2] = v->data[8][0];


    //fill offset into transform matrix
    for (i = 0; i < 3; i++) {
        T->data[3][i] = offset->data[i][0];
    }

}

void create_auxiliary_matrix_B3_from_B4(mf16 *B3, mf16 *B4) {
    B3->rows = 3;
    B3->columns = 3;

    mf16_fill(B3, F16(0.0));
    mf16_add(B3, B3, B4);
    mf16_div_s(B3, B3, -B4->data[3][3]);
}

void calculate_eigenvalues(mf16 *eigenvalues, mf16 *luB3, uint8_t *p, mf16 *B3) {
    fix16_t a, b, c, Q, Qpower3, R, S, x1, x2, x3, phi = F16(0.0);

    //B3 characteristic equation's coefficents
    mf16_trace(&a, B3);
    mf16_lup_determinant(&c, luB3, p);
    mf16_mul(B3, B3, B3);
    mf16_trace(&b, B3);
    b = b - fix16_mul(fix16_sq(a), a);
    b = fix16_div(b, F16(2.0));

    fix16_t theta = fix16_mul(fix16_div(F16(2.0), F16(3.0)), a);

    //Viet's trigonometric formula
    Q = fix16_div(fix16_sub(fix16_sq(a), fix16_mul(b, F16(3.0))), F16(9.0));
    R = fix16_mul(fix16_mul(fix16_sq(a), a), F16(2.0));
    R = R - fix16_mul(fix16_mul(a, b), F16(2.0));
    R = R + fix16_mul(c, F16(27.0));
    R = fix16_div(R, F16(54.0));

    Qpower3 = fix16_mul(fix16_sq(Q), Q);
    S = Qpower3 - fix16_sq(R);

    if (S > 0) {
        phi = fix16_div(fix16_acos(fix16_div(R, fix16_sqrt(Qpower3))), F16(3.0));
        x1 = fix16_mul(F16(-2.0), fix16_sqrt(Q));
        x2 = x3 = x1;
        x1 = fix16_mul(x1, fix16_cos(phi)) - fix16_div(a, F16(3.0));
        x2 = fix16_mul(x2, fix16_cos(phi + theta)) - fix16_div(a, F16(3.0));
        x3 = fix16_mul(x3, fix16_cos(phi - theta)) - fix16_div(a, F16(3.0));
    }

    eigenvalues->data[0][0] = x1;
    eigenvalues->data[0][1] = x2;
    eigenvalues->data[0][2] = x3;
}

void calculate_eigenvectors(mf16 *R, mf16 *eigenvalues, mf16 *luB3, uint8_t *p) {
    mf16 tmp, tempSystem;

    fix16_t *b = (fix16_t *) malloc(sizeof(fix16_t) * 3);

    memset(b, F16(0.0), sizeof(fix16_t));

    //initialize temporary matrixes
    tmp.columns = 3;
    tmp.rows = 3;

    tempSystem.columns = 3;
    tempSystem.rows = 3;

    R->columns = 3;
    R->rows = 3;

    mf16_fill(&tmp, F16(0.0));
    mf16_fill(&tempSystem, F16(0.0));
    mf16_fill(R, F16(0.0));

    for (uint8_t i = 0; i < 3; i++) {
        mf16_fill(&tmp, F16(0.0));
        mf16_fill_diagonal(&tmp, -(eigenvalues->data[0][i]));
        mf16_add(&tempSystem, luB3, &tmp);

        mf16_lup_solve(&(R->data[i][0]), &tempSystem, b, p);
    }

}

void calculate_axis_gains(mf16 *gains, mf16 *eigenvalues) {

    for (uint8_t i = 0; i < 3; i++) {
        gains->data[0][i] = fix16_div(F16(1.0), fix16_sqrt(eigenvalues->data[0][i]));
    }

}

void copy_matrix_to_array(fix16_t *vHIcomp, fix16_t (*mSIcomp)[3], mf16 *offset, mf16 *R) {

    for (uint8_t i = 0; i < 3; i++) {
        vHIcomp[i] = offset->data[0][i];
    }

    for (uint8_t i = 0; i < 3; i++) //row
    {
        for (uint8_t j = 0; j < 3; j++) //column
        {
            mSIcomp[i][j] = R->data[i][j];
        }
    }
}

void R_mpy(fix16_t *dest, const fix16_t *a, const fix16_t (*b)[3]) {

    for (uint8_t row = 0; row < 3; row++) {
        dest[row] = fa16_dot(a, 1, &b[row][0], 1, 3);
    }
}

void correct_offset(fix16_t *dest, const fix16_t *offset) {
    for (uint8_t row = 0; row < 3; row++) {
        dest[row] -= offset[row];
    }
}

uvector_t getBf(void) {
    return vBf;
}

void setBp_setBf(uint16_t x, uint16_t y, uint16_t z) {

    vBp.x = fix16_mul(fix16_from_int(x), MAX_MAG_VALUE_XY);
    vBp.y = fix16_mul(fix16_from_int(y), MAX_MAG_VALUE_XY);
    vBp.z = fix16_mul(fix16_from_int(z), MAX_MAG_VALUE_Z);

    memcpy(vBf.arr, vBp.arr, sizeof(fix16_t) * 3);

    correct_offset(vBf.arr, vHIcomp);
    R_mpy(vBf.arr, vBf.arr, mSIcomp);

}
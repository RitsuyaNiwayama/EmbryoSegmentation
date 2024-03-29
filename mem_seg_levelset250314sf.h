#ifndef __MEM_SEG_LEVELSET120214SF_H
#define __MEM_SEG_LEVELSET120214SF_H
#include <iostream>
#include "tiffio.h"
#include "tiff.h"
#include "tiffvers.h"
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <cstdio>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <time.h>
//#include <matrix.h>
#include "llistRN.h"
#include "lsops3cRN.h"
#include "energy3cRN.h"
#include "types_sf.h"
#include "parameters-tools_sf.h"
#include "utility_sf.h"

using namespace std;
//for nuclei file io
//mainly for tiff io
#define TOTFILENAME (3000)
#define MAX_CHAR_NUM (256)
#define OTHER_FILE_NUM (2)
#define WBAND (10)
#define WRESET (3)
#define numdims (3)
#define totalpixels (260*260*270+100)
#define maxpixels (130000)
#define DIMX dims[1]
#define DIMY dims[0]
#define DIMZ dims[2]
#define DIMXY dims[3]
#define NUMEL dims[4]

#define OFFX dims[0]
#define OFFY 1
#define OFFZ dims[3]

//mainly for solving equation
#define ILIM 25000
//#define ILIM 5

#define OUTPUTINVAL 500
#define EPSILON 0.001
#define DT 0.15
#define KURILIM 1
#define EPS_ITEL 0.000000001
#define wa 1.0
#define wd 0.2
#define DELTA 0.35
#define LARGE_K 1000.0
#define MAXCELLNUM 10000
#define wa_G 1.0
#define wd_G 0.2
#define DELTA_G -0.05
//SOR
#define om 1.5
#define nend 3000
#define convp 0.00000000000001
//not used
//for InnerCell
#define DWBAND 1.0
#define LIM_RAD 6
#define z_factor 1.00
#define fin_valu 0.002
#define INITIALTIME 1
#define ENDTIME 60
//for Globalembryo
#define OFS_MIN_Z 5
#define OFS_MIN_Y 5
#define OFS_MIN_X 5
#define OFS_MAX_Z (255)
#define OFS_MAX_Y (255)
#define OFS_MAX_X (255)
#define CENT_X 130.0
#define CENT_Y 130.0
#define CENT_Z 130.0
#define RADMAX 110.0
/*
#define OFS_MIN_Z 5
#define OFS_MIN_Y 5
#define OFS_MIN_X 5
#define OFS_MAX_Z (96-5)
#define OFS_MAX_Y (97-5)
#define OFS_MAX_X (99-5)
#define CENT_X 50.0
#define CENT_Y 50.0
#define CENT_Z 50.0
#define RADMAX 100.0
*/
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

#define DO_GLOBAL_SEGMENTATION (0)

extern char FILENAME[TOTFILENAME][MAX_CHAR_NUM];

extern double dd_out[totalpixels];
extern double d_phi[totalpixels];
extern double dd_g[totalpixels];
extern double d_mask[totalpixels];
extern double d_label[totalpixels];
//newly added 020714
extern int count_file2;
extern DIR* dir2;
extern struct dirent* dp2;
extern char FILENAME2[TOTFILENAME][MAX_CHAR_NUM];
extern char ddirname2[MAX_CHAR_NUM];
//end newly added


extern DIR* dir;
extern struct dirent* dp;
extern DIR* dir_g;
extern struct dirent* dp_g;

extern int count_file;

extern char ddirname[MAX_CHAR_NUM];
extern char diroutname[MAX_CHAR_NUM];
extern char dirlogname[MAX_CHAR_NUM];
extern char dnucdirname[MAX_CHAR_NUM];
extern char headername[MAX_CHAR_NUM];
extern char globalname[MAX_CHAR_NUM];
//for membrane and NarrowBand
extern int NumDetectedNuc;
extern int CircleMap[WBAND+1][300*(WBAND+1)][3];
extern int NCircleMap[WBAND+1];
extern int i_ResetInval;
extern int i_NucX[MAXCELLNUM];
extern int i_NucY[MAXCELLNUM];
extern int i_NucZ[MAXCELLNUM];
extern int i_NucT[MAXCELLNUM];
extern int i_NucC[MAXCELLNUM];
extern int GRIDSIZE_X;
extern int GRIDSIZE_Y;
extern int GRIDSIZE_Z;
extern long dims[5];
extern long mdims[3];
extern TOOL_t Tools; 
extern double new_count;
extern double old_count;
extern FILE* FileLog;
//for nuclei-io
//mainfunction
double ReInitialization(LL *Lz,double setd);
int Indix(int ii, int kk, int jj);
void Calcg();
void CalcRight_SF(LL *Lz);
void CalcQs_SF(LL *Lz);
void SolveUsingSOR_SF(int my_rank,LL *Lz);
void SolveUsingSOR_SF(int my_rank);
double Check_Variation(LL *Lz);
//tiff-io
void SetFilename();
void SetFilename2();
void SetFilename_glob();
void ReadTIFFile2forSF();
void ReadTIFFileMaskforSF();
//Inner cell
void OptimizeNucleiPosition(int my_rank,TOOL_t * const tools);
double Inten(const TOOL_t * const tools, int xx, int yy, int zz,int rad);
void InitializeFrontPositionCellbyCellMPIforSF(int ca,int my_rank);
//Global embryo segmentation
void InitializeFrontPositionGMPIforSF(int ca,int my_rank);
void CalcRight_SF_G(LL *Lz);
void CalcQs_SF_G(LL *Lz);
void SolveUsingSOR_SF_G(int my_rank,LL *Lz);
void WriteFrontPixelsforSF_G(int img_num,LL *Lz,int step_num);
//MPI-io
bool startsWith(const char *pre, const char *str);
void Writeparameters(int myrank,FILE* filp);
void ScanDataALL(int  my_rank);
void ScanData(FILE *fp, int  my_rank);
void WriteFrontPixelsforSF(int my_rank,int img_num, LL *Lz,int step_num);


///Shawn Lankton functions
void lrbac_vessel_yz(double *img, double *phi, double *label, long *dims,
                  LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,int iter, double rad, double lambda, double dthresh,int display);

// "Soft Plaque Detection and Vessel Segmentation" 
// - Lankton, et. al.
void lrbac_vessel_cv(double *img, double *phi, double *label, long *dims,
		     LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,int iter, double rad, double lambda, double dthresh,int display);

// "Localizing Region-based Active Contours" 
// - Lankton, et. al.
void lrbac_chanvese(double *img, double *phi, double *label, long *dims,
                    LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,int iter, double rad, double lambda,int display);

// "Image Segmentation Using Active Contours Driven by the Bhattacharyya
// Gradient Flow" - Michailovich, et al.
void bhattacharyya(double *img, double *phi, double *label, long *dims,
                   LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,int iter,double lambda,int display);

// "Active Contours Without Edges" - Chan and Vese
void chanvese(double *img, double *phi, double *label, long *dims,
              LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,int iter,double lambda,int display);

// "A Variational Framework for Active and Adaptive Segmentation of Vector
// Valued Images" - Rousson and Deriche
void meanvar(double *img, double *phi, double *label, long *dims,
             LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,int iter,double lambda,int display);

// "A Fully Global Approach to Image Segmentation via Coupled Curve Evolution
// Equations" - Yezzi, et. al.
void yezzi(double *img, double *phi, double *label, long *dims,
           LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
           int iter,double lambda,int display);

void grow(double *img, double *phi, double *label, long *dims,
           LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
           int iter,double lambda,int display);

// "Soft Plaque Detection and Vessel Segmentation:
// - Lankton et. al.
void shrink(double *img, double *phi, double *label, long *dims,
           LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
            int iter,double rad,double lambda,int display);

void kappa(double *img, double *phi, double *label, long *dims,
           LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
           int iter,double lambda,int display);
//End of ShawnLankton functions
#endif

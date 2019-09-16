#include "mem_seg_levelset250314sf.h"
void InitializeFrontPositionGMPIforSF(int ca,int my_rank){
  //sphere
  
  int n=0;
  switch (ca) {
  case 1:

    break;
  case 2:
    for(int ii=0;ii<GRIDSIZE_Z;ii++){
      for(int jj=0;jj<GRIDSIZE_Y;jj++){
	for(int kk=0;kk<GRIDSIZE_X;kk++){
	  double dick = (double)(ii-CENT_Z)*(double)(ii-CENT_Z);
	  dick += (double)(jj-CENT_Y)*(double)(jj-CENT_Y);
	  dick += (double)(kk-CENT_X)*(double)(kk-CENT_X);
	  if(ii>OFS_MIN_Z&&ii<OFS_MAX_Z&&
	     jj>OFS_MIN_Y&&jj<OFS_MAX_Y&&
	     kk>OFS_MIN_X&&kk<OFS_MAX_X&&dick<RADMAX*RADMAX){
	    //  d_mask[Indix(ii,jj,kk)]=0;
	  }
	  else{
	    //d_mask[Indix(ii,jj,kk)]=1;
	  }
	}
      }
    }
    break;
  }
}



void CalcRight_SF_G(LL *Lz){
  int x,y,z,idx;
  ll_init(Lz); //begining of list;
  while(Lz->curr != NULL){ 
    x = Lz->curr->x;
    y = Lz->curr->y; z = Lz->curr->z; idx = Lz->curr->idx;
    int h=z;
    int i=y;
    int j=x;
    Lz->curr->d_Right=d_phi[Indix(h,i,j)];
    Lz->curr->d_oldphi=d_phi[Indix(h,i,j)];
    if(h<1||h>GRIDSIZE_Z-3||i<1||i>GRIDSIZE_Y-3||j<1||j>GRIDSIZE_X-3){
    }
    else{
      double d_Grax=(dd_g[Indix(h,i,j+1)]-dd_g[Indix(h,i,j-1)])/2.0;
      double d_Gray=(dd_g[Indix(h,i+1,j)]-dd_g[Indix(h,i-1,j)])/2.0;
      double d_Graz=(dd_g[Indix(h+1,i,j)]-dd_g[Indix(h-1,i,j)])/2.0;
      Lz->curr->d_Right+=-wa_G*DT*max(-d_Grax,0.0)*(d_phi[Indix(h,i,j)]-d_phi[Indix(h,i,j-1)]);
      Lz->curr->d_Right+=-wa_G*DT*min(-d_Grax,0.0)*(d_phi[Indix(h,i,j+1)]-d_phi[Indix(h,i,j)]);	
      Lz->curr->d_Right+=-wa_G*DT*max(-d_Gray,0.0)*(d_phi[Indix(h,i,j)]-d_phi[Indix(h,i-1,j)]);	
      Lz->curr->d_Right+=-wa_G*DT*min(-d_Gray,0.0)*(d_phi[Indix(h,i+1,j)]-d_phi[Indix(h,i,j)]);
      Lz->curr->d_Right+=-wa_G*DT*max(-d_Graz,0.0)*(d_phi[Indix(h,i,j)]-d_phi[Indix(h-1,i,j)]);
      Lz->curr->d_Right+=-wa_G*DT*min(-d_Graz,0.0)*(d_phi[Indix(h+1,i,j)]-d_phi[Indix(h,i,j)]);
    }
    ll_step(Lz);
  }
}

void CalcQs_SF_G(LL *Lz){
  int x,y,z,idx;
  ll_init(Lz); //begining of list;
  while(Lz->curr != NULL){
    x = Lz->curr->x;
    y = Lz->curr->y;
    z = Lz->curr->z;
    idx = Lz->curr->idx;
    int h=z;
    int i=y;
    int j=x;
    double d_Uppo, d_Upmo, d_Umpo, d_Ummo;
    double d_Upop, d_Upom, d_Umop, d_Umom;
    double d_Uopp, d_Uopm, d_Uomp, d_Uomm;
    if(h>0&&h<GRIDSIZE_Z-1&&i>0&&i<GRIDSIZE_Y-1){
      d_Uppo=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i+1,j)]+d_phi[Indix(h+1,i,j)]+d_phi[Indix(h+1,i+1,j)])/4.0;
      d_Upmo=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i-1,j)]+d_phi[Indix(h+1,i,j)]+d_phi[Indix(h+1,i-1,j)])/4.0;
      d_Umpo=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i+1,j)]+d_phi[Indix(h-1,i,j)]+d_phi[Indix(h-1,i+1,j)])/4.0;
      d_Ummo=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i-1,j)]+d_phi[Indix(h-1,i,j)]+d_phi[Indix(h-1,i-1,j)])/4.0;
    }
    else{
      d_Uppo=(d_phi[Indix(h,i,j)]);
      d_Upmo=(d_phi[Indix(h,i,j)]);
      d_Umpo=(d_phi[Indix(h,i,j)]);
      d_Ummo=(d_phi[Indix(h,i,j)]);
    }
    if(h>0&&h<GRIDSIZE_Z-1&&j>0&&j<GRIDSIZE_X-1){
      d_Upop=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i,j+1)]+d_phi[Indix(h+1,i,j)]+d_phi[Indix(h+1,i,j+1)])/4.0;
      d_Upom=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i,j-1)]+d_phi[Indix(h+1,i,j)]+d_phi[Indix(h+1,i,j-1)])/4.0;
      d_Umop=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i,j+1)]+d_phi[Indix(h-1,i,j)]+d_phi[Indix(h-1,i,j+1)])/4.0;
      d_Umom=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i,j-1)]+d_phi[Indix(h-1,i,j)]+d_phi[Indix(h-1,i,j-1)])/4.0;
    }
    else{
      d_Upop=(d_phi[Indix(h,i,j)]);
      d_Upom=(d_phi[Indix(h,i,j)]);
      d_Umop=(d_phi[Indix(h,i,j)]);
      d_Umom=(d_phi[Indix(h,i,j)]);
    }
    if(i>0&&i<GRIDSIZE_Y-1&&j>0&&j<GRIDSIZE_X-1){
      d_Uopp=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i,j+1)]+d_phi[Indix(h,i+1,j)]+d_phi[Indix(h,i+1,j+1)])/4.0;
      d_Uopm=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i,j-1)]+d_phi[Indix(h,i+1,j)]+d_phi[Indix(h,i+1,j-1)])/4.0;
      d_Uomp=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i,j+1)]+d_phi[Indix(h,i-1,j)]+d_phi[Indix(h,i-1,j+1)])/4.0;
      d_Uomm=(d_phi[Indix(h,i,j)]+d_phi[Indix(h,i,j-1)]+d_phi[Indix(h,i-1,j)]+d_phi[Indix(h,i-1,j-1)])/4.0;
    }
    else{
      d_Uopp=(d_phi[Indix(h,i,j)]);
      d_Uopm=(d_phi[Indix(h,i,j)]);
      d_Uomp=(d_phi[Indix(h,i,j)]);
      d_Uomm=(d_phi[Indix(h,i,j)]);
    }
    double d_delpooz,d_delpooy,d_delpoox,d_delmooz,d_delmooy,d_delmoox;
    double d_delopoz,d_delopoy,d_delopox,d_delomoz,d_delomoy,d_delomox;
    double d_deloopz,d_deloopy,d_deloopx,d_deloomz,d_deloomy,d_deloomx;	
    d_delpooz=d_phi[Indix(h+1,i,j)]-d_phi[Indix(h,i,j)];
    d_delmooz=d_phi[Indix(h,i,j)]-d_phi[Indix(h-1,i,j)];
    d_delpooy=d_Uppo-d_Upmo;
    d_delmooy=d_Umpo-d_Ummo;
    d_delpoox=d_Upop-d_Upom;
    d_delmoox=d_Umop-d_Umom;
    d_delopoz=d_Uppo-d_Umpo;
    d_delomoz=d_Upmo-d_Ummo;
    d_delopoy=d_phi[Indix(h,i+1,j)]-d_phi[Indix(h,i,j)];
    d_delomoy=d_phi[Indix(h,i,j)]-d_phi[Indix(h,i-1,j)];
    d_delopox=d_Uopp-d_Uopm;
    d_delomox=d_Uomp-d_Uomm;
    d_deloopz=d_Upop-d_Umop;
    d_deloomz=d_Upom-d_Umom;
    d_deloopy=d_Uopp-d_Uomp;
    d_deloomy=d_Uopm-d_Uomm;
    d_deloopx=d_phi[Indix(h,i,j+1)]-d_phi[Indix(h,i,j)];
    d_deloomx=d_phi[Indix(h,i,j)]-d_phi[Indix(h,i,j-1)];
    double d_Qpoo,d_Qmoo,d_Qopo,d_Qomo,d_Qoop,d_Qoom;
    d_Qpoo=sqrt(EPSILON*EPSILON+d_delpooz*d_delpooz+d_delpooy*d_delpooy+d_delpoox*d_delpoox);
    d_Qmoo=sqrt(EPSILON*EPSILON+d_delmooz*d_delmooz+d_delmooy*d_delmooy+d_delmoox*d_delmoox);
    d_Qopo=sqrt(EPSILON*EPSILON+d_delopoz*d_delopoz+d_delopoy*d_delopoy+d_delopox*d_delopox);
    d_Qomo=sqrt(EPSILON*EPSILON+d_delomoz*d_delomoz+d_delomoy*d_delomoy+d_delomox*d_delomox);
    d_Qoop=sqrt(EPSILON*EPSILON+d_deloopz*d_deloopz+d_deloopy*d_deloopy+d_deloopx*d_deloopx);
    d_Qoom=sqrt(EPSILON*EPSILON+d_deloomz*d_deloomz+d_deloomy*d_deloomy+d_deloomx*d_deloomx);
    double d_Qbarkarisum=0.0;
    double d_Qbarkari=d_delpooz*d_delpooz+d_delpooy*d_delpooy+d_delpoox*d_delpoox;
    if(d_Qbarkari<0.0){
      d_Qbarkari=0.0;
    }
    d_Qbarkarisum+=sqrt(d_Qbarkari);
    
    d_Qbarkari=d_delmooz*d_delmooz+d_delmooy*d_delmooy+d_delmoox*d_delmoox;
    if(d_Qbarkari<0.0){
      d_Qbarkari=0.0;
    }
    d_Qbarkarisum+=sqrt(d_Qbarkari);
    
    d_Qbarkari=d_delopoz*d_delopoz+d_delopoy*d_delopoy+d_delopox*d_delopox;
    if(d_Qbarkari<0.0){
      d_Qbarkari=0.0;
    }
    d_Qbarkarisum+=sqrt(d_Qbarkari);
    
    d_Qbarkari=d_delomoz*d_delomoz+d_delomoy*d_delomoy+d_delomox*d_delomox;
    if(d_Qbarkari<0.0){
      d_Qbarkari=0.0;
    }
    d_Qbarkarisum+=sqrt(d_Qbarkari);
    
    d_Qbarkari=d_deloopz*d_deloopz+d_deloopy*d_deloopy+d_deloopx*d_deloopx;
    if(d_Qbarkari<0.0){
      d_Qbarkari=0.0;
    }
    d_Qbarkarisum+=sqrt(d_Qbarkari);
    
    d_Qbarkari=d_deloomz*d_deloomz+d_deloomy*d_deloomy+d_deloomx*d_deloomx;
    if(d_Qbarkari<0.0){
      d_Qbarkari=0.0;
    }
    d_Qbarkarisum+=sqrt(d_Qbarkari);
    d_Qbarkarisum/=6.0;
    double d_Qbar;
    d_Qbar=sqrt(EPSILON*EPSILON+d_Qbarkarisum*d_Qbarkarisum);
    Lz->curr->d_czp =    -wd_G*dd_g[Indix(h,i,j)]*DT*(d_Qbar/d_Qpoo);
    Lz->curr->d_czm =    -wd_G*dd_g[Indix(h,i,j)]*DT*(d_Qbar/d_Qmoo);
    Lz->curr->d_cyp =    -wd_G*dd_g[Indix(h,i,j)]*DT*(d_Qbar/d_Qopo);
    Lz->curr->d_cym =    -wd_G*dd_g[Indix(h,i,j)]*DT*(d_Qbar/d_Qomo);
    Lz->curr->d_cxp =    -wd_G*dd_g[Indix(h,i,j)]*DT*(d_Qbar/d_Qoop);
    Lz->curr->d_cxm =    -wd_G*dd_g[Indix(h,i,j)]*DT*(d_Qbar/d_Qoom);
    Lz->curr->d_cce = 1.0+wd_G*dd_g[Indix(h,i,j)]*DT*(d_Qbar/d_Qpoo+d_Qbar/d_Qmoo+d_Qbar/d_Qopo+d_Qbar/d_Qomo+d_Qbar/d_Qoop+d_Qbar/d_Qoom);
    double d_Dpoo,d_Dmoo,d_Dopo,d_Domo,d_Doop,d_Doom;//zyx
    d_Dpoo =min(d_phi[Indix(h+1,i,j)]-d_phi[Indix(h,i,j)],0);
    d_Dpoo*=min(d_phi[Indix(h+1,i,j)]-d_phi[Indix(h,i,j)],0);
    d_Dmoo =min(d_phi[Indix(h-1,i,j)]-d_phi[Indix(h,i,j)],0);
    d_Dmoo*=min(d_phi[Indix(h-1,i,j)]-d_phi[Indix(h,i,j)],0);
    d_Dopo =min(d_phi[Indix(h,i+1,j)]-d_phi[Indix(h,i,j)],0);
    d_Dopo*=min(d_phi[Indix(h,i+1,j)]-d_phi[Indix(h,i,j)],0);
    d_Domo =min(d_phi[Indix(h,i-1,j)]-d_phi[Indix(h,i,j)],0);
    d_Domo*=min(d_phi[Indix(h,i-1,j)]-d_phi[Indix(h,i,j)],0);
    d_Doop =min(d_phi[Indix(h,i,j+1)]-d_phi[Indix(h,i,j)],0);
    d_Doop*=min(d_phi[Indix(h,i,j+1)]-d_phi[Indix(h,i,j)],0);
    d_Doom =min(d_phi[Indix(h,i,j-1)]-d_phi[Indix(h,i,j)],0);
    d_Doom*=min(d_phi[Indix(h,i,j-1)]-d_phi[Indix(h,i,j)],0);
    double d_Mioo = max(d_Dpoo,d_Dmoo);
    double d_Moio = max(d_Dopo,d_Domo);
    double d_Mooi = max(d_Doop,d_Doom);
    double ssq = d_Mioo+d_Moio+d_Mooi;
    if(ssq<0.0){
      ssq = 0.0;
    }
    else{
      ssq = sqrt(ssq);
    }
    Lz->curr->d_Right+=DELTA_G*DT*dd_g[Indix(h,i,j)]*ssq;
    ll_step(Lz);
  }
}


	

void SolveUsingSOR_SF_G(int my_rank,LL *Lz){
  
  int enu;
  int i,j,k;
  double delp,pdel;
  double omg,t1;
  int flagg = 1;
  int i1,j1,k1;
  int ie,je,ke;
  int ii,jj,kk;
  double pdelmax;
  // I have to solve from 0toGRIDSIZE-1
  // in this program I set 0 and GRIDSIZE-1 Dirichlet boundary
  i1 = GRIDSIZE_Z-1-1;
  j1 = GRIDSIZE_Y-1-1;
  k1 = GRIDSIZE_X-1-1;
  for(enu=1; enu<=nend;enu++){
    delp = 0.0;
    pdelmax =-1.0;    
    if(flagg==1){
      int x,y,z,idx;
      ll_init(Lz); //begining of list;
      while(Lz->curr != NULL){
	x = Lz->curr->x;
	y = Lz->curr->y;
	z = Lz->curr->z;
	idx = Lz->curr->idx;
	int i=z;
	int j=y;
	int k=x;
	if(i<1||i>GRIDSIZE_Z-3||j<1||j>GRIDSIZE_Y-3||k<1||k>GRIDSIZE_X-3)continue;
	omg = om/Lz->curr->d_cce;
	t1= (Lz->curr->d_cce*d_phi[Indix(i,j,k)]
	     +Lz->curr->d_czp*d_phi[Indix(i+1,j,k)]
	     +Lz->curr->d_czm*d_phi[Indix(i-1,j,k)]
	     +Lz->curr->d_cyp*d_phi[Indix(i,j+1,k)]
	     +Lz->curr->d_cym*d_phi[Indix(i,j-1,k)]
	     +Lz->curr->d_cxp*d_phi[Indix(i,j,k+1)]
	     +Lz->curr->d_cxm*d_phi[Indix(i,j,k-1)]-Lz->curr->d_Right);
	pdel = -omg*t1;
	if(abs(t1)>pdelmax){
	  pdelmax=t1;
	}
	delp = delp + pdel*pdel;
	d_phi[Indix(i,j,k)] =d_phi[Indix(i,j,k)] +pdel;
	ll_step(Lz);
      }
      delp = sqrt(delp/(double)((GRIDSIZE_X-2)*(GRIDSIZE_Y-2)*(GRIDSIZE_Z-2)));
      //rintf("res %f\n",delp);
      if(pdelmax<convp){
	//	printf("fin\n");
	flagg=0;
      }
      
    }
  }
  if(flagg==1){
    //  printf("notfin %d\n",my_rank);
  }
}
void WriteFrontPixelsforSF_G(int img_num,LL *Lz, int step_num){
  char infile[MAX_CHAR_NUM];
  infile[0] = '\0';
  strcpy(infile,diroutname);
  char infile2[MAX_CHAR_NUM];
  sprintf(infile2,"out_t%04d_globalFR%04d.txt",img_num+1,step_num);
  strcat(infile, infile2);
  if(step_num%5==0){
    FILE* fq;
    fq=fopen(infile,"w");
    ll_init(Lz);
    int count = 0;
    if(Lz!=NULL){
      while(Lz->curr != NULL){
	fprintf(fq,"%ld\t%ld\t%ld\n", Lz->curr->z,Lz->curr->y,Lz->curr->x);
	ll_step(Lz);
	// printf("bullsshit\n");
	count++;
      }
    }
    fclose(fq);
  }
  ll_init(Lz);
  int countn = 0;
  if(Lz!=NULL){
    while(Lz->curr != NULL){
      ll_step(Lz);
      // printf("bullsshit\n");
      countn++;
    }
  }
  new_count = countn;
  
}

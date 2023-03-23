#include <stdio.h>
#include <math.h>
#include <stdlib.h> 
#include <string.h> 
#include <cpgplot.h>
#include "/Users/stephencurran/C/functions/mean_func.h"
#include "smooth_func.h"
// #define npts 1000000000

// ON MAC
// gcc -c interspec.c -I/usr/local/pgplot -o interspec.o; gfortran -o interspec interspec.o -I/opt/local/include -L/opt/local/lib -lcpgplot -lpgplot -lpng -lz;./interspec 

int main()
{
  int i,j,k,q,n,numChannels,npts=99999;
  float *xline,*yline,*x,*y,*X,*Y, text = 1.5, c = 2.99792458e5; // km/s;
  char infile[200],outfile[200],psfile[100],systemCall[300],lineInput[npts],format[10],label[200];
  FILE *f1,*f2,*f3,*f4,*f5,*finfo, *fsize;;
  //////////////// MEAN /////////////////////////
  float VAL[99999], *val, V, *mean, prob, E, e[99999],expec, v[99999], var, med,ave, sigma, Ave;
  char sigma_string[200]; 
  /////////////// SPECIFIC /////////////////////////
  int colour, style, poly_order;
  double coeff[10];
  float obsfreq, xmin, xmax, ymin,ymax,chave, pubflux, fracerr, tau_g, delta_v, delta_f, d_f_ave, delta_z,d_z_ave;
  float *cont, *noise, RMS1, RMS2, RMS, *a,*area,A,  v_ave, ypeak, vpeak, v_1, v_2, FWHM,optdepth, inttau, *s, *prod;
  char units[10],smooth[10],intact[10],emabs[10],dir[300], trans_text[100],zoom[10];
  float G[50],H[50],L[50],K[50], *vx1, *vy1, *vx2, *vy2, VX1, VY1, VX2, VY2; //cursor stuff */
  float *x1, *y1,*x2, *y2, X1,X2,Y1,Y2,*xstart,*xend,XSTART,XEND,*ystart,*yend,YSTART,YEND; 
  float res, sf, *ex, *why; //,*EX, *WHY; // changed from arrays to pointers as seg fault */
  FILE *fspec;
  system("ls *.dat");
  system("rm cursor.dat line.dat vertices.dat whole.dat bpass.dat zoom.dat zero.dat");

  //strcpy(infile, "1555-140.dat"); // for testing
  printf("\nEnter the input infile: "); scanf("%s",infile);
  while ( (fspec=fopen(infile, "r") ) == NULL)
    {
      printf("Error opening file %s\n", infile);system("ls *.txt");printf("Try again: ");scanf("%s",infile);
    }
  ////////////// remove the gunk at the top of files put out by GBTIDL ///////////////////////////
  sprintf(systemCall,"sed 's/E/e/g' %s > temp.dat", infile); system(systemCall); // changing E to e
  sprintf(systemCall,"awk 'match($0,/[A-Z]/) == 0 {print $0}' temp.dat > temp2.dat"); system(systemCall);
  sprintf(systemCall,"mv temp2.dat %s", infile); system(systemCall);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  sprintf(systemCall,"wc -l %s > /tmp/filesize",infile); system(systemCall); fsize = fopen("/tmp/filesize","r");
  fscanf(fsize,"%d",&numChannels);fclose(fsize);//  printf("numChannels = %d\n", numChannels);

  // x = (float *)calloc(npts,sizeof(float)); y = (float *)calloc(npts,sizeof(float));
  //xline = (float *)calloc(npts,sizeof(float)); yline = (float *)calloc(npts,sizeof(float));
  x= (float*)malloc(npts*sizeof(float)); y = (float*)malloc(npts*sizeof(float)); 
  xline = (float*)malloc(npts*sizeof(float)); yline = (float*)malloc(npts*sizeof(float)); 

  ////////////////////////////// X-AXIS//////////////////////////////////
  sprintf(systemCall,"head -2 %s", infile);system(systemCall);
  sprintf(systemCall,"tail -2 %s", infile);system(systemCall);
  printf("\n==========================================================");
  printf("\n If from GBTIDL stick with frequency version as velocities");
  printf("\n  being doubled screwing up FHWM and redshifts");
  printf("\n==========================================================\n");
  printf("\nVelocity in km/s [k],\nm/s (as output by kvis) [m],\nthat stupid cz [c] \nMHz or GHz [h]: "); scanf("%s",units);
  //strcpy(units, "m"); // for testing

  if (strcmp(units, "c") == 0) {printf("\nInput reference redshift: "); scanf("%f",&obsfreq);} 

  if (strcmp(units, "h") == 0) {
    printf("\nInput reference frequency, where line is expected [MHz or GHz - however given above]: "); scanf("%f",&obsfreq);}
  
  printf("Is flux intact [n - if already in tau / normalised to zero / emission]? "); scanf("%s", intact);
 // strcpy(intact,"y"); // for testing

  if (strcmp(intact, "n") == 0){
    printf("Emission [e] or absorption [any]? "); scanf("%s", emabs);
    if (strcmp(emabs, "e") == 0) pubflux=0;
    else {
      printf("Input flux given in paper [same units as file - input 0 if already in tau]: ");
      scanf("%f",&pubflux);
    }
  }

  ////////////////////////////// Y-AXIS//////////////////////////////////

  i=0;
  while(fgets(lineInput, npts, fspec) > 0){
    if(lineInput[0]!='#'){
      sscanf(lineInput,"%e%e", &x[i],&y[i]); {
	if  (strcmp(units, "k") == 0) x[i] = x[i];
	else if (strcmp(units, "m") == 0) x[i] = x[i]/1000;
	else if (strcmp(units, "c") == 0) x[i] = x[i] - obsfreq*c ; // actually redshift
	else if (strcmp(units, "h") == 0) x[i] = c*((pow(obsfreq,2) - pow(x[i],2))/(pow(obsfreq,2) + pow(x[i],2)));

	if (strcmp(emabs, "e") == 0) y[i] = -1*y[i]; // inverting may be easiest way
	if (strcmp(intact, "n") == 0 && pubflux >0) y[i] = y[i] + pubflux;

	if(i==0) xmin = xmax = x[i];
	else{
	  if(xmin>x[i]) xmin=x[i];
	  if(xmax<x[i]) xmax=x[i];
	}
	if(i==0) ymin = ymax = y[i];
	else{
	  if(ymin>y[i]) ymin=y[i];
	  if(ymax<y[i]) ymax=y[i];
	}
	f3 = fopen ("whole.dat", "a");fprintf(f3,"%1.3f %1.4f\n", x[i], y[i]);fclose(f3); // whole spectrum
	i++;
      }
    }
  } fclose(fspec);
  ///////////////////////////// ON SCREEN PLOT /////////////////////
  cpgopen("/xs"); cpgpap(8.,0.72); cpgslw(3); cpgsch(text);
  cpgenv(xmin,xmax,ymin-((ymax-ymin)/8),ymax+((ymax-ymin)/8),0,0);
  if (strcmp(intact, "y") == 0) sprintf(label,"Flux density [Jy]");
  else{
    if (pubflux==0){
      if (strcmp(emabs, "e") == 0) sprintf(label,"Inverted flux density [Jy]");
      else sprintf(label,"Optical depth, \\gt\\dobs");
    }
    else sprintf(label,"Flux density [mJy/Jy]");
  }
  cpgmtxt("B",2.5,0.5,0.5,"Velocity [km/s]"); cpgmtxt("L",2.5,0.5,0.5,label); cpgmtxt("T",2.0,0.5,0.5,infile);
  

  colour =7;
  cpgsci(colour); cpgline(i,x,y); // draw spectrum

  chave = (xmax - xmin)/(numChannels-1); // -1 puts closer to mean
  printf("No. of channels = %d, mean width =  %1.2f km/s mean\n",numChannels,chave);
 
  ///////////////////////////////// ZOOM ////////////////////////////
  
   printf("Zoom in [y/n]? "); scanf("%s",zoom);
   //strcpy(zoom,"n"); // for testing
  do{
    system("rm zoom.dat");
    if (strcmp(zoom, "y") == 0) {
      system("rm vertices.dat");
      printf("Mark 2 points for top left and bottom right of range");
      cpgsci(1);
      n = 0; // THIS IS CRUCIAL OR K[0],L[0] JUST GET REPEATED!!!
      cpgncur(50,&n,K,L,5); // max has to be more than 2 as repeating loop
      f1 = fopen ("vertices.dat", "a+");
      for (q=0; q<50; ++q)fprintf(f1,"%1.3f %1.3f\n", K[q],L[q]); fclose(f1);
             
      vx1 = &K[0]; vy2 = &L[0]; vx2 = &K[1]; vy1 = &L[1];  // vertex 1 & 2 - 2 first in y as upper level
      VX1 = *vx1;  VY1 = *vy1; VX2 = *vx2; VY2 = *vy2; //converting from pointers to numbers
 
      cpgask(0); // switch off page prompting
      
      if (VY2 > VY1) cpgenv(VX1,VX2,VY1,VY2,0,0);
      else cpgenv(VX1,VX2,VY2,VY1,0,0);
      cpgmtxt("B",2.5,0.5,0.5,"Velocity [km/s]"); cpgmtxt("L",2.5,0.5,0.5,label);
      colour = 7; cpgsci(colour);cpgline(i,x,y); cpgsci(1);

      // printf("\nVX1 = %1.2f VX2 = %1.2f\n", VX1, VX2);
      f3 = fopen ("zoom.dat", "a");
      for (i=0;i<numChannels;i++) if (x[i]>VX1 && x[i]<VX2) fprintf(f3,"%1.3f %1.4f\n", x[i], y[i]);
      fclose(f3);

      printf("\nZoom in [hard copy range will be determined by cursor points] [y/n]: "); scanf("%s",zoom);
    }
    else{ 
      system("cp whole.dat zoom.dat"); 
      VX1 = xmin,VX2 = xmax,VY1 = ymin-((ymax-ymin)/8),VY2 = ymax+((ymax-ymin)/8); // as originally
    }

  } while(strcmp(zoom, "y") == 0);
  
  f1 = fopen ("zoom.dat", "r"); for (i = 1; fscanf(f1,"%f %f",  &x[i], &y[i]) !=EOF; ++i) {} //re-reading
  
  sprintf(systemCall,"wc -l zoom.dat > /tmp/filesize"); system(systemCall); fsize = fopen("/tmp/filesize","r");
  fscanf(fsize,"%d",&numChannels);fclose(fsize);//  printf("numChannels = %d\n", numChannels);
  
  /////////////// SMOOTHING DATA - moving average ///////////////////////////
  // now after zoom becasue of memoery when dealing with huge GBT spectra - e.g.  0108+388_HI.dat
  int start,count,old_count,inc = 100; // inc IS THE FACTOR EACH PAIR OF POINTS IS BEING SPLIT INTO - 10 NOT ACCURATE ENOUGH
     
  // X= (float *)calloc(npts,sizeof(float)); Y = (float *)calloc(npts,sizeof(float));
  //ex = (float *)calloc(npts,sizeof(float));  why = (float *)calloc(npts,sizeof(float));
   X= (float*)malloc(npts*sizeof(float)); Y = (float*)malloc(npts*sizeof(float)); 
   ex= (float*)malloc(npts*sizeof(float)); why = (float*)malloc(npts*sizeof(float)); 

  if (strcmp(emabs, "e") != 0){ // 18423+5938_CI.dat a mess
     printf("Smooth data [y/n]?: "); scanf("%s",smooth); // can read smoothed to file if want rms
     if (strcmp(smooth, "y") == 0) {// chave was chwidth
       system("cp zoom.dat zoom.txt");
      printf("New resolution, must be > %1.2f km/s: ", chave); scanf("%f",&res); // the smooth factor is therefore
      while (res < chave){
	printf("Error as %1.2f < %1.2f, try again: \n", res,chave); scanf("%f",&res);
      }
      smooth_func(i,j,k,n,numChannels,q,inc,sf,res,chave,start,count,old_count,x,y,X,Y,ex,why,systemCall,lineInput,f1);
      printf("\n======== AWKING HERE ==========\n");
      system("awk < smooth_out.dat '{print \" \" $2 \" \" $3}' > zoom.dat"); // need escape characters round quotes
    }
    else res = chave;
  }
  else res = chave;

  ///////////////////////////// MARKING POINTS ////////////////////////////////////
  colour = 2, style = 1;
 
  do{
    system("rm cursor.dat line.dat gauss.dat vel.dat");
    printf("\nMark 2 points to remove continuum [EVEN IF PUTTING IN FLUX BY HAND, THIS IS IMPORTANT IN DETERMING TAU AND FWHM]");
    printf("\nMark start and end points of profile - \n A [left mouse] to add, DO ONCE FOR LEFT AND ONCE FOR RIGHT THEN, \n X [right mouse] to complete [exit] for baseline [red points],\n then repeat procedure for velocity limits of profile [green]\n use D to delete any points ");
    n = 0;
    cpgsci(colour); cpgncur(4,&n,K,L,17); //need to do a linear fit between these 2 points
    for (k=0; k<n; ++k){
      f1 = fopen ("line.dat", "a+"); fprintf(f1,"%1.3f %1.3f\n", K[k],L[k]); fclose(f1);
      x1 = &K[0]; y1 = &L[0]; x2 = &K[1]; y2 = &L[1];
    }
    X1 = *x1; X2 = *x2; Y1 = *y1; Y2 = *y2;
    //----------------------------------------//
    k = 0;
    cpgsci(colour+1); cpgncur(2,&k,K,L,17);
    for (q=0; q<k; ++q){
      f1 = fopen ("cursor.dat", "a+"); if (K[q]!=0)fprintf(f1,"%1.3f %1.3f\n", K[q],L[q]); fclose(f1);
    }
    f1 = fopen ("cursor.dat", "r"); //open file and read
    for (q=0; q<k; ++q){
      xstart = &K[0]; ystart = &L[0];xend = &K[1]; yend = &L[1];
    }
    XSTART = *xstart; YSTART = *ystart; XEND = *xend; YEND = *yend;//want to convert from pointer t
    printf("\nBetween X1 = %1.2f and XSTART = %1.2f /  XEND = %1.2f  X2 =  %1.2f\n", X1, XSTART, XEND, X2);
  
    //---- Polynomail fitting and bandpass -----//
     
    f1 = fopen ("zoom.dat", "r");
    for (i = 1; fscanf(f1,"%f %f",  &x[i], &y[i]) !=EOF; ++i) {
      f2 = fopen ("bpass.dat", "a"); if ((x[i] > X1 && x[i] < XSTART) || (x[i] > XEND && x[i] < X2)) fprintf(f2,"%1.3f %1.3f\n", x[i], y[i]); fclose(f2);
      f3 = fopen ("mean_out.txt", "a"); if ((x[i] > X1 && x[i] < XSTART) || (x[i] > XEND && x[i] < X2)) fprintf(f3,"%1.3f\n",y[i]); fclose(f3);
    } fclose(f1);                                  // over defined range avoiding line
    system("cp mean_out.txt mean_check.txt");

    mean_func(numChannels, i, systemCall, fsize, mean, VAL, val, V, f1, prob, E, e,expec,v,var,sigma_string,med);
    f1= fopen("mean_values.txt","r");while (fscanf(f1,"%f %f %d", &ave, &sigma, &n) !=EOF);
    fclose(f1);

    if  (strcmp(intact, "n") == 0 && pubflux >0) ave = pubflux;

    ////////////////////////////
    
    system("cp bpass.dat xrange.dat"); // for poly-fit_interspec.py
    printf("\nFit polynomial order [0-9, zero to not bother for RFI dominated]: ");  scanf("%d", &poly_order);
    sprintf(systemCall,"pwd > /tmp/filesize"); system(systemCall);
    f2 = fopen("/tmp/filesize","r"); fscanf(f2,"%s",dir); fclose(f2);
    sprintf(systemCall,"/Users/stephencurran/python/poly-fit_interspec.py %s %d", dir, poly_order); system(systemCall);
    
    cont = (float *)calloc(npts,sizeof(float));  a = (float *)calloc(npts,sizeof(float));
    s = (float *)calloc(npts,sizeof(float));  prod = (float *)calloc(npts,sizeof(float));

    f1 = fopen ("poly-fit_coeff.txt", "r");for (j = poly_order; fscanf(f1,"%lf", &coeff[j]) !=EOF; --j){} fclose(f1); //
    
    f1 = fopen ("zoom.dat", "r");
    for (i = 1; fscanf(f1,"%f %f",  &x[i], &y[i]) !=EOF; ++i) {
      cont[i] = coeff[9]*pow(x[i],9) +  coeff[8]*pow(x[i],8) +  coeff[7]*pow(x[i],7) +  coeff[6]*pow(x[i],6) +  coeff[5]*pow(x[i],5) + coeff[4]*pow(x[i],4) +  coeff[3]*pow(x[i],3) +  coeff[2]*pow(x[i],2) + coeff[1]*pow(x[i],1) + coeff[0]*pow(x[i],0);
      cpgsci(colour);cpgsls(style); cpgsch(0.7*text); cpgpt1(x[i],cont[i],17);
    
      y[i] = y[i]-cont[i]; //subtracting continuum fit   ****
      // Y[i] = y[i]+ave;cpgsci(2);cpgsls(3); cpgslw(2); cpgline(i,x,cont); cpgsls(1); cpgslw(3);
      //cpgsci(2); cpgslw(2);cpgsls(3); cpgmove(X1,Y1); cpgdraw(X2,Y2);
     
      //f2 = fopen ("zero.dat", "a");fprintf(f2,"%1.3f %1.4f\n", x[i], y[i]); fclose(f2); // save zeroed to file

      if (x[i] > X1 &&  x[i] <= XSTART){
	Y[i]= pow(y[i],2); Y[i]=Y[i]+Y[i-1]; noise =  &Y[i]; RMS = *noise;
	n = (XSTART-X1)/res; RMS1 =sqrt(RMS/n);
	//	printf("x[%d] = %1.1f y[%d] = %1.3e Y[%d]= %1.3e Y[%d]= %1.3e RMS1 = %1.3e\n", i, x[i], i, y[i],i,Y[i], i-1, Y[i-1], RMS1);
      }
      if (x[i] < XEND + res && x[i-1] > XEND-res){ //printf("Counting from i = %d\n", i);
	start = (float)i;}
      if (i > start  &&  x[i] < X2){
	Y[start] = 0; // need this
	Y[i]= pow(y[i],2); Y[i]=Y[i]+Y[i-1];
	noise =  &Y[i]; RMS = *noise;
	n = (X2-XEND)/res; RMS2 = sqrt(RMS/n);
	// printf("n = %d x[%d] = %1.1f y[%d] = %1.3e Y[%d]= %1.3e Y[%d]= %1.3e RMS2 = %1.3e\n", n, i, x[i], i, y[i],i,Y[i], i-1, Y[i-1], RMS2);
      }
      RMS = (RMS1 + RMS2)/2;
           
      if (x[i] < XSTART + res && x[i-1] > XSTART-res){ //printf("Abs starts at i = %d\n", i);
	start = (float)i;}

      if (x[i] > XSTART && x[i] < XEND){
       	a[i]=(x[i]-x[i-1])*(y[i] - 0.5*(y[i]-y[i-1])); // see PHS343 notes
	a[i]= a[i] +a[i-1];
	//	printf("x[%d]-x[%d] = %1.2f, y[%d] = %1.3e, y[%d] = %1.3e, a[%d]= %1.3f a[%d]= %1.3f\n", i, i-1, x[i]-x[i-1], i, y[i],i-1, y[i-1],  i, a[i], i-1, a[i-1]);
      	area = &a[i];A = *area;   // SUBTRACT BASELINE AND FIT GAUSSIAN - for area, vel, etc
	
	prod[i] = x[i]*y[i]; // GETTING WEIGHTED VELOCITY OFFSET
	prod[i] = prod[i] + prod[i-1];
	s[i] = y[i];// +y[i-1];
	s[i] = s[i]+s[i-1];
	v_ave = prod[i]/s[i]; // NOT FWHM BUT OFFSET

	//y[i] = y[i]+cont[i]; // to test
	if(i==start+1) {
	  ymin = y[i]; vpeak = x[i];
	}
	else{
	  if(ymin>y[i]) {ymin=y[i]; vpeak = x[i];}
	}
	ypeak = ymin;
	if (strcmp(intact, "y") == 0 || pubflux >0){optdepth = -1*log(1 - ypeak/ave); inttau = A/ave; fracerr= sigma/ave;}
	else {optdepth = ypeak; inttau = A; fracerr = sigma/ypeak;}

	if ((y[i] <= ypeak/2) && y[i-1] > ypeak/2) v_1 = ((ypeak/2 - y[i-1])/((y[i] - y[i-1])/(x[i] - x[i-1]))) + x[i-1];
	if ((y[i] > ypeak/2) && y[i-1] <= ypeak/2) v_2 = ((ypeak/2 - y[i])/((y[i] - y[i-1])/(x[i] - x[i-1]))) + x[i];
      }
      f2 = fopen ("gauss.dat", "a"); fprintf(f2,"%1.3f %1.5e\n", x[i], y[i]); fclose(f2);
      FWHM = v_2 -v_1;
    } fclose(f1);
  
    float ynorm, centre, sig, A_gauss;
    sprintf(systemCall,"./gauss_fit-input.py %1.3e %1.3f %1.3f", ypeak, vpeak, FWHM/2.35); system(systemCall);
    // printf("\nTO TEST TRY %s\n",systemCall);
    //https://ned.ipac.caltech.edu/level5/Leo/Stats2_3.html
    system("/Users/stephencurran/bin/./col2row.csh gauss_fit");
    f3= fopen("gauss_fit_cols.txt","r"); while (fscanf(f3,"%f %f %f", &ynorm, &centre, &sig) !=EOF);
    fclose(f3);

    //if (ynorm < 0) ynorm =-1*ynorm;
    A_gauss = ynorm*sig/0.3989;
    // https://www.physicsforums.com/threads/area-under-gaussian-peak-by-easy-measurements.419285/
    if (strcmp(intact, "y") == 0) tau_g =  -1*log(1 - ynorm/ave);
    else tau_g = ynorm;

    for (q=0 ; q<100 ; q++) {
      xline[q]= 0.9*XSTART + ((1.1*XEND - 0.9*XSTART)/100)*q;
      yline[q] = ave + (ynorm*exp(-1*pow(xline[q] - centre,2)/(2*sig*sig)));
    }
    cpgsci(12);cpgsls(1); cpgslw(3); cpgline(100,xline,yline);cpgslw(3);
    
    // printf("%1.2f to  %1.2f km/s RMS1 = %1.2e and %1.2f to  %1.2f km/s RMS2 = %1.2e => RMS = %1.2e\n", X1, XSTART, RMS1, XEND, X2, RMS2, RMS);
    //printf("FWZM of line is %1.1f km/s, continuum = %1.3f rms = %1.2e Jy area = %1.3e inttau = %1.3f cf. %1.3f\n     peak = %1.3e (%1.3f)  optdepth = %1.3e at %1.1f v = %1.2f (peak) cf. %1.2f (FWHM) km/s\n", XEND - XSTART, ave, RMS, A, inttau, A/ave, ypeak, ave+ypeak,optdepth, vpeak, v_ave, v_2 - v_1);
    //printf("\nFrom Gaussian fit, peak = %1.3e (%1.3f) tau = %1.3e at %1.1f and FWHM = %1.1f km/s\n", ynorm, ynorm+ave, -1*log(1 - ynorm/ave), centre, 2.35*sig);

    cpgsci(2); cpgsls(3); cpgmove(v_ave,(v_ave*-999));cpgdraw(v_ave,(v_ave*999)); // vertical line for offset
 
    cpgsci(3); cpgslw(5); cpgsls(4); cpgmove(v_ave-0.5*A/ypeak,  ypeak/2 +ave); cpgdraw(v_ave+0.5*A/ypeak,  ypeak/2 +ave);
    //cpgsci(8);  cpgsls(4); cpgmove(vpeak - v_ave, ave+ ypeak/2);cpgdraw(vpeak + v_ave, ave+ ypeak/2);
    cpgsci(8);  cpgslw(3); cpgsls(2); cpgmove(v_1, ave+ ypeak/2);cpgdraw(v_2, ave+ ypeak/2);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    if (inttau < 0) inttau = -1*inttau;
    if (ypeak < 0) ypeak =-1*ypeak;

    printf("====================================================================================================\n");
    printf("RMS = %1.3e over %1.0f channels and RMS = %1.3e over %1.0f channels => %1.2e Jy per %1.2f km/s\n", RMS1, (XSTART-X1)/res, RMS2, (X2-XEND)/res, RMS, res);
    printf("\nContinuum level is %1.3g +/- %1.3g, area = %1.3f +/- %1.3f cf. %1.3f (Gaussian)\n", ave, sigma, A, A*fracerr, A_gauss);
 
    if (strcmp(emabs, "e") != 0)   printf("\n      => inttau = %1.3f  +/- %1.3f km/s N = %1.3e +/- %1.3e logN = %1.2f [for HI]\n", inttau, inttau*fracerr, 1.823e18*inttau, 1.823e18*inttau*fracerr, log10(1.823e18*inttau)); //frac error
   
    printf("\nFWZM = %1.2f, 'FWHM' = %1.2f (half-peak width - orange), %1.2f (area/peak - green), %1.2f (Gaussian) km/s\n", XEND - XSTART, v_2 - v_1,inttau/optdepth,  2.35*sig);

    printf("\nPeak = %1.3e (y = %1.3f) tau = %1.3e at %1.1f cf. tau = %1.3e at %1.1f km/s (Gaussian)\n", ypeak, ave+ypeak, optdepth, vpeak, tau_g, centre);
    printf("       \nWeighted offset at %1.2f +/- %1.2f km/s [vertical red line]\n", v_ave, v_ave*fracerr);
    printf("====================================================================================================\n");
 
    printf("Would you like to redefine the points/polynomial fit [y/n]? "); scanf("%s", smooth);
    colour = colour + 1; style = style + 1;
    coeff[9]=coeff[8]=coeff[7]=coeff[6]=coeff[5]=coeff[4]=coeff[3]=coeff[2]=coeff[1]=coeff[0]=0;
  } while(strcmp(smooth,"y")==0);
 
   ////////////////////////////////////// REF FREQ STUFF HERE //////////////////////////////
  float ref_f_z, ref_freq, f_b, f_t, f_peak, f_ave, z_peak, z_ave; 
  char trans[20];

   if (strcmp(units, "h") == 0) { 
      if (obsfreq > 10) ref_f_z = obsfreq;  // MHz already have - may as well use
      else ref_f_z = obsfreq*1000; // for GHz
    }
   else if (strcmp(units, "c") == 0) ref_f_z = obsfreq; 
   else{
      printf("Enter reference frequency [MHz] or redshift [z] (usually optical) [0 if velocity in z * c]\nNOTE! Frequency the safer option if detection offset from band centre,\n e.g. 1301 MHz, rather than z = 0.0971, for 1555-140: "); scanf("%f", &ref_f_z);
    }
      printf("\nRedshift relative to HI [any], OH (1665+1667) main [m], 1612, [12], 1665 [65], 1667 [67], 1720 [20]:"); scanf("%s",trans);

  //ref_f_z = 1301; strcpy(trans,"a"); // for testing

  if (strcmp(trans, "m") == 0 ) ref_freq = 1667.3590; // SET TO 1667 MHZ LINE 1666.3804;  // mean of main lines
  else if (strcmp(trans, "12") == 0 ) ref_freq = 1612.2310;
  else if (strcmp(trans, "65") == 0 ) ref_freq = 1665.40180;
  else if (strcmp(trans, "67") == 0 ) ref_freq = 1667.3590;
  else if (strcmp(trans, "20") == 0 ) ref_freq = 1720.5300;
  else ref_freq = 1420.405752; //obsfreq = 1420.405752/(zopt + 1);

  delta_v = v_ave*fracerr;

  if (ref_f_z != 0){
    if (ref_f_z < 10 && ref_f_z > 0) ref_f_z = ref_freq/(ref_f_z + 1); // redshift to freq
    else  ref_f_z = ref_f_z;  // if > 10 probably already in MHz   GHZ TOTALLY CONFUSES THIS, THINKS IT'S THE REDSHIFT!
    
    f_b =  ref_f_z * pow(((1 - X2/c)/(1 + X2/c)),0.5); 
    f_t =  ref_f_z * pow(((1 - X1/c)/(1 + X1/c)),0.5);
    f_ave = ref_f_z * pow(((1 - v_ave/c)/(1 + v_ave/c)),0.5); 
    d_f_ave = ref_f_z * (1 - pow(((1 - delta_v/c)/(1 + delta_v/c)),0.5));
    f_peak =  ref_f_z * pow(((1 - vpeak/c)/(1 + vpeak/c)),0.5);
    z_ave = ref_freq/f_ave - 1; 
    delta_f = f_ave -d_f_ave;
    d_z_ave = z_ave  - (ref_freq/delta_f  - 1);
    z_peak = ref_freq/f_peak -1;
  }
  else {z_ave = v_ave/c; z_peak = vpeak/c;} // for that stupid v/c definition */
  printf("\n--------------------------------------------------------------------------------------\n");
  printf("For a reference frequency of * %1.4f MHz *  freq_mean =  %1.5f +/- %1.5f\n", ref_f_z,f_ave, d_f_ave);
  printf("                          for redshifted %s  => z_mean =  %1.6f +/- %1.6f\n",trans_text, z_ave, d_z_ave);
  printf("For a reference frequency of * %1.4f MHz *  freq_peak =  %1.5f\n", ref_f_z,  f_peak);
  printf("                          for redshifted %s  => z_peak =  %1.6f\n", trans_text, z_peak);  

  printf("Velocity range of %1.1f to %1.1f km/s is %1.2f to %1.2f MHz (z = %1.5f - %1.5f) \n if defined by frequency or redshift",  X1, X2, f_b, f_t, (ref_freq/f_t) - 1, (ref_freq/f_b) - 1);
  printf("                         or * %1.1f to %1.1f MHz * if defined by velocity\n",  ref_freq/(X2/2.99792458e5 + 1), ref_freq/(X1/2.99792458e5 + 1) );
  printf("--------------------------------------------------------------------------------------------\n");  
  colour = style = 1;
  ///////////////////////////////////  NOW ONTO HARD COPY ////////////////////////////////////
  char plot[10], y_axis[10], plot_label[20],name[50],ps_out[200], ps_out_dev[200],freq_vel_label[50], freq_vel[10]; 
  float YMID, *tau, *z, *x_peak;

  tau = (float *)calloc(npts,sizeof(float)); z = (float *)calloc(npts,sizeof(float)); x_peak = (float *)calloc(npts,sizeof(float));

  if  (strcmp(intact, "y") == 0) YMID = ave;
  else YMID = pubflux;

  X1 = VX1; X2 = VX2; Y1 = VY1; Y2 = VY2; 
  // printf("\n\nYMID = %1.3f, X1 = %1.2f X2 = %1.2f Y1 = %1.2f Y2 = %1.2f\n\n", YMID, X1, X2, Y1, Y2);
 
  char  mantran[20], mol[20],axls[10], ztop[10],centre_v[10], opt[10], axl[10];
  float axl_size;
  int line_thick;
 
   printf("Do you want to produce a plot [y/n?] "); scanf("%s", plot);
   //strcpy(plot,"y"); // for testing
  if (strcmp(plot, "y") == 0){

    if (strcmp(intact, "y") == 0){
      printf("Plot in observed [o], actual [t] optical depth or flux [a]: "); scanf("%s", y_axis);
      if (strcmp(y_axis, "o") == 0) strcpy(plot_label,"tau_obs");
      else if (strcmp(y_axis, "t") == 0) strcpy(plot_label,"tau_act");
      else strcpy(plot_label,"flux");
  }
    else strcpy(plot_label,"fixed");

    printf("Want plot in freq [f] or velocity [v]?: "); scanf("%s",freq_vel);
    if (strcmp(freq_vel, "f") == 0) strcpy(freq_vel_label,"freq");
    else strcpy(freq_vel_label,"vel");

    sprintf(ps_out, "%s-%s_poly%d_%s_%1.0fkms.eps",infile,freq_vel_label,poly_order,plot_label,res);
    sprintf(ps_out_dev, "%s-%s_poly%d_%s_%1.0fkms.eps/cps",infile,freq_vel_label,poly_order,plot_label,res);
    cpgopen(ps_out_dev);

    printf("Text sizes? Large, medium or by hand \nLarge [l], medium [m] or input cpgsch by hand [h]: "); scanf("%s",axls);
    if (strcmp(axls, "l") == 0) axl_size = 1.3;
    else if (strcmp(axls, "m") == 0) axl_size = 1;
    else{printf("Factor of current text size? "); scanf("%f",&axl_size);}
    //if (axl_size > 2)
    if (strcmp(axls, "l") == 0) line_thick = 4;
    else line_thick = 3;

    cpgpap(8.,0.75); cpgslw(3); cpgsch(axl_size*text);
     
    /*     j = 1; // first character */
    /*     while (infile[j]){ */
    /*       if (infile[j] == '_') strncpy(name,infile,j);  // cutting file name at first underscore */
    /*       ++j; // need this or hangs */
    /*     } */
     
       printf("Colour of trace? Black [1], red [2], green [3], blue [4], cyan [5],\nmagenta [6], yellow [7], orange [8], lime [9], turquoise [10],\nviolet [11], purple [12], pink [13], dark grey [14], light grey [15]: "); scanf("%d",&colour);
    
   
    if (strcmp(freq_vel, "f") != 0){
      printf("Use the central defined points as spectrum centre (only if detection) [y/n]: "); scanf("%s",centre_v);
    }

    //////////// REOPEN AND PLOT SPECTRUM ////////////////////////////////////////     
    fspec = fopen ("gauss.dat", "r"); for (i = 0; fscanf(fspec,"%f %f",  &x[i], &y[i]) !=EOF; ++i) {
      if (strcmp(freq_vel, "f") == 0 ) { 
	x[i] = ref_f_z*pow(((1 - x[i]/c)/(1 + x[i]/c)),0.5); 
	z[i] = (ref_freq/x[i]) -1;
      }
      else{
	x[i] = x[i]; //velocity option
	z[i] = (ref_freq/f_peak -1) + x[i]/c; // re-centering reference to peak
	if (strcmp(centre_v, "y") == 0) x[i] = x[i] - v_ave; // centering on zero from mean velocity
      }
  
      x_peak[i] = x[i] - vpeak; // from peak velocity
   
      if (strcmp(intact, "y") == 0){
	if (strcmp(y_axis, "a") == 0) tau[i] = y[i]+YMID;
	if (strcmp(y_axis, "o") == 0){
	  if (pubflux == 0) tau[i] = y[i];
	  else  tau[i] =  y[i]/YMID;
	}
      }
      
      else if (strcmp(y_axis, "t") == 0)  tau[i] =  -1*log(1 - y[i]/YMID); 
      else if (strcmp(emabs, "e") == 0) tau[i] =  y[i]; // now working but have to scales y-axis accordingly
      else tau[i] =  y[i] + YMID;
      // printf("X1 = %1.1f X2 = %1.2f i = %d, x[i] = %1.2f, y[i] = %1.3f, YMID = %1.3f  \n", X1, X2, i, x[i], y[i], YMID); 
    } fclose(fspec);
    //////////////// AXIS RANGES ////////////////
    
    if (strcmp(freq_vel, "f") == 0) {
      X1= ref_f_z*sqrt((1-X1/c)/(1+X1/c)); // relativsitic conversion of velocity to freq
      X2 = ref_f_z*sqrt((1-X2/c)/(1+X2/c)); 
      v_ave = ref_f_z*sqrt((1-v_ave/c)/(1+v_ave/c));
    }
    else{
      if (strcmp(centre_v, "y") == 0) {X1 = X1; X2 = X2; v_ave = v_ave; z_ave = z_ave;}
      else  z_ave = (ref_freq/ref_f_z) -1;
    }
    
    float VEL_OFF; 
    if (X2 - v_ave < v_ave - X1) VEL_OFF = X2 - v_ave; 
    else VEL_OFF = v_ave - X1; 

    if (strcmp(centre_v, "y") == 0) { 
      if (strcmp(freq_vel, "f") != 0) {X1 = -1*VEL_OFF; X2 = VEL_OFF;}
      //  else // NOT WORKING FOR FREQ, BUT IS THERE ANY POINT?
    }
    
    if (strcmp(y_axis, "o") == 0) {Y1 = 1.2*ypeak; Y2 = -0.4*ypeak;}
    else if (strcmp(y_axis, "t") == 0) {Y1 = 1.2*optdepth; Y2 = -0.4*optdepth;}
    printf("tau = %1.4f Y1 = %1.2f Y2 = %1.2f\n", optdepth, Y1, Y2);
      
    ///////////////////// LABELS AND SIZES /////////////////////////////////////
    printf("Enter source name to appear in plot (one word) [FILENAME IS %s]: ", infile); scanf("%s",name); 
    printf("Secondary label - manual, e.g. 'WSRT' [m] or transition, based on above [t]: ");
   
    scanf("%s",mantran);
    if (strcmp(mantran, "m") == 0){printf("Enter string with no breaks [e.g. WSRT]"); scanf("%s",label);}
    else {
      if (strcmp(trans, "m") == 0 ) strcpy(mol,"OH\\d1665,1667");
      else if (strcmp(trans, "12") == 0 ) strcpy(mol,"OH\\d1612");
      else if (strcmp(trans, "65") == 0 ) strcpy(mol,"OH\\d1665");
      else if (strcmp(trans, "67") == 0 ) strcpy(mol,"OH\\d1667");
      else if (strcmp(trans, "20") == 0 ) strcpy(mol,"OH\\d1720");
      else strcpy(mol,"HI\\(2748)21-cm"); //  \\(2748)   is a blank
    }
    
    cpgsch(axl_size*text); cpgslw(line_thick);
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    printf("Labels - labels only [l], axes only [a], both [b] or none[n]?: "); scanf("%s",axl);
    //strcpy(axl,"b"); // for testing
       
    ///////////////////// ADDITIONAL SCALES AND ARROWS ///////////////////////////////////////////////
    float Z1, Z2,x1, x2, y1, y2;
    
    cpgqvp(0,&x1,&x2,&y1,&y2);
    cpgsvp(x1+0.05,x2,y1+0.08,y2); // need these for swin
    cpgsch(axl_size*text);

    printf("Show redshift along top axis [y/n]: "); scanf("%s",ztop);
      
    if (strcmp(ztop, "y") == 0) {
      if (strcmp(freq_vel, "f") == 0) {Z1 = (ref_freq/X1) -1; Z2 = (ref_freq/X2) -1;  cpgswin(Z2,Z1, X2,X1);}
      else  {Z1 = (ref_freq/f_ave -1) - VEL_OFF/c ; Z2 =  (ref_freq/f_ave -1) + VEL_OFF/c; cpgswin(Z1, Z2, X1, X2);}
      cpgbox("cmts",0,0.,"",0.,0.);  // m - n above viewport
  }
    cpgsch(axl_size*text); cpgsci(1); cpgslw(line_thick); 

  ///////////////////////////////////////////////
    
  float temp1, temp2; 
    
  //====================================================================================================================//
  if (strcmp(freq_vel, "f") == 0)  {temp1 = X1; temp2 = X2; X1 = temp2; X2 = temp1;} // easier for labels just to swap?
   cpgswin(X1,X2,Y1,Y2);
     
  if (strcmp(ztop, "y") != 0) cpgbox("bcnts",0,0.,"bcnts",0.,0.);
  else cpgbox("bnts",0,0.,"bcnts",0.,0.);   

  //=====================================================================================================================//

  float fopt, arrow_end, bar_length, fopt_65, velrange, velrange1, velrange2, velrange65_1, velrange65_2;
  
  if (strcmp(freq_vel, "f") == 0) {printf("Show optical redshift arrow [y/n]: "); scanf("%s",opt);} // ahh only for freq

  arrow_end = Y1+(Y2-Y1)/6; bar_length = (arrow_end - Y2)/80;
  
  if (strcmp(opt, "y") == 0) {   
    fopt = ref_f_z; //does them all - for main line want to show secondary line, which is offset at (1667.35900 - 1665.40180)/z+1
    fopt_65 = 1665.40180*fopt/1667.3590; // doing without direct z - fopt should be that for 1667MHz
    
    printf("Velocity range [+/- km/s] - 0 not to show: ");scanf("%f",&velrange);

    velrange1 = fopt*(1-(velrange/c)); velrange2 = fopt*(1+(velrange/c));
    velrange65_1 = fopt_65*(1-(velrange/c)); velrange65_2 = fopt_65*(1+(velrange/c)); // for main option
    cpgslw(4); cpgsci(14); cpgsch(text); cpgsah(1,45,0.5); 
    
    cpgarro(fopt,arrow_end,fopt,Y1);//  HAS TO BE AFTER SETTING OF AXES (cpgswin/cpgenv)
    printf("\nfopt = %1.2f ,arrow_end = %1.4f Y1 = %1.4f  Y2 = %1.4f velrange1 = %1.0f \n", fopt, arrow_end,Y1, Y2, velrange1);
    
    if (velrange > 0){
      cpgmove(velrange1,arrow_end); cpgdraw(velrange2,arrow_end);
      cpgmove(velrange1,arrow_end - bar_length); cpgdraw(velrange1,arrow_end + bar_length); 
      cpgmove(velrange2,arrow_end - bar_length); cpgdraw(velrange2,arrow_end + bar_length);

      if (strcmp(trans, "m") == 0) {
	cpgarro(fopt_65,arrow_end,fopt_65,Y2); // also 1665 MHz for main lines option
	cpgmove(velrange65_1,arrow_end); cpgdraw(velrange65_2,arrow_end); //bar and then  wee error bar bits 
	cpgmove(velrange65_1,arrow_end - bar_length); cpgdraw(velrange65_1,arrow_end + bar_length);
	cpgmove(velrange65_2,arrow_end - bar_length); cpgdraw(velrange65_2,arrow_end + bar_length);
      }
    }
  }
 cpgsch(axl_size*text); cpgsci(1); cpgslw(line_thick);cpgsch(text);

 if (strcmp(y_axis, "o") == 0) sprintf(label,"Optical depth, \\gt\\dobs\\u"); // has to be after env
 else if (strcmp(y_axis, "t") == 0) 	sprintf(label,"Optical depth, \\gt");
 else if (strcmp(y_axis, "a") == 0) sprintf(label,"Flux density, S\\dobs\\u [Jy]");
 else sprintf(label,"Optical depth, \\gt\\dobs\\u");

 if (strcmp(axl, "a") == 0 || strcmp(axl, "b") == 0) cpgmtxt("L",2.5,0.5,0.5, label);
  
 if (strcmp(freq_vel, "f") == 0) sprintf(label,"Observed frequency, \\gn\\dobs\\u [MHz]");
 else sprintf(label,"Velocity offset relative to z = %1.5f [km s\\u-1\\d]",z_ave);
 if (strcmp(axl, "a") == 0 || strcmp(axl, "b") == 0)  cpgmtxt("B",2.5,0.5,0.5, label);
   
 float L_label = X1 + (X2 - X1)/16, R_label = X2 - (X2 - X1)/4, B_label = Y1 + (Y2 - Y1)/3.5, T_label = Y2 - (Y2 - Y1)/12;
    
 if (strcmp(axl, "l") == 0 || strcmp(axl, "b") == 0)  cpgtext(L_label, T_label, name);  cpgtext(R_label, T_label, mol);
    
  //////////////////////////////////////////////////////////////////////////////////////////////////////
 cpgslw(line_thick); cpgsci(colour); cpgline(i,x,tau); // TRACE
 // for (i = 0; i <=numChannels; i++) printf("i = %d, x[%d] = %1.3f y[%d] = %1.3f\n", i, i, x[i], i, tau[i]); //to check

  printf("Postscipt written to %s\n", ps_out);
  sprintf(systemCall,"gv %s &", ps_out); system(systemCall);

}
cpgend();
}

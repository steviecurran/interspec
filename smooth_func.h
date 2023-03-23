void smooth_func(int i, int j, int k, int n, int numChannels,int q, int inc, float sf, float res, float chave, int start, int count, int old_count,float *x, float *y, float *X, float *Y, float *ex, float *why, char systemCall[],char lineInput[], FILE *f1){

  // MOVING AVERAGE, BUT SO SMOOTH BY NON-INTEGER FACTOR
  // the way to do this would be to oversample the data by drawing a line between each point and placing
  // extra points at equal spacings along each of these and then smoothing by  sf*no_points 

  system("rm smooth_out.dat");
   X[j]= x[i] = x[i-1]= Y[j] = y[i] = y[i-1] = X[n] = Y[n] = 0; 

  sf = res/chave;
  if (sf < 0) sf = -1*sf;
  k = inc*sf;   // so k is an integer value to which to smooth 
  if (sf >= k + 0.5) k = k+1;
  else if (sf <= sf - 0.5) k = k -1;
  else if ((sf <= k + 0.5) && (sf >= k - 0.5)) k = k;

  printf("Smoothing by factor of %1.2f\n", (float)k/inc); 
 
  start = 2; // first point not always trustworthy
  for(i = start; i <= numChannels-start+1; i++){
    //  printf("i = %d  x[i] = %1.3f,  y[i] = %1.3f\n", i, x[i],y[i]);
    for (j = 0; j<inc; j++){ // will do as per rekuired resolution
      X[j] =  x[i-1]  + chave*((float)j+1)/inc; 
      Y[j] = (y[i] - y[i-1])*(X[j] - x[i-1])/(x[i] - x[i-1]) + y[i-1];// ekuation of straight line 
      X[n] = X[j]; Y[n] = Y[j];  // NEED TO SUM THE INDEX
      n = (inc*(i-start)) + j +1; 
     
      // printf("n = %d j = %d, X[j] = %1.3f,  Y[j] = %1.3f\n", n, j, X[j],Y[j]);
      old_count = count;
      count = n/k + 1;
      // printf("k = %d i = %d, j = %d n = %d count = %d old_count = %d n/k = %d, X[n] = %1.2f Y[n] = %1.3e\n", k, i, j, n, count, old_count, n/k, X[n-1], Y[n-1]);
      if (count > old_count) q=0;
      else {q++;
	ex[n] = X[n-1]; // USE n RATHER THAN q OTHERWISE SMOOTHING DOESN'T WORK  
	ex[count] = X[n-1];// - res/inc; // LAST BIT CORRECTS IN x BUT y OUT SO JUST USE inc = 100 TO GET NEAR 
	why[n] = Y[n-1]; why[n] = why[n] + why[n-1]; why[count] = why[n]/q; // summing up
        
	if (q==1) {
	  f1 = fopen ("smooth_out.dat", "a"); if (ex[count] == ex[count] && why[count] == why[count]) fprintf(f1, "%d %1.3f %1.5f\n",  count, ex[count], why[count]); fclose(f1);
	}                                          // to get rid of nan
      }
    }
  }
  f1 = fopen ("smooth_out.dat", "r"); // reading in smoothed data 
  for (i = 1; fscanf(f1,"%d %f %f", &i, &x[i-1], &y[i-1]) !=EOF; ++i){
    // printf("%1.2f, %1.2f\n", x[i-1], y[i-1]);
  } fclose(f1);
  cpgslw(3);cpgsci(5); cpgline(i-1,x,y); cpgslw(3);
} 

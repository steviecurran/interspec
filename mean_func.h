void mean_func (int numChannels, int i, char systemCall[], FILE *fsize, float *mean , float VAL[], float *val, float V, FILE *f1, float prob , float E, float e[], float expec, float v[], float var, char sigma_string[], float med){

  ///////////////// determines size of array /////////////////////////////////////////////
  system("wc -l mean_out.txt > /tmp/filesize");  // will create and remove  mean_out.txt each time
  fsize = fopen("/tmp/filesize","r");
  fscanf(fsize,"%d",&numChannels);
  fclose(fsize);
   /////////////////////////////////////////////////////////////////////////////////////////
  
   f1= fopen("mean_out.txt","r"); 
   i=0; // THIS IS VERY IMPORTANT!!!! 
    while(fscanf(f1,"%f",&VAL[i])!=EOF){ 
     prob = 1/(float)numChannels;
      e[i] = (VAL[i]*prob);
      i++;
      // printf("Val = %1.2f  numChannels = %d prob = %1.2e e[%d] = %1.2e\n", VAL[i], numChannels, prob, i, e[i]);
   }fclose(f1); // avoid memory leaks
   expec = 0;  
   for(i=0; i<numChannels; i++){expec += e[i];
  }
   sprintf(sigma_string,"%s","sigma");

f1= fopen("mean_out.txt","r"); 
 i=0; // THIS IS VERY IMPORTANT!!!! 
    while(fscanf(f1,"%f",&VAL[i])!=EOF){ 
     prob = 1/(float)numChannels;
      v[i] =  (VAL[i] - expec)*(VAL[i] - expec)*prob;
      i++;
    } fclose(f1); // avoid memory leaks
var = 0;  
for(i=0; i<numChannels; i++) var += v[i];
//printf("Mean is %1.3f +/- %1.3g [where sigma = %1.3g] for n = %d\n", expec, sqrt(var/numChannels), sqrt(var), numChannels);
 system("rm mean_out.txt"); // for next time
 system("rm mean_values.txt");
 f1 = fopen ("mean_values.txt", "a"); fprintf(f1,"%1.3f %1.3f %d\n", expec, sqrt(var),numChannels); fclose(f1);
 //system("rm mean_values.txt");
 }

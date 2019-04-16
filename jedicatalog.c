#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "fitsio.h"

/*jedicatalog
 *This program creates a realistic catalog of galaxies.
 */

#define BANDHEIGHT 2048
#define BINS_PER_MAG 100
#define DELTA_MAG 0.01
#define DATABASE_PIXSCALE 0.065	/*orginal value 0.03*/
/////////////////////////////////////////////////////////
//nfw number distribution for galaxies
#define PI 3.14159265
#define NUM_INTV 8000	// For NFW (distance) bin: ~ r200 pix
#define MAXLEN 2048    //number of chars in array
/////////////////////////////////////////////////////////

char *help[] = {
    "Creates a realistic catalog of galaxies.",
    "Usage: jedicatalog config_file",
    "Arguments: config_file - a text file with the configuration settings for this program.",
    0};
/*2d array i.e. pointer array*/
/*help? global variable*/

typedef struct {
    float   radius; //r50 radius of the stamp in arcseconds     /*maybe from SExtractor?*/
    float   pixscale;   //pixel scale of the stamp in arcseconds per pixel
    float   magnitude;  //magnitude of the stamp
    char    *name;   //the file name for the postage stamp
} Stamp;


typedef struct {
    char    name[MAXLEN]; //the name of the galaxy
    float   x;  //the x position of the galaxy
    float   y;  //the y position of the galaxy
    float   angle;  //the angle to rotate the galaxy (in degrees)
    float   redshift;   //the redshift of the galaxy
    float   pixscale;   //the pixelscale of the galaxy in arcseconds per pixel
    float   old_mag;    //the original magnitude of the galaxy
    float   old_rad;    //the original radius (pixels) of the galaxy //IMPORTANT!!!!
    float   new_mag;    //the final magnitude of the galaxy
    float   new_rad;    //the final radius of the galaxy //NEW DISTRIBUTION!!!!
    char    stamp_name[MAXLEN]; //the filepath for the postage stamp of this galaxy after its parameters are set
    char    dis_name[MAXLEN];   //the filepath for the postage stamp of this galaxy after it has been distorted
} Galaxy;


int compareInt(const void* p1, const void* p2);	
int compareStamp(const void* p1, const void* p2);
void print_galaxy(Galaxy *gal, FILE *fptr);
float rand_float();	/*float random!*/
float get_mag(float *CDF, int min_mag, int num_bins);
float original_function_nfw(float x);
float rand_nfw(float ofn_xmax, float ofn_ymax, float ofn_xy[NUM_INTV+1][2]);

int main(int argc, char *argv[]){	/*argc number of char, argv 2D array i.e. some sentences*/
    FILE        *config_fp;      //filepointer for the configuration file
    char        buffer[MAXLEN], buffer3[MAXLEN];   //string buffers to read in the config file
    char        *buffer2;        //another string buffer

    //physics settings
    int         single_redshift = 0;    //0 -> get redshift from db; 1-> fixed redshift for all galaxies
    float       fixed_redshift;     //the fixed redshift to use if single_redshift=1
    float       pix_scale;      //pixel scale to use for the simulation in arcseconds per pixel
    long int    nx, ny;         //image size to simulate
    long int    xborder, yborder;   //borders on the image so we don't have galaxies going off the edges
    long int    ngalaxies;      //the total number of galaxies to simulate
    int         min_mag, max_mag;   //the minimum and maximum magnitude
    char        radius_db_folder[MAXLEN], red_db_folder[MAXLEN];    //the folders for the redshift and radius databases
    float       mag_power;      //the power for the power law distribution of galactic magnitudes //ATTENTION!!!!!!!!!!!!!!!!!!!!
    
    //output settings
    char        output_folder[MAXLEN];    //the output folder for postage stamps
    char        prefix[MAXLEN];           //the prefix for all filenames associated with this trial

    //catalog file settings
    char        temp_catalog_file[MAXLEN], catalog_file[MAXLEN];//catalog of galaxy parameters
    FILE        *catalog_fp;
    char        temp_convlist_file[MAXLEN], convlist_file[MAXLEN];//list of postage stamps to be convolved
    FILE        *convlist_fp;
    char        temp_distortedlist_file[MAXLEN], distortedlist_file[MAXLEN];//list of distorted postage stamps
    FILE        *distortedlist_fp;
    char        temp_convolvedlist_file[MAXLEN], convolvedlist_file[MAXLEN];//list of convolved postage stamps
    FILE        *convolvedlist_fp;

    //source image settings
    long int    num_source_images;  //the number of source images
    Stamp       *source_images = NULL;  //array of source images//ATTENTION!
    long int    nimage=0;       //counts how many image filenames have been parsed so far /*parse=analyse?*/

    //cosmic reality databases//IN DATABASE FOLDER
    float       **radius_db;      //radius database: an array of radii for each integer bin of magnitudes, from min_mag to max_mag//ATTENTION!!! THAT IS WHY NO CHANGE???
    float       **red_db;         //redshift database: same format as radius_db
    long int    *radius_db_bin_sizes;   //array that lists the size of each magnitude bin for the radius database
    long int    *red_db_bin_sizes;   //array that lists the size of each magnitude bin for the redshift database
    float       *CDF;   //cumulative distribution function for the power law dist. of magnitudes  //NONONO new distribution!!!!!
    float x1,x2,y1,y2,k,b;      // For linear interpolation


    //parse command line input
    if(argc != 2){
        int line;
        for(line = 0; help[line] != 0; line++){
            fprintf(stderr, "%s\n", help[line]);
        }
        exit(1);
    }
    config_fp = fopen(argv[1],"r");
    if(config_fp == NULL){
        fprintf(stderr,"Error: cannot open config file.");
        exit(1);
    }

    //parse config file
    while(fgets(buffer, MAXLEN, config_fp) != NULL){
        //read in whatever isn't a comment
        if(buffer[0] != '#'){
            buffer2 = strtok(buffer,"#\n");	//breaks string str into a series of tokens using the delimiter delim
            sscanf(buffer2,"%[a-z_]=",buffer3);	/*regular expression*/
            //physics settings	/*can we use switch here?*/
            if(strcmp(buffer3,"pix_scale")==0) sscanf(buffer2,"pix_scale=%f",&pix_scale);
            else if(strcmp(buffer3,"nx")==0) sscanf(buffer2,"nx=%li", &nx);
            else if(strcmp(buffer3,"ny")==0) sscanf(buffer2,"ny=%li", &ny);
            else if(strcmp(buffer3,"x_border")==0) sscanf(buffer2,"x_border=%li", &xborder);
            else if(strcmp(buffer3,"y_border")==0) sscanf(buffer2,"y_border=%li", &yborder);
            else if(strcmp(buffer3,"num_galaxies")==0) sscanf(buffer2,"num_galaxies=%li", &ngalaxies);
            else if(strcmp(buffer3,"min_mag")==0) sscanf(buffer2,"min_mag=%i",&min_mag);
            else if(strcmp(buffer3,"max_mag")==0) sscanf(buffer2,"max_mag=%i",&max_mag);
            else if(strcmp(buffer3,"power")==0) sscanf(buffer2,"power=%f", &mag_power);
            else if(strcmp(buffer3,"single_redshift")==0) sscanf(buffer2,"single_redshift=%i", &single_redshift);
            else if(strcmp(buffer3,"fixed_redshift")==0) sscanf(buffer2,"fixed_redshift=%f", &fixed_redshift);

            //output settings
            //^ means 'not'
            //%[chars] format
            else if(strcmp(buffer3,"output_folder")==0) sscanf(buffer2, "output_folder=\"%[^\"]", output_folder);/*brilliant regular expression!*/
            else if(strcmp(buffer3,"prefix")==0) sscanf(buffer2, "prefix=\"%[^\"]", prefix);

            //catalog file settings
            else if(strcmp(buffer3,"catalog_file")==0) sscanf(buffer2, "catalog_file=\"%[^\"]", temp_catalog_file);

            else if(strcmp(buffer3,"convlist_file")==0) sscanf(buffer2, "convlist_file=\"%[^\"]", temp_convlist_file);
            else if(strcmp(buffer3,"distortedlist_file")==0) sscanf(buffer2, "distortedlist_file=\"%[^\"]", temp_distortedlist_file);
            else if(strcmp(buffer3,"convolvedlist_file")==0) sscanf(buffer2, "convolvedlist_file=\"%[^\"]", temp_convolvedlist_file);

            //cosmic reality databases
            else if(strcmp(buffer3,"radius_db_folder")==0) sscanf(buffer2, "radius_db_folder=\"%[^\"]", radius_db_folder);
            else if(strcmp(buffer3,"red_db_folder")==0) sscanf(buffer2, "red_db_folder=\"%[^\"]", red_db_folder);
            
            //image settings
            else if(strcmp(buffer3,"num_source_images")==0){
                sscanf(buffer2,"num_source_images=%li",&num_source_images);
                source_images = (Stamp*) calloc(num_source_images, sizeof(Stamp));	/*calloc() zero-initializes the buffer, while malloc() leaves the memory uninitialized to save time.*/
                if(source_images == NULL){
                    fprintf(stderr,"Error: could not allocate list of source images.\n");
                    exit(1);
                }
            } 
            //read stamp path and save into source_image[] (in memory rather than disk)
            else if(strcmp(buffer3,"image")==0 && nimage < num_source_images){
                if(source_images == NULL){
                    fprintf(stderr, "Error: the parameter 'num_galaxies' must come before any 'image' parameters\n");
                    exit(1);
                }
                char temp[MAXLEN];//new temp in each loop? occupying more space?
                sscanf(buffer2, "image=\"%[^\"]", temp);
                source_images[nimage].name = (char*) calloc(strlen(temp)+1, sizeof(char));
                if(source_images[nimage].name == NULL){
                    fprintf(stderr,"Error: could not allocate source image struct %li.\n", nimage);
                    exit(1);
                }
                strcpy(source_images[nimage].name,temp);//name is path!
                nimage++;
            }

        }
    }
    
    sprintf(catalog_file, "%s%s%s", output_folder, prefix, temp_catalog_file);
    sprintf(convlist_file, "%s%s%s", output_folder, prefix, temp_convlist_file);
    sprintf(distortedlist_file, "%s%s%s", output_folder, prefix, temp_distortedlist_file);
    sprintf(convolvedlist_file, "%s%s%s", output_folder, prefix, temp_convolvedlist_file);

/*
    //print out what was just read in
    //The printf function is equivalent to fprintf with the argument stdout interposed before the arguments to printf.
    fprintf(stdout,"pixscale: %f\n", pix_scale);
    fprintf(stdout,"nx: %li\n", nx);
    fprintf(stdout,"ny: %li\n", ny);
    fprintf(stdout,"x_border: %li\n", xborder);
    fprintf(stdout,"y_border: %li\n", yborder);
    fprintf(stdout,"num_galaxies: %li\n",ngalaxies);
    fprintf(stdout,"min_mag: %li\n", min_mag);
    fprintf(stdout,"max_mag: %li\n", max_mag);
    fprintf(stdout,"radius_db_folder: %s\n", radius_db_folder);
    fprintf(stdout,"red_db_folder: %s\n", red_db_folder);
    fprintf(stdout,"output_folder: %s\n", output_folder);
    fprintf(stdout,"power: %f\n", mag_power);
    fprintf(stdout,"catalog_file: %s\n", catalog_file);
    fprintf(stdout,"convlist_file: %s\n", convlist_file);
    fprintf(stdout,"distortedlist_file: %s\n", distortedlist_file);
    fprintf(stdout,"convolvedlist_file: %s\n", convolvedlist_file);
    fprintf(stdout,"num_source_images: %li\n", num_source_images);
*/

    //get magnitude and radius of each of the source galaxies with cfitsio
    //READ FROM STAMP FITS FILE
    {
        fitsfile    *sfptr;     //source image fits file pointer
        int         status;     //fits statusAN
        for(nimage = 0; nimage < num_source_images; nimage++){
            if(source_images[nimage].name == NULL){//name is path!
                fprintf(stderr,"Error: image file path %li is null.\n", nimage);
                exit(1);
            }
            char temp[MAXLEN];
            fits_open_file(&sfptr, source_images[nimage].name, READONLY, &status);/*open fits file*/
            
            fits_read_key_str(sfptr, "MAG", temp, NULL, &status);
            sscanf(temp,"%f", &source_images[nimage].magnitude);

            fits_read_key_str(sfptr,"RADIUS", temp, NULL, &status);
            sscanf(temp,"%f", &source_images[nimage].radius);

            fits_read_key_str(sfptr,"PIXSCALE", temp, NULL, &status);
            sscanf(temp, "%f", &source_images[nimage].pixscale);
            source_images[nimage].radius *= source_images[nimage].pixscale;//ATTENTION:RADIUS TO SEC
            //printf("source_images[nimage=%d].radius=%f\n",nimage,source_images[nimage].radius);
            
            fits_close_file(sfptr,&status);
            if(status){
                fits_report_error(stderr,status);
                exit(1);
            }
        }
        qsort(source_images, num_source_images, sizeof(Stamp), compareStamp);
    }
//compareStamp: diff = ((Stamp*)p1)->radius - ((Stamp*)p2)->radius;
//void qsort(void *base, size_t nitems, size_t size, int (*compar)(const void *, const void*)) sorts an array.
/*

    base -- This is the pointer to the first element of the array to be sorted.

    nitems -- This is the number of elements in the array pointed by base.

    size -- This is the size in bytes of each element in the array.

    compar -- This is the function that compares two elements.
*/


    //print out the source images
    //for(nimage = 0; nimage < num_source_images; nimage++)
      //  fprintf(stdout,"image %i:\n\tName:%s\n\tMagnitude: %f\n\tRadius: %f (arcseconds)\n", nimage, source_images[nimage].name, source_images[nimage].magnitude, source_images[nimage].radius);
/*

    //read in the radius and redshift databases
    {
        radius_db = (float **) calloc((max_mag-min_mag+1), sizeof(float*));//2d array: mags-radii
        radius_db_bin_sizes = (long int *) calloc((max_mag-min_mag+1), sizeof(long int));
        red_db_bin_sizes = (long int *) calloc((max_mag-min_mag+1), sizeof(long int));
        if(radius_db == NULL){
            fprintf(stderr,"Error: could not allocate radius database.\n");
            exit(1);
        } 
//ATTENTION!!!
        fprintf(stdout,"Allocated radius db with %i bins.\n", (max_mag-min_mag+1));

        red_db = (float **) calloc((max_mag-min_mag+1), sizeof(float*));
        if(red_db == NULL){
            fprintf(stderr,"Error: could not allocate redshift database.\n");
            exit(1);
        } 
//ATTENTION!!!
        fprintf(stdout,"Allocated redshift db with %i bins.\n", (max_mag-min_mag+1));

        int mag_bin;
        FILE    *radiusfp, *redfp;  //radius and redshift FILEpointers
        for(mag_bin=0; mag_bin <= (max_mag-min_mag); mag_bin++){//start from 0, so amount add one
            char radius_filename[MAXLEN], redshift_filename[MAXLEN];    //names of the redshift and radius data files
            //using convention that the file names are radius_db_folder/n.dat where n from min_mag to max_mag
            sprintf(radius_filename,"%s%i.dat",radius_db_folder,mag_bin+min_mag);//add path
            sprintf(redshift_filename,"%s%i.dat",red_db_folder,mag_bin+min_mag);

            //import RADIUS bin
            radiusfp = fopen(radius_filename,"r");
            long int     nradius = 0;    //radius counter
            while(fgets(buffer, MAXLEN, radiusfp) != NULL){
                nradius++;                
            }
            radius_db[mag_bin] = (float *) calloc(nradius, sizeof(float));
            radius_db_bin_sizes[mag_bin] = nradius;
            if(radius_db[mag_bin]==NULL){
                fprintf(stderr,"Error: could not allocate radius_db bin %i.\n", mag_bin);
                exit(1);
            }
//ATTENTION!!!
            fprintf(stdout,"Allocated radius db bin %i of size %i.\n", mag_bin, nradius);
            rewind(radiusfp);
            nradius = 0;
            while(fgets(buffer, MAXLEN, radiusfp) != NULL){
                if(sscanf(buffer,"%f",&radius_db[mag_bin][nradius]) != 1){
                    fprintf(stderr,"Error");
                    exit(1);
                }
                radius_db[mag_bin][nradius] *= DATABASE_PIXSCALE;//pay attention here!
                //fprintf(stdout,"bin: %i\tnradius: %li\tsize: %li\tvalue: %.3f\tbuffer: %s\n", mag_bin+min_mag, nradius, size, radius_db[mag_bin][nradius], buffer);
                nradius++;
            }
            (radiusfp);
            //sort the radii
            qsort(radius_db[mag_bin], nradius, sizeof(float), compareInt);



            //import redshift bin

            redfp = fopen(redshift_filename,"r");
            long int     nred = 0;    //radius counter
            while(fgets(buffer, MAXLEN, redfp) != NULL){
                nred++;                
            }
            red_db[mag_bin] = (float *) calloc(nred, sizeof(float));
            red_db_bin_sizes[mag_bin] = nred;
            if(radius_db[mag_bin]==NULL){
                fprintf(stderr,"Error: could not allocate radius_db bin %i.\n", mag_bin);
                exit(1);
            }
            fprintf(stdout,"Allocated red db bin %i of size %i.\n", mag_bin, nred);
            rewind(redfp);
            nred = 0;
            while(fgets(buffer, MAXLEN, redfp) != NULL){
                if(sscanf(buffer,"%f",&red_db[mag_bin][nred]) != 1){
                    fprintf(stderr,"Error");
                    exit(1);
                }
                //fprintf(stdout,"bin: %i\tnradius: %li\tvalue: %.3f\tbuffer: %s\n", mag_bin+min_mag, nred, red_db[mag_bin][nred], buffer);
                nred++;
            }
            (redfp);

        }

        for(mag_bin = 0; mag_bin <=(max_mag-min_mag); mag_bin++)
           fprintf(stdout,"mag bin: %i\t redshift size: %li\n", mag_bin, red_db_bin_sizes[mag_bin]); 
    
    }

*/

    //make and print the catalog
    long int g;  //galaxy counter
    fprintf(stdout,"catalog_file: %s\n", catalog_file);
    catalog_fp = fopen(catalog_file, "w");
    if(catalog_fp == NULL){
        fprintf(stderr, "Error: could not open catalog file.\n");
        exit(1);
    }
    distortedlist_fp = fopen(distortedlist_file, "w");


// NFW number distribution
    
    float rou, theta;    
    
    float ofn_xmax = 2.0;    // This is for c=2; if c=4, this would be 4.0 
    float ofn_ymax = original_function_nfw(ofn_xmax);
    float ofn_xy[NUM_INTV+1][2]={0};//original function of nfw function x-y pairs (bins); x=0, y=0
    int counter=0;
    
    for(counter=1;counter<=NUM_INTV;counter++){
	ofn_xy[counter][0] = counter*1.0/NUM_INTV*ofn_xmax;
        ofn_xy[counter][1] = original_function_nfw(ofn_xy[counter][0]);
    }

// Radius
    FILE *rd_datafile;//,*fp;
    int counter_new = 0;//counter
    int index;
    //char oneline[MAXLEN] = {'\0'};//MAXLEN microsoft notepad limit
    float **cdf_mat = NULL;//cdf matrix
    int row,col;//row,col of cdf matrix
    int mat_rows,mat_cols;//number of rows & cols
    float r_len,d_len;//lengths of r, d's interval
    //float *radius_sec = NULL;
    //float *radius_density = NULL;
    float unif_rand;//uniform in [0,1]
    //float radius;
    char *pstr;//pointer to string 
    
    printf("read in cdf_matrix from file\n");
    rd_datafile = fopen("./simdatabase/clash/mat_dr_cdf_whole.txt","r");
    if(rd_datafile==NULL){
        printf("ERROR cannot open ./simdatabase/clash/mat_dr_cdf_whole.txt.\n");
        exit(1);
    }
    //read in parameters in header
    mat_rows = 0;
    mat_cols = 0;
    while(fgets(buffer,MAXLEN,rd_datafile)!=NULL){
        //printf("buffer=%s\n",buffer);
        if(buffer[0]!='#'){
            buffer2 = strtok(buffer,"#\n");//printf("buffer2=%s\n",buffer2);
            sscanf(buffer2,"%[a-z_]=",buffer3);//printf("buffer3=%s\n",buffer3);
            if(strcmp(buffer3,"mat_rows")==0) sscanf(buffer2,"mat_rows=%i",&mat_rows);
            else if(strcmp(buffer3,"mat_cols")==0) sscanf(buffer2,"mat_cols=%i",&mat_cols);
            else if(strcmp(buffer3,"r_len")==0) sscanf(buffer2,"r_len=%f",&r_len);
            else if(strcmp(buffer3,"d_len")==0) sscanf(buffer2,"d_len=%f",&d_len);
            else if(strcmp(buffer3,"cdf_matrix")==0){
                printf("we need to stop right now. break.\n");
                break;//normally it has started to read the first line of matrix
            }
        }
    } 
    printf("mat_rows=%d,mat_cols=%d,r_len=%f,d_len=%f\n",mat_rows,mat_cols,r_len,d_len);
    //read in cdf matrix
    //first build up a matrix (2d pointer)
    //rows
    cdf_mat = (float **)malloc(sizeof(float*)*mat_rows);
    if(cdf_mat==NULL||mat_rows==0){
        printf("ERROR cannot malloc cdf_mat.\n");
        exit(1);
    }
    //columns
    for(row=0;row<mat_rows;row++){
        cdf_mat[row] = (float *)malloc(sizeof(float)*mat_cols); 
        if(cdf_mat[row]==NULL||mat_cols==0){
            printf("ERROR cannot malloc cdf_mat[row=%d].\n",row);
            exit(1);
        }
    }
    //then put data from file to here (memory)
    //put first line of matrix into memory first
    for(row=0;row<mat_rows;row++){
        fgets(buffer,MAXLEN,rd_datafile);
        pstr = strtok(buffer," \n");
        //printf("read #%d row of cdf matrix.\n",row);
        for(col=0;col<mat_cols;col++){
            cdf_mat[row][col] = atof(pstr);
            //printf("read in cdf_mat[%d][%d] = %f\n",row,col,cdf_mat[row][col]);
            pstr = strtok(NULL, " \n");
        }
    }
    int distance_maxrow = mat_rows;    

    //now matrix is in memory
/*    
    //check it!
    for(row=0;row<mat_rows;row++){
        for(col=0;col<mat_cols;col++){
            printf("%f ",cdf_mat[row][col]);
        }
        printf("\n");
    }
*/    
/*
    rewind(rd_datafile);//after counting the amount of bins, go back and start reading data
    counter_new = 0;
    while(1){
        fgets(oneline,MAXLEN,rd_datafile);//A terminating null character is automatically appended after the characters copied to str
        //printf("%s",oneline);
        sscanf(oneline,"%f,%f",&radius_sec[counter_new],&radius_density[counter_new]);//i.e. x and f(x
        radius_sec[counter_new]*=DATABASE_PIXSCALE;//ATTENTION:RADIUS PIX TO SEC
        //printf("radius_sec=%f   radius_density=%f\n",radius_sec[counter_new],radius_density[counter_new]);
        counter_new++;
        if(feof(rd_datafile)) break;
    }
*/
/*
    for(counter_new=0;counter_new<3;counter_new++){
        printf("%f    %f\n",radius_sec[counter_new],radius_density[counter_new]);
    }
*/
/*
    srand(time(NULL));
    //fp = fopen("output.txt","w");
    for(counter_new=0;counter_new<1000;counter_new++){
        unif_rand = rand()*1.0/RAND_MAX;
        index = 0;
        while(unif_rand>radius_density[index]){
            index++;
        }
        radius = radius_sec[index];//output: radius_sec[index]
        //fprintf(fp,"%f\n",radius_sec[index]);
    }
*/
/*
    free(radius_sec);
    free(radius_density);
*/
    fclose(rd_datafile);



// Magnitude
    FILE *rm_datafile;//,*fp;
    //int counter_new = 0;//counter
    //int index;
    //char oneline[MAXLEN] = {'\0'};//MAXLEN microsoft notepad limit
    float **cdf_mag_mat = NULL;//cdf matrix
    //int row,col;//row,col of cdf matrix
    //int mat_rows,mat_cols;//number of rows & cols
    //float r_len,d_len;//lengths of r, d's interval
    float mag_len,mag_zero;
    //float *radius_sec = NULL;
    //float *radius_density = NULL;
    //float unif_rand;//uniform in [0,1]
    //float radius;
    //char *pstr;//pointer to string 
    
    printf("read in cdf_mag_matrix from file\n");
    rm_datafile = fopen("./simdatabase/clash/mat_rm_cdf.txt","r");
    if(rm_datafile==NULL){
        printf("ERROR cannot open ./simdatabase/clash/mat_rm_cdf.txt.\n");
        exit(1);
    }
    //read in parameters in header
    mat_rows = 0;
    mat_cols = 0;
    while(fgets(buffer,MAXLEN,rm_datafile)!=NULL){
        //printf("buffer=%s\n",buffer);
        if(buffer[0]!='#'){
            buffer2 = strtok(buffer,"#\n");//printf("buffer2=%s\n",buffer2);
            sscanf(buffer2,"%[a-z_]=",buffer3);//printf("buffer3=%s\n",buffer3);
            if(strcmp(buffer3,"mat_rows")==0) sscanf(buffer2,"mat_rows=%i",&mat_rows);
            else if(strcmp(buffer3,"mat_cols")==0) sscanf(buffer2,"mat_cols=%i",&mat_cols);
            //else if(strcmp(buffer3,"r_len")==0) sscanf(buffer2,"r_len=%f",&r_len);
            else if(strcmp(buffer3,"mag_len")==0) sscanf(buffer2,"mag_len=%f",&mag_len);
            else if(strcmp(buffer3,"mag_zero")==0) sscanf(buffer2,"mag_zero=%f",&mag_zero);
            else if(strcmp(buffer3,"cdf_matrix")==0){
                printf("we need to stop reading header right now. break. we will read matrix next.\n");
                break;//normally it has started to read the first line of matrix
            }
        }
    } 
    printf("mat_rows=%d,mat_cols=%d,r_len=%f, mag_len=%f\n",mat_rows,mat_cols,r_len,mag_len);
    //read in cdf matrix
    //first build up a matrix (2d pointer)
    //rows
    cdf_mag_mat = (float **)malloc(sizeof(float*)*mat_rows);
    if(cdf_mag_mat==NULL||mat_rows==0){
        printf("ERROR cannot malloc cdf_mag_mat.\n");
        exit(1);
    }
    //columns
    for(row=0;row<mat_rows;row++){
        cdf_mag_mat[row] = (float *)malloc(sizeof(float)*mat_cols); 
        if(cdf_mag_mat[row]==NULL||mat_cols==0){
            printf("ERROR cannot malloc cdf_mag_mat[row=%d].\n",row);
            exit(1);
        }
    }
    //then put data from file to here (memory)
    //put first line of matrix into memory first
    for(row=0;row<mat_rows;row++){
        fgets(buffer,MAXLEN,rm_datafile);
        pstr = strtok(buffer," \n");//split with blank or return;get first part
        //printf("read #%d row of cdf mag matrix.\n",row);
        for(col=0;col<mat_cols;col++){
            cdf_mag_mat[row][col] = atof(pstr);
            //printf("read in cdf_mag_mat[%d][%d] = %f\n",row,col,cdf_mag_mat[row][col]);
            pstr = strtok(NULL, " \n");
        }
    }

    int radius_maxrow = mat_rows;
/*
    //check it!
    for(row=0;row<mat_rows;row++){
        for(col=0;col<mat_cols;col++){
            printf("%f ",cdf_mag_mat[row][col]);
        }
        printf("\n");
    }
*/    
    
    fclose(rm_datafile);



// Galaxy    
    for(g = 0; g < ngalaxies; g++){
        srand(time(NULL)+rand());   //the code runs too fast to use the time in miliseconds as the seed
        Galaxy gal;
        
        //set galaxy position and angle 
	//distance to center satisfying NFW distribution
	
        rou = rand_nfw(ofn_xmax,ofn_ymax,ofn_xy)*4006.644477; // In pixel; this is for c=2; for c=4 it's 2003.322238
        //actually nx (also xborder) here should be the smaller one between nx and ny
	theta = rand_float()*2*PI;        
        gal.x = rou * cos(theta) + nx/2;//in pix
        gal.y = rou * sin(theta) + ny/2;
//        gal.x = xborder + rand_float()*(nx-2*xborder);
//        gal.y = yborder + rand_float()*(ny-2*yborder);
        //if they go outside... checked still same distribution inside image's square region
        if(gal.x<=xborder||gal.x>=(nx-xborder)||gal.y<=yborder||gal.y>=(ny-yborder)){
//            g--; 
            continue;       // It won't write into catalog
        }

//        rou = sqrt((gal.x-nx/2)*(gal.x-nx/2) + (gal.y-ny/2)*(gal.y-ny/2));
        gal.angle = rand_float()*360;
////////////////////////////////////////////////////////////////////////////////////////////////////////



//NEW METHOD
//printf("largest_radius=%f\n",source_images[num_source_images-1].radius);
//choose a stamp
//while(1){}//stop point
        printf("-------catalog start g=%li---------\n",g);
        //generate a valid (not too large) radius obeying distribution
        //first find the distance of this gal: rou
        //coarse-grain it too index/row in cdf_matrix: int(rou/d_len), which starts from 0
        //its corresponding cdf should be cdf_mat[int(rou/d_len)][index]
        //then gal.new_rad = (index+1)*r_len //add one to avoid radius=0
        row = (int)(rou/d_len);//distance
        if(row>=distance_maxrow){ 
            row = distance_maxrow-1;
            printf("The row is larger than max row in matrix. Use value at max row\n"); 
        }
        printf("rou=%f, choose distance row=%d\n", rou, row);
        while(1){
            unif_rand = rand_float();
            index = 0;
            while(unif_rand>=cdf_mat[row][index]){
                index++;
            }
            if (index==0){
                x1 = 0;
                y1 = 0;                
            }
            else{
                x1 = (index-1)+0.5;
                y1 = cdf_mat[row][index-1];
            }
            x2 = index+0.5;
            y2 = cdf_mat[row][index];
            if (x1==x2) {
                printf("Error: x1==x2 at mat dr sampling!\n");
                exit(1);
            }
            k = (y2-y1)/(x2-x1);
            b = (y1*x2-x1*y2)/(x2-x1);
            gal.new_rad = (unif_rand-b)/k*r_len*DATABASE_PIXSCALE;//printf("index=%d,radius_sec(new_radius)=%f,largest=%f\n",index,radius_sec[index],source_images[num_source_images-1].radius);
            //IN SEC!!!
            if(gal.new_rad<source_images[num_source_images-1].radius){ 
                //smaller than "the largest old radius (qsorted)"
                break;
            }
        }
        printf("index=%d, x1=%f, x2=%f, y1=%f, y2=%f, k=%f, b=%f\n",index,x1,x2,y1,y2,k,b);
        printf("unif_rand=%f, (unif_rand-b)/k=%f\n", unif_rand, (unif_rand-b)/k);
        printf("gal.new_rad=%f\n",gal.new_rad);




        //after we know the radius, choose a stamp (with proper galaxy size)
        int image_num;
        while(1){
            image_num = rand() % (num_source_images);//randomly pick one stamp            
            gal.old_rad = source_images[image_num].radius;
            //if old radius is smaller than new radius, which means expand the stamp, 
            //we should give up this stamp and look for another stamp
            //still same radius
            if(gal.old_rad>gal.new_rad){    
                printf("gal.old_rad=%f>gal.new_rad=%f   find it!\n",gal.old_rad,gal.new_rad);
                sprintf(gal.name,"%s", source_images[image_num].name);
                gal.pixscale = source_images[image_num].pixscale;
                gal.old_mag = source_images[image_num].magnitude;
                break;
            }
        }
        //choose mag
        //the radius we need is (index+1)*r_len
        row = (int)(gal.new_rad/(r_len*DATABASE_PIXSCALE)); //this corresponds to the radius we need 
        printf("choose mag row=%d\n", row);
        unif_rand = rand_float();
        index = 0;
        while(unif_rand>=cdf_mag_mat[row][index]){
            index++;
        }
        if(index==0){
            x1 = 0;
            y1 = 0;
        }
        else{
            x1 = (index-1)+0.5;
            y1 = cdf_mag_mat[row][index-1];
        }
        x2 = index+0.5;
        y2 = cdf_mag_mat[row][index];
        if (x1==x2){ 
            printf("Error: x1==x2 at mat rm sampling!\n");
            exit(1);
        }
        k = (y2-y1)/(x2-x1);
        b = (y1*x2-x1*y2)/(x2-x1);
        gal.new_mag = (unif_rand-b)/k*mag_len+mag_zero;

        printf("index=%d, x1=%f, x2=%f, y1=%f, y2=%f, k=%f, b=%f\n",index,x1,x2,y1,y2,k,b);
        printf("unif_rand=%f, (unif_rand-b)/k=%f\n", unif_rand, (unif_rand-b)/k);
        printf("gal.old_mag=%f;gal.new_mag=%f\n",gal.old_mag,gal.new_mag);
        printf("new index %d (might not be same due to position)\n", (int)((gal.new_mag-mag_zero)/mag_len));


       
printf("\nCHECK THIS!!!------------------1\n");
///////////////////////////////////////////////////////////////////////////////////////////////////////
/*
        //double-check, in case there's a bug
        if(gal.new_rad > gal.old_rad){
            fprintf(stderr,"Error: selected illegal new radius.");
            exit(1);
        }
*/
        //get a redshift from the database
        //WE DON'T NEED TO CONSIDER THIS
        gal.redshift = fixed_redshift;
/*
        if(single_redshift == 0){
            int redshift_num = rand() % red_db_bin_sizes[mag-min_mag];
            gal.redshift = red_db[mag-min_mag][redshift_num];
        } else if(single_redshift == 1){
            gal.redshift = fixed_redshift;
        } else {
            fprintf(stderr, "Error: single_redshift must be either 0 or 1.");
        }
*/


        //set galaxy filepaths
        sprintf(gal.stamp_name, "%sstamp_%li/stamp_%li.fits.gz", output_folder, g/1000, g);

		//can't write to .gz files with FITSIO, so this is disabled so jedigrid can work
        //sprintf(gal.dis_name, "%sdistorted_%i/distorted_%i.fits.gz", output_folder, g/1000, g);
        sprintf(gal.dis_name, "%sdistorted_%li/distorted_%li.fits", output_folder, g/1000, g);
        

        
        //print_galaxy(&gal, stdout);
        print_galaxy(&gal, catalog_fp);//WRITE this galaxy to our catalog, still in "for loop"

        fprintf(distortedlist_fp, "%s\n", gal.dis_name);
    }

printf("\nCHECK THIS!!!-----------------2\n");
    convolvedlist_fp = fopen(convolvedlist_file, "w");
    int band, num_bands = ny/BANDHEIGHT;
    char conv_name[MAXLEN];
    for(band = 0; band < num_bands; band++){   
        sprintf(conv_name, "%sconvolved/convolved_band_%i.fits", output_folder, band);
        fprintf(convolvedlist_fp, "%s\n", conv_name);
    }

    fclose(config_fp);
    fclose(catalog_fp);
    fclose(distortedlist_fp);
    fclose(convolvedlist_fp);
printf("\nCHECK THIS!!!-----------------3\n");
//free pointers
//    free(radius_sec);
//    free(radius_density);

    for(row=0;row<distance_maxrow;row++){
        free(cdf_mat[row]);
    }
    for(row=0;row<radius_maxrow;row++){
        free(cdf_mag_mat[row]);
    }
    free(cdf_mat);
    free(cdf_mag_mat);


    for(nimage = 0; nimage < num_source_images; nimage++){
        free(source_images[nimage].name);
    }
    free(source_images);
/*
    free(radius_db);
    free(radius_db_bin_sizes);
    free(red_db_bin_sizes);
    free(red_db);
*/
    //mag-cdf(from old code, my radius cdf is a function)
printf("\nCHECK THIS!!!-----------------4\n");

    //fprintf(stdout,"Success! Catalog made with %i entries.\n", ngalaxies);
    return 0;
}

//comparator for sorting a list of integers
int compareInt(const void* p1, const void* p2){
    return (*(int*)p1 - *(int*)p2); 
}

//comparator for sorting a list of stamps
int compareStamp(const void* p1, const void* p2){
    float diff = ((Stamp*)p1)->radius - ((Stamp*)p2)->radius;
    if(diff > 0) return 1;
    else if(diff <0) return -1;
    else return 0; 
}

//prints out the given galaxy to the given file stream
void print_galaxy(Galaxy *gal, FILE *fptr){
    fprintf(fptr, "%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s\n",gal->name, gal->x, gal->y, gal->angle, gal->redshift, gal->pixscale, gal->old_mag, gal->old_rad, gal->new_mag, gal->new_rad, gal->stamp_name, gal->dis_name);
    //fprintf(fptr,"%f\t%f\n",gal->x, gal->y);
}

//returns a random float in [0,1)
float rand_float(){
    return ((float) rand())/RAND_MAX;
}

float get_mag(float *CDF, int min_mag, int num_bins){/*random mag following distribution*/
    float r = rand_float();    
    int j=0;
    while(CDF[j]<= r && r < num_bins) 
        j++;
    //fprintf(stdout,"min_mag: %i\tnum_bins: %i\trand: %f\tj: %li\tCDF[j]: %f\n",min_mag, num_bins, r, j, CDF[j]); 
    return (float) min_mag + j*DELTA_MAG;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//nfw distribution definition (number distribution)
float original_function_nfw(float x){
    float y = 0;
    //the original function of nfw function (without const factor delta_c, rou_c, r_s; x=R/r_s)
    if (x>=0 && x<1) y = log(x/2) + 2/(sqrt(1-x*x))*atanh(sqrt((1-x)/(1+x)));
    else if (x==1) y = log(1.0/2.0) + 1;
    else if (x>1) y = log(x/2) + 2/(sqrt(x*x-1))*atan(sqrt((x-1)/(1+x)));
    else y = 0;//or printf("error!\n");exit(1);
    return y;
}

//F^-1(x)
float rand_nfw(float ofn_xmax, float ofn_ymax, float ofn_xy[NUM_INTV+1][2]) {
    float unf_rand_y; //= rand_float()*ofn_ymax;
    float output = 1;
    int counter = 0;
    float xx1,xx2,yy1,yy2,kk,bb;      // For linear interpolation
    // original_function_nfw*normalize<=1
    unf_rand_y = rand_float()/1.653987; //this is for c=2; if c=4, it's 0.967602
    for(counter=0;counter<(NUM_INTV+1);counter++){
	if(unf_rand_y>=ofn_xy[counter][1]&&unf_rand_y<=ofn_xy[counter+1][1]) break;
    }
    // Here we suppose the the value at slot is the mid value
    // We also suppose x=0 is at the middle of 0th slot
    yy2 = ofn_xy[counter+1][1];
    yy1 = ofn_xy[counter][1];
    xx2 = ofn_xy[counter+1][0];
    xx1 = ofn_xy[counter][0];
    kk = (yy2-yy1)/(xx2-xx1);
    bb = (yy1*xx2-xx1*yy2)/(xx2-xx1);
    output = (unf_rand_y-bb)/kk;
    //output = (ofn_xy[counter][0]+ofn_xy[counter+1][0])/2.0;//tested!
    return output;	//rand_nfw() generates a number satisfying nfw proj-distribution with specific c
}

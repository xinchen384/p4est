#include <stdio.h>
#include <stdlib.h>

/*

  int *** test_array = new int **[len];
        for (int i = 0; i < len; i++){
                test_array[i] = new int *[len];
                for (int j = 0; j < len; j++)
                        test_array[i][j] = new int[len];
        }


        std::ifstream infile("./stl-files/myvolume");
        std::string line;
        for (int i = 0; i < len; i++)
        for (int j = 0; j < len; j++){
                if (std::getline(infile, line)){
                        if (line.size() != len){
                                cout << "reading file to 3D array, length error!!!" << endl;
                                break;
                        }
                        for (int k = 0; k != line.size(); k++){
                                test_array[i][j][k] = (int)line[k]-48;

                        }
                }
        }

*/

int*** test_array;

void main(){
  printf(" hello world !!!\n" );
  int i, j, k;
  int len;
  FILE * fp;
  char *line = NULL;
  ssize_t read;
  size_t slen = 0;
  int data;
  
  fp = fopen("/home/xin/Dropbox/3d-printing-paper/vs-projects/Project1/Project1/stl-files/myvolume", "r");
  if (fp == NULL){
    printf(" can not open file ! \n");
    exit(1);
  }
  read = getline(&line, &slen, fp);
  len = (int)read - 2;
  test_array = (int***) malloc(len * sizeof(int **));
  for (i=0; i<len; i++){
    test_array[i] = (int **) malloc ( len* sizeof(int*) );
    for (j=0; j<len; j++)
      test_array[i][j] = (int *) malloc ( len* sizeof(int) );
  }
  printf ("length %d line: %s \n", len, line);
  fseek(fp, 0, SEEK_SET);

  for ( i=0; i<len; i++ )
  for ( j=0; j<len; j++ ) {
    if ((read = getline(&line, &slen, fp)) != -1){
      //printf ("length %zu line: %s \n", read, line);
      for ( k=0; k<len; k++ ){
        test_array[i][j][k] = (int)line[k] - 48;
      }
    } 
  }
  fclose(fp);

  i = 0;
  for (j=0; j<len; j++){
  for (k=0; k<len; k++)
    printf("%d", test_array[i][j][k]);
  printf("\n");
  }

}

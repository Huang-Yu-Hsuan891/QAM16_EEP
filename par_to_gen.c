#include <stdio.h>
#include <stdlib.h>

int main()  {
    int n, rc;
    int e;
    int n1, rc1;
    int i, j, m, k;
    int **Cmask;
    int Cmaskrow;
    int Cmaskcolumn;
    int **H;
    int Hrow;
    int Hcolumn;
    int **H1;
    int H1row,H1column;

    FILE *fpr;
    fpr = fopen("paritymatrixCmask.txt","r");
    fscanf(fpr,"%d",&e);
    fscanf(fpr, "%d",&n);
    fscanf(fpr, "%d",&rc);
    fscanf(fpr, "%d",&rc1);
    fscanf(fpr, "%d",&n1);
    
    n1 = n1 * 2;

    Cmaskrow = rc1;
    Cmaskcolumn = n1;
    Cmask = (int **)malloc(Cmaskrow * sizeof (int *));
    for (i = 0; i < Cmaskrow; i++) Cmask[i] = (int *)malloc(Cmaskcolumn * sizeof(int));

    for (i = 0; i < rc1; i++) {
        for (j = 0; j < n1; j++) {
            fscanf(fpr, "%d", &Cmask[i][j]);
            //printf("%d ",Cmask[i][j]);
        }
        //printf("\n");
    }
    fclose(fpr);

    Hrow = rc;
    Hcolumn = n;
    H = (int **)malloc(Hrow * sizeof(int *));
    for (i = 0; i < Hrow; i++) H[i] = (int *)malloc(Hcolumn * sizeof(int));
    
    H1row = rc;
    H1column = n;
    H1 = (int **)malloc(H1row * sizeof(int *));
    for (i = 0; i < H1row; i++) H1[i] = (int *)malloc(H1column * sizeof(int));

    for (i = 0; i < Hrow; i++) {
        for (j  = 0; j < Hcolumn; j++) {
            H[i][j] = 0;
        }
    }
    n1 = n1 / 2;
    int temp;

    printf("n1 = %d\n", n1);
    for (i = 0; i < rc1; i++) {
        for (j = 0; j < n1; j++) {
            if (Cmask[i][2 * j] == 0) {
                printf("no\n");
                continue;
            }
            if (Cmask[i][2 * j] == 1) {
                printf("yes\n");
                for (m = 0; m < e; m++) {
                   H[e * i + m][e * j + (m + Cmask[i][2 * j + 1]) % e] = 1; 
                }    
            }
               
        }
    }
    int num = 0;
    printf("Hrow = %d; Hcolumn = %d\n",  Hrow, Hcolumn);
    for (i = 0; i < Hrow; i++) {
        for (j  = 0; j < Hcolumn; j++) {
            printf("%d ", H[i][j]);
            H1[i][j] = H[i][j];
            if (H[i][j] == 1) num++;
        }
        printf("\n");
    }
    printf("num = %d\n", num);

    // GAUSS JORDAN METHOD
    int temprow = 0;    // for store pivot row 
    int stop = 0;
    int temptrans;
    int s = 0;
    printf("gauss jordan form!\n");
    for (i = 0; i < rc + s&& temprow < rc/*rc - s*//*n*/; i++) {                        // find pivot
        stop = 0;
        //printf("i = %d\n", i);
        for (j = temprow/*0*/; j < rc && stop == 0; j++) {      // where the pivot in the row
            if (H[j][i]  ==  1 && j == temprow) {
                temprow = j;
                //printf("j = %d temprow = %d \n", j, temprow);
                stop = 1;
            }
            else if (H[j][i]  ==  1 && j != temprow) {
                stop = 1;
                //printf("j = %d temprow = %d\n", j, temprow);
                for (k = 0; k < n; k++) {
                    temptrans = H[j][k];
                    H[j][k] = H[temprow][k];
                    H[temprow][k] = temptrans;
                }
            }
        }
        if(j == rc && stop == 0) {
            //i++;
            s++;
            printf("s = %d\n",s);
            continue;
        }
        // Eliminate the column have 1
        printf(" i = %d\n", i);
        for (j = 0; j < rc; j++) {
            if (j == temprow) continue;                    // skip the original row because the row no need Elimination
            else {
                //printf("h\n");
                if (H[j][i] == 1) {
                    //printf("h\n");
                    for (k = 0; k < n; k++) {
                        H[j][k] = (H[j][k] + H[temprow][k]) % 2;
                    }
                }
            }
        }
        temprow++;
        //printf("temprow = %d\n", temprow);
    }
    // bug is need to replace the position

    for (i = 0; i < rc; i++) {
        for (j = 0; j < n; j++) {
            printf("%d ",H[i][j]); 
            //if (H[i][j] == 1) printf("%d ",j + 1);//;printf("%d ",H[i][j]); 
        }
        printf("\n");
    }
    FILE *outfp1;
    outfp1 = fopen("paritycheckmatrix.txt","w");
    fprintf(outfp1,"%d ",n);
    fprintf(outfp1,"%d ",rc);
    fprintf(outfp1,"\n");
    fprintf(outfp1,"3 ");
    fprintf(outfp1,"6 ");
    fprintf(outfp1,"\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < rc; j++) {
            //fprintf(outfp1,"%d ",H[i][j]);
            if (H1[j][i] == 1) fprintf(outfp1,"%d ",j+1);
        }
        fprintf(outfp1,"\n");
    }
    for (int i = 0; i < rc; i++) {
        for (int j = 0; j < n; j++) {
            //fprintf(outfp1,"%d ",H[i][j]);
            if (H1[i][j] == 1) fprintf(outfp1,"%d ",j+1);
        }
        fprintf(outfp1,"\n");
    }

    
    fclose(outfp1);

    int a = 0;
    int a1 = 0;
    int c[1320];
    int c1[1320];
    m = 0;
    int **Hsyst;
    int Hsystrow = 1320;
    int Hsystcolumn = 1320;
    int **Hsyst1;
    int Hsyst1row = 1320;
    int Hsyst1column = 1320;
    int **Gsys;
    int Gsysrow = 1320;
    int Gsyscolumn = 2640;
    int **G;
    int Grow = 1320;
    int Gcolumn = 2640;

    Hsyst = (int **)malloc(Hsystrow * sizeof(int *));
    for (i = 0; i < Hsystrow; i++) Hsyst[i] = (int *)malloc(Hsystcolumn * sizeof(int));
    Hsyst1 = (int **)malloc(Hsyst1row * sizeof(int *));
    for (i = 0; i < Hsyst1row; i++) Hsyst1[i] = (int *)malloc(Hsyst1column * sizeof(int));
    Gsys = (int **) malloc (Gsysrow * sizeof(int *));
    for (i = 0; i < Gsysrow; i++) Gsys[i] = (int *) malloc(Gsyscolumn *sizeof(int));
    G = (int **)malloc(Grow * sizeof(int *));
    for (i = 0; i < Grow; i++) G[i] = (int *)malloc(Gcolumn * sizeof(int));

    for (j = 0; j < n; j++) {
        m = 0;
        for (i = 0; i < rc; i++) {
            if (H[i][j] == 0) m = m + 1;
        }
        if (m == rc - 1) {
            for (i = 0; i < rc; i++) {
                Hsyst[i][a] = H[i][j]; 
            }
            c[a] = j;
            a++;
        }
        else {
            for (i = 0; i < rc; i++) {
                Hsyst1[i][a1] = H[i][j];
            }
            c1[a1] = j;
            a1++;
        }
    }
    for (i = 0; i < 1320; i++) printf("%d ", c[i]);
    printf("\n");
    for (i = 0; i < 1320; i++) printf("%d ",c1[i]);
    printf("\n");
    m = 0;
    for (j = 0; j < n; j++) {
        //m = 0;
        if (j < 1320) {
            for (i = 0; i < rc; i++) {
                Gsys[i][j] = Hsyst[i][j];
            }
        } else {
            /*for (i = 0; i < rc; i++) {*/
                for (k = rc; k < n; k++) {
                    Gsys[m][k] = Hsyst1[k-1320][j-1320];
                }  
                m++;
            //}
            
        }
        
    }
    for (i = 0; i < rc; i++) {
        for (j = 0; j < n; j++) {
            printf("%d ", Gsys[i][j]);
        }
    }
    //printf("nn\n");
    a = 0;
    a1 = 0;
    for (j = 0; j < n; j++) {
        if (c[a] == j) {
            printf("nna\n");
            for (i = 0; i < rc; i++) {
                G[i][j] = Gsys[i][1320 + a];
            }
            a++;
        }
        if (c1[a1] == j) {
            printf("nna1\n");
            for (i = 0; i < rc; i++) {
                G[i][j] = Gsys[i][a1];
            }
            a1++;
        }

    }
    printf("y\n");
    for (i = 0; i < rc; i++) {
        for (j = 0; j < n; j++) {
            if (G[i][j] == 1) printf("%d ", j + 1);
        }
        printf("\n");
    }

    FILE *outfp;
    outfp = fopen("generator1.txt","w");
    for (int i = 0; i < rc; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(outfp,"%d ",G[i][j]);
        }
        fprintf(outfp,"\n");
    }
    fclose(outfp);
    printf("\n\n");
    int temp1;
    for (i = 0; i < rc; i++) {
        temp1 = 0;
        for (j = 0; j < n; j++) {
            temp1 += (G[i][j] * H1[i][j]);
            temp1 = temp1 % 2;   
        }
        printf("%d ", temp1);
    }

    for (i = 0; i < H1row; i++) free(H1[i]);
    free(H1);
    for (i = 0; i < Hrow; i++) free(H[i]);
    free(H);
    for (i = 0; i < Cmaskrow; i++) free(Cmask[i]);
    free(Cmask);
    for (i = 0; i < Hsystrow; i++) free(Hsyst[i]);
    free(Hsyst);
    for (i = 0; i < Hsyst1row; i++) free(Hsyst1[i]);
    free(Hsyst1);
    for (i = 0; i < Gsysrow; i++) free(Gsys[i]);
    free(Gsys);
    for (i = 0; i < Grow; i++) free(G[i]);
    free(G);


    return 0;
}
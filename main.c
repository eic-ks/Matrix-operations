#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<cblas.h>
#include<lapacke.h>
#define N 10
typedef struct{
    int len;
    double *p;

}vec;
typedef struct{
    int row;
    int col;
    double *p;
}mat;
typedef struct{
    int num_matrix;
    mat matrix_set[N];
}Data;
Data Read(char *file);
void Case_sum(Data *data);
void Case_product(Data *data);
void Case_eigen(Data *data);
void Case_show(Data *data);

//コマンドラインでデータを渡す
int main(int argc, char *argv[]){
    char operation[100];
    Data data = Read(argv[1]);
    while( 1 ){
        printf("####################\n"
               "# choose operation #\n"
               "####################\n");
        printf("operation:【sum,product,eigen,show,end】\n");

        scanf(" %s",operation);
        if( strcmp(operation,"sum") == 0 ){
            Case_sum(&data);
            continue;

        }else if( strcmp(operation,"product") == 0){
            Case_product(&data);
            continue;

        }else if( strcmp(operation,"eigen") == 0 ){
            Case_eigen(&data);
            continue;
        }else if( strcmp(operation,"show") == 0 ){
            Case_show(&data);
            continue;
        }else if( strcmp(operation,"end") == 0 ){
            break;
        }else{
            printf("[%s] is Not available\n\n",operation);
            continue;
        }
    }

    for(int i = 0; i < data.num_matrix; i++){
        mat *matrix = &data.matrix_set[i];
        free(matrix->p);
    }
    return 0;
}
int Checkdata(char *file,int *row_data,int *col_data){
    //データ行列の数、それぞれの行数と列数を数える。
    FILE *fp;
    char c;
    fp = fopen(file,"r");
    if(fp == NULL){
        printf("input error(Checkdata)\n");
        exit(1);
    }
    int linebreak,count_row,count_col;
    int num_matrix = 0;
    while( 1 ){
        count_row = 0;
        count_col = 0;
        linebreak = 0;
        while( 1 ){
            c = fgetc(fp);
            if( c == ' ' ){
                count_col++;
            }
            if( c == '\n' ){
                linebreak++;
                if( linebreak > 1 ){
                    if(  count_col % count_row == 0 ){
                        count_col /= count_row;
                        row_data[num_matrix] = count_row;
                        col_data[num_matrix] = count_col;
                        num_matrix++;
                        break;
                    }else{
                        printf("%d番目の行列の入力形式に誤りがあります。\n",num_matrix+1);
                        exit(1);
                    }
                }
                count_row++;
                count_col++;
            }else if( c != '\r' ){
                linebreak = 0;
            }
            if( c == EOF ){
                if(  count_col % count_row == 0 ){
                    count_col /= count_row;
                    row_data[num_matrix] = count_row;
                    col_data[num_matrix] = count_col;
                    num_matrix++;
                    return num_matrix;
                }else{
                    printf("%d番目の行列の入力形式に誤りがあります。\n",num_matrix+1);
                    printf("%d %d\n",count_row,count_col);
                    exit(1);
                }
            } 
        }
        //制限数以上の読み取りは行わない
        if( num_matrix == N ){
            return num_matrix;
        }
    }
    fclose(fp);
}
Data Read(char *file){
    int i,j,k;
    int row_data[N];
    int col_data[N];
    int num_matrix = Checkdata(file,row_data,col_data);
    Data data = {num_matrix};

    for(i = 0; i < num_matrix; i++){
        mat *matrix = &data.matrix_set[i];
        matrix->row = row_data[i];
        matrix->col = col_data[i];
        matrix->p = (double*)calloc(matrix->row * matrix->col,sizeof(double));
    }
    //行列の読み込み
    FILE *fp;
    fp = fopen(file,"r");
    if(fp == NULL){
        printf("input error(Read)\n");
        exit(1);
    }
    for(k = 0; k < num_matrix; k++){
        mat *matrix = &data.matrix_set[k];
        for(i = 0; i < matrix->row; i++){
            for(j = 0; j < matrix->col; j++){
                fscanf(fp,"%lf",&matrix->p[i + matrix->row*j]);
            }
        }
    }
    fclose(fp);

    return data;
}
bool Existmatrix(Data *data,int num){
        if( num > data->num_matrix || num < 0 ){
            printf("番号%dに対応する行列がありません\n",num);
            return false;
        }else{
            return true;
        }
}
void Show(Data *data,int n){
    if( !Existmatrix(data,n) ){
        return;
    }
    int i,j;
    //n番目の行列(1~)
    mat *matrix = &data->matrix_set[n-1];
    for(i = 0; i < matrix->row; i++){
        for(j = 0; j < matrix->col; j++){
            printf("%10.3lf ",matrix->p[i + matrix->row*j]);
        }
        printf("\n");
    }
    return;
}
void Summatrix(Data *data,int numA,int numB,double a,double b){
    if( !(Existmatrix(data,numA) * Existmatrix(data,numB)) ){
        return;
    }
    mat *A = &data->matrix_set[numA-1];
    mat *B = &data->matrix_set[numB-1];
    printf("行列A\n");
    Show(data,numA);
    printf("行列B\n");
    Show(data,numB);
    if( A->row != B->row || A->col != B->col ){
        printf("行列の和A + Bは定義できません\n");
        return;
    }
    mat C;
    C.row = A->col;
    C.col = A->col;
    C.p = (double*)calloc(C.row * C.col,sizeof(double));

    int i,j;
    for(i = 0; i < C.row; i++){
        for(j = 0; j < C.col; j++){
            if( i == j ){
                C.p[i + C.row*j] = 1;
            }
        }
    }
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
                A->row,C.col,A->col,a,A->p,A->row,C.p,C.row,b,B->p,B->row);

    printf("%.3lfA + %.3lfB = \n",a,b);
    for(i = 0; i < B->row; i++){
        for(j = 0; j < B->col; j++){
            printf("%10.3lf ",B->p[i + B->row*j]);
        }
        printf("\n");
    }
    return;
}
void Productmatrix(Data *data,int numA,int numB,double a){
    if( !(Existmatrix(data,numA) * Existmatrix(data,numB)) ){
        return;
    }
    mat *A = &data->matrix_set[numA-1];
    mat *B = &data->matrix_set[numB-1];
    printf("行列A\n");
    Show(data,numA);
    printf("行列B\n");
    Show(data,numB);
    if( A->col != B->row ){
        printf("行列の積A * Bは定義できません\n");
        return;
    }
    mat C;
    C.row = A->row;
    C.col = B->col;
    C.p = (double*)calloc(C.row * C.col,sizeof(double));
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
                A->row,B->col,A->col,a,A->p,A->row,B->p,B->row,0,C.p,C.row);

    int i,j;
    printf("%.3lfA * B = \n",a);
    for(i = 0; i < C.row; i++){
        for(j = 0; j < C.col; j++){
            printf("%10.3lf ",C.p[i + C.row*j]);
        }
        printf("\n");
    }
    return;
}
void Showeigenvalue(Data *data,int num,char mode){
    if( !Existmatrix(data,num) ){
        return;
    }
    mat *original = &data->matrix_set[num-1];
    mat copymat;
    mat *A = &copymat;
    A->row = original->row;
    A->col = original->col;
    A->p = (double*)calloc(A->row * A->col,sizeof(double));
    int i,j;
    for(i = 0; i < A->row; i++){
        for(j = 0; j < A->col; j++){
            A->p[i + A->row*j] = original->p[i + original->row*j];
        }
    }
    printf("行列\n");
    Show(data,num);
    if( A->row != A->col ){
        printf("正方行列ではありません\n");
        return;
    }else {
        bool is_symmetric = true;
        for(i = 0; i < A->row; i++){
            for(j = 0; j <= i; j++){
                if( A->p[i + A->row*j] != A->p[j + A->row*i] ){
                    is_symmetric = false;
                    break;
                }
            }
            if( !is_symmetric ){
                printf("対称行列ではありません\n");
                return;
            }
        }
    }
    double w[A->row];
    LAPACKE_dsyev(LAPACK_COL_MAJOR,mode,'U',A->row,A->p,A->row,w);
    printf("固有値:\n");
    for(i = 0; i < A->row; i++){
        printf("%.3lf\n",w[i]);
    }
    if( mode == 'V' ){
        printf("\n固有ベクトル:\n");
        for(i = 0; i < A->row; i++){
            for(j = 0; j < A->col; j++){
                printf("%10.3lf ",A->p[i + A->row*j]);
            }
            printf("\n");
        }
    }
    return;
}
void Case_sum(Data *data){
    int numA,numB;
    double alpha,beta;
    printf("choose matrices in αA + βB\n");
    printf("matrix:【1~%d】\n",data->num_matrix);
    printf("行列A:");
    scanf("%d",&numA);
    printf("行列B:");
    scanf("%d",&numB);
    printf("スカラーα:");
    scanf("%lf",&alpha);  
    printf("スカラーβ:");
    scanf("%lf",&beta); 
    Summatrix(data,numA,numB,alpha,beta);
    printf("\n");
}
void Case_product(Data *data){
    int numA,numB;
    double alpha;
    printf("choose matrices in αA * B\n");
    printf("matrix:【1~%d】\n",data->num_matrix);
    printf("行列A:");
    scanf("%d",&numA);
    printf("行列B:");
    scanf("%d",&numB);
    printf("スカラーα:");
    scanf("%lf",&alpha);  
    Productmatrix(data,numA,numB,alpha);
    printf("\n");
}
void Case_eigen(Data *data){
    int num;
    char mode;
    printf("choose matrices\n");
    printf("matrix:【1~%d】\n",data->num_matrix);
    printf("行列:");
    scanf("%d",&num);
    printf("固有ベクトルの表示【y/n】:\n");
    scanf(" %c",&mode);
    if( mode == 'y' ){
        Showeigenvalue(data,num,'V');

    }else{
        Showeigenvalue(data,num,'N');
    }
    printf("\n");
}
void Case_show(Data *data){
    int num;
    printf("choose matrix:\n");
    printf("matrices:【1~%d】\n",data->num_matrix);
    scanf("%d",&num);
    printf("%d番目の行列\n",num);
    Show(data,num);
    printf("\n");
}

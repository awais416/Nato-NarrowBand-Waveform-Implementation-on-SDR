#include<stdio.h>
#include <math.h>
#include<conio.h>
int *decimal_to_bin(int decimalNumber[],int binary_data[],int binaryNumber[])
{
    int b,i=0,a;
    int quotient;
    int counter;
    int counter1=0;
    int counter2=0;
    for (b=0;b<4;b++)
    {
        quotient = decimalNumber[b];
        counter=0;
        while(quotient!=0){

         binaryNumber[i]= quotient % 2;
         i++;
         quotient = quotient / 2;
         counter++;
        }
        counter2=i-1;
        for (a=0;a<counter;a++)
        {
            binary_data[counter1]=binaryNumber[counter2];
            counter1++;
            counter2--;
        }
        printf("\n");

    }
    return binary_data;
}
int *threshhold(int binary_data[],int binary_data1[],int final_size,int size_binary)
{
    int j;
    for (j=0;j<final_size;j++)
    {
        binary_data1[j]=binary_data[j];
        if(j>=size_binary)
        {

            binary_data1[j]=0;
        }
        printf("%d",binary_data1[j]);
    }
    return binary_data1;
}
int *overlapping(int grouped_bits[][3],int encoded_bits[],int nn[],int newary[],int groupeds,int maryc[])
{
    int i1,i5,x,a9,j1=0;
    for (i1=0;i1<groupeds;i1++)
    {
        encoded_bits[i1]=grouped_bits[i1][j1]*100+grouped_bits[i1][j1+1]*10+grouped_bits[i1][j1+2]*1;

        x=0;
        for(i5=0;encoded_bits[i1]!=0;++i5)
        {
            a9=encoded_bits[i1]%10;
            x=(a9)*(pow(2,i5))+x;
            encoded_bits[i1]=encoded_bits[i1]/10;
        }
            nn[i1]=x;
            newary[i1]=maryc[nn[i1]];

    }
    return newary;
}
int *upsample(int upsampled_bits[],int newary[],int groupeds,int oversampled)
{
    int b4,a5,a6=0;
    for (b4=0;b4<groupeds;b4++)
    {
        for (a5=0;a5<oversampled;a5++)
        {
            if (a5==0)
            upsampled_bits[a6]=newary[b4];
            else
            upsampled_bits[a6]=0;

            a6=a6+1;
        }

    }
    return upsampled_bits;
}
int *upsample1(int upsampled_bits1[],int upsampled_bits[],int upsampling_size1)
{
    int a4,a6=0,a8=0;
     for (a4=0;a4<upsampling_size1;a4++)
    {
        if (a6<16)
        {
            upsampled_bits1[a6]=0;
        }
        else
        {
            upsampled_bits1[a6]=upsampled_bits[a8];
            a8=a8+1;
        }
        a6=a6+1;
    }
    return upsampled_bits1;
}
int *upsample2(int upsampled_bits2[],int upsampled_bits[],int upsampling_size1,int upsampling_size)
{
    int a4,a6=0,a8=0;
    for (a4=0;a4<upsampling_size1;a4++)
    {
        if (a6>=upsampling_size)
        {
            upsampled_bits2[a6]=0;
        }
        else
        {
            upsampled_bits2[a6]=upsampled_bits[a8];
            a8=a8+1;
        }
        a6=a6+1;
    }
    return upsampled_bits2;
}
int *final_phase(int overlapped_bits[],int phase_smoother[],int pulse_shaper[],int final_phase[],int upsampling_size1)
{
final_phase[0]=0;
    int a4,k2=0,k3=1;

     for (a4=0;a4<upsampling_size1;a4++)
    {
        pulse_shaper[a4]=overlapped_bits[a4]*phase_smoother[k2];
        if(k2==32)
        {
            k2=-1;
        }
        k2=k2+1;
        final_phase[k3]=pulse_shaper[a4]+final_phase[k3-1];
        k3=k3+1;
    }
    return final_phase;
}
int *final_block(int upsampling_size1,int final_phase [],float final_phase1[])
{
    int a4;
    for (a4=0;a4<upsampling_size1;a4++)
    {
    final_phase1[a4]=2*3.14*(.0909)*final_phase[a4];
    }
    return final_phase1;
}
float *final_block_real(int upsampling_size1,float final_phase1[],float real_bits[])
{
    int a4;
    for (a4=0;a4<upsampling_size1;a4++)
    {
            real_bits[a4]=cos(final_phase1[a4]);
    }
    return real_bits;
}
float *final_block_complex(int upsampling_size1,float final_phase1 [],float complex_bits[])
{
    int a4;
    for (a4=0;a4<upsampling_size1;a4++)
    {
            complex_bits[a4]=sin(final_phase1[a4]);
    }
    return complex_bits;
}
int main(){
    int size=5;
    int *binary_data_p;
    int size_binary=12;
    int binary_data[size_binary];
    int *decimalto_bin_p;
    int decimalNumber[4]={ 9, 7, 4, 2};
    int *binaryNumber_p;
    int binaryNumber[25];
    int mod_three=0;
    mod_three=size_binary%3;
    int final_size= (size_binary )+ (3-mod_three);
    int *binary_data1_p;
    int binary_data1[final_size];
    int *maryc_p;
    int maryc[8]={7,5,3,1,-1,-3,-5,-7};
    int groupingsize=3;
    int groupeds=final_size/groupingsize;
    int *grouped_bits_p;
    int grouped_bits[groupeds][groupingsize];
    int encoded_bits_p;
    int encoded_bits[groupeds];
    int phase_smoother[48]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int *nn_p;
    int nn[groupeds];
    int dec;
    int i2=0;
    int d;
    int bn,dn,j4=1,remainder;
    int n;
    int *newary_p;
    int newary[groupeds];
    //....................................//
    int oversampled = 16;
    int upsampling_size=oversampled*groupeds;
    int upsampling_size1=upsampling_size+16;
    int upsampled_bits[upsampling_size];
    int *upsampled_bits_p;
    int *upsampled_bits1_p;
    int upsampled_bits1[upsampling_size1];
    int *upsampled_bits2_p;
    int upsampled_bits2[upsampling_size1];
    int a8=0,a5=0;
    int *overlapped_bits_p;
    int overlapped_bits[upsampling_size1];
    int phase_length=32;
    int *final_phase_p;
    int final_phase[upsampling_size1];
    int *pulse_shaper_p;
    int pulse_shaper[upsampling_size1];
    float *final_phase1_p;
    float final_phase1[upsampling_size1];
    float *real_bits_p;
    float real_bits[upsampling_size1];
    float *complex_bits_p;
    float complex_bits[upsampling_size1];

    decimalto_bin_p=decimal_to_bin(decimalNumber,binary_data,binaryNumber);
    binary_data1_p=threshhold(binary_data,binary_data1,final_size,size_binary);
    int a1,a2,a3;

    a3=0;
    for(a1=0;a1<groupeds;a1++)
    {
        for(a2=0;a2<groupingsize;a2++)
        {
            grouped_bits[a1][a2]=binary_data1[a3];
            a3=a3+1;
        }
    }
    newary_p=overlapping(grouped_bits,encoded_bits,nn,newary,groupeds,maryc);

    printf("\n");
    int a4;
    for (a4=0;a4<groupeds;a4++)
    {
           printf("%d",newary[a4]);
    }


        int a6;
    a6=0;
    a8=0;
    upsampled_bits_p=upsample( upsampled_bits,newary,groupeds,oversampled);
    upsampled_bits1_p=upsample1( upsampled_bits1, upsampled_bits, upsampling_size1);
    upsampled_bits2_p=upsample2( upsampled_bits2,upsampled_bits,upsampling_size1,upsampling_size);
    printf("\n upsampling");
    for (a4=0;a4<upsampling_size;a4++)
    {
        printf("%d",upsampled_bits[a4]);
    }

    printf("\n upsampled1 \n");
    for (a4=0;a4<upsampling_size1;a4++)
    {
        printf("%d",upsampled_bits1[a4]);
    }

    printf("\n upsampled2 \n");
    for (a4=0;a4<upsampling_size1;a4++)
    {
        printf("%d",upsampled_bits2[a4]);
    }
    for (a4=0;a4<upsampling_size1;a4++)
    {
        overlapped_bits[a4]=upsampled_bits1[a4]+upsampled_bits2[a4];
    }
    printf("\n overlapped \n");
    for (a4=0;a4<upsampling_size1;a4++)
    {
        printf("%d",overlapped_bits[a4]);
    }

     printf("\n phasesmoothing \n");
    for (a4=0;a4<48;a4++)
    {
        printf("%d",phase_smoother[a4]);
    }
    final_phase[0]=0;

        final_phase1_p=final_block( upsampling_size1, final_phase , final_phase1);
        real_bits_p=final_block_real(upsampling_size1, final_phase1, real_bits);
        complex_bits_p=final_block_complex( upsampling_size1, final_phase1 , complex_bits);

    printf("\n final_phase1");
    for (a4=0;a4<upsampling_size1;a4++)
    {
        printf("\n%d",final_phase[a4]);
    }

    return 0;
}

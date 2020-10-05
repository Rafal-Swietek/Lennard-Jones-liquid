#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <unistd.h>

using namespace std;

const double PI = 3.141592653;
const double r_c = 2.5;
const double x_step = 0.01;
const double y_step = 0.01;
const double r_step = 0.05;
const double dr = 0.1;
const double ro = 0.4;
const double T = 0.7;
const unsigned int L = 40;
const unsigned int mcs = 230000;
const unsigned int sim_time = 30000;
const unsigned int check_time = 100;
constexpr unsigned int N = static_cast<int>(ro*L*L);
const double Rmax = static_cast<double>(L/2.-1.0);
const unsigned int r_size = static_cast<int>( Rmax/r_step + 1 );
const unsigned int mcs_size = static_cast<int>( (mcs-sim_time)/check_time +1);
unsigned int counter;

fstream matrix;
fstream Prob_liq;

void showmatrix2(double* x, double* y)
{
    matrix.open("matrix.txt", ios::out);
    for(unsigned int i=1;i<=N; i++){
        matrix<<x[i]<<" "<<y[i]<<endl;
    }
    matrix.close();
}

void Liquid2(double* x,double* y)
{
    double d;
    unsigned int M, i;
    unsigned int k = static_cast<int>(N/2.);
    x[1] = 0.5; y[1] = 0.5; M = (L-1)*(L-1);
    if(N>=M){
       d = 1.0;
    }
    else{
        d = 1.0 + 1./(ro*L);
    }
    cout<<"d = "<<d<<endl;
    for(i=2; i<k;i++){
        x[i] = x[i-1] + d; y[i] = y[i-1];
        if(x[i]>=L-0.5){
            y[i] += d;
            x[i] = 0.5;
        }
    }
    x[k] = L - 0.5;
    y[k] = L - 0.5;
    for(i=k+1; i<=N; i++){
        x[i] = x[i-1] - d; y[i] = y[i-1];
        if(x[i]<=0.5){
            y[i] -= d;
            x[i] = L - 0.5;
        }
    }
}

int main()
{

    //srand(time(NULL));
    cout<<"N = "<<N<<endl;
    double* x;
    double* y;
    x = new double[N+1];
    y = new double[N+1];
    Liquid2(x,y);
    showmatrix2(x,y);
    unsigned int i,j,k;
    double xnew,ynew,dx,dy,r_ij2,Rnew2;
    double dU,dx_new,dy_new;
    for(k=1;k<=mcs;k++){
        for(i=1;i<=N;i++){
            dU = 0;
            xnew = x[i] + (rand()/RAND_MAX-0.5)*x_step;
            ynew = y[i] + (rand()/RAND_MAX-0.5)*y_step;
            if(xnew > L) xnew = xnew - L;
            if(ynew > L) ynew = ynew - L;
            if(xnew < 0.0) xnew += L;
            if(ynew < 0.0) ynew += L;
            for(j=1;j<=N;j++){
                if(i!=j){
                  dx = fabs(x[i] - x[j]);
                  dy = fabs(y[i] - y[j]);
                  if(dx > L/2.) dx = L - dx;
                  if(dy > L/2.) dy = L - dy;
                  r_ij2 = dx*dx + dy*dy;

                  dx_new = fabs( xnew-x[j] );
                  dy_new = fabs( ynew-y[j] );
                  if(dx_new > L/2.) dx_new = L - dx_new;
                  if(dy_new > L/2.) dy_new = L - dy_new;
                  Rnew2 = dx_new*dx_new + dy_new*dy_new;
                  if(Rnew2<=1.0) break;
                  if(r_ij2 <= r_c && Rnew2 <= r_c){
                    dU += 4.*(pow(Rnew2,-6) - pow(r_ij2,-6));
                  }
                }
            }
            if(Rnew2>1.0 && rand()/RAND_MAX <= min(1.0,exp(-dU/T))){
                x[i] = xnew; y[i] = ynew;
            }
        }
        showmatrix2(x,y);
        cout<<k<<endl;
    }
    showmatrix2(x,y);


    return 0;
}



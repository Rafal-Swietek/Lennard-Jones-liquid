#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

const double PI = 3.141592653;
const double r_c = 2.5;
const double x_step = 0.01;
const double y_step = 0.01;
const double r_step = 0.05;
const double dr = 0.1;
const double ro = 0.4;
const double T = 0.7;
const unsigned int L = 35;
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
fstream time_log;

void showmatrix(double* x, double* y)
{
    for(unsigned int i=1;i<=N; i++){
        matrix<<x[i]<<" "<<y[i]<<endl;
    }
}

void Liquid(double* x,double* y)
{
    double d;
    unsigned int M, i;
    unsigned int k = static_cast<int>(N/2.);
    x[1] = 0.5; y[1] = 0.5; M = (L-1)*(L-1);
    if(N>=M){
       d = 1.0;
    }
    else{
        d = 1.0 + 1./L;
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

double Probability(double* x, double* y,double r)
{
    unsigned int i,j,M;
    double deltaR2;
    double P = 0;
    if(r<=dr) return 0;
    for(i=1;i<=N;i++){
        M = 0;
        for(j=i+1;j<=N;j++){
            deltaR2 = sqrt( pow(x[i]-x[j],2) + pow(y[i]-y[j],2) );
            if(deltaR2>=r-dr && deltaR2<= r) M += 1;
        }
        P += M/(2.*PI*r*dr*ro);
    }
    P = P/(2*N +0.0);
    return P;
}

void save_prob(double tab[r_size][mcs_size])
{
    double P;
    for(unsigned int i=0; i<r_size; i++){
        P = 0.0;
        for(unsigned int j=0; j<mcs_size; j++){
            P += tab[i][j];
        }
        P = P/mcs_size;
        cout<<P<<" "<<i*r_step<<endl;
        Prob_liq<<P<<" "<<i*r_step<<endl;
    }
}

int main()
{
    long double time_start = time (NULL);
    Prob_liq.open("prob_gas.txt", ios::out);
    matrix.open("matrix.txt", ios::out);
    time_log.open("time_log.txt", ios::out);
    cout<<"N = "<<N<<endl;
    double* x;
    double* y;
    x = new double[N+1];
    y = new double[N+1];
    cout<<"$";
    Liquid(x,y);
    showmatrix(x,y);
    srand(time(NULL));
    cout<<"Probability"<<"             "<<"  r"<<endl;
    unsigned int i,j,k,m;
    double xnew,ynew,dx,dy,r_ij2,Rnew2;
    double dU,dx_new,dy_new;
    double P_r[r_size][mcs_size];
    counter = 0;
    for(k=1;k<=mcs;k++){
        for(i=1;i<=N;i++){
            dU = 0;
            cout<<"#";
            xnew = x[i] + (rand()/RAND_MAX-0.5)*x_step;
            ynew = y[i] + (rand()/RAND_MAX-0.5)*y_step;
            if(xnew > L-0.5) xnew = xnew - L;
            if(ynew > L-0.5) ynew = ynew - L;
            if(xnew < 0.5) xnew += L-1.0;
            if(ynew < 0.5) ynew += L-1.0;
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
        //showmatrix();
        if(k>=sim_time && k%check_time==0){
            double r = 0.0; m = 0;
            while(r<=Rmax){
                P_r[m][counter] = Probability(x,y,r);
                r += r_step;
            }
            counter ++;
            cout<<k<<endl;
        }
    }
    save_prob(P_r);
    time_log << "Czas dzialania programu: " << static_cast<long double>(time (NULL)) - time_start << " s" << endl;

    time_log.close();
    Prob_liq.close();
    matrix.close();

    return 0;
}


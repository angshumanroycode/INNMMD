#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int index(int i,int j,int k,int l,int a,int b,int c,int d){
  int kx=d;
  int jx=c*kx;
  int ix=b*jx;
  int x=i*ix+j*jx+k*kx+l;
  return x;
}

// [[Rcpp::export]]
NumericVector PARAMETERS(int n,int d,int iteration){
  Function medianC("median");
  NumericVector parameters(n-2);
  for(int m=2;m<n;m++){
    double meanmeddist=0;
    for(int itr=0;itr<iteration;itr++){
      NumericVector distvec(n*(n-1)/2);
      for(int k=0;k<d-1;k++){
        int l=0;
        NumericVector rnd=runif(n,0.0,m+0.0);
        NumericVector ind=ceil(rnd);
        for(int i=0;i<n-1;i++){
          for(int j=i+1;j<n;j++){
            int dif=ind[i]-ind[j];
            distvec[l]=distvec[l]+dif*dif;
            l++;
          }
        }
      }
      NumericVector meddist=medianC(distvec);
      meanmeddist=meanmeddist+sqrt(meddist[0]);
    }
    parameters[m-2]=meanmeddist/(iteration+0.0);
  }
  parameters=parameters/sqrt(2);
  return parameters;
}

// [[Rcpp::export]]
List INNMMDPRESET(int n,int d,NumericVector parameters){
  List preset;
  for(int i=0;i<n-2;i++){
    int m=i+2;
    List subpreset;
    NumericVector s1v(m);
    double mp=1/parameters[i];
    for(int j=0;j<m;j++){
      double x=j*mp;
      double y=-x*x/2;
      s1v[j]=exp(y);
    }
    subpreset.push_back(s1v);
    NumericVector s1vc=cumsum(s1v);
    NumericVector s2v(m);
    for(int j=0;j<m;j++)
      s2v[j]=(s1vc[j]+s1vc[m-1-j]-1)/(m+0.0);
    subpreset.push_back(s2v);
    double s3=0;
    for(int j=0;j<m;j++){
      s3=s3+(m-j)*s1v[j];
    }
    s3=(2*s3-m)/(m*m+0.0);
    s3=pow(s3,d-1);
    subpreset.push_back(s3);
    preset.push_back(subpreset);
  }
  return preset;
}

// [[Rcpp::export]]
double MMD(NumericMatrix X,int m,List subpreset){
  int n=X.nrow();
  int d=X.ncol();
  NumericVector s1v=subpreset[0];
  NumericVector s2v=subpreset[1];
  double s3=subpreset[2];
  IntegerVector distvec0=rep(1,n*(n-1)/2);
  NumericVector distvec=as<NumericVector>(distvec0);
  for(int k=0;k<d;k++){
    NumericVector vec=X(_,k);
    int l=0;
    for(int i=0;i<n-1;i++){
      for(int j=i+1;j<n;j++){
        int dif=vec[i]-vec[j];
        if(dif<0)
          dif=-dif;
        distvec[l]=distvec[l]*s1v[dif];
        l++;
      }
    }
  }
  double s1=sum(distvec);
  s1=(2*s1+n)/(n*n+0.0);
  double s2=0;
  for(int i=0;i<n;i++){
    double pd=1;
    for(int j=0;j<d;j++){
      pd=pd*s2v[X(i,j)-1];
    }
    s2=s2+pd;
  }
  s2=s2/(n+0.0);
  double mmdsq=s1-2*s2+s3;
  double mmd=0;
  if(mmdsq>0)
    mmd=sqrt(mmdsq);
  return mmd;
}

// [[Rcpp::export]]
NumericVector INNMMD(List dlist,int n,int d,List preset){
  Function orderC("order");
  NumericVector mats(n*(n-2)*d*(d-1));
  for(int u=0;u<n;u++){
    NumericMatrix R(n-1,d);
    NumericMatrix O(n-1,d);
    for(int v=0;v<d;v++){
      NumericMatrix dlistv=dlist[v];
      NumericVector vec=dlistv(_,u);
      vec.erase(u);
      NumericVector vec2=orderC(vec);
      O(_,v)=vec2-1;
      NumericVector vec3=orderC(vec2);
      R(_,v)=vec3;
    }
    NumericMatrix M(d,d-1);
    for(int v=0;v<d;v++){
      NumericMatrix U(n-1,d-1);
      int wv=0;
      NumericVector vec=O(_,v);
      for(int w=0;w<d;w++)
        if(w!=v){
          NumericVector vec2=R(_,w);
          NumericVector vec3=vec2[vec];
          U(_,wv)=vec3;
          wv++;
        }
        NumericMatrix V=U(Range(0,n-3),_);
        for(int t=1;t<n-2;t++)
          for(int a=t;a<n-2;a++)
            for(int b=0;b<d-1;b++)
              if(U(a,b)>U(t-1,b))
                V(a,b)=V(a,b)-1;
        for(int a=0;a<n-2;a++)
          for(int b=0;b<d-1;b++){
            int ind=index(u,v,a,b,n,d,n-2,d-1);
            mats[ind]=V(a,b);
          }
    }
  }
  NumericVector Ts(d);
  for(int v=0;v<d;v++){
    double Ta=0;
    for(int a=0;a<n-2;a++){
      NumericMatrix Y(n,d-1);
      for(int u=0;u<n;u++)
        for(int b=0;b<d-1;b++){
          int ind=index(u,v,a,b,n,d,n-2,d-1);
          Y(u,b)=mats[ind];
        }
      List subpreset=preset[n-3-a];
      double mmd=MMD(Y,n-a-1,subpreset);
      Ta=Ta+mmd;
    }
    Ts[v]=Ta;
  }
  NumericVector T(2);
  T[0]=sum(Ts);
  T[1]=max(Ts);
  return T;
}

// [[Rcpp::export]]
List INNMMDTEST(List D,int n,int d,List preset,int iteration){
  Function sampleC("sample.int");
  NumericVector tstat=INNMMD(D,n,d,preset);
  NumericMatrix ND(2,iteration);
  for(int itr=0;itr<iteration;itr++){
    List E;
    for(int k=0;k<d;k++){
      NumericMatrix Dk=D[k];
      NumericMatrix Ek(n,n);
      IntegerVector ind=sampleC(n);
      ind=ind-1;
      for(int i=0;i<n-1;i++){
        for(int j=i+1;j<n;j++){
          Ek(i,j)=Dk(ind(i),ind(j));
          Ek(j,i)=Ek(i,j);
        }
      }
      E.push_back(Ek);
    }
    NumericVector ndist=INNMMD(E,n,d,preset);
    ND(_,itr)=ndist;
  }
  List result=List::create(tstat,ND);
  return result;
}
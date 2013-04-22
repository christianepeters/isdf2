\\ author: Christiane Peters, http://christianepeters.wordpress.com
\\ April 2013
\\ run using 'gp isdf2.gp'
log2(x)=log(x*1.)/log(2.)

\\ checking parameters from Section 7 in
\\ Daniel J. Bernstein, Tanja Lange, Christiane Peters.
\\ "Attacking and defending the McEliece cryptosystem"
\\ published at PQCrypto'08

n_=[2046,1632,2048,2960,6624,1744, 2480, 3408, 4624, 6960]
t_=[32,33,27,56,115,35, 45, 67, 95, 119]
w_=[32,34,27,57,117,36, 46, 68, 97, 121]

{for(i=1,10,
  n=n_[i]; k=n-ceil(log2(n))*t_[i]; w=w_[i]; x=floor(k/2);
  mincost=10000000; bestp=0; bestl=0;
  for(p=1,6,\
    Anum=binomial(x,p); Bnum=binomial(k-x,p);
    for(l=1,floor(log(Anum))+10,\
      ops=0.5*(n-k)^2*(n+k)\
          + ((0.5*k-p+1)+(Anum+Bnum))*l\
          + 2*(w-2*p+1)*2*p*Anum*Bnum/2^l;
      prob=Anum*Bnum*binomial(n-k-l,w-2*p)/binomial(n,w);
      cost=log2(ops)-log2(prob);
      if(cost<mincost,
        mincost=cost;
      bestp=p; bestl=l;
    )));
  print("\n["n,",",k,",",w,"]");
  printf("sec level 2^%03.2f\n",mincost);
  printf("%03.2f kilobytes\n",k*(n-k)/1024./8)
  );
}
